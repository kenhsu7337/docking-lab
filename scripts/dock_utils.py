"""Small helpers used by the docking notebook.

Every function does one step of the workflow and prints what it did, so
students can read the code and see how each step maps to the handout.
Tested with Vina 1.2.7 (Python bindings), Meeko 0.8 and RDKit 2024+.
"""
import subprocess
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

RCSB_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_LIGAND_URLS = [
    "https://files.rcsb.org/ligands/download/{name}_ideal.sdf",
    "https://files.rcsb.org/ligands/view/{name}_ideal.sdf",
]


# ---------------------------------------------------------------- download
def _download(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": "docking-lab"})
    with urllib.request.urlopen(req, timeout=60) as r:
        Path(dest).write_bytes(r.read())


def fetch_pdb(pdb_id, dest=None):
    """Download a PDB entry from RCSB and return the local path."""
    pdb_id = pdb_id.upper()
    dest = Path(dest or f"{pdb_id}.pdb")
    _download(RCSB_PDB.format(pdb_id=pdb_id), dest)
    print(f"Downloaded {pdb_id} -> {dest}")
    return dest


def fetch_ligand_sdf(resname, dest=None):
    """Download the idealized 3D SDF for a ligand residue name (e.g. MK1)."""
    resname = resname.upper()
    dest = Path(dest or f"{resname}_ideal.sdf")
    for url in RCSB_LIGAND_URLS:
        try:
            _download(url.format(name=resname), dest)
            print(f"Downloaded ideal coordinates for {resname} -> {dest}")
            return dest
        except Exception:
            continue
    raise RuntimeError(
        f"Could not download {resname}_ideal.sdf automatically. "
        f"Download it from the RCSB ligand page for {resname} and save it as {dest}."
    )


# ------------------------------------------------------- split the complex
def split_complex(pdb_path, ligand_resname, receptor_out="receptor.pdb",
                  ligand_out="ligand_crystal.pdb", chains=None):
    """Write the protein and the crystal ligand to separate files.

    Keeps protein ATOM records (alternate location A only), keeps the first
    copy of the ligand, and drops waters, ions and other heteroatoms.
    """
    receptor_lines, ligand_lines, dropped = [], [], Counter()
    ligand_key = None
    for line in Path(pdb_path).read_text().splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        line = line[:16] + " " + line[17:]
        resname, chain = line[17:20].strip(), line[21]
        if chains and chain not in chains:
            continue
        if rec == "ATOM":
            receptor_lines.append(line)
        elif resname == ligand_resname.upper():
            key = (chain, line[22:27])
            ligand_key = ligand_key or key
            if key == ligand_key:
                ligand_lines.append(line)
        else:
            dropped[resname] += 1
    if not ligand_lines:
        raise ValueError(f"No ligand named {ligand_resname} found in {pdb_path}")
    Path(receptor_out).write_text("\n".join(receptor_lines) + "\nEND\n")
    Path(ligand_out).write_text("\n".join(ligand_lines) + "\nEND\n")
    print(f"Receptor: {len(receptor_lines)} atoms -> {receptor_out}")
    print(f"Ligand {ligand_resname}: {len(ligand_lines)} atoms -> {ligand_out}")
    if dropped:
        summary = ", ".join(f"{k} x{v}" for k, v in dropped.most_common())
        print(f"Dropped heteroatoms (atom counts): {summary}")
        print("Check this list: a metal or cofactor near the site may need to stay.")
    return Path(receptor_out), Path(ligand_out)


# ------------------------------------------------------------------ prep
def _run(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def prepare_receptor(receptor_pdb, basename="receptor"):
    """Convert a receptor PDB to PDBQT with Meeko (adds hydrogens and charges)."""
    _run(["mk_prepare_receptor.py", "--read_pdb", str(receptor_pdb),
          "-o", basename, "-p"])
    out = Path(f"{basename}.pdbqt")
    print(f"Receptor PDBQT -> {out}")
    return out


def ligand_from_smiles(smiles, out_sdf="ligand.sdf", seed=42):
    """Make a 3D ligand from SMILES. Write the SMILES in the protonation
    state you want at pH 7.4 (for example, charged amines and carboxylates)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise RuntimeError("RDKit could not generate 3D coordinates for this SMILES")
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, str(out_sdf))
    print(f"3D ligand -> {out_sdf}")
    return Path(out_sdf)


def prepare_ligand(ligand_sdf, out_pdbqt="ligand.pdbqt"):
    """Convert a 3D ligand SDF to PDBQT with Meeko and report rotatable bonds."""
    _run(["mk_prepare_ligand.py", "-i", str(ligand_sdf), "-o", str(out_pdbqt)])
    lines = Path(out_pdbqt).read_text().splitlines()
    n_branch = sum(1 for l in lines if l.startswith("BRANCH"))
    print(f"Ligand PDBQT -> {out_pdbqt} ({n_branch} rotatable bonds)")
    return Path(out_pdbqt)


# ------------------------------------------------------------------- box
def read_coords(path):
    """Read x, y, z of every atom from a PDB, PDBQT or SDF/MOL file."""
    path = Path(path)
    if path.suffix.lower() in (".sdf", ".mol"):
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
        return mol.GetConformer().GetPositions()
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def box_from_ligand(path, padding=8.0):
    """Box centered on the ligand, padded by `padding` Å on every side."""
    xyz = read_coords(path)
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    center = ((lo + hi) / 2).round(2)
    size = (hi - lo + 2 * padding).round(1)
    print(f"center = {tuple(float(c) for c in center)}  "
          f"size = {tuple(float(s) for s in size)}")
    if size.max() > 30:
        print("Note: a side is longer than 30 Å. Keep exhaustiveness at 32 or "
              "higher, or reduce the padding.")
    return tuple(float(c) for c in center), tuple(float(s) for s in size)


# ------------------------------------------------------------------ dock
def dock(receptor_pdbqt, ligand_pdbqt, center, size, out_pdbqt,
         exhaustiveness=32, n_poses=9, seed=42, cpu=0):
    """Run AutoDock Vina through its Python bindings. Returns affinities."""
    from vina import Vina
    v = Vina(sf_name="vina", cpu=cpu, seed=seed, verbosity=1)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(str(out_pdbqt), n_poses=n_poses, overwrite=True)
    affinities = read_affinities(out_pdbqt)
    print(f"{len(affinities)} poses written. Affinities (kcal/mol):", affinities)
    return affinities


def read_affinities(poses_pdbqt):
    """Affinity (kcal/mol) of every pose in a Vina output file, best first.
    Vina only writes poses within 3 kcal/mol of the best one by default."""
    values = []
    for line in Path(poses_pdbqt).read_text().splitlines():
        if line.startswith("REMARK VINA RESULT"):
            values.append(float(line.split()[3]))
    return values


def poses_to_sdf(poses_pdbqt, out_sdf):
    """Convert Vina output to SDF so bond orders are kept."""
    _run(["mk_export.py", str(poses_pdbqt), "-s", str(out_sdf)])
    return Path(out_sdf)


# ---------------------------------------------------------------- RMSD
def rmsd_to_crystal(poses_sdf, crystal_file, template_sdf=None):
    """Heavy-atom RMSD of each pose to the crystal ligand, no re-alignment,
    symmetry aware. If the crystal ligand is a PDB file (no bond orders),
    pass the ideal SDF as `template_sdf`."""
    crystal_file = Path(crystal_file)
    if crystal_file.suffix.lower() in (".sdf", ".mol"):
        crystal = Chem.MolFromMolFile(str(crystal_file), removeHs=False)
    else:
        crystal = Chem.MolFromPDBFile(str(crystal_file))
        template = Chem.RemoveHs(Chem.MolFromMolFile(str(template_sdf)))
        crystal = AllChem.AssignBondOrdersFromTemplate(template, crystal)
    crystal = Chem.RemoveHs(crystal)
    values = []
    for pose in Chem.SDMolSupplier(str(poses_sdf), removeHs=False):
        values.append(round(rdMolAlign.CalcRMS(Chem.RemoveHs(pose), crystal), 2))
    for i, v in enumerate(values, 1):
        print(f"pose {i}: RMSD = {v:.2f} Å")
    return values


# ------------------------------------------------------------- contacts
def _protein_atoms(receptor_pdb):
    atoms = []
    for line in Path(receptor_pdb).read_text().splitlines():
        if line.startswith("ATOM"):
            atoms.append({
                "res": f"{line[17:20].strip()}{line[22:26].strip()}:{line[21]}",
                "name": line[12:16].strip(),
                "elem": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
                "xyz": np.array([float(line[30:38]), float(line[38:46]),
                                 float(line[46:54])]),
            })
    return atoms


def contacts(receptor_pdb, poses_sdf, pose_index=0, cutoff=4.0, hbond=3.3):
    """Residues within `cutoff` Å of one pose, plus polar (N/O) contacts
    shorter than `hbond` Å. Distances are geometric only: they flag
    candidate hydrogen bonds, they do not prove them."""
    pose = Chem.RemoveHs(Chem.SDMolSupplier(str(poses_sdf), removeHs=False)[pose_index])
    lig_xyz = pose.GetConformer().GetPositions()
    lig_polar = [a.GetIdx() for a in pose.GetAtoms() if a.GetSymbol() in ("N", "O")]
    near, polar = {}, []
    for atom in _protein_atoms(receptor_pdb):
        d = np.linalg.norm(lig_xyz - atom["xyz"], axis=1)
        if d.min() <= cutoff:
            near[atom["res"]] = min(near.get(atom["res"], 99.0), float(d.min()))
        if atom["elem"] in ("N", "O"):
            for i in lig_polar:
                if d[i] <= hbond:
                    polar.append((atom["res"], atom["name"],
                                  pose.GetAtomWithIdx(i).GetSymbol(), round(float(d[i]), 2)))
    near = dict(sorted(near.items(), key=lambda kv: kv[1]))
    return near, sorted(polar, key=lambda t: t[3])
