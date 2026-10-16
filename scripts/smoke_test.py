"""Check that the whole pipeline works, using the example files that ship
with AutoDock Vina (Abl kinase with imatinib, PDB 1IEP).

Run:  python scripts/smoke_test.py
Takes about a minute. It uses a low exhaustiveness, so it tests that every
step runs, not that the docking is accurate.
"""
import sys
import tempfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dock_utils as du

BASE = ("https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/"
        "develop/example/basic_docking/data/")


def main():
    work = Path(tempfile.mkdtemp(prefix="dock_smoke_"))
    for name in ("1iep_receptorH.pdb", "1iep_ligand.sdf"):
        urllib.request.urlretrieve(BASE + name, work / name)

    # The example receptor starts with one stray atom line; drop it.
    atoms = [l for l in (work / "1iep_receptorH.pdb").read_text().splitlines()
             if l.startswith("ATOM")][1:]
    (work / "receptor.pdb").write_text("\n".join(atoms) + "\nEND\n")

    import os
    os.chdir(work)
    receptor = du.prepare_receptor("receptor.pdb", "receptor")
    ligand = du.prepare_ligand("1iep_ligand.sdf", "ligand.pdbqt")
    center, size = du.box_from_ligand("1iep_ligand.sdf", padding=4)
    affinities = du.dock(receptor, ligand, center, size, "out.pdbqt",
                         exhaustiveness=2)
    sdf = du.poses_to_sdf("out.pdbqt", "poses.sdf")
    rmsd = du.rmsd_to_crystal(sdf, "1iep_ligand.sdf")
    near, polar = du.contacts("receptor.pdb", sdf, 0)

    assert affinities and affinities[0] < -5, "no sensible docking score"
    assert len(rmsd) == len(affinities), "pose count mismatch"
    assert near, "no contacts found"
    print(f"\nSmoke test passed in {work}")


if __name__ == "__main__":
    main()
