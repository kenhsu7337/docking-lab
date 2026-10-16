# Docking lab: AutoDock Vina and PyMOL

A guided lab for generating preliminary docking data for a proposal. Students prepare a protein, dock a known drug back into it, check the result against the crystal structure (RMSD), then dock their own compound the same way.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/OWNER/REPO)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/OWNER/REPO/HEAD?urlpath=lab/tree/notebooks/docking_exercise.ipynb)

Replace `OWNER/REPO` in the two links above with your own (for example `jsmith/docking-lab`). The course page links work without editing.

## What is in this repository

| Path | What it is |
| --- | --- |
| `docs/index.html` | The course page (GitHub Pages): launch buttons, box and config builder, step checklist, results log, methods paragraph generator |
| `notebooks/docking_exercise.ipynb` | The lab itself, in five parts, with a 3D viewer |
| `scripts/dock_utils.py` | Helper functions the notebook calls (one per workflow step) |
| `scripts/smoke_test.py` | One-minute check that the whole pipeline runs |
| `.devcontainer/` | Codespaces setup: conda environment, Vina program, and a PyMOL desktop in the browser |
| `binder/requirements.txt` | Binder setup (pip only, no PyMOL) |
| `environment.yml` | The conda environment used by Codespaces |

## Set it up (about 10 minutes)

1. Create a new repository on GitHub and upload everything in this folder, keeping the structure. Use the **Add file > Upload files** button, or `git push` if you prefer.
2. Turn on the course page: **Settings > Pages > Build and deployment**, set **Source** to **Deploy from a branch**, choose your main branch and the **/docs** folder, and save. The page appears at `https://OWNER.github.io/REPO/` after a minute or two. Its launch buttons find the repository from the address, so they need no editing.
3. Open the page, press **Open in Codespaces**, and run the notebook once yourself before class (see below).
4. Give students the page address. If they use Binder, open the Binder link once first so the image is built before class, otherwise the first student waits for the build.

## Two ways for students to run it

| | Codespaces | Binder |
| --- | --- | --- |
| Account | GitHub account | None |
| PyMOL | Yes, in a browser desktop (open port 6080, password `vscode`, then run `DISPLAY=:1 pymol &`) | No, the notebook uses a built-in viewer |
| Speed | Dedicated 4 cores | Shared, so slower |
| Files | Kept in the codespace | Deleted when the session times out |

Codespaces has a free monthly allowance, and GitHub Education accounts get more. Check GitHub's current limits before class. A student who runs out can switch to Binder.

## Test it before class

In Codespaces, open a terminal and run:

```bash
python scripts/smoke_test.py
```

It downloads the example files that ship with AutoDock Vina (Abl kinase with imatinib), prepares them, docks at low exhaustiveness, and checks the scoring and RMSD code. Then run the notebook top to bottom with the defaults (HIV-1 protease, PDB 1HSG, ligand MK1) and note the affinity and RMSD for the control so you know what students should see. The same check runs on GitHub under **Actions > Smoke test**.

## Changing the teaching system

Edit the settings cell at the top of the notebook: `PDB_ID`, `LIGAND_RESNAME`, `EXHAUSTIVENESS`, and `PADDING`. Any PDB entry with a drug-like ligand works. For a very flexible ligand, keep exhaustiveness at 32 or higher. For a classroom without internet access to RCSB, put the files in `notebooks/work/` and set `LOCAL_COMPLEX` and `LOCAL_LIGAND_SDF`.

## Known limits

- The pipeline was tested with Vina 1.2.7, Meeko 0.8 and RDKit 2024+. The Meeko command-line options have changed between releases, so pin versions in `environment.yml` and `binder/requirements.txt` if you want a fixed behavior for the term.
- Docking scores are estimates used to prioritize hypotheses, not measured affinities. The course page and notebook say this to students.
- The Codespaces desktop and the ideal-ligand download from RCSB depend on services outside this repository. If the ligand download fails, the notebook says which file to download by hand.
