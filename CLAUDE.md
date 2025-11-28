# CLAUDE.md

## Purpose

This project evaluates how different "coding personalities" (A/B/C) influence the quality of automatically generated figure captions from plotting code.

Your tasks as Claude are:

1. **Understand the dummy data and figure structure** by reading `script/s1_make_dummy_datasets.R`.
2. **Generate personality-specific data copies** under `data/generated_*/run_{timestamp}/figure/`, including renaming files and variables according to each personality.
3. **Generate personality-specific plotting scripts** under `data/generated_*/run_{timestamp}/figure_script/` that:
   - Import the personality-specific data files you created in step 2.
   - Recreate the same figure structure as `script/s1_figure.R` (same panels and layout) while reflecting the personality’s code style.

You must keep the *semantic structure* of the figure identical across personalities. Only naming style, code clarity, and annotation quality should differ.


## Repository structure (reference)

You can assume the repository structure is:

- `CLAUDE.md` (this file)
- `LICENSE`
- `README.md`
- `data/`
  - `captions/`
    - ...
  - `figure/`
    - `s1_cells.csv`
    - `s1_cluster_counts.csv`
    - `s1_dev_events.csv`
    - `s1_taxonomy.csv`
    - `s1_tree_edges.csv`
    - `s1_tree_nodes.csv`
  - `generated_A/`
  - `generated_B/`
  - `generated_C/`
  - `papers/`
    - ...
  - `personality/`
    - `A.md`
    - `B.md`
    - `C.md`
  - `prompts/`
    - ...
- `figure/`
  - ...
- `logs/`
  - ...
- `script/`
  - `s1_figure.R`
  - `s1_make_dummy_datasets.R`


## Personality descriptions

The personalities for A/B/C are described in:

- `data/personality/A.md`
- `data/personality/B.md`
- `data/personality/C.md`

Each file describes *how that person tends to code and name things*, including:

- Style of file/variable names,
- Degree of clarity and documentation,
- How carefully they label axes, legends, etc.

You must:

- Read the relevant personality file before generating any data or scripts for that type.
- Use the personality description to decide:
  - How to rename data files,
  - How to rename variables and objects inside the plotting script,
  - How verbose or clear the plotting code and labels should be.