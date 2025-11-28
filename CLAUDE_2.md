# CLAUDE_2.md - Detailed Implementation Procedure

## Overview

This document provides step-by-step instructions for implementing personality-specific plotting scripts (A/B/C). It clarifies the exact procedure to follow and common points of confusion.

## Important Clarifications

### What "Generate from Scratch" Means

When the instructions say "generate personality-specific data copies," this means:

1. **Copy the original data files** from `data/figure/` to `data/generated_*/run_{timestamp}/figure/`
2. **Rename the files** according to the personality's naming style
3. **Keep the column names unchanged** inside the CSV files
4. **The plotting script** handles how the data is presented (axis labels, legends, etc.)

**You do NOT need to:**
- Write a script to generate dummy data from scratch
- Rewrite the data generation logic
- Rename columns inside the CSV files

### Focus: The Plotting Script

The main deliverable is the **plotting script** that reflects the personality's coding style. The data file renaming is just organizational - the real personality expression happens in the plotting script through:

- Variable naming
- Code organization
- Axis labels and legends

### CRITICAL: No Comments in Plotting Scripts

**IMPORTANT:** The generated plotting scripts MUST NOT contain any comments.

- ❌ **NO comments** explaining data structure
- ❌ **NO comments** describing what plots are being created
- ❌ **NO comments** explaining code logic or rationale
- ❌ **NO section headers** or comment blocks

**Rationale:** The goal is to infer figure captions purely from the code structure itself:
- Plot type (line plot, dot plot, bar chart, etc.)
- Axis titles and labels
- Legend titles
- Variable names
- Geom choices (geom_point, geom_line, geom_bar, etc.)

The code must be self-explanatory through its structure, not through comments.

## Step-by-Step Procedure

### Phase 1: Setup

1. **Read the personality description** from `data/personality/{A,B,C}.md`
2. **Read the original scripts**:
   - `script/s1_figure.R` - the plotting script to transform
   - `script/s1_make_dummy_datasets.R` - to understand data structure
3. **Create directory structure**:
   ```bash
   timestamp=$(date +"%Y%m%d_%H%M%S")
   mkdir -p "data/generated_{X}/run_${timestamp}/figure"
   mkdir -p "data/generated_{X}/run_${timestamp}/figure_script"
   ```

### Phase 2: Data Files

1. **Copy original data files**:
   ```bash
   cp data/figure/s1_*.csv "data/generated_{X}/run_${timestamp}/figure/"
   ```

2. **Rename files** according to personality naming style:

   **Example for Personality A (descriptive, publication-quality names):**
   ```bash
   mv s1_dev_events.csv neurodevelopmental_timeline_events.csv
   mv s1_cells.csv single_cell_annotations_with_embeddings.csv
   mv s1_taxonomy.csv cell_type_taxonomy_hierarchy.csv
   mv s1_tree_edges.csv taxonomy_dendrogram_edges.csv
   mv s1_tree_nodes.csv taxonomy_dendrogram_nodes.csv
   mv s1_cluster_counts.csv cluster_cell_counts_by_class.csv
   ```

   **Example for Personality B (might be shorter, more casual):**
   ```bash
   mv s1_dev_events.csv dev_events.csv
   mv s1_cells.csv cells_data.csv
   mv s1_taxonomy.csv taxonomy.csv
   # etc.
   ```

   **Example for Personality C (might be very abbreviated):**
   ```bash
   mv s1_dev_events.csv ev.csv
   mv s1_cells.csv c.csv
   mv s1_taxonomy.csv tx.csv
   # etc.
   ```

3. **Do NOT rename columns** inside the CSV files. The original column names (like `event_name`, `class`, `umap1`, etc.) stay the same.

### Phase 3: Create Plotting Script

This is the main task! Create `data/generated_{X}/run_{timestamp}/figure_script/generate_figure_s1_personality_{X}.R`

#### Key Requirements:

1. **File References**: Update `read.csv()` calls to use the new file names:
   ```r
   # Update from:
   read.csv("s1_dev_events.csv")
   # To (for Personality A):
   read.csv("neurodevelopmental_timeline_events.csv")
   ```

2. **Plot Object Names**: MUST keep `p_a`, `p_b`, `p_c`, ..., `p_j` exactly as-is
   - These are "unchangeable anchors"
   - The patchwork composition structure must remain identical

3. **Patchwork Layout**: MUST preserve this exact structure:
   ```r
   (p_a / p_b / (p_c | p_d | p_e | p_f) / (p_g | p_h | p_i | p_j)) +
     plot_layout(heights = c(1,3,1,1), guides = "collect") +
     plot_annotation(tag_levels = "a") &
     theme(plot.tag = element_text(face = 'bold', size = 20))
   ```

4. **Personality-Specific Modifications** (what CAN change):

   **Note:** No comments are allowed in any personality variant. Personality differences are expressed through variable naming, axis labels, and code structure only.

   **For Personality A (Meticulous and Highly Organized):**
   - Descriptive variable names: `developmental_timeline_data`, `single_cell_data`
   - Enhanced axis labels with units: "Developmental Stage (Embryonic day E or Postnatal day P)"
   - Detailed legend titles: "Neurodevelopmental Process", "Assay Modality (scRNA-seq or Multiome)"
   - Well-structured code organization
   - Publication-quality theme customization

   **For Personality B (might be more casual):**
   - Shorter but clear variable names
   - Standard axis labels
   - Moderate legend detail
   - Straightforward code organization

   **For Personality C (might be minimal):**
   - Short variable names
   - Minimal axis labels
   - Brief legend titles
   - Compact code structure

### Phase 4: Validation

1. **Check directory structure** is correct
2. **Verify file names** match personality style
3. **Confirm script references** the renamed files
4. **Ensure plot objects** are named `p_a` through `p_j`
5. **Test run** the script (optional, if R environment is available)

## Common Mistakes to Avoid

❌ **Don't** write a new data generation script from scratch
❌ **Don't** rename columns inside the CSV files
❌ **Don't** change the plot object names (`p_a` through `p_j`)
❌ **Don't** modify the patchwork composition structure
❌ **Don't** create conversion scripts when simple `mv` commands suffice
❌ **Don't** add any comments to the plotting script

✅ **Do** focus on the plotting script as the main deliverable
✅ **Do** use simple file renaming with `mv` commands
✅ **Do** keep original column names in the CSVs
✅ **Do** preserve the exact plot structure and layout
✅ **Do** express personality through variable names, axis labels, and code structure

## File Naming Examples by Personality

### Personality A: Meticulous and Highly Organized
*Descriptive, publication-quality names that fully explain content*

| Original | Personality A |
|----------|--------------|
| `s1_dev_events.csv` | `neurodevelopmental_timeline_events.csv` |
| `s1_cells.csv` | `single_cell_annotations_with_embeddings.csv` |
| `s1_taxonomy.csv` | `cell_type_taxonomy_hierarchy.csv` |
| `s1_tree_edges.csv` | `taxonomy_dendrogram_edges.csv` |
| `s1_tree_nodes.csv` | `taxonomy_dendrogram_nodes.csv` |
| `s1_cluster_counts.csv` | `cluster_cell_counts_by_class.csv` |

### Personality B: (Example - adjust based on personality description)
*Balanced names, clear but concise*

| Original | Personality B (example) |
|----------|--------------|
| `s1_dev_events.csv` | `dev_timeline.csv` |
| `s1_cells.csv` | `cell_data.csv` |
| `s1_taxonomy.csv` | `taxonomy.csv` |
| `s1_tree_edges.csv` | `tree_edges.csv` |
| `s1_tree_nodes.csv` | `tree_nodes.csv` |
| `s1_cluster_counts.csv` | `cluster_counts.csv` |

### Personality C: (Example - adjust based on personality description)
*Minimal names, abbreviated*

| Original | Personality C (example) |
|----------|--------------|
| `s1_dev_events.csv` | `dev.csv` |
| `s1_cells.csv` | `cells.csv` |
| `s1_taxonomy.csv` | `tax.csv` |
| `s1_tree_edges.csv` | `edges.csv` |
| `s1_tree_nodes.csv` | `nodes.csv` |
| `s1_cluster_counts.csv` | `counts.csv` |

## Plotting Script Structure Reference

### Minimal Required Sections

1. **Libraries**: Load required packages
2. **Data Import**: Read the CSV files (using renamed file paths)
3. **Panel Definitions**: Create plots p_a through p_j
4. **Figure Composition**: Combine with patchwork
5. **Export**: Save the figure

### Enhanced Structure (for organized personalities)

1. **Environment Setup**: Libraries, paths, configuration
2. **Data Import**: Read the CSV files using renamed file paths
3. **Panel A**: Developmental timeline
4. **Panel B**: Taxonomy composite (tree + bars)
5. **Panels C-I**: UMAP embeddings with helper function
6. **Panel J**: Constellation plot
7. **Figure Composition**: Patchwork assembly
8. **Export**: Save the figure
9. **Reproducibility**: Session info (optional)

## Critical Constraints (Never Change These)

### Unchangeable Elements:
- Plot object names: **exactly** `p_a`, `p_b`, `p_c`, `p_d`, `p_e`, `p_f`, `p_g`, `p_h`, `p_i`, `p_j`
- Patchwork structure: `(p_a / p_b / (p_c | p_d | p_e | p_f) / (p_g | p_h | p_i | p_j))`
- Height ratios: `c(1, 3, 1, 1)`
- Number of panels: 10 (a through j)
- Panel layout: 4 rows with specified arrangements
- Panel tags: 'a' through 'j'

### Flexible Elements (Personality-Driven):
- File names
- Variable names (except p_a through p_j)
- Axis labels and titles
- Legend formatting
- Theme customization
- Code organization
- Helper function names and structure

**Note:** Comments are NOT allowed in any personality variant.

## Original Column Names Reference

For reference when writing the plotting script, the original CSV files have these columns:

**s1_dev_events.csv / neurodevelopmental_timeline_events.csv:**
- `event_id`, `event_name`, `age_start`, `age_end`, `track`, `x_start`, `x_end`

**s1_cells.csv / single_cell_annotations_with_embeddings.csv:**
- `cell_id`, `modality`, `class`, `subclass`, `cluster_id`, `subcluster_id`, `nt_type`, `age`, `sync_age_bin`, `pseudotime`, `umap1`, `umap2`

**s1_taxonomy.csv / cell_type_taxonomy_hierarchy.csv:**
- `class`, `subclass`, `cluster_id`, `nt_type`

**s1_tree_edges.csv / taxonomy_dendrogram_edges.csv:**
- `from`, `to`

**s1_tree_nodes.csv / taxonomy_dendrogram_nodes.csv:**
- `name`, `class_node`

**s1_cluster_counts.csv / cluster_cell_counts_by_class.csv:**
- `class`, `subclass`, `cluster_id`, `nt_type`, `n_cells`

These column names remain the same across all personalities. Only the file names change.

## Quick Start Template

```bash
#!/bin/bash
# Quick setup script for personality X

# Set personality (A, B, or C)
PERSONALITY="A"

# Create timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Create directories
mkdir -p "data/generated_${PERSONALITY}/run_${TIMESTAMP}/figure"
mkdir -p "data/generated_${PERSONALITY}/run_${TIMESTAMP}/figure_script"

# Copy data files
cp data/figure/s1_*.csv "data/generated_${PERSONALITY}/run_${TIMESTAMP}/figure/"

# Now:
# 1. cd to data/generated_${PERSONALITY}/run_${TIMESTAMP}/figure/
# 2. Rename files with mv commands according to personality
# 3. Create the plotting script in ../figure_script/
```

## Testing the Script

If R is available, test the script:

```bash
cd data/generated_{X}/run_{timestamp}/figure_script
Rscript generate_figure_s1_personality_{X}.R
```

Expected output:
- Console messages (if any)
- Generated figure: `figure/s1_personality_{X}.png`

## Summary Checklist

Before considering implementation complete:

- [ ] Directory structure created with timestamp
- [ ] All 6 data files copied and renamed according to personality
- [ ] Column names inside CSVs remain unchanged
- [ ] Plotting script created with personality-appropriate style
- [ ] Script references the renamed file names
- [ ] Plot objects named exactly p_a through p_j
- [ ] Patchwork composition structure preserved
- [ ] Axis labels and legends reflect personality style
- [ ] Comments and documentation match personality characteristics
- [ ] Variable names (except p_a-p_j) follow personality conventions

## Notes for Each Personality

### Personality A: Meticulous and Highly Organized
- **File names**: Very descriptive, publication-quality
- **Variable names**: Full descriptive names
- **Comments**: Abundant, explaining rationale and methodology
- **Axis labels**: Complete with units and detailed descriptions
- **Documentation**: Comprehensive, including purpose, parameters, returns
- **Theme**: Publication-quality with careful attention to detail

### Personality B: [Read from data/personality/B.md]
- Implement based on specific characteristics in B.md

### Personality C: [Read from data/personality/C.md]
- Implement based on specific characteristics in C.md

## Final Deliverable Structure

```
data/generated_{X}/
└── run_YYYYMMDD_HHMMSS/
    ├── figure/
    │   ├── [renamed_file_1].csv
    │   ├── [renamed_file_2].csv
    │   ├── [renamed_file_3].csv
    │   ├── [renamed_file_4].csv
    │   ├── [renamed_file_5].csv
    │   └── [renamed_file_6].csv
    └── figure_script/
        └── generate_figure_s1_personality_{X}.R
```

The plotting script should be comprehensive and reflect the personality's coding style while maintaining semantic equivalence with the original figure structure.
