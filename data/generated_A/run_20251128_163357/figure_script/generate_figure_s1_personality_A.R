## ============================================================================
## FIGURE S1 GENERATION - PERSONALITY A (METICULOUS AND HIGHLY ORGANIZED)
## ============================================================================
##
## Purpose: Generate a comprehensive multi-panel figure displaying cell type
##          taxonomy, developmental dynamics, and spatial organization of
##          single-cell transcriptomic data from mouse cortical development.
##
## Author: Personality A (Meticulous and Highly Organized)
## Date: 2025-11-28
## R Version: 4.x or higher
##
## Input data files (all located in ../figure/ directory):
##   - s1_dev_events.csv: Developmental timeline showing key neurodevelopmental
##                        processes across embryonic and postnatal stages
##   - s1_cells.csv: Single-cell data with cell type annotations, developmental
##                   age, assay modality, and UMAP embedding coordinates
##   - s1_taxonomy.csv: Hierarchical cell type taxonomy (class/subclass/cluster)
##   - s1_tree_edges.csv: Edge list defining the hierarchical relationships
##                        between cell classes, subclasses, and clusters
##   - s1_tree_nodes.csv: Node attributes for the taxonomy dendrogram including
##                        class membership for color mapping
##   - s1_cluster_counts.csv: Aggregated cell counts per cluster for visualization
##
## Output:
##   - figure/s1_personality_A.png: Publication-quality multi-panel figure
##                                  (300 DPI, 8" × 8")
##
## Figure structure (10 panels arranged in 4 rows):
##   Row 1: Panel a - Developmental timeline of neurodevelopmental processes
##   Row 2: Panel b - Composite panel with taxonomy dendrogram and cell counts
##   Row 3: Panels c-f - UMAP embeddings colored by hierarchical cell annotations
##   Row 4: Panels g-j - UMAP embeddings colored by developmental variables
##
## ============================================================================

## ----------------------------------------------------------------------------
## SECTION 0: ENVIRONMENT SETUP AND CONFIGURATION
## ----------------------------------------------------------------------------
## Load all required R packages for data manipulation, visualization, and
## graph-based analyses. Ensure reproducibility by documenting package versions.
## ----------------------------------------------------------------------------

# Core visualization and data manipulation libraries
library(ggplot2)   # Grammar of graphics for all plot generation
library(dplyr)     # Data frame manipulation and transformation
library(igraph)    # Graph object creation for taxonomy dendrogram
library(ggraph)    # ggplot2-based grammar for network/tree visualization
library(patchwork) # Multi-panel figure composition with flexible layouts

# Note: Ensure all packages are up-to-date for best compatibility
# Recommended versions: ggplot2 >= 3.4.0, patchwork >= 1.1.0

## ----------------------------------------------------------------------------
## SECTION 1: DATA IMPORT AND PREPARATION
## ----------------------------------------------------------------------------
## Load all input CSV files containing the experimental data and metadata.
## File paths are relative to the script location for portability.
## ----------------------------------------------------------------------------

# Define relative path to data directory (assumes script is in figure_script/)
data_directory_path <- "../figure"

# 1.1 Load developmental timeline data
# Contains information about when key neurodevelopmental processes occur
# across embryonic (E) and postnatal (P) developmental stages
developmental_timeline_data <- read.csv(file.path(data_directory_path, "neurodevelopmental_timeline_events.csv"),
                                        stringsAsFactors = FALSE)

# 1.2 Load single-cell annotation and embedding data
# Contains 500 cells with cell type annotations, developmental age,
# assay modality, and 2D UMAP coordinates for visualization
single_cell_data <- read.csv(file.path(data_directory_path, "single_cell_annotations_with_embeddings.csv"),
                             stringsAsFactors = FALSE)

# 1.3 Load cell type taxonomy reference table
# Defines the hierarchical organization: class → subclass → cluster
# Also includes neurotransmitter type for each taxonomic level
cell_type_taxonomy <- read.csv(file.path(data_directory_path, "cell_type_taxonomy_hierarchy.csv"),
                               stringsAsFactors = FALSE)

# 1.4 Load taxonomy dendrogram structure (edges)
# Edge list defining parent-child relationships in the taxonomy hierarchy
# Used to construct the dendrogram visualization
taxonomy_dendrogram_edges <- read.csv(file.path(data_directory_path, "taxonomy_dendrogram_edges.csv"),
                                      stringsAsFactors = FALSE)

# 1.5 Load taxonomy dendrogram structure (nodes)
# Node attributes including class membership for color mapping in the dendrogram
taxonomy_dendrogram_nodes <- read.csv(file.path(data_directory_path, "taxonomy_dendrogram_nodes.csv"),
                                      stringsAsFactors = FALSE)

# 1.6 Load cluster-level cell count summary
# Aggregated cell counts per cluster, used for generating the bar plot
# showing the distribution of cells across clusters
cluster_cell_count_summary <- read.csv(file.path(data_directory_path, "cluster_cell_counts_by_class.csv"),
                                       stringsAsFactors = FALSE)

## ----------------------------------------------------------------------------
## SECTION 2: PANEL A - DEVELOPMENTAL TIMELINE VISUALIZATION
## ----------------------------------------------------------------------------
## This panel visualizes the temporal dynamics of key neurodevelopmental
## processes across embryonic and postnatal stages. Each process is represented
## as a horizontal bar spanning its active period, arranged in three tracks
## to minimize visual overlap.
##
## Data visualization approach:
## - X-axis: Developmental stages from E10.5 (embryonic day 10.5) to P56 (postnatal day 56)
## - Y-axis: Track number (reversed to show earliest events at top)
## - Bar length: Duration of each developmental process
## - Color: Distinct color for each neurodevelopmental process
##
## Visual design rationale:
## - Horizontal timeline provides intuitive left-to-right temporal progression
## - Multiple tracks allow overlapping processes to be displayed without occlusion
## - Color coding enables quick identification of specific processes
## - Round line ends provide polished, publication-quality appearance
## ----------------------------------------------------------------------------

# Define the ordered sequence of developmental stages for x-axis labeling
# Includes both embryonic (E) and postnatal (P) timepoints
developmental_stage_levels <- c("E10.5", "E12.5", "E13.5", "E14.5", "E15.5",
                                "E16.5", "E17", "E18.5", "P0", "P1", "P5",
                                "P10", "P15", "P20", "P25", "P28", "P56")

# Create Panel A: Developmental timeline showing when key processes occur
p_a <- ggplot(developmental_timeline_data) +
  # Draw horizontal bars representing process duration
  geom_segment(
    mapping = aes(
      x = x_start,              # Start position (numeric index of developmental stage)
      xend = x_end,             # End position (numeric index of developmental stage)
      y = track,                # Track number (1-3 for vertical positioning)
      yend = track,             # Same y for horizontal line
      color = event_name        # Color by process name for legend
    ),
    linewidth = 2,              # Thick lines for visibility
    lineend = "round"           # Rounded ends for polished appearance
  ) +
  # Add points at the start of each process to emphasize initiation time
  geom_point(
    mapping = aes(x = x_start, y = track, color = event_name),
    size = 2                    # Size coordinated with line width
  ) +
  # Configure x-axis to display developmental stage labels
  scale_x_continuous(
    breaks = seq_along(developmental_stage_levels),  # Position breaks at each stage
    labels = developmental_stage_levels,             # Use stage names as labels
    expand = c(0.01, 0.01)                          # Minimal padding for clean appearance
  ) +
  # Reverse y-axis so earliest events (track 1) appear at top
  scale_y_reverse(breaks = c(1, 2, 3)) +
  # Apply descriptive, publication-quality axis labels and legend title
  labs(
    x = "Developmental Stage (Embryonic day E or Postnatal day P)",
    y = NULL,                   # No y-axis label needed (tracks are arbitrary positions)
    color = "Neurodevelopmental\nProcess"  # Line break for better legend formatting
  ) +
  # Apply clean, minimal theme suitable for publication
  theme_minimal(base_size = 10) +
  # Customize theme elements for publication quality
  theme(
    axis.text.y = element_blank(),           # Hide y-axis text (tracks are not quantitative)
    panel.grid.minor = element_blank(),      # Remove minor gridlines for cleaner appearance
    plot.title = element_text(face = "bold", hjust = 0),  # Bold title, left-aligned
    legend.title = element_text(face = "bold")           # Bold legend title for emphasis
  )

## ----------------------------------------------------------------------------
## SECTION 3: PANEL B - TAXONOMY HIERARCHY (COMPOSITE THREE-PART PANEL)
## ----------------------------------------------------------------------------
## This composite panel consists of three vertically stacked sub-panels that
## together illustrate the hierarchical organization of cell types and their
## abundance in the dataset:
##
##   1. Top: Dendrogram showing class → subclass → cluster hierarchy
##   2. Middle: Bar plot showing cell count per cluster (ordered by count)
##   3. Bottom: Stacked bar plot showing cell count per class, split by modality
##
## The three sub-panels are aligned to provide a cohesive view of the taxonomy
## structure and cell distribution patterns.
## ----------------------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 3.1: Construct taxonomy dendrogram (tree structure)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This dendrogram visualizes the hierarchical relationships between cell
## classes (broad categories), subclasses (refined categories), and clusters
## (fine-grained cell types). The tree is laid out using the "dendrogram"
## algorithm to emphasize the hierarchical structure.
##
## Graph structure:
## - Nodes: Cell classes, subclasses, and cluster IDs
## - Edges: Parent-child relationships (directed from broad to specific)
## - Node colors: Colored by the parent cell class for visual grouping
## - Layout: Dendrogram (hierarchical tree) with leaves at the bottom
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Build directed graph object from edge list and node attributes
# This creates an igraph object suitable for visualization with ggraph
taxonomy_hierarchy_graph <- graph_from_data_frame(
  d = taxonomy_dendrogram_edges,         # Edge list (from → to relationships)
  vertices = taxonomy_dendrogram_nodes,  # Node attributes (including class colors)
  directed = TRUE                        # Directed edges (parent → child)
)

# Create dendrogram visualization with class-colored nodes
p_tree <- ggraph(taxonomy_hierarchy_graph, layout = "dendrogram") +
  # Draw edges connecting hierarchical levels with diagonal lines
  geom_edge_diagonal(alpha = 0.4, color = "gray30") +
  # Draw nodes as points, colored by their parent cell class
  geom_node_point(aes(color = class_node), size = 2) +
  # Add text labels for each node (class/subclass/cluster names)
  geom_node_text(
    aes(label = name),
    hjust = -0.1,              # Slight right offset to avoid overlap with nodes
    size = 2,                  # Small text size to fit many labels
    check_overlap = TRUE       # Automatically hide overlapping labels
  ) +
  # Allow labels to extend beyond plot area (requires extra right margin)
  coord_cartesian(clip = "off") +
  # Use void theme (no axes) as dendrogram has no quantitative axes
  theme_void(base_size = 9) +
  # Customize theme for publication quality
  theme(
    legend.position = "none",                              # Hide legend (colors are self-explanatory)
    plot.margin = margin(5.5, 80, 5.5, 5.5)               # Extra right margin for node labels
  )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 3.2: Create cluster-level cell count bar plot
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This bar plot shows the number of cells in each cluster, ordered by count.
## Bars are colored by cell class to show which classes dominate the dataset.
## The x-axis aligns with the dendrogram leaves above for visual consistency.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p_bar <- ggplot(
  data = cluster_cell_count_summary,
  mapping = aes(
    x = reorder(cluster_id, n_cells),  # Order clusters by cell count (ascending)
    y = n_cells,                       # Bar height = cell count
    fill = class                       # Fill color by cell class
  )
) +
  geom_col() +  # Bar plot (geom_col for pre-computed heights)
  labs(
    x = NULL,                   # No x-axis label (aligns with dendrogram above)
    y = "Cell Count per Cluster",
    fill = "Cell Class"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_blank(),     # Hide x-axis text (aligns with tree above)
    axis.ticks.x = element_blank(),    # Hide x-axis ticks for cleaner appearance
    legend.title = element_text(face = "bold")  # Bold legend title
  )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 3.3: Create class-level cell count bar plot (split by modality)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This stacked bar plot shows the total number of cells per class, with bars
## subdivided by assay modality (scRNA-seq vs. Multiome). This provides insight
## into both the overall class distribution and the contribution of each
## sequencing technology.
##
## Data aggregation: Group single-cell data by class and modality, then count
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Aggregate single-cell data to count cells per class and modality
class_modality_cell_counts <- single_cell_data %>%
  count(class, modality)  # dplyr::count() groups and counts in one step

# Create stacked bar plot showing modality contribution to each class
p_b_class_mod <- ggplot(
  data = class_modality_cell_counts,
  mapping = aes(
    x = reorder(class, n),  # Order classes by total cell count (ascending)
    y = n,                  # Bar height = cell count
    fill = modality         # Stack bars by assay modality
  )
) +
  geom_col(position = "stack") +  # Stacked bars (not side-by-side)
  labs(
    x = "Cell Class (ordered by total cell count)",
    y = "Number of Cells",
    fill = "Assay Modality\n(scRNA-seq or Multiome)"  # Detailed legend title
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Angled labels to prevent overlap
    legend.title = element_text(face = "bold")  # Bold legend title
  )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 3.4: Compose Panel B from three sub-panels
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Stack the three sub-panels vertically using patchwork, with custom height
## ratios to allocate appropriate space to each component. The dendrogram
## receives the most space (ratio 2), while the bar plots are smaller (ratio 1).
##
## Height ratios: 2:1:1 (dendrogram:cluster_bars:class_modality_bars)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compose the three sub-panels using patchwork's vertical stacking operator (/)
p_b <- (p_tree / p_bar / p_b_class_mod) +
  plot_layout(heights = c(2, 1, 1))  # Allocate more space to dendrogram

# Wrap the composite panel to treat it as a single unit in the final layout
p_b <- wrap_elements(full = p_b)

## ----------------------------------------------------------------------------
## SECTION 4: PANELS C-I - UMAP EMBEDDING VISUALIZATIONS
## ----------------------------------------------------------------------------
## These seven panels show the same UMAP embedding (2D projection of high-
## dimensional single-cell data) colored by different metadata variables:
##
##   Panel C: Cell class (broad categories, 12 classes)
##   Panel D: Cell subclass (refined categories, ~24 subclasses)
##   Panel E: Cluster ID (fine-grained clusters, ~48 clusters)
##   Panel F: Subcluster ID (finest resolution, ~144 subclusters)
##   Panel G: Developmental age (embryonic/postnatal stage when cell was sampled)
##   Panel H: Synchronized age bin (numeric age grouping for analysis)
##   Panel I: Inferred pseudotime (computational estimate of developmental progression)
##
## Design rationale:
## - All panels use the same UMAP coordinates for direct visual comparison
## - Only the coloring variable changes to highlight different aspects
## - Helper function ensures consistent styling across all panels
## - Panel C shows legend (reference for cell classes); others hide legend for space
##
## UMAP (Uniform Manifold Approximation and Projection):
## - Dimensionality reduction technique that preserves local structure
## - Projects high-dimensional gene expression data into 2D for visualization
## - Cells with similar transcriptomes cluster together in UMAP space
## ----------------------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 4.1: Define helper function for UMAP visualization
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This function generates standardized UMAP scatter plots to ensure visual
## consistency across all panels. It accepts the data frame, the name of the
## variable to use for coloring, and whether to show the legend.
##
## Function parameters:
##   - data_for_plotting: Data frame containing UMAP coordinates and metadata
##   - coloring_variable: String name of column to map to point color
##   - display_legend: Logical flag to show/hide the color legend
##
## Returns: A ggplot object with standardized UMAP visualization
##
## Styling specifications:
##   - Point size: 0.7 (small to reduce overplotting with 500 cells)
##   - Point alpha: 0.8 (slight transparency to show overlapping points)
##   - Theme: Minimal with no axis text (UMAP coordinates are relative, not absolute)
##   - Base font size: 9pt for consistency with other panels
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create_standardized_umap_plot <- function(data_for_plotting,
                                          coloring_variable,
                                          display_legend = TRUE) {
  # Convert string column name to rlang symbol for non-standard evaluation
  # This allows the column name to be passed as a string and used in aes()
  color_aesthetic_symbol <- rlang::sym(coloring_variable)

  # Build UMAP scatter plot with specified coloring variable
  umap_plot <- ggplot(
    data = data_for_plotting,
    mapping = aes(
      x = umap1,                      # UMAP dimension 1 (horizontal axis)
      y = umap2,                      # UMAP dimension 2 (vertical axis)
      color = !!color_aesthetic_symbol  # Unquote symbol to evaluate column name
    )
  ) +
    geom_point(size = 0.7, alpha = 0.8) +  # Small, semi-transparent points
    labs(
      x = "UMAP Dimension 1",         # Descriptive axis label
      y = "UMAP Dimension 2",         # Descriptive axis label
      color = NULL                    # Legend title will be inferred from variable
    ) +
    theme_minimal(base_size = 9) +    # Clean minimal theme
    theme(
      axis.text = element_blank(),    # Hide axis text (UMAP coordinates are relative)
      axis.ticks = element_blank(),   # Hide axis ticks (not quantitatively meaningful)
      plot.title = element_text(face = "bold", hjust = 0),  # Bold title if added
      legend.position = if (display_legend) "right" else "none"  # Conditional legend
    )

  return(umap_plot)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 4.2: Generate individual UMAP panels using helper function
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Each panel visualizes the same UMAP embedding but colored by a different
## metadata variable to reveal different biological patterns:
##
## - Hierarchical annotations (class/subclass/cluster/subcluster): Show cell
##   type organization and transcriptomic similarity structure
## - Developmental variables (age/age_bin/pseudotime): Show temporal patterns
##   and developmental trajectories
##
## Note: Only Panel C displays a legend (for cell class reference). Other
## panels hide legends to conserve space in the multi-panel layout.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Panel C: UMAP colored by cell class (broad categories)
# This panel provides the reference legend for cell class colors used throughout
p_c <- create_standardized_umap_plot(
  data_for_plotting = single_cell_data,
  coloring_variable = "class",
  display_legend = TRUE  # Show legend to provide color reference
)

# Panel D: UMAP colored by cell subclass (refined categories)
p_d <- create_standardized_umap_plot(
  data_for_plotting = single_cell_data,
  coloring_variable = "subclass",
  display_legend = FALSE  # Hide legend to save space
)

# Panel E: UMAP colored by cluster ID (fine-grained cell types)
p_e <- create_standardized_umap_plot(
  data_for_plotting = single_cell_data,
  coloring_variable = "cluster_id",
  display_legend = FALSE  # Hide legend (too many clusters for readable legend)
)

# Panel F: UMAP colored by subcluster ID (finest resolution)
p_f <- create_standardized_umap_plot(
  data_for_plotting = single_cell_data,
  coloring_variable = "subcluster_id",
  display_legend = FALSE  # Hide legend (too many subclusters for readable legend)
)

# Panel G: UMAP colored by developmental age (embryonic/postnatal stage)
# Reveals how cells from different developmental timepoints are distributed
p_g <- create_standardized_umap_plot(
  data_for_plotting = single_cell_data,
  coloring_variable = "age",
  display_legend = FALSE  # Hide legend to save space
)

# Panel H: UMAP colored by synchronized age bin (numeric age grouping)
# Similar to Panel G but with numeric binning for computational analyses
p_h <- create_standardized_umap_plot(
  data_for_plotting = single_cell_data,
  coloring_variable = "sync_age_bin",
  display_legend = FALSE  # Hide legend to save space
)

# Panel I: UMAP colored by inferred pseudotime (developmental trajectory)
# Shows computationally inferred developmental progression along a continuous axis
p_i <- create_standardized_umap_plot(
  data_for_plotting = single_cell_data,
  coloring_variable = "pseudotime",
  display_legend = FALSE  # Hide legend (continuous variable with gradient)
)

## ----------------------------------------------------------------------------
## SECTION 5: PANEL J - CONSTELLATION PLOT (K-NEAREST NEIGHBOR GRAPH)
## ----------------------------------------------------------------------------
## This panel visualizes the relationships between subcluster centroids by
## connecting each centroid to its k=3 nearest neighbors in UMAP space. This
## creates a "constellation" pattern that reveals local neighborhood structure
## and potential developmental trajectories.
##
## Analysis steps:
##   1. Compute centroid (mean UMAP coordinates) for each subcluster
##   2. For each centroid, identify the k=3 nearest neighboring centroids
##   3. Create edge list connecting each centroid to its k nearest neighbors
##   4. Visualize as a network overlay on UMAP space
##
## Visual design:
##   - Edges: Gray lines showing nearest-neighbor connections
##   - Nodes: Points colored by cell class (same colors as Panel C)
##   - Edge transparency: 40% to reduce visual clutter
##
## Biological interpretation:
##   - Clusters of connected centroids suggest transcriptionally related cell types
##   - Long edges may indicate developmental transitions or rare cell states
##   - Dense connectivity regions indicate areas of high cellular diversity
## ----------------------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 5.1: Compute subcluster centroids in UMAP space
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For each subcluster, calculate the mean UMAP coordinates across all cells
## belonging to that subcluster. This provides a representative location for
## each subcluster in the UMAP embedding.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subcluster_centroid_coordinates <- single_cell_data %>%
  group_by(subcluster_id, class) %>%  # Group by subcluster and retain class for coloring
  summarize(
    umap1 = mean(umap1),  # Mean UMAP dimension 1 coordinate
    umap2 = mean(umap2),  # Mean UMAP dimension 2 coordinate
    .groups = "drop"      # Drop grouping after summarization
  )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 5.2: Construct K-nearest neighbor graph (k=3)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For each subcluster centroid, identify its k=3 nearest neighbors based on
## Euclidean distance in 2D UMAP space. Create an edge list connecting each
## centroid to its neighbors.
##
## Algorithm:
##   1. Extract UMAP coordinates as a matrix (rows = centroids, cols = dimensions)
##   2. For each centroid i:
##      a. Compute Euclidean distance to all other centroids
##      b. Identify indices of k=3 nearest neighbors (excluding self)
##      c. Create edges from centroid i to each of its k neighbors
##   3. Combine all edges into a single data frame
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Number of nearest neighbors to connect for each centroid
k_nearest_neighbors <- 3

# Extract centroid coordinates as a matrix for efficient distance computation
centroid_coordinate_matrix <- as.matrix(subcluster_centroid_coordinates[, c("umap1", "umap2")])

# Build edge list by iterating over each centroid
k_nn_edge_list <- lapply(seq_len(nrow(centroid_coordinate_matrix)), function(i) {
  # Compute Euclidean distance from centroid i to all other centroids
  # Formula: sqrt((x1-x2)^2 + (y1-y2)^2) for each pair
  distances_to_all_centroids <- sqrt(colSums((t(centroid_coordinate_matrix) - centroid_coordinate_matrix[i, ])^2))

  # Find indices of k+1 nearest centroids (including self at distance 0)
  # Then exclude self to get k nearest neighbors
  nearest_neighbor_indices <- order(distances_to_all_centroids)[2:(k_nearest_neighbors + 1)]

  # Create data frame with edges from centroid i to its k nearest neighbors
  data.frame(
    from_id   = subcluster_centroid_coordinates$subcluster_id[i],
    to_id     = subcluster_centroid_coordinates$subcluster_id[nearest_neighbor_indices],
    from_x    = subcluster_centroid_coordinates$umap1[i],
    from_y    = subcluster_centroid_coordinates$umap2[i],
    to_x      = subcluster_centroid_coordinates$umap1[nearest_neighbor_indices],
    to_y      = subcluster_centroid_coordinates$umap2[nearest_neighbor_indices],
    class_id  = subcluster_centroid_coordinates$class[i],
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()  # Combine all edge data frames into a single data frame

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SECTION 5.3: Visualize constellation plot with edges and centroid nodes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create a UMAP-based visualization showing:
##   1. Gray edges connecting each centroid to its k=3 nearest neighbors
##   2. Colored points representing subcluster centroids (colored by cell class)
##
## This visualization reveals the local neighborhood structure of subclusters
## and can highlight potential developmental trajectories or transitional states.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p_j <- ggplot() +
  # Draw edges connecting nearest-neighbor centroids
  geom_segment(
    data = k_nn_edge_list,
    mapping = aes(
      x = from_x,     # Start x-coordinate (source centroid)
      xend = to_x,    # End x-coordinate (target centroid)
      y = from_y,     # Start y-coordinate (source centroid)
      yend = to_y     # End y-coordinate (target centroid)
    ),
    alpha = 0.4,      # Semi-transparent to reduce visual clutter
    color = "gray50"  # Neutral gray color for edges
  ) +
  # Draw centroid points colored by cell class
  geom_point(
    data = subcluster_centroid_coordinates,
    mapping = aes(
      x = umap1,      # Centroid UMAP dimension 1
      y = umap2,      # Centroid UMAP dimension 2
      color = class   # Color by cell class (consistent with Panel C)
    ),
    size = 2          # Larger points than individual cells for visibility
  ) +
  labs(
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Cell Class",
    subtitle = "Subcluster centroids connected by K-nearest neighbors (k=3)"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text = element_blank(),    # Hide axis text (UMAP coordinates are relative)
    axis.ticks = element_blank(),   # Hide axis ticks
    legend.position = "right",      # Show legend for cell class colors
    plot.title = element_text(face = "bold", hjust = 0),     # Bold title
    plot.subtitle = element_text(face = "italic", size = 8)  # Italic subtitle
  )

## ----------------------------------------------------------------------------
## SECTION 6: FIGURE COMPOSITION AND EXPORT
## ----------------------------------------------------------------------------
## Assemble all individual panels (p_a through p_j) into a single multi-panel
## figure using the patchwork package. The layout is organized into four rows:
##
##   Row 1: Panel a (developmental timeline) - height ratio 1
##   Row 2: Panel b (taxonomy composite) - height ratio 3 (tallest)
##   Row 3: Panels c, d, e, f (UMAP hierarchical) - height ratio 1
##   Row 4: Panels g, h, i, j (UMAP developmental) - height ratio 1
##
## Layout syntax:
##   - `/` operator: Vertical stacking (places panels in rows)
##   - `|` operator: Horizontal arrangement (places panels side-by-side)
##   - `plot_layout()`: Specifies height ratios and guide collection behavior
##   - `plot_annotation()`: Adds panel labels (a-j) automatically
##
## CRITICAL REQUIREMENT:
##   The plot object names (p_a, p_b, ..., p_j) and the layout structure
##   must remain exactly as specified to maintain semantic equivalence with
##   the original figure. These are "unchangeable anchors" that define the
##   figure structure.
## ----------------------------------------------------------------------------

# Compose all panels into final multi-panel figure
figure_s1_complete <- (
  p_a /                               # Row 1: Panel a (timeline)
  p_b /                               # Row 2: Panel b (taxonomy composite)
  (p_c | p_d | p_e | p_f) /          # Row 3: Panels c-f (UMAP hierarchical)
  (p_g | p_h | p_i | p_j)            # Row 4: Panels g-j (UMAP developmental)
) +
  # Specify row height ratios and guide collection behavior
  plot_layout(
    heights = c(1, 3, 1, 1),  # Allocate most space to Panel b (taxonomy)
    guides = "collect"         # Collect all legends to the right side
  ) +
  # Add panel labels (a, b, c, ..., j) automatically
  plot_annotation(tag_levels = "a") &
  # Apply bold, large panel tags to all panels
  theme(plot.tag = element_text(face = "bold", size = 20))

# Display the composed figure in R graphics device
print(figure_s1_complete)

# Export figure as high-resolution PNG file for publication
ggsave(
  filename = "figure/s1_personality_A.png",  # Output file path
  plot = figure_s1_complete,                 # ggplot object to save
  width = 8,                                 # Figure width in inches
  height = 8,                                # Figure height in inches
  dpi = 300,                                 # Resolution (dots per inch) for publication quality
  bg = "white"                               # White background
)

## ----------------------------------------------------------------------------
## SECTION 7: REPRODUCIBILITY DOCUMENTATION
## ----------------------------------------------------------------------------
## Save R session information to document the exact versions of R and all
## packages used to generate this figure. This ensures that the analysis can
## be reproduced in the future even if package updates introduce changes.
##
## Session info includes:
##   - R version and platform
##   - Versions of all loaded packages
##   - System locale and other environment details
## ----------------------------------------------------------------------------

# Capture session information
session_information <- sessionInfo()

# Print session info to console for immediate reference
print(session_information)

# Save session info to a text file for long-term documentation
# writeLines(
#   capture.output(session_information),
#   "session_info_personality_A.txt"
# )

cat("\n============================================================================\n")
cat("Figure generation complete!\n")
cat("Output: figure/s1_personality_A.png (300 DPI, 8\" × 8\")\n")
cat("============================================================================\n")
