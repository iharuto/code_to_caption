set.seed(123)

library(dplyr)
library(tibble)

## ------------------------------------------------------------
## 1. Developmental timeline data (Fig. 1a)
## ------------------------------------------------------------

age_levels <- c("E10.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5", "E17", "E18.5",
                "P0", "P1", "P5", "P10", "P15", "P20", "P25", "P28", "P56")

dev_events <- data.frame(
  event_id   = paste0("ev", seq_len(10)),
  event_name = c(
    "Neurogenesis",
    "Tangential migration (GABA)",
    "Radial migration (GABA)",
    "Corticothalamic projections",
    "Thalamocortical projections",
    "Oligodendrogenesis",
    "Astrogenesis",
    "Synapse formation",
    "Synaptic pruning",
    "Blood–brain barrier"
  ),
  age_start  = factor(
    c("E10.5", "E13.5", "E15.5", "E15.5", "E16.5",
      "P0",    "P0",    "P5",    "P15",   "P5"),
    levels = age_levels, ordered = TRUE
  ),
  age_end    = factor(
    c("P1",   "P5",    "P10",   "P10",   "P15",
      "P20",  "P20",   "P15",   "P28",   "P15"),
    levels = age_levels, ordered = TRUE
  ),
  track      = c(1,1,1,2,2,2,2,3,3,3),
  stringsAsFactors = FALSE
)

age_to_num <- function(x) as.numeric(factor(x, levels = age_levels, ordered = TRUE))

dev_events <- dev_events %>%
  mutate(
    x_start = age_to_num(age_start),
    x_end   = age_to_num(age_end)
  )

## ------------------------------------------------------------
## 2. Taxonomy and cell-level data (for UMAPs etc.)
## ------------------------------------------------------------

classes <- c(
  "NEC",
  "CR Glut",
  "RG",
  "IP",
  "IMN",
  "NonIT Glut",
  "IT Glut",
  "GABA",
  "Glioblast",
  "OPC-Oligo",
  "Vascular",
  "Immune"
)

subclasses_list <- list(
  "NEC"         = c("NEC"),
  "CR Glut"     = c("CR Glut"),
  "RG"          = c("RG"),
  "IP"          = c("IP nonIT", "IP IT"),
  "IMN"         = c("IMN nonIT", "IMN IT"),
  "NonIT Glut"  = c("L6 CT CTX Glut", "L5 ET CTX Glut", "L6b CTX Glut"),
  "IT Glut"     = c("L2/3 IT CTX Glut", "L4/5 IT CTX Glut",
                    "L5 IT CTX Glut", "L6 IT CTX Glut"),
  "GABA"        = c("Vip Gaba", "Sst Gaba", "Pvalb Gaba", "Lamp5 Gaba"),
  "Glioblast"   = c("Glioblast"),
  "OPC-Oligo"   = c("OPC NN", "Oligo NN"),
  "Vascular"    = c("Endo NN", "Peri NN", "VLMC NN"),
  "Immune"      = c("Microglia NN", "BAM NN")
)

nt_map <- c(
  "NEC" = "Other",
  "CR Glut" = "Glutamate",
  "RG"  = "Other",
  "IP"  = "Other",
  "IMN" = "Other",
  "NonIT Glut" = "Glutamate",
  "IT Glut"    = "Glutamate",
  "GABA"       = "GABA",
  "Glioblast"  = "Other",
  "OPC-Oligo"  = "Other",
  "Vascular"   = "Other",
  "Immune"     = "Other"
)

## taxonomy: class / subclass / cluster
taxonomy_list <- lapply(names(subclasses_list), function(cl){
  subs <- subclasses_list[[cl]]
  do.call(rbind, lapply(subs, function(sc){
    n_clusters <- 2
    data.frame(
      class      = cl,
      subclass   = sc,
      cluster_id = paste0("cl_", gsub(" ", "_", cl), "_",
                          gsub(" ", "_", sc), "_", seq_len(n_clusters)),
      stringsAsFactors = FALSE
    )
  }))
})
taxonomy <- do.call(rbind, taxonomy_list)

taxonomy$nt_type <- unname(nt_map[taxonomy$class])

## subclusters: 3 per cluster
subclusters_list2 <- lapply(seq_len(nrow(taxonomy)), function(i){
  cl_row <- taxonomy[i, ]
  n_sub <- 3
  data.frame(
    class        = cl_row$class,
    subclass     = cl_row$subclass,
    cluster_id   = cl_row$cluster_id,
    subcluster_id = paste0(cl_row$cluster_id, "_sub", seq_len(n_sub)),
    nt_type      = cl_row$nt_type,
    stringsAsFactors = FALSE
  )
})
subclusters <- do.call(rbind, subclusters_list2)

## cell-level table
n_cells <- 500

sub_idx <- sample(seq_len(nrow(subclusters)), size = n_cells, replace = TRUE)
sub_info <- subclusters[sub_idx, ]

age_vec <- sample(age_levels, size = n_cells, replace = TRUE)
sync_age_bin <- as.integer(factor(age_vec, levels = age_levels, ordered = TRUE))

modality <- sample(c("scRNA", "Multiome"), size = n_cells, replace = TRUE,
                   prob = c(0.7, 0.3))

pseudotime <- runif(n_cells)
umap1 <- rnorm(n_cells)
umap2 <- rnorm(n_cells)

cells <- data.frame(
  cell_id        = paste0("cell_", seq_len(n_cells)),
  modality       = modality,
  class          = sub_info$class,
  subclass       = sub_info$subclass,
  cluster_id     = sub_info$cluster_id,
  subcluster_id  = sub_info$subcluster_id,
  nt_type        = sub_info$nt_type,
  age            = factor(age_vec, levels = age_levels, ordered = TRUE),
  sync_age_bin   = sync_age_bin,
  pseudotime     = pseudotime,
  umap1          = umap1,
  umap2          = umap2,
  stringsAsFactors = FALSE
)

## ------------------------------------------------------------
## 3. Extra data for Fig. 1b: tree + bar
## ------------------------------------------------------------

## Edges for tree (class -> subclass -> cluster)
edges_class_subclass <- taxonomy %>%
  distinct(class, subclass) %>%
  transmute(from = class, to = subclass)

edges_subclass_cluster <- taxonomy %>%
  distinct(subclass, cluster_id) %>%
  transmute(from = subclass, to = cluster_id)

edges <- bind_rows(edges_class_subclass, edges_subclass_cluster)

## Node table with class annotation (for tree colors)
nodes <- tibble(name = unique(c(edges$from, edges$to)))

nodes <- nodes %>%
  left_join(
    taxonomy %>%
      distinct(subclass, class) %>%
      rename(node = subclass),
    by = c("name" = "node")
  ) %>%
  left_join(
    taxonomy %>%
      distinct(cluster_id, class) %>%
      rename(node = cluster_id),
    by = c("name" = "node"),
    suffix = c("_from_subclass", "_from_cluster")
  ) %>%
  mutate(
    class_from_subclass = ifelse(is.na(class_from_subclass),
                                 class_from_cluster, class_from_subclass),
    class_node = dplyr::case_when(
      name %in% classes ~ name,
      !is.na(class_from_subclass) ~ class_from_subclass,
      TRUE ~ NA_character_
    )
  ) %>%
  select(name, class_node)

## Dummy cell counts per cluster (for bar plot)
set.seed(123)
cluster_counts <- taxonomy %>%
  mutate(n_cells = sample(50:300, size = n(), replace = TRUE))

## ------------------------------------------------------------
## 4. Save everything to a single RDS
## ------------------------------------------------------------

# export CSV
write.csv(dev_events, "data/figure/s1_dev_events.csv", row.names = FALSE)
write.csv(cells,      "data/figure/s1_cells.csv",      row.names = FALSE)
write.csv(taxonomy,   "data/figure/s1_taxonomy.csv",   row.names = FALSE)
write.csv(edges,      "data/figure/s1_tree_edges.csv", row.names = FALSE)
write.csv(nodes,      "data/figure/s1_tree_nodes.csv", row.names = FALSE)
write.csv(cluster_counts, "data/figure/s1_cluster_counts.csv", row.names = FALSE)
