## ---------------------------------------------
## Dummy data to mimic Fig. 1 (panels a–j)
## ---------------------------------------------

set.seed(123)

# Developmental ages
age_levels <- c("E10.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5", "E17", "E18.5",
                "P0", "P1", "P5", "P10", "P15", "P20", "P25", "P28", "P56")

age_factor <- factor(age_levels, levels = age_levels, ordered = TRUE)

# --- Fig. 1a dummy events -----------------------------------
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
  track      = c(1,1,1,2,2,2,2,3,3,3)
)

# Helper numeric encoding for plotting segments on x-axis
age_to_num <- function(x) as.numeric(factor(x, levels = age_levels, ordered = TRUE))
dev_events$x_start <- age_to_num(dev_events$age_start)
dev_events$x_end   <- age_to_num(dev_events$age_end)

# --- Cell-level data (Fig. 1b–j) ----------------------------
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

taxonomy <- lapply(names(subclasses_list), function(cl){
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
}) |> do.call(rbind, args = _)

taxonomy$nt_type <- unname(nt_map[taxonomy$class])

subclusters <- lapply(seq_len(nrow(taxonomy)), function(i){
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
}) |> do.call(rbind, args = _)

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
