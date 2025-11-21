library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(patchwork)

## ------------------------------------------------------------
## 1. Load data
## ------------------------------------------------------------

dev_events <- read.csv("data/figure/s1_dev_events.csv")
cells <- read.csv("data/figure/s1_cells.csv")
taxonomy <- read.csv("data/figure/s1_taxonomy.csv")
edges <- read.csv("data/figure/s1_tree_edges.csv")
nodes <- read.csv("data/figure/s1_tree_nodes.csv")
cluster_counts <- read.csv("data/figure/s1_cluster_counts.csv")

## ------------------------------------------------------------
## 2. Fig. 1a – developmental timeline
## ------------------------------------------------------------

p_a <- ggplot(dev_events) +
  geom_segment(aes(x = x_start, xend = x_end,
                   y = track,   yend = track,
                   color = event_name),
               linewidth = 2, lineend = "round") +
  geom_point(aes(x = x_start, y = track, color = event_name),
             size = 2) +
  scale_x_continuous(
    breaks = seq_along(age_levels),
    labels = age_levels,
    expand = c(0.01, 0.01)
  ) +
  scale_y_reverse(breaks = c(1, 2, 3)) +
  labs(x = "Developmental age", y = NULL, color = "Event") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0)
  )

## ------------------------------------------------------------
## 3. Fig. 1b – tree + bar plot
## ------------------------------------------------------------

## 3.1 Tree (class -> subclass -> cluster)

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

p_tree <- ggraph(g, layout = "dendrogram") +
  geom_edge_diagonal(alpha = 0.4) +
  geom_node_point(aes(color = class_node), size = 2) +
  geom_node_text(
    aes(label = name),
    hjust = -0.1,
    size = 2,
    check_overlap = TRUE
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 9) +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 80, 5.5, 5.5) # space for labels on the right
  )

## 3.2 Bar plot (cluster cell counts)

p_bar <- ggplot(cluster_counts,
                aes(x = reorder(cluster_id, n_cells),
                    y = n_cells,
                    fill = class)) +
  geom_col() +
  labs(
    x = NULL,
    y = "Number of cells",
    fill = "Class"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_blank(),  # visually lines up with tree
    axis.ticks.x = element_blank()
  )


class_counts <- cells %>%
  count(class, modality)

p_b_class_mod <- ggplot(class_counts,
                        aes(x = reorder(class, n), y = n, fill = modality)) +
  geom_col(position = "stack") +
  labs(x = NULL, y = "Number of cells", fill = "Modality") +
  theme_minimal(base_size = 10)


p_b <- (p_tree / p_bar / p_b_class_mod) +
  plot_layout(heights = c(2, 1, 1))
p_b <- wrap_elements(full = p_b)



## 4.2 Helper for UMAP panels
make_umap_plot <- function(data, color_var, legend = TRUE) {
  aes_col <- rlang::sym(color_var)
  p <- ggplot(data, aes(x = umap1, y = umap2, color = !!aes_col)) +
    geom_point(size = 0.7, alpha = 0.8) +
    labs(x = "UMAP1", y = "UMAP2", color = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = if (legend) "right" else "none"
    )
  p
}

p_c <- make_umap_plot(cells, "class")
p_d <- make_umap_plot(cells, "subclass", legend = FALSE)
p_e <- make_umap_plot(cells, "cluster_id", legend = FALSE)
p_f <- make_umap_plot(cells, "subcluster_id", legend = FALSE)
p_g <- make_umap_plot(cells, "age", legend = FALSE)
p_h <- make_umap_plot(cells, "sync_age_bin", legend = FALSE)
p_i <- make_umap_plot(cells, "pseudotime", legend = FALSE)

## 4.3 Simple constellation-like plot (Fig. 1j)

subcluster_centroids <- cells %>%
  group_by(subcluster_id, class) %>%
  summarize(
    umap1 = mean(umap1),
    umap2 = mean(umap2),
    .groups = "drop"
  )

k <- 3
coords <- as.matrix(subcluster_centroids[, c("umap1", "umap2")])

edge_list <- lapply(seq_len(nrow(coords)), function(i) {
  d <- sqrt(colSums((t(coords) - coords[i, ])^2))
  nn <- order(d)[2:(k + 1)]
  data.frame(
    from   = subcluster_centroids$subcluster_id[i],
    to     = subcluster_centroids$subcluster_id[nn],
    from_x = subcluster_centroids$umap1[i],
    from_y = subcluster_centroids$umap2[i],
    to_x   = subcluster_centroids$umap1[nn],
    to_y   = subcluster_centroids$umap2[nn],
    class  = subcluster_centroids$class[i]
  )
}) %>%
  bind_rows()

p_j <- ggplot() +
  geom_segment(data = edge_list,
               aes(x = from_x, xend = to_x,
                   y = from_y, yend = to_y),
               alpha = 0.4) +
  geom_point(data = subcluster_centroids,
             aes(x = umap1, y = umap2, color = class),
             size = 2) +
  labs(x = "UMAP1", y = "UMAP2",
       color = "Class") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0)
  )

## 4.4 Example full Fig.1 layout (if you want it)

fig1_full <- (p_a /
  p_b /
  (p_c | p_d | p_e | p_f) /
  (p_g | p_h | p_i | p_j)) +
  plot_layout(heights = c(1,3,1,1), guides = "collect") + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 20))

fig1_full
ggsave("figure/s1.png", width = 8, height = 8, dpi = 300)
