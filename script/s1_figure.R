library(ggplot2)
library(patchwork)

## --- 1a: developmental timeline -----------------------------

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
  labs(x = "Developmental age", y = NULL, color = "Event",
       title = "a  Developmental events") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0)
  )

## --- 1b: composition barplot --------------------------------
# Example: cells per class by modality
class_counts <- cells |>
  count(class, modality)

p_b <- ggplot(class_counts,
              aes(x = reorder(class, n), y = n, fill = modality)) +
  geom_col(position = "stack") +
  coord_flip() +
  labs(x = NULL, y = "Number of cells", fill = "Modality",
       title = "b  Class × modality") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0)
  )

## Helper to make UMAP plots with different color encodings ----
make_umap_plot <- function(data, color_var, label, legend = TRUE) {
  aes_col <- rlang::sym(color_var)
  p <- ggplot(data, aes(x = umap1, y = umap2, color = !!aes_col)) +
    geom_point(size = 0.7, alpha = 0.8) +
    coord_equal() +
    labs(x = "UMAP1", y = "UMAP2",
         title = label, color = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = if (legend) "right" else "none"
    )
  p
}

## 1c–1i: UMAPs with different colorings -----------------------

p_c <- make_umap_plot(cells, "class",        "c  UMAP by class")
p_d <- make_umap_plot(cells, "subclass",     "d  UMAP by subclass",  legend = FALSE)
p_e <- make_umap_plot(cells, "cluster_id",   "e  UMAP by cluster",   legend = FALSE)
p_f <- make_umap_plot(cells, "subcluster_id","f  UMAP by subcluster",legend = FALSE)
p_g <- make_umap_plot(cells, "age",          "g  UMAP by age",       legend = FALSE)
p_h <- make_umap_plot(cells, "sync_age_bin", "h  UMAP by sync age",  legend = FALSE)
p_i <- make_umap_plot(cells, "pseudotime",   "i  UMAP by pseudotime",legend = FALSE)

## --- 1j: simple constellation-like plot ---------------------
# Subcluster centroids in UMAP space, with dummy edges between nearest neighbors

subcluster_centroids <- cells |>
  group_by(subcluster_id, class) |>
  summarize(
    umap1 = mean(umap1),
    umap2 = mean(umap2),
    .groups = "drop"
  )

# Create edges: link each subcluster to its 3 nearest neighbors
# (very crude graph, but enough to illustrate)
k <- 3
coords <- as.matrix(subcluster_centroids[, c("umap1", "umap2")])

edge_list <- lapply(seq_len(nrow(coords)), function(i) {
  d <- sqrt(colSums((t(coords) - coords[i, ])^2))
  nn <- order(d)[2:(k + 1)]  # skip itself
  data.frame(
    from = subcluster_centroids$subcluster_id[i],
    to   = subcluster_centroids$subcluster_id[nn],
    from_x = subcluster_centroids$umap1[i],
    from_y = subcluster_centroids$umap2[i],
    to_x   = subcluster_centroids$umap1[nn],
    to_y   = subcluster_centroids$umap2[nn],
    class  = subcluster_centroids$class[i]
  )
}) |> bind_rows()

p_j <- ggplot() +
  geom_segment(data = edge_list,
               aes(x = from_x, xend = to_x,
                   y = from_y, yend = to_y),
               alpha = 0.4) +
  geom_point(data = subcluster_centroids,
             aes(x = umap1, y = umap2, color = class),
             size = 2) +
  coord_equal() +
  labs(x = "UMAP1", y = "UMAP2",
       title = "j  Constellation-like plot",
       color = "Class") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0)
  )

## ============================================================
## 3. Patchwork: assemble into a single Fig. 1-like layout
## ============================================================

# Example layout:
# Row 1: [a] spanning full width
# Row 2: [b][c][d]
# Row 3: [e][f][g]
# Row 4: [h][i][j]

row2 <- p_b | p_c | p_d
row3 <- p_e | p_f | p_g
row4 <- p_h | p_i | p_j

fig1 <- p_a /
  row2 /
  row3 /
  row4 +
  plot_layout(heights = c(1.2, 1, 1, 1)) &
  theme(plot.margin = margin(5, 5, 5, 5))

# Print full figure
fig1