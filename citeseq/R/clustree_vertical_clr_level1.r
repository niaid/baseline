library(Seurat)
fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

library(clustree)
library(ggraph)

dir.create("clustree_sig", showWarnings = F)

  # Get the graph layout
  layout <- clustree(h1@meta.data %>% dplyr::select(K1:K3), prefix = "K", 
                     return = "layout")
  df.nodes = h1@meta.data %>% dplyr::select(K1:K3) %>% 
    gather("K","cluster") %>% 
    mutate(node=paste0(K,"C",cluster)) %>% 
    distinct() %>% 
    mutate(clr = str_extract(cluster, "C\\d"))
  inode = match(layout$node, df.nodes$node)

  # Add the statistic column (or whatever you want to show)
  layout$clr <- df.nodes$clr[inode]

  # Plot the graph (modified from clustree.R)
  gg <- ggraph::ggraph(layout)
  
  # Add the edges
  gg <- gg + geom_edge_link(#arrow = arrow(length = unit(1.5 * 5, "points"), ends = "last"),
                            end_cap = circle(9.5 * 1, "points"),
                            start_cap = circle(9.5 * 1, "points"),
                            aes_(
                                 alpha = ~in_prop,
                                 edge_width = ~is_core) ) +
    scale_edge_width_manual(values = c(1.5, 1.5),
                            guide = "none") +
    scale_edge_colour_gradientn(colours = viridis::viridis(256)) +
    scale_edge_alpha(limits = c(0, 1))
  
  # Add the node points (replace "statistic" with the column you want to show)
  gg <- gg + clustree:::add_node_points("clr", "size", 1, colnames(layout))
  
  # Add the node text

  # Plot theme
  gg <- gg + #scale_size(range = c(4, 15)) +
    ggraph::theme_graph(base_family = "",
                        plot_margin = ggplot2::margin(2, 2, 2, 2))

  cm = pals::glasbey(13)[c(1:3,5:9,13,11)]# %>% as.character()
  gg + 
    scale_color_manual(values = cm) +
    scale_size_continuous(range = c(5, 30)) + 
    coord_flip() +
    scale_x_reverse() +
    scale_y_reverse() +
    guides(size="none", edge_alpha="none", edge_color="none", color="none")
  
  dn.fig = "figures/cluster_annotation"
  dir.create(dn.fig, showWarnings = F, recursive = T)
  ggsave(glue::glue("{dn.fig}/clustree_clr_level1.png"), w=3, h=15)
  ggsave(glue::glue("{dn.fig}/clustree_clr_level1.pdf"), w=3, h=15, useDingbats=F, colormodel="rgb")

