library(Seurat)
library(clustree)
library(ggraph)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

sig.list = readRDS("sig/sig.list.RDS")

res = fread("results/test_sig.genes_in_clusters_high_vs_low_responders.txt")

res = res %>% dplyr::filter(clustering!="K0", test.dir=="IN")
df.res = res %>% 
  dplyr::mutate(signif = case_when(
    w.pv<=0.001 ~ "***",
    w.pv<=0.01 ~ "**",
    w.pv<=0.05 ~ "*",
    TRUE ~ ""
  ))

dn.fig = "figures/sig_test_clustree"
dir.create(dn.fig, showWarnings = F, recursive = T)

test.name = "SigScore"
for(s in names(sig.list)) {
  cat (s, "\n")
  df.nodes = df.res %>% 
    dplyr::filter(test==test.name, sig==s) %>% 
    dplyr::mutate(HL.score = -log10(w.pv), node = paste0(clustering,"C",cluster)) %>% 
    dplyr::mutate(label_signif = paste0(cluster, signif)) %>% 
    dplyr::select(node, t.stat, w.pv, signif, cluster, label_signif)
    
  
  # Get the graph layout
  layout <- clustree(h1@meta.data %>% dplyr::select(K1:K3), prefix = "K", 
                     return = "layout")
  identical(layout$node, df.nodes$node)
  inode = match(layout$node, df.nodes$node)

  # Add the statistic column (or whatever you want to show)
  layout$clr <- df.nodes$t.stat[inode]
  layout$label <- df.nodes$signif[inode] # or label_signif
  
  # Plot the graph (modified from clustree.R)
  gg <- ggraph::ggraph(layout)
  
  # Add the edges
  gg <- gg + geom_edge_link(#arrow = arrow(length = unit(1.5 * 5, "points"), ends = "last"),
                            end_cap = circle(9.5 * 1, "points"),
                            start_cap = circle(9.5 * 1, "points"),
                            # aes_(colour = ~count,
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
  gg <- gg + geom_node_text(aes_(label = ~label), size = 12,
                            colour = "black")
  
  # Plot theme
  gg <- gg + #scale_size(range = c(4, 15)) +
    ggraph::theme_graph(base_family = "",
                        plot_margin = ggplot2::margin(2, 2, 2, 2))
  gg + 
    scale_color_gradient2(low="blue", mid="grey90", high="red", midpoint = 0, 
                          breaks=c(-3,0,3), limits=c(-3,3), oob=scales::squish, name="T stat") +
    scale_size_continuous(range = c(5, 30)) + 
    coord_flip() +
    scale_x_reverse() +
    scale_y_reverse() +
    guides(size="none", edge_alpha="none", edge_color="none")
  
  ggsave(glue::glue("{dn.fig}/clustree_{test.name}_P_high_vs_low_{s}.png"), w=3, h=15)
  ggsave(glue::glue("{dn.fig}/clustree_{test.name}_P_high_vs_low_{s}.pdf"), w=3, h=15, useDingbats=F)
}

