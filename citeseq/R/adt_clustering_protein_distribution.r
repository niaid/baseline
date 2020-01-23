suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))

library(ggridges)
library(viridis)

fn = file.path("data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds")
h1 = readRDS(fn)

# hclust average protein data by celltype 
prot_use = c("CD2","CD3","CD5","CD7","CD4","CD8","CD62L","CD45RA","CD45RO","CD27","CD28","CD278 ", "CD25", 
             "CD127","CD161", "KLRG1","CD195","CD314 ", "CD194","CD103", "CD56","CD57","CD244",
             "CD16","CD14", "CD11b","CD11c", "CD1d", "CD33", "CD13", "CD31", "CD64", "CD163","CD86", "HLA-DR", 
             "CD123", "CD141", "CD71", "CD303", "CD117","CD34", "CD38","CD39", "CD1c", "CD32", 
             "IgM","IgD","IgA","CD19","CD20","CD21","CD24","CD40","CD185","CD196","BTLA")


adt.l = h1@assay$CITE@data %>% t %>% as.data.frame() %>% 
  setNames(str_remove(rownames(h1@assay$CITE@data), "_PROT")) 
celltypes = h1@meta.data %>% dplyr::select(K3)
adt.l$celltype = celltypes$K3
adt.l = adt.l %>% 
  gather(key = prot, value = dsb_count, AnnexinV:CD20) %>% 
  dplyr::filter(prot %in% prot_use)

adt.l$prot = factor(adt.l$prot, levels = rev(prot_use))
adt.l$celltype = factor(adt.l$celltype)#, levels=celltype_order)


cu = pals::alphabet(n = 26)
cu = colorRampPalette(cu)(length(prot_use))
p = ggplot(adt.l, aes(x = dsb_count, y = prot, fill = prot, color = prot, alpha = 1)) + 
  geom_density_ridges2(show.legend = F, inherit.aes = T, scale = 2) + 
  theme_ridges(font_size = 10, center_axis_labels = T, grid = FALSE) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
  scale_color_manual(values = cu) + 
  scale_fill_manual(values = cu) +
  facet_wrap(~celltype, ncol = 24, scales = "free_x") + 
  ylab("protein") +
  theme(strip.background =element_blank())+
  theme(strip.text = element_text(colour = 'black', size = 6, face = "bold")) + 
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.y = element_text(size = 6, face = "bold")) 
dn.fig = "figures/cluster_annotation/"
ggsave(p, filename = file.path(dn.fig, "prot_dist_sub.png"), width = 12.5, height = 8)
ggsave(p, filename = file.path(dn.fig, "prot_dist_sub.pdf"), width = 12.5, height = 8)

