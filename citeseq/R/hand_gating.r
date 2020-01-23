library(Seurat)

source("R/functions/GenePlot2.r")

dn.fig = "figures/hand_gating"
dir.create(dn.fig, showWarnings = F, recursive = T)

fn = "data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"
h1 = readRDS(fn)

bc = SubsetData(h1, subset.name = "K1", accept.value = "C3")

GenePlot2(bc, "CD20_PROT","CD38_PROT")
ggsave(file.path(dn.fig, "CD20_CD38_biplot_Bc.png"), w=6, h=4)

dat = h1@assay$CITE@data
df = dat %>% as.data.frame() %>% 
  tibble::rownames_to_column("protein") %>% 
  gather("cell","value", -protein)
df = df %>% 
  left_join(h1@meta.data %>% dplyr::select("barcode_check",starts_with("K")), 
            by=c("cell"="barcode_check"))
df.bc = df %>% 
  dplyr::filter(K1=="C3")

prot = c("CD3_PROT", "CD56_PROT", "CD14_PROT", "CD19_PROT", "CD20_PROT", "CD38_PROT")

ggplot(df %>% dplyr::filter(protein %in% prot), aes(value)) +
  geom_density(fill = "grey80", alpha = 0.5) +
  geom_density(data = df.bc %>% dplyr::filter(protein %in% prot), fill = "lightgreen", alpha = 0.5) +
  facet_wrap(~protein) +
  scale_x_continuous(breaks = seq(-5, 15, 5), minor_breaks = seq(-5 , 16, 1)) +
  coord_cartesian(xlim = c(-3,16)) +
  theme_bw()
ggsave(file.path(dn.fig, "Density_Bc_markers.png"))

  # GenePlot2(h1, "CD3_PROT","CD56_PROT")

CBSig = function(SeuratObject, scale.norm = F, 
                 cd3.th = 3, cd56.th = 3, cd14.th = 2, cd19.th = 3, cd20.th = 3.75, cd38.th = 7) {
  cat(paste(cd3.th, cd56.th, cd14.th, cd19.th, cd20.th, cd38.th,"\n"))
  gated.idx <- SeuratObject@assay$CITE@data["CD3_PROT",  ] < cd3.th &
    SeuratObject@assay$CITE@data["CD14_PROT", ] < cd14.th &
    SeuratObject@assay$CITE@data["CD56_PROT", ] < cd56.th &
    SeuratObject@assay$CITE@data["CD19_PROT", ] > cd19.th &
    SeuratObject@assay$CITE@data["CD20_PROT", ] > cd20.th &
    SeuratObject@assay$CITE@data["CD38_PROT", ] > cd38.th
  gated.cells = SeuratObject@cell.names[gated.idx]
  cat(length(gated.cells), " gated cells\n")
  gated <- SubsetData(SeuratObject, cells.use = gated.cells, do.clean = F)
  if (scale.norm == TRUE) {
    gated <- NormalizeData(gated)
    gated <- ScaleData(gated) 
  }
  else  {
    return(gated)
  }
  return(gated)
}

bc38 = CBSig(h1)
bc.hg = CBSig(h1, cd20.th = -Inf, cd38.th = -Inf)
# GenePlot2(bc.hg, "CD20_PROT","CD38_PROT")
# ggsave("CD20_CD38_biplot_Bc_hg.png", w=6, h=4)

df.b = data.frame(x = -Inf, y = -Inf, label=glue::glue("{length(bc@cell.names)} cells"))
df.b.sel = data.frame(x = 4, y = 10.5, label=glue::glue("{length(bc38@cell.names)} cells"))
GenePlot2(bc.hg, "CD20_PROT","CD38_PROT") +
  geom_rect(xmin=3.75, xmax=12.1, ymin=6.5, ymax=10.6, col="red", fill=NA, size=0.5) +
  geom_text(data=df.b, aes(x, y, label=label), hjust = -0.1, vjust = -0.5, col="black", inherit.aes = F) +
  geom_text(data=df.b.sel, aes(x, y, label=label), hjust = 0, vjust = 1, col="red", inherit.aes = F) +
  xlab("CD20") + ylab("CD38")
ggsave(file.path(dn.fig, "CD20_CD38_biplot_Bc_hg_gate.png"), w=4, h=3)
ggsave(file.path(dn.fig, "CD20_CD38_biplot_Bc_hg_gate.pdf"), w=4, h=3)

df.subj = h1@meta.data %>% dplyr::mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response) %>% 
  dplyr::mutate(subject = factor(subject)) %>% 
  distinct() %>% 
  dplyr::arrange(subject)

df.bc = bc@meta.data %>% 
  dplyr::mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response)
df.bc38 = bc38@meta.data %>% 
  dplyr::mutate(response = str_remove(adjmfc.time, "d0 ")) %>% 
  dplyr::select(subject=sampleid, response)

bc.sum = df.bc %>% 
  group_by(subject, response) %>% 
  summarise(nBC = n())
cd38.sum = df.bc38 %>% 
  group_by(subject, response) %>% 
  summarise(nCD38 = n())
DF = left_join(bc.sum, cd38.sum, by=c("subject","response")) %>% 
  dplyr::mutate(nCD38 = ifelse(is.na(nCD38), 0, nCD38)) %>% 
  dplyr::mutate(CD38.freq = nCD38/nBC) %>% 
  dplyr::mutate(response = factor(response, levels=c("low","high")))
fwrite(DF, file="results/CITEseq_CD38hi_cell_data.txt", sep="\t")

# Compare low and high responders ----------------------

X = DF$CD38.freq
Y = DF$response
w.pv = wilcox.test(X[Y=="high"], X[Y=="low"], exact=F, alternative="greater")$p.value

df.text = data.frame(label.w.pv=sprintf("p = %.2g",w.pv))

df.text$x = 0.25
df.text$y = 0.1

ggplot(DF, aes(x=response, y=CD38.freq, group=response)) +
  geom_boxplot(alpha=1, outlier.colour = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=response)) +
  scale_fill_manual(values=c("white", "black"), name="Response",
                    breaks=c(0,2), labels=c("low","high")) +
  geom_text(data = df.text, aes(label=label.w.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
  xlab("Response") + ylab("CD38++ cells (% of alive B cells)") +
  # coord_cartesian(xlim=c(0.8, 2.2)) +
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())

fn.fig = file.path(dn.fig, "CITEseq_CBSig_vs_response")
ggsave(paste0(fn.fig, ".png"), w=3, h=3.5)
ggsave(paste0(fn.fig, ".pdf"), w=3, h=3.5, useDingbats=F)

