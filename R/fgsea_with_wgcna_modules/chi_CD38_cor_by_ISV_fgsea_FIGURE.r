library(tmod)

dn.out = file.path(PROJECT_DIR, "generated_data/fgsea_with_wgcna_modules/")
dir.create(dn.out, showWarnings = F)
dn.fig = file.path(PROJECT_DIR, "figure_generation/SLE-Sig")
dir.create(dn.fig, showWarnings = F)


isv.seq = c(0.5,0.75)

# load FGSEA resuls for BTM modules
mset = "LI"
fn.res = file.path(PROJECT_DIR, "generated_data/fgsea_with_btm_modules", 
                   sprintf("BTM.%s_genes_fgsea_results.rds", mset))
res = readRDS(fn.res)

iuse = names(res) %in% paste0("ISV.",isv.seq)
res = res[iuse]

# extract p-values and AUC as effect size
df.pv = data.frame(pathway = res[[1]]$pathway, stringsAsFactors = F)
df.padj = df.pv
for(i in 1:length(res)) {
  tmp = res[[i]] %>% dplyr::select(pathway, pval)
  names(tmp)[2] = names(res)[i]
  df.pv = full_join(df.pv, tmp, by="pathway")
  tmp = res[[i]] %>% dplyr::select(pathway, padj)
  names(tmp)[2] = names(res)[i]
  df.padj = full_join(df.padj, tmp, by="pathway")
}
df.btm = df.pv %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("pathway") %>% 
  data.matrix()
df.btm.padj = df.padj %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("pathway") %>% 
  data.matrix()

# load FGSEA resuls for WGCNA
fn.res = file.path(PROJECT_DIR, "generated_data", "fgsea_with_wgcna_modules", 
                   "WGCNA_genes_fgsea_results.rds")
res = readRDS(fn.res)

# select ISV threshold to plot and filter the loaded results
iuse = names(res) %in% paste0("ISV.",isv.seq)
res = res[iuse]

# extract p-values and AUC as effect size
df.pv = data.frame(pathway = res[[1]]$pathway, stringsAsFactors = F)
df.padj = df.pv
for(i in 1:length(res)) {
  tmp = res[[i]] %>% dplyr::select(pathway, pval)
  names(tmp)[2] = names(res)[i]
  df.pv = full_join(df.pv, tmp, by="pathway")
  tmp = res[[i]] %>% dplyr::select(pathway, padj)
  names(tmp)[2] = names(res)[i]
  df.padj = full_join(df.padj, tmp, by="pathway")
}
df.wgcna = df.pv %>% 
  dplyr::filter(pathway %in% "brown") %>%
  tibble::remove_rownames() %>% tibble::column_to_rownames("pathway") %>% 
  data.matrix()
df.wgcna.padj = df.padj %>% 
  dplyr::filter(pathway %in% "brown") %>%
  tibble::remove_rownames() %>% tibble::column_to_rownames("pathway") %>% 
  data.matrix()

df.comb = rbind(df.wgcna, df.btm)
df.comb.padj = rbind(df.wgcna.padj, df.btm.padj)
df.comb = apply(df.comb, 2, p.adjust, method="BH")

# filter modules
pval.row.in = 0.05 # p threshold for modules to be included
pval.col.in = 0.05 # p threshold for columns to be included
pval.min = 0.1 # p threshold for an enrichment to be shown

i.row.in = apply(df.comb, 1, function(x) any(x <= pval.row.in, na.rm = T))

df.comb = df.comb[i.row.in, ]

rownm = rownames(df.comb)
ibrown = rownm %in% "brown"
data(tmod)
rowttl = tmod$MODULES$Title[match(rownm, tmod$MODULES$ID)]
rownm2 = sprintf("%s (%s)", rownm, rowttl)
rownm2[ibrown] = "WGCNA brown module"

i.hide = df.comb > pval.min
df.comb[i.hide] = NA

df.comb = df.comb %>% as.data.frame() %>% 
  mutate(module = rownm2) %>% 
  gather("comparison","P.Value", -module)

df.mod = df.comb %>% 
  dplyr::filter(comparison == "ISV.0.5", !grepl("TBA", module), !is.na(P.Value))

iord = order(-log10(df.mod$P.Value[-1]))
iord = c(iord+1, 1)
df.mod = df.mod %>% 
  mutate(module = factor(module, levels=module[iord]))

p = ggplot(df.mod, aes(module, -log10(P.Value))) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(c(0.05, 0.01)), col="red", lty=2) +
  geom_text(data=data.frame(y=-log10(c(0.05,0.01)), x=nrow(df.mod), label=c("5% FDR", "1% FDR")), 
            aes(x=x,y=y,label=label), col="red", size=2, vjust=-2, hjust=0.5) +
  xlab("") + ylab("adjusted p-value,-log10") +
  coord_flip() +
  theme_bw() +
  theme(
    # panel.border = element_blank(), 
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5,5,0,0,"mm"))
# Disable clip-area.
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

library(grid)
library(gridExtra)
grid.draw(gt)

fn.fig = file.path(dn.fig, "brown_BTM_modules_enrichment_CD38.cor.ranked_ISV0.5_barplot")
ggsave(paste0(fn.fig, ".png"), gt, w=5.5, h=3)
ggsave(paste0(fn.fig, ".pdf"), gt, w=5.5, h=3)

