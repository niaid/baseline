library(pROC)

fn.poB = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_percent_of_B_filtered_day0.txt")
flow.poB = fread(fn.poB) %>% 
  select(sample, matches("Gate")) %>% 
  dplyr::rename(CD38high = Gate3, 
                CD38high.CD10pos = Gate2,
                CD38high.CD10neg = Gate1,
                CD38pos = Gate4)


fn.info = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_sample_info_filtered_day0.txt")
flow.info = fread(fn.info)

fn.titer = file.path(PROJECT_DIR, "data", "CHI", "phenotypes", "titer_processed.txt")
df.titer = fread(fn.titer) %>% 
  mutate(Subject = as.character(Subject)) %>%
  mutate(Response = ifelse(adjMFC_class==0, "low",
                           ifelse(adjMFC_class==2, "high", "middle")))

flow.titer = flow.info %>%
  left_join(df.titer %>% dplyr::select(Subject, Response), by=c("subject"="Subject")) %>% 
  dplyr::filter(Response %in% c("low","high"),
         time %in% c(-7,0,70)) %>% 
  mutate(Response = factor(Response,levels=c("low","high")))

# get originally gated flow data
fn.old = file.path(PROJECT_DIR, "data", "CHI/flow/original_gates", "day0.log10.txt")
flow.old = read.table(fn.old, sep="\t", header=T, row.names=1, stringsAsFactors=F)
flow.old = as.data.frame(t(flow.old)) %>% 
  dplyr::select(ID87, ID91, ID96, ID103, ID108) %>% 
  tibble::rownames_to_column("subject") %>% 
  mutate(subject=sub("X","",subject) %>% as.numeric())

flow = flow.poB %>% 
  inner_join(flow.titer, by="sample") %>% 
  inner_join(flow.old, by="subject")

pops = names(flow)[c(rev(2:5),12)]
pops.name = pops
pops.name[1:4] = 1:4
pops.group = (1:3)[c(2,2,2,2,3)] %>% factor()

df.auc = data.frame(Population=factor(pops, levels=rev(pops)), group=pops.group, AUC=NA)
for (p in seq_along(pops)) {
  df.auc$AUC[p] = roc(flow[["Response"]],flow[[pops[p]]], direction="<", quiet = T)$auc
}

ggplot(df.auc, aes(Population, AUC, fill=group)) + geom_bar(stat="identity") +
  geom_hline(yintercept = 0.5, col="black", lty=1, size=1) +
  scale_x_discrete(labels=rev(pops.name)) +
  xlab("") +
  coord_flip() + 
  theme_bw() + theme(legend.position="none")

fn.fig = file.path(PROJECT_DIR, "figure_generation", sprintf("CHI_flow_selected_gates_AUC"))
ggsave(paste0(fn.fig, ".png"), w=4,h=3)
ggsave(paste0(fn.fig, ".pdf"), w=4,h=3)
