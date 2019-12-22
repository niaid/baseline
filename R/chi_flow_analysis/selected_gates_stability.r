fn.poB = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_percent_of_B_filtered.txt")
flow.poB = fread(fn.poB) %>% 
  select(sample, matches("Gate")) %>% 
  dplyr::rename(CD38high = Gate3, 
                CD38high.CD10pos = Gate2,
                CD38high.CD10neg = Gate1,
                CD38pos = Gate4)

fn.info = file.path(PROJECT_DIR, "generated_data", "CHI", "flow_sample_info_filtered.txt")
flow.info = fread(fn.info)


# get originally gated flow data
flow.old = data.frame()
for (tp in c("day0", "pre7","day70")) {
  fn.old = file.path(PROJECT_DIR, "data", "CHI/flow/original_gates", sprintf("%s.log10.txt", tp))
  df.tmp = read.table(fn.old, sep="\t", header=T, row.names=1, stringsAsFactors=F)
  df.tmp = as.data.frame(t(df.tmp)) %>% 
    dplyr::select(ID87) %>% 
    tibble::rownames_to_column("subject") %>% 
    mutate(subject=sub("X","",subject), time=tp)
  flow.old = rbind(flow.old, df.tmp)
}
flow.old = flow.old %>% 
  mutate(sample = paste(subject, time, sep="_")) %>% 
  dplyr::select(-time, -subject)

flow = inner_join(flow.poB, flow.info, by="sample") %>% 
  dplyr::filter(time %in% c(0, -7, 70)) %>%
  mutate(subject = as.character(subject)) %>% 
  inner_join(flow.old, by="sample")


pops = names(flow)[c(rev(2:5),11)]
pops.name = c(
  "CD20+CD38+ of total B cells",
  "CD20+CD38++ of total B cells",
  "CD20+CD38++CD10+ of total B cells\n(transitional-like)",
  "CD20+CD38++CD10- of total B cells\n(memory-like)",
  "CD20-CD38++CD27++ of total B cells\n(plasmablasts)"
)
pops.group = (1:3)[c(2,2,2,2,3)] %>% factor()

df.stab = data.frame(Population=factor(pops, levels=rev(pops)), group=pops.group, ISV=NA)
for (i in seq_along(pops)) {
  form = as.formula(sprintf("%s ~ subject", pops[i]))
  fit = aov(form, data=flow)
  ss = summary(fit)[[1]]["Sum Sq"][[1]]
  ssn = ss / sum(ss)
  df.stab$ISV[i] = ssn[1]
}

ggplot(df.stab, aes(Population, ISV, fill=group)) + geom_bar(stat="identity") +
  scale_x_discrete(labels=rev(pops.name), position="top") +
  xlab("") +
  ylab("ISV") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position="none", axis.ticks.y=element_blank())

fn.fig = file.path(PROJECT_DIR, "figure_generation", sprintf("CHI_flow_selected_gates_ISV"))
ggsave(paste0(fn.fig, ".png"), w=6,h=3)
ggsave(paste0(fn.fig, ".pdf"), w=6,h=3)
