fn.si = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_sample_info_2_sle.txt")
info = fread(fn.si, data.table=F)

fn.ge = file.path(PROJECT_DIR, "generated_data", "SLE", "SLE_ge_matrix_gene_sle.txt")
dat = fread(fn.ge, data.table = F) %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% 
  data.matrix()


info = info %>% 
  mutate(DA = factor(DA, levels=c("high","mid","low")), 
         SUBJECT = factor(SUBJECT, levels=
                            unique(SUBJECT[order(as.numeric(sub("SLE-","",SUBJECT)))]))) #%>%
  # group_by(SUBJECT) %>%
  # mutate(VISIT.INCLUDE = VISIT>=min(VISIT[DA=="low"],na.rm=T), # 48 or 35 subj
  #        # mutate(VISIT.INCLUDE = !(CUMULATIVE_TIME==0 & DA=="high"), # 42 subj
  #        # INCLUDE = sum(SLEDAI<4)>=2 & max(SLEDAI[VISIT.INCLUDE],na.rm=T)>=4 #35 subj
  #        INCLUDE = sum(SLEDAI<4)>=2 & max(SLEDAI)>=4 #48 subj
  # ) %>%
  # ungroup() %>%
  # arrange(SUBJECT,VISIT) %>%
  # filter(INCLUDE)# & VISIT.INCLUDE)

yvar="SLEDAI"
subj = sort(unique(info$SUBJECT))
s1 = 113;s2 = 158 #up to 158
# ggplot(info %>% subset(SUBJECT %in% subj[s1:s2]), 
#        aes_string(x="CUMULATIVE_TIME", y=yvar, group="SUBJECT", col="DA")) + 
#   geom_line(col="gray20",size=0.5) + geom_point(size=1) + facet_wrap(~SUBJECT) +
#   scale_color_manual(values=c("green","blue","red")) +
#   theme_bw()
# ggsave(sprintf("SLE_%s_TIME_subj.%d.%d.png",yvar,s1,s2),w=10,h=6)

# s = c("SLE-55", "SLE-65", "SLE-169", "SLE-201")
s = c("SLE-201")
ggplot(info %>% subset(SUBJECT %in% s), 
       aes_string(x="CUMULATIVE_TIME", y=yvar, group="SUBJECT", fill="DA")) + 
  # geom_rect(xmin=-Inf, xmax=Inf, ymin=3, ymax=7, fill="gray90", col=NA) +
  geom_hline(yintercept = c(3, 7), col="black", lty=2) +
  # geom_rect(xmin=90, xmax=1120, ymin=1.5, ymax=2.5, col="blue", fill=NA) +
  # geom_rect(xmin=-30, xmax=700, ymin=9.5, ymax=16.5, col="red", fill=NA) +
  geom_line(col="gray20",size=0.5) + geom_point(size=3, pch=21) + facet_wrap(~SUBJECT) +
  scale_fill_manual(values=c("black","grey","white")) +
  xlab("Time, days") +
  coord_cartesian(ylim=c(0,17)) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=12),
        legend.position = c(0.9, 0.7), 
        legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"))
# ggsave(sprintf("SLE_%s_TIME_subj.%s.png",yvar,s), w=4, h=2)
fn.fig = file.path(PROJECT_DIR, "figure_generation", 
                   sprintf("SLE_%s_TIME_subj.%s", yvar, paste(sub("SLE-","",s),collapse=".")))
ggsave(paste0(fn.fig, ".png"), w=5, h=3)
ggsave(paste0(fn.fig, ".pdf"), w=5, h=3)
