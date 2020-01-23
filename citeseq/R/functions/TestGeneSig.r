TestGeneSig <- function(tobj, sig.list, subj.col, df.subj, 
                         test.col, test.col.levels, test.to.print=NA,
                         out.dir = ".") {
  
  library(Seurat)
  library(pROC)
  
  # Calculate score from genes x samples matrix as average z-score
  get_score <- function(x) {
    x = t(scale(t(x)))
    return (colMeans(x, na.rm=T))
  }
  
  dat = tobj@data
  meta = tobj@meta.data %>% 
    dplyr::rename_(subject = subj.col) %>% 
    dplyr::mutate(subject = factor(subject, levels=df.subj$subject))

  df.test = data.frame()
  
  for(s in seq_along(sig.list)) {
    sig = sig.list[[s]]
    sig.name = names(sig.list)[s]
    cat(sig.name, " ")
    gi = rownames(dat) %in% sig
  
    dat2 = aggregate(t(as.matrix(dat[gi,,drop=F])), list(subject = meta$subject), mean, drop=F)
    dat2[is.na(dat2)] = 0
    df = dat2 %>% 
      left_join(df.subj %>% dplyr::select_("subject",response=test.col), by="subject")
    
    mat = df %>% 
      tibble::column_to_rownames("subject") %>% 
      dplyr::select(-response) %>% 
      data.matrix()
    
    if(sum(mat)==0) next
  
    for(test in c("Ngene","SigScore")) {
      if(test == "Ngene") {
        X = rowSums(mat>0)
      } else if(test == "SigScore") {
        X = get_score(t(mat))
      }
      Y = df$response %>% factor(levels=test.col.levels)
      r = roc(Y,X, direction="<", quiet = T)
      w.pv = wilcox.test(X[Y==test.col.levels[2]], X[Y==test.col.levels[1]], alternative = "greater", exact=F)$p.value
      mean.diff = mean(X[Y==test.col.levels[2]], na.rm=T) - mean(X[Y==test.col.levels[1]], na.rm=T)
      t.stat = ifelse(length(unique(X)) == 1, 0,
                      t.test(X[Y==test.col.levels[2]], X[Y==test.col.levels[1]], alternative = "greater")$statistic)
      
      df.test = rbind(df.test, data.frame(sig = sig.name, test, auc=r$auc, mean.diff, t.stat, w.pv, N = length(tobj@cell.names)))
      
      if(test %in% test.to.print) {
        r.df = data.frame(Specificity = r$specificities, Sensitivity=r$sensitivities) %>%
          dplyr::arrange(Sensitivity)
        
        df.fig = data.frame(subject=df$subject, score=X, Response=factor(Y))
        df.text = data.frame(label.auc=sprintf("AUC = %.2f", r$auc), 
                             label.w.pv=sprintf("p = %.2g",w.pv))
        df.text$x = 0.25
        df.text$y = 0.1
        
        p1 = ggplot(df.fig, aes(x=Response, y=score, group=Response)) +
          geom_boxplot(alpha=1, outlier.colour = NA) +
          geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Response)) +
          scale_fill_manual(values=c("white", "black"), name="Response",
                            breaks=c(0,2), labels=c("low","high")) +
          geom_text(data = df.text, aes(label=label.w.pv), x=1.5, y=Inf,vjust=1.1, hjust=0.5, size=4, inherit.aes = F) +
          xlab("Response") + ylab("Score") +
          theme_bw() + theme(legend.position="none") +
          theme(panel.border = element_blank(), strip.background = element_blank(), 
                strip.text.x = element_text(size=12),
                panel.spacing = unit(0,"mm"), panel.grid.major.x = element_blank(),
                axis.ticks = element_blank())
        
        p2 = ggplot(r.df, aes(x=Specificity, y=Sensitivity)) + 
          geom_line(size=1) +
          geom_abline(intercept = 1, slope = 1, lty=2) +
          scale_x_reverse() + 
          geom_text(data=df.text, aes(x=x,y=y,label=label.auc)) +
          coord_fixed() + 
          theme_bw() + theme(panel.grid = element_blank())
        
        pp = plot_grid(p1, p2, ncol=1, rel_heights = c(1, 1))
        
        fn.fig = file.path(out.dir, glue::glue("{sig.name}_score_vs_response"))
        # ggsave(paste0(fn.fig, ".png"), plot=pp, w=3,h=7)
        ggsave(paste0(fn.fig, ".pdf"), plot=pp, w=3,h=7)
      }
    }
  }
  return(df.test)
}
