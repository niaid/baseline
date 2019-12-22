for (study in c("SDY212", "SDY400", "SDY404")) {
  
  fn.ge = file.path(PROJECT_DIR, "generated_data", "HIPC", 
                    paste0(study, "_GE_matrix_gene.txt"))
  dat = fread(fn.ge, data.table = F)
  
  fn.si = file.path(PROJECT_DIR, "generated_data", "HIPC",
                    paste0(study, "_sample_info.txt"))
  info = fread(fn.si)
  
  si.d0 = info$time == "d0"
  
  dat.d0 = dat[,c(1,which(si.d0)+1)]
  info.d0 = info[si.d0,]
  
  fn.ge.d0 = sub(".txt", "_day0.txt", fn.ge)
  fn.si.d0 = sub(".txt", "_day0.txt", fn.si)
  
  fwrite(dat.d0, fn.ge.d0, sep="\t", quote=T)
  fwrite(info.d0, fn.si.d0, sep="\t", quote=T)
  
  si.hl = info.d0$Response %in% c("low","high")
  dat.hl = dat.d0[,c(1,which(si.hl)+1)]
  info.hl = info.d0[si.hl,]
  
  fn.ge.hl = sub(".txt", "_ResponseLoHi.txt", fn.ge.d0)
  fn.si.hl = sub(".txt", "_ResponseLoHi.txt", fn.si.d0)
  
  fwrite(dat.hl, fn.ge.hl, sep="\t", quote=T)
  fwrite(info.hl, fn.si.hl, sep="\t", quote=T)
}
