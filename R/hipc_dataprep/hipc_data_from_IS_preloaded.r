# copy preloaded data from ImmuneSpace using hipc_data_from_IS.r in case ImmuneSpace won't wok in future
# data were loaded 10/11/2017

dn.from = file.path("./data/HIPC/IS_preloaded/")
dn.to = file.path("./generated_data/HIPC")
dir.create(dn.to, showWarnings = F)
file.copy(file.path(dn.from, "HIPC_IS_esets.rds"), dn.to)
file.copy(file.path(dn.from, "HIPC_IS_esets_d0.rds"), dn.to)
