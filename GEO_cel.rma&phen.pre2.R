library(data.table)
library(tidyverse)
library(sva)
source("./0.fun/lzXENA.R")
source("./0.fun/lzGEO.R")
source("./0.fun/lzPR.R")

## 1. cleadata phen数据
phens.f <- list.files("./data/", pattern = "cleandata.csv$", full.names = T) 
phens <- lapply(phens.f, function(x)  fread(x, data.table=F) )
phens8 <- lapply(phens, function(x) x[,1:8])
names(phens8[[1]])
names(phens8[[2]])
names(phens8[[3]])
names(phens8[[4]])
phens8[[1]]$OS.time %>% range()
phens8[[1]]$OS.time <- round(phens8[[1]]$OS.time * 365)
phens8[[2]]$OS.time %>% range()
phens8[[2]] %>% view()
phens8[[3]] %>% view()
phens8[[4]] %>% view()
phens8[[4]]$OS.time %>% range()
phens8[[4]]$OS.time <- round(phens8[[4]]$OS.time * 365)
phen.cbind <- do.call(rbind, phens8)       

## 2. 4set rma数据
eset <- fread("./data/0_cbind.eset.rma_GSE29013&GSE31210&GSE37745&GSE50081.csv.gz",
              data.table = F)
data <- as.matrix(eset[, 2:ncol(eset)])
rownames(data) <- eset[,1]
data[1:4,1:4]
data <- data[, intersect(phen.cbind$sample_id, colnames(data))]
phen.cbindz <- phen.cbind[phen.cbind$sample_id %in% colnames(data), ]
dim(data)
all(colnames(data) == phen.cbindz$sample_id)

###
batch <- phen.cbindz$platform %>% as.factor()
mod <- NULL
combat_data <- ComBat(dat=data, batch=batch, mod=mod, par.prior=TRUE, prior.plot=FALSE)
combat_data[1:4,1:4]
write.csv(combat_data, file = './data/0_cbind.eset.rma_4data.combat.csv', 
          row.names = T, quote = F)


#annot
annot <- getGEO(GEO = "GPL570", AnnotGPL = T, getGPL = T, 
                destdir = "E:/OneDrive/Desktop/DM.R1/geo")
annot <- annot@dataTable@table
annot <- annot %>% dplyr::select(ID, `Gene symbol`)
annot %>% head
names(annot) <- c('probe_id', 'symbol')

annot.addr <- "E:/BaiduSyncdisk/Public_TCGA/UCSCxena/annot.2/gencode.v22.annot.RData"
self <- new.env(parent = emptyenv())
load(annot.addr, envir = self)
####
annot_s <- GEO_split_annot(annot_df = annot, taget_col = 'symbol', p = "///")
annot_s %>% head
annot_s$type <- self$all_anot$gencodev22.gtf[match(
  annot_s$symbol, self$all_anot$gencodev22.gtf$symbol), 3]
annot_s %>% head
# write.csv(annot_s, file = 'E:/Public_data/GEO.annot/GPL570.clean.annot.csv')


## 3. luad数据 p1数据
phen.cbindz$Histology %>% table
phen.p1 <- phen.cbindz %>% filter(Histology == "LUAD")
data.p1 <- combat_data[, phen.p1$sample_id] %>% data.frame(check.names = F)
data.p1[1:4,1:4]
data.p1$id <- rownames(data.p1)
data.p1 <- merge(annot_s, data.p1, by='id', all=T)
data.p1[1:5, 1:5]
data.p1$id <- NULL
eset_m <- data.p1 %>% filter(type == "protein_coding")
eset_m[1:5, 1:5]
eset_m <- GEO_UniqueDF(eset_m)
# 
phen.p1[1:5,]
eset_clin <- GEO_MExprClin(eset_m, phen.p1)
checkna(eset_clin[,1:12])
eset_clin[1:5,1:12]
nrow(eset_clin)
write.csv(eset_clin, file = './data/1_cbind.geoluad.rma_4data.combat.csv')
