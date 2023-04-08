rm(list = ls());gc()
library(GEOquery);library(affyPLM)
library(affy);library(hgu133plus2cdf)
library(dplyr);library(RColorBrewer);library(tidyverse)
source('0.fun/lzGEO.R')
source('0.fun/p1.fun.R')
setproxy()
# 1. 设置数据编号和目录
# -----------------------
# 需要工作目录有这些文件夹存在， data/st res/
# 设置数据路径
p.raw.dir <- list.dirs('E:/Public_data/GEO/GPL570.raw/meta') # g
# #
raw.dir <- p.raw.dir %>% str_subset("GSE[:digit:]+_RAW$")            # h
gse <- raw.dir %>% str_extract("GSE[:digit:]+")
gse

# 2. 读取数据、编号重命名
# -----------------------
# 读取数据
data.raw <- lapply(raw.dir, function(x) {
  ReadAffy(celfile.path = x) #去读当前目录下所有的CEL文件
})
names(data.raw) <- gse
sampleNames(data.raw[[1]])
for (i in 1:length(data.raw)) {
  sampleNames(data.raw[[i]]) <- sampleNames(data.raw[[i]]) %>% 
    str_extract("GSM[:digit:]+")
}

# 3. 标准化及输出
# -----------------------
# rma标准化
data.rma <- lapply(data.raw, rma)
# 提取矩阵
eset.rma <- lapply(data.rma, exprs)
# 合并数据
eset.rmadf <- do.call(cbind, eset.rma)
eset.rmadf[1:5,1:5]
dim(eset.rmadf)
out.name.rma <- paste0('data/0_geo.eset.cbind4set.rma_', 
                       paste(names(eset.rma), collapse = '&'),'.csv')
out.name.rma
# 保存行名，但读取时不可设置为行名(lz格式)
write.csv(eset.rmadf, file = out.name.rma, quote = F, row.names = T)
#::: 保存原始的cel文件的rma后的原始表达数据矩阵表 (剔除了几个不合格cel文件)


# 4. 提取分组信息,并导出指定的数据矩阵
# -----------------------
gses_l <- lapply(gse, function(id) {
  getGEO(id, getGPL = F, destdir = "E:/Public_data/GEO/GPL570.raw/meta")
})
for (i in 1:4) {
  pd <- pData(gses_l[[i]][[1]])
  pd_id <- gse[i]
  write.csv(pd, file = paste0("data/0_geo.phen_", pd_id, '.csv'), 
            quote = F, row.names = T)
}
#pData(gses_l[[4]][[1]])[1:4,1:4]
i = 2
pd <- pData(gses_l[[i]][[1]])
pd_id <- gse[i]
pd_id;dim(pd)
pd[1:4,1:4]
pd %>% View()
pd <- pd[, 48:ncol(pd)]
paste0("data/0_geo.phen_", pd_id, '_cleandata.csv')
write.csv(pd, file = paste0("data/0_geo.phen_", pd_id, '_cleandata.csv'), 
          quote = F, row.names = T)
i = 3
pd <- pData(gses_l[[i]][[1]])
pd_id <- gse[i]
pd_id;dim(pd)
pd[1:4,1:4]
pd %>% View()
pd <- pd[, 41:ncol(pd)]
paste0("data/0_geo.phen_", pd_id, '_cleandata.csv')
write.csv(pd, file = paste0("data/0_geo.phen_", pd_id, '_cleandata.csv'), 
          quote = F, row.names = T)
i = 4
pd <- pData(gses_l[[i]][[1]])
pd_id <- gse[i]
pd_id;dim(pd)
pd[1:4,1:4]
pd %>% View()
pd <- pd[, 44:ncol(pd)]
paste0("data/0_geo.phen_", pd_id, '_cleandata.csv')
write.csv(pd, file = paste0("data/0_geo.phen_", pd_id, '_cleandata.csv'), 
          quote = F, row.names = T)
