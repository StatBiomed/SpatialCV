## More details
# https://github.com/xzhoulab/SPARK-X-Analysis
# http://htmlpreview.github.io/?https://github.com/xzhoulab/SPARK-X-Analysis/blob/main/analysis/SV2_SPARKX.html

## TO INSTALL spark / sparkx
## https://xzhoulab.github.io/SPARK/04_installation/
# devtools::install_github('xzhoulab/SPARK')

library(SPARK)
library(data.table)


## TO DOWNLOAD the data
# download.file(
#   "https://github.com/xzhoulab/SPARK/raw/master/data/Layer2_BC_Count.rds",
#   "Layer2_BC_Count.rds"
#   )


## Load data
load("./Layer2_BC_Count.rds")

rawcount = Matrix::Matrix(rawcount, sparse = TRUE)

spot_info <- cbind.data.frame(
  x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
  y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)),
  total_counts=apply(rawcount,2,sum)
)
rownames(spot_info) <- colnames(rawcount)
location <- spot_info[, 1:2]

## Run SPATRKX
sparkx_res <- sparkx(rawcount, location, numCores=1, option="mixture")

save(sparkx_res, file="./Layer2_BC_Count_sparkx.rds")


## Save data
Matrix::writeMM(rawcount, './Layer2_BC_Count.mtx')
write.table(spot_info, "./Layer2_BC_Count_spots.tsv", quote = FALSE, sep='\t')
write.table(sparkx_res$res_mtest, "./Layer2_BC_Count_sparkx.tsv", quote = FALSE)

