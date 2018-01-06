library("dplyr")
library("data.table")
library("stringr")
itemrv <- fread("../data/itemsrv.csv")
itemrv <- itemrv[,c("item_nbr","perishable","family_nbr","class_nbr")]
flist <- list.files("../trainitem",full.names=TRUE)
inbr <- str_extract(flist,"[0-9]+")

len <- length(flist)
datlist <- list()
timeline <- as.Date("2017-08-15")
nwd <- 7
frac <- 0.08
for(i in 1:len)
{
    dat <- fread(flist[i])
    dat[,date := as.Date(date)]
    dat <- dat[date <= timeline,]
    nobs <- nrow(dat)
    nstore <- nrow(dat[!duplicated(store_nbr),])
    nsp <- ceiling(nobs*nfrac/nstore/nwd)
    spdat <- dat %>% group_by(dow,store_nbr) %>% sample_n(nsp) %>% data.table()
    spdat[,item_nbr := as.integer(inbr)]
    spdat <- merge(spdat,itemrv,by="item_nbr",all.x=TRUE)
    datlist[[i]] <- copy(spdat)
}
dfsp <- rbindlist(datlist)
nfile <- length(list.files("../sample"))
fname <- paste("../sample/sample",1+nfile,".csv",sep="")
fwrite(dfsp,fname)
