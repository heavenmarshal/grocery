trans <- fread("transactions.csv/transactions.csv")
sumprom <- fread("sumprom.csv")
transal <- alignmts(dseq,trans,field="transactions")
sumpromal <- alignmts(dseq,sumprom,field="sumprom")
sumpromalt <- melt(sumpromal,value.names="sumprom",id.var="store_nbr",variable.name="date",
                   variable.factor=FALSE,value.factor=FALSE)
sumpromalt <- dcast(sumpromalt,date~store_nbr)
res <- melt(transal,value.name="transactions",id.var="store_nbr",
            variable.name="date",variable.factor=FALSE,value.factor=FALSE)
res[,date:=as.Date(date)]
transalt <- dcast(res,date~store_nbr)
mcity <- summarizeCate(transal,"store_nbr","city_nbr",stores,mean)
mstate <- summarizeCate(transal,"store_nbr","state_nbr",stores,mean)
mtype <- summarizeCate(transal,"store_nbr","store_type_nbr",stores,mean)
mcluster <- summarizeCate(transal,"store_nbr","cluster",stores,mean)
lags <- c(1,2,7,14,21,28)
malags <- c(2,3,7,14,21,28,35,70,140,280)
pmlags <- c(70,140,210,280)

substore <- stores[,c("store_nbr","city_nbr","state_nbr","store_type_nbr","cluster")]
res <- merge(res,substore,by="store_nbr",all.x=TRUE)
holcity[,date := as.Date(date)]
res1 <- merge(res,holcity,by=c("date","city_nbr"),all.x=TRUE)
mnames <- setdiff(names(holcity),names(res))
fill0 <- function(x) ifelse(is.na(x),0,x)
res1[,mnames] <- res1[,lapply(.SD,fill0),.SDcols=mnames]
sumprom[,date := as.Date(date)]
res1 <- merge(res1,sumprom,by=c("date","store_nbr"),all.x=TRUE)
res1[,sumprom := fill0(sumprom)]

res1[,year := year(date)-2013]
res1[,month := month(date)]
res1[,doy := yday(date)]
res1[,dow := wday(date)]

res1 <- proclagma(transalt,"store_nbr","translag",lags,res1,shift)
res1 <- proclagma(transalt,"store_nbr","transma",malags,res1,lagma)
res1 <- proclagma(transalt,"store_nbr","transmasd",malags,res1,lagsd)
res1 <- procperiod(transalt,"store_nbr","transpm",pmlags,res1,periodmean)
res1 <- procperiod(transalt,"store_nbr","transpmsd",pmlags,res1,periodsd)

res1 <- proclagma(hcityalt,"city_nbr","hollag",lags,res1,shift)
res1 <- proclagma(hcityalt,"city_nbr","holcnt",malags,res1,lagsum)
res1 <- procgap(hcityalt,"city_nbr","holpregap",res1,lasttruegap)
res1 <- procgap(hcityalt,"city_nbr","holnextgap",res1,nexttruegap)

res1 <- proclagma(mcity,"city_nbr","mcitylag",lags,res1,shift)
res1 <- proclagma(mcity,"city_nbr","mcityma",malags,res1,lagma)
res1 <- procperiod(mcity,"city_nbr","mcitypm",pmlags,res1,periodmean)

res1 <- proclagma(mstate,"state_nbr","mstatelag",lags,res1,shift)
res1 <- proclagma(mstate,"state_nbr","mstatema",malags,res1,lagma)
res1 <- procperiod(mstate,"state_nbr","mstatepm",pmlags,res1,periodmean)

res1 <- proclagma(mtype,"store_type_nbr","mtypelag",lags,res1,shift)
res1 <- proclagma(mtype,"store_type_nbr","mtypema",malags,res1,lagma)
res1 <- procperiod(mtype,"store_type_nbr","mtypepm",pmlags,res1,periodmean)

res1 <- proclagma(mcluster,"cluster","mclusterlag",lags,res1,shift)
res1 <- proclagma(mcluster,"cluster","mclusterma",malags,res1,lagma)
res1 <- procperiod(mcluster,"cluster","mclusterpm",pmlags,res1,periodmean)

res1 <- proclagma(sumpromalt,"store_nbr","spromlag",lags,res1,shift)
res1 <- proclagma(sumpromalt,"store_nbr","spromma",malags,res1,lagma)
res1 <- procperiod(sumpromalt,"store_nbr","sprompm",pmlags,res1,periodmean)

cutdate <- as.Date("2014-08-01")
res1 <- res1[date >= cutdate,]
fwrite(res1,"transfeature,csv")
