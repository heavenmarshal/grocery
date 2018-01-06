procitem <- function(itemno,con,dseq,transalt,stores,hcityalt,lags,malags,pmlags,folder)
{
    sqlfmt <- "select * from train where item_nbr=%d"
    sqlquery <- sprintf(sqlfmt,itemno)
    rs <- dbSendQuery(con,sqlquery)
    item <- data.table(fetch(rs,n=Inf))
    item[,onpromotion := procpromotion(onpromotion)]
    itemalsale <- alignmts(dseq,item)
    itemalprom <- alignmts(dseq,item,field="onpromotion")
    itemprom <- melt(itemalprom,value.name="onpromotion",id.vars="store_nbr",
                     variable.name="date",variable.factor=FALSE,value.factor=FALSE)

    res <- melt(itemalsale,value.name="unit_sales",id.vars="store_nbr",
                variable.name="date",variable.factor=FALSE,value.factor=FALSE)
    itemalsalet <- dcast(res,date~store_nbr,value.var="unit_sales")
    res <- merge(res,itemprom,by=c("date","store_nbr"),all.x=TRUE)
    res[,date := as.Date(date)]
    res[,year:=year(date)-2013]
    res[,month:=month(date)]
    res[,doy:=yday(date)]
    res[,dow:=wday(date)]
    redstore <- stores[,.(store_nbr,city_nbr,state_nbr,store_type_nbr,cluster)]
    res <- merge(res,redstore,by="store_nbr",all.x=TRUE)

    itemalpromt <- dcast(itemprom,date~store_nbr,value.var="onpromotion")
    mcitysale <- summarizeCate(itemalsale,"store_nbr","city_nbr",stores,mean)
    mstatesale <- summarizeCate(itemalsale,"store_nbr","state_nbr",stores,mean)
    mtypesale <- summarizeCate(itemalsale,"store_nbr","store_type_nbr",stores,mean)
    mclussale <- summarizeCate(itemalsale,"store_nbr","cluster",stores,mean)
    res1 <- proclagma(itemalsalet,"store_nbr","salelag",lags,res,shift)
    res1 <- proclagma(itemalsalet,"store_nbr","salema",malags,res1,lagma)
    res1 <- proclagma(itemalsalet,"store_nbr","salemasd",malags,res1,lagsd)
    res1 <- procperiod(itemalsalet,"store_nbr","salepm",pmlags,res1,periodmean)
    res1 <- procperiod(itemalsalet,"store_nbr","salepmsd",pmlags,res1,periodsd)
    ## process holiday
    res1 <- proclagma(hcityalt,"city_nbr","hollag",lags,res1,shift)
    res1 <- proclagma(hcityalt,"city_nbr","holcnt",malags,res1,lagsum)
    res1 <- procgap(hcityalt,"city_nbr","holpregap",res1,lasttruegap)
    res1 <- procgap(hcityalt,"city_nbr","holnextgap",res1,nexttruegap)
    ## process promotion
    res1 <- proclagma(itemalpromt,"store_nbr","promlag",lags,res1,shift)
    res1 <- proclagma(itemalpromt,"store_nbr","promcnt",malags,res1,lagsum)
    res1 <- procgap(itemalpromt,"store_nbr","prompregap",res1,lasttruegap)
    res1 <- procgap(itemalpromt,"store_nbr","prompnextgap",res1,nexttruegap)
    ## process city
    res1 <- proclagma(mcitysale,"city_nbr","mcitylag",lags,res1,shift)
    res1 <- proclagma(mcitysale,"city_nbr","mcityma",malags,res1,lagma)
    res1 <- procperiod(mcitysale,"city_nbr","mcitypm",pmlags,res1,periodmean)
    ## process state
    res1 <- proclagma(mstatesale,"state_nbr","mstatelag",lags,res1,shift)
    res1 <- proclagma(mstatesale,"state_nbr","mstatema",malags,res1,lagma)
    res1 <- procperiod(mstatesale,"state_nbr","mstatepm",pmlags,res1,periodmean)
    ## process type
    res1 <- proclagma(mtypesale,"store_type_nbr","mtypelag",lags,res1,shift)
    res1 <- proclagma(mtypesale,"store_type_nbr","mtypema",malags,res1,lagma)
    res1 <- procperiod(mtypesale,"store_type_nbr","mtypepm",pmlags,res1,periodmean)
    ## process cluster
    res1 <- proclagma(mclussale,"cluster","mclusterlag",lags,res1,shift)
    res1 <- proclagma(mclussale,"cluster","mclusterma",malags,res1,lagma)
    res1 <- procperiod(mclussale,"cluster","mclusterpm",pmlags,res1,periodmean)
    ## process transactions
    res1 <- proclagma(transalt,"store_nbr","translag",lags,res1,shift)
    res1 <- proclagma(transalt,"store_nbr","transma",malags,res1,lagma)
    res1 <- procperiod(transalt,"store_nbr","transpm",pmlags,res1,periodmean)
    fname <- paste("inbr",itemno,".csv",sep="")
    fname <- paste(folder,fname,sep="/")
    fwrite(res1,fname)
}
