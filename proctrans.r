wmean <- function(vec,weight)
{
    wmean <- sum(vec*weight)/sum(weight)
}
ensemblepred <- function(modellist,scores,predset)
{
    date <- predset[["date"]]
    store_nbr <- predset[["store_nbr"]]
    predset <- predset[,-c("date","transactions")]
    py <- lapply(modellist,predict,as.matrix(predset))
    py <- matrix(unlist(py),ncol=length(modellist))
    df <- data.table(file=names(modellist))
    df <- merge(df,scores,all.x=TRUE,by="file")
    weight <- 1/df[["score"]]
    py <- apply(py,1,wmean,weight)
    pydf <- data.table(date=date,store_nbr=store_nbr,transactions=py)
}
predtrans <- function(initft,preddate,modellist,scores,dseq,origdf,
                      ftdf,stores,lags,malags,pmlags)
{
    predft <- initft
    py <- ensemblepred(modellist,scores,predft)
    pred <- copy(py)
    len <- length(preddate)
    idx <- 2
    while(idx <= len)
    {
        newdate <- preddate[idx]
        predft <- gennewfeature(dseq,origdf,py,newdate,ftdf,stores,lags,malags,pmlags)
        py <- ensemblepred(modellist,scores,predft)
        pred <- rbind(pred,py)
        idx <- idx+1
    }
    return(pred)
}
gennewfeature <- function(dseq,origdf,pred,newdate,ftdf,stores,lags,malags,pmlags)
{
    origdf[,date := as.character(date)]
    pred[,date := as.character(date)]
    newdf <- rbind(origdf,pred)
    transal <- alignmts(dseq,newdf,field="transactions")
    res <- melt(transal,value.name="transactions",id.var="store_nbr",
                variable.name="date",variable.factor=FALSE,value.factor=FALSE)
    res[,date:=as.Date(date)]
    transalt <- dcast(res,date~store_nbr)
    mcity <- summarizeCate(transal,"store_nbr","city_nbr",stores,mean)
    mstate <- summarizeCate(transal,"store_nbr","state_nbr",stores,mean)
    mtype <- summarizeCate(transal,"store_nbr","store_type_nbr",stores,mean)
    mcluster <- summarizeCate(transal,"store_nbr","cluster",stores,mean)
    substore <- stores[,c("store_nbr","city_nbr","state_nbr","store_type_nbr","cluster")]
    res <- merge(res,substore,by="store_nbr",all.x=TRUE)
    res <- proclagma(transalt,"store_nbr","translag",lags,res,shift)
    res <- proclagma(transalt,"store_nbr","transma",malags,res,lagma)
    res <- proclagma(transalt,"store_nbr","transmasd",malags,res,lagsd)
    res <- procperiod(transalt,"store_nbr","transpm",pmlags,res,periodmean)
    res <- procperiod(transalt,"store_nbr","transpmsd",pmlags,res,periodsd)

    res <- proclagma(mcity,"city_nbr","mcitylag",lags,res,shift)
    res <- proclagma(mcity,"city_nbr","mcityma",malags,res,lagma)
    res <- procperiod(mcity,"city_nbr","mcitypm",pmlags,res,periodmean)
    res <- proclagma(mstate,"state_nbr","mstatelag",lags,res,shift)
    res <- proclagma(mstate,"state_nbr","mstatema",malags,res,lagma)
    res <- procperiod(mstate,"state_nbr","mstatepm",pmlags,res,periodmean)

    res <- proclagma(mtype,"store_type_nbr","mtypelag",lags,res,shift)
    res <- proclagma(mtype,"store_type_nbr","mtypema",malags,res,lagma)
    res <- procperiod(mtype,"store_type_nbr","mtypepm",pmlags,res,periodmean)

    res <- proclagma(mcluster,"cluster","mclusterlag",lags,res,shift)
    res <- proclagma(mcluster,"cluster","mclusterma",malags,res,lagma)
    res <- procperiod(mcluster,"cluster","mclusterpm",pmlags,res,periodmean)
    res <- res[date==as.Date(newdate),]
    remnames <- c("date","store_nbr",setdiff(names(ftdf),names(res)))
    resrem <- ftdf[date==as.Date(newdate),remnames,with=FALSE]
    res <- merge(res,resrem,by=c("date","store_nbr"),all.x=TRUE)
    setcolorder(res,names(ftdf))
    return(res)
}
