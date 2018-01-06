lagma <- function(vec,lag)
{
    len <- length(vec)
    out <- .C("lagsum",as.double(vec),as.integer(len),
              as.integer(lag),ans=double(len))
    ans <- out$ans/lag
    ans[1:lag] <- NA
    return(ans)
}
lagmax <- function(vec,lag)
{
    lmx <- rollmax(vec,lag,na.pad=TRUE,align="right")
    lmx <- shift(lmx,1)
    return(lmx)
}
lagsd <- function(vec,lag)
{
    len <- length(vec)
    ## length must be > 1, for speed consideration, do not check it.
    out <- .C("lagvar",as.double(vec),as.integer(len),
              as.integer(lag), ans = double(len))
    ans <- sqrt(out$ans)
    ans[1:lag] <- NA
    return(ans)
}
lagsum <- function(vec,lag)
{
    len <- length(vec)
    out <- .C("lagsum",as.double(vec),as.integer(len),
              as.integer(lag), ans = double(len))
    ans <- out$ans
    ans[1:lag] <- NA
    return(ans)
}
periodmean <- function(vec,start,end,by)
{
    len <- length(vec);
    out <- .C("periodsum",as.double(vec),as.integer(len),
              as.integer(start),as.integer(end),
              as.integer(by), rend=integer(1),
              slen = integer(1), ans = double(len))
    ans <- out$ans
    ans <- ans/out$slen
    ans[1:out$rend] <- NA
    return(ans)
}
periodsd <- function(vec,start,end,by)
{
    len <- length(vec)
    out <- .C("periodvar",as.double(vec),as.integer(len),
              as.integer(start),as.integer(end),
              as.integer(by),rend=integer(1),
              slen = integer(1), ans = double(len))
    ans <- sqrt(out$ans)
    ans[1:out$rend] <- NA
    return(ans)
}
periodapply <- function(vec,start,end,by,fun,...)
{
    len <- length(vec)
    lags <- seq(start,end,by)
    shifv <- unlist(shift(vec,lags))
    shifv <- matrix(shifv,nrow=len)
    redperiod <- apply(shifv,1,fun,...)
    return(redperiod)
}
summarytsmat <- function(itemts,stores,byvar,fun)
{
    stores <- as.data.frame(stores)
    storeidx <- itemts$stores
    ctype <- stores[storeidx,byvar]
    df <- data.table(t(itemts$tsmat))
    df[,byvar] <- ctype
    dat <- df[,lapply(.SD,mean),by=byvar]
    return(dat)
}
procpromotion <- function(vprom)
{
    oprom <- ifelse(vprom=="True",1,0)
}
alignmts <- function(dseq,tabitem,id="store_nbr",field="unit_sales")
{
    dummy <- data.table(date=as.character(dseq))
    dummy[,(id):=-1]
    form <- formula(paste(id,"~","date",sep=""))
    itemrv <- rbind(tabitem,dummy,use.names=TRUE,fill=TRUE)
    itemrv <- dcast(itemrv,formula=form,value.var=field,fill=0)
    residx <- drop(itemrv[,id,with=FALSE] != -1)
    itemrv <- itemrv[residx,]
    return(itemrv)
}
alignmtshol <- function(dseq,tabitem,id="city_nbr",field="isholiday")
{
    dummy <- data.table(date=as.character(dseq))
    dummy[,(id):=-1]
    form <- formula(paste(id,"~","date",sep=""))
    itemrv <- rbind(tabitem,dummy,use.names=TRUE,fill=TRUE)
    itemrv <- dcast(itemrv,formula=form,fun.aggregate=max,value.var=field,fill=0)
    residx <- drop(itemrv[,id,with=FALSE] != -1)
    itemrv <- itemrv[residx,]
    return(itemrv)
}

summarizeCate <- function(alignmat,id,cate,dict,fun)
{
    idx <- match(alignmat[[id]],dict[[id]])
    mtcate <- dict[idx,cate,with=FALSE]
    dat <- cbind(mtcate,alignmat[,-id,with=FALSE])
    sumdat <- dat[,lapply(.SD,fun),by=cate]
    sumdat <- melt(sumdat,id.vars=cate,variable.name="date")
    form <- formula(paste("date~",cate,sep=""))
    sumdat <- dcast(sumdat,formula=form)
    return(sumdat)
}
applytocate <- function(itemalt,valname,varname,vartype,fun,...)
{
    cols <- setdiff(names(itemalt),"date")
    res <- itemalt[,lapply(.SD,fun,...),.SDcols=cols]
    res[,date := itemalt[,date]]
    transfun <- get(paste("as.",vartype,sep=""))
    res <- melt(res,value.name=valname,variable.name=varname,id.vars="date",
                variable.factor=FALSE,value.factor=FALSE)
    res[,date := as.Date(date)]
    res[,(varname) := lapply(.SD,transfun),.SDcols=varname]
    return(res)
}
holidaycitytable <- function(holidays,citytab)
{
    vars <- citytab[["city_nbr"]]
    setkey(citytab,"city_nbr")
    ret <- holidays[,c("date","type","locale","locale_name","transferred")]
    for(vari in vars)
    {
        lset <- c(as.character(citytab[.(vari),.(city,state)]),"Ecuador")
        ret[,(as.character(vari)):= as.integer(locale_name %in% lset)]
    }
    ret[,locale_name := NULL]
    idvar <- setdiff(names(ret),as.character(vars))
    ret <- melt(ret,value.name="isholiday",id.vars=idvar,variable.name="city_nbr",
                variable.factor=FALSE,value.factor=FALSE)
    ret[,type:=ifelse(isholiday,type,"not_in_city")]
    ret[,locale:=ifelse(isholiday,locale,"not_in_city")]
    ret[,transferred:=ifelse(isholiday,transferred,FALSE)]
    return(ret)
}
lasttruegap <- function(vec)
{
    len <- length(vec)
    out <- .C("lasttruegap",as.integer(vec),as.integer(len),
              ans=integer(len))
    ans <- out$ans
    ans <- ifelse(ans>=0,ans,NA)
    return(ans)
}
nexttruegap <- function(vec)
{
    rvec <- rev(vec)
    ans <- rev(lasttruegap(rvec))
    return(ans)
}
proclagma <- function(df,id,opname,lags,res,fun)
{
    df[,date := as.Date(date)]
    lenlag <- length(lags)
    cols <- setdiff(names(df),"date")
    for(i in 1:lenlag)
    {
        lag <- lags[i]
        newdf <- df[,lapply(.SD,fun,lag),.SDcols=cols]
        newdf[,date:=df[,date]]
        vname <- paste(opname,lag,sep="_")
        newdf <- melt(newdf,value.name=vname,id.vars="date",variable.name=id,
                      variable.factor=FALSE,value.factor=FALSE)
        newdf[,(id):= lapply(.SD,as.integer),.SDcols=id]
        newdf[,date:=as.Date(date)]
        res <- merge(res,newdf,by=c("date",id),all.x=TRUE)
    }
    return(res)
}
## period mean or period sd
procperiod <- function(df,id,opname,lags,res,fun)
{
    df[,date := as.Date(date)]
    lenlag <- length(lags)
    cols <- setdiff(names(df),"date")
    for(i in 1:lenlag)
    {
        lag <- lags[i]
        newdf <- df[,lapply(.SD,fun,7,lag,7),.SDcols=cols]
        newdf[,date:=df[,date]]
        vname <- paste(opname,lag,sep="_")
        newdf <- melt(newdf,value.name=vname,id.vars="date",variable.name=id,
                      variable.factor=FALSE,value.factor=FALSE)
        newdf[,(id):= lapply(.SD,as.integer),.SDcols=id]
        newdf[,date:=as.Date(date)]
        res <- merge(res,newdf,by=c("date",id),all.x=TRUE)
    }
    return(res)
}
procgap <- function(df,id,opname,res,fun)
{
    cols <- setdiff(names(df),"date")
    newdf <- df[,lapply(.SD,fun),.SDcols=cols]
    newdf[,date:=df[,date]]
    newdf <- melt(newdf,value.name=opname,id.vars="date",variable.name=id,
                  variable.factor=FALSE,value.factor=FALSE)
    newdf[,(id):= lapply(.SD,as.integer),.SDcols=id]
    newdf[,date:=as.Date(date)]
    newdf[,(id):= lapply(.SD,as.integer),.SDcols=id]
    res <- merge(res,newdf,by=c("date",id),all.x=TRUE)
    return(res)
}
