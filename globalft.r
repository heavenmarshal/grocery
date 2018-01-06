lasttrueidx <- function(idx,logvec)
{
    len <- length(logvec)
    vidx <- 1:len
    ltidx <- max(which(vidx<=idx & logvec))
    return(ltidx)
}
nexttrueidx <- function(idx,logvec)
{
    len <- length(logvec)
    vidx <- 1:len
    ntidx <- min(which(vidx>=idx & logvec))
    return(ntidx)
}
procdate <- function(dseq)
{
    ynum <- year(dseq)-2013
    mon <- month(dseq)
    days <- mday(dseq)
    doy <- yday(dseq)
    dow <- weekdays(dseq)
    dateft <- data.frame(date=dseq,year=ynum,month=mon,
                         day=days,dayofyear=doy,dayofweek=dow)
    return(dateft)
}
## build lookup table for holidays
## processing store information and holiday, since holiday is regional
## holiday type: 0: national 1: state 2: city -1 no holiday

procholiday <- function(citystate,dseq,extdate,holidays)
{
    city <- citystate[1]
    state <- citystate[2]
    listdate <- as.Date(holidays$date)
    inlist <- extdate %in% listdate
    holidx <- match(extdate,listdate)
    locale <- as.character(holidays$locale_name[holidx])
    isnation <- locale == "Ecuador"
    iscity <- locale == city
    isstate <- locale == state
    isholiday <- inlist & (isnation | isstate | iscity)
    isholiday[is.na(isholiday)] <- FALSE
    hollocale <- ifelse(isnation, 1,
                        ifelse(isstate,2,
                               ifelse(iscity,3,0)))
    hollocale[is.na(hollocale)] <- 0
    holtype <- as.integer(holidays$type[holidx])
    holtype[!isholiday] = 0
    transferred <- as.integer(holidays$transferred[holidx])-1
    transferred[!isholiday] = 0

    seqidx <- match(dseq,extdate)
    islag1 <- isholiday[seqidx-1]
    islag2 <- isholiday[seqidx-2]
    islag3 <- isholiday[seqidx-3]
    islag7 <- isholiday[seqidx-7]
    islag14 <- isholiday[seqidx-14]
    islag21 <- isholiday[seqidx-21]
    lidx <- sapply(seqidx,lasttrueidx,isholiday)
    nidx <- sapply(seqidx,nexttrueidx,isholiday)
    last_holiday <- seqidx-lidx
    next_holiday <- nidx - seqidx
    dat <- data.frame(date=dseq,isholiday=isholiday[seqidx],
                      locale=hollocale[seqidx],type=holtype[seqidx],
                      transferred=transferred[seqidx],islag1=islag1,
                      islag2=islag2,islag3=islag3,islag7=islag7,
                      islag14=islag14,islag21=islag21,last_holiday=last_holiday,
                      next_holiday=next_holiday)
    return(dat)
}
