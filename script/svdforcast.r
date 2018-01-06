svdforcast <- function(tsmat,nsv,npred,start=c(2013,1,1),freq=365)
{
    svobj <- svd(tsmat)
    umat <- svobj$u[,1:nsv]
    dvec <- svobj$d[1:nsv]
    vmat <- svobj$v[,1:nsv]
    pumat <- matrix(nrow=npred,ncol=nsv)
    for(i in 1:nsv)
    {
        tsu <- ts(umat[,i],start=start,frequency=freq)
        armatsu <- auto.arima(tsu)
        fortsu <- forecast(armatsu,h=npred)
        pumat[,i] <- fortsu$mean
    }
    predmat <- pumat%*%diag(dvec)%*%t(vmat)
    return(predmat)
}
