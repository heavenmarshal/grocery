alignfun <- function(date,sales,dseq)
{
    date <- as.Date(date)
    isin <- dseq %in% date
    tlen <- length(dseq)
    vec <- rep(0,tlen)
    vec[isin] <- sales
    return(vec)
}
