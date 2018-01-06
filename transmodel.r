library("data.table")
library("lightgbm")
trans <- fread("transfeature.csv")
trans[,date := as.Date(date)]
cutdate <- as.Date("2017-08-15")
ttrain <- trans[date<= cutdate,]
cutblow <- as.Date("2015-08-01")
ttrain <- ttrain[date>=cutblow,]
lab <- ttrain[["transactions"]]
dropname <- c("date","transactions")
ttrain <- ttrain[,-dropname,with=FALSE]
ntrain <- nrow(ttrain)
cate_vars <- c("store_nbr","cluster","store_type_nbr","state_nbr","city_nbr",
               "htype_nbr","locale_nbr","dow")
param <- list(learning_rate=.05,objective="regression",metric="l2_root")

nsample <- 10
nlev <- c(31,15)
lambdal1 <- c(.5,1.5)
folder <- "transmodel/"
nmd <- nsample*length(nlev)*length(lambdal1)
scoredf <- data.frame(file=rep("",nmd),score=rep(0,nmd),stringsAsFactors=FALSE)
l <- 1
for(i in 1:nsample)
{
    valididx <- sample(1:ntrain,ceiling(.2*ntrain))
    train <- ttrain[-valididx,]
    valid <- ttrain[valididx,]
    dtrain <- lgb.Dataset(data=as.matrix(train),label=lab[-valididx],categorical_feature=cate_vars)
    dvalid <- lgb.Dataset(data=as.matrix(valid),label=lab[valididx],categorical_feature=cate_vars)
    for(j in nlev)
    {
        for(k in lambdal1)
        {
            param$num_leaves=j
            param$lambda_l1=k
            lgbmd <- lgb.train(param,dtrain,nrounds=5000,valids=list(test=dvalid),early_stopping_round=40,verbose=0)
            fname <- paste(folder,"sam",i,"lev",j,"lamb",k,".lgb",sep="")
            scoredf[l,] <- c(fname,lgbmd$eval_valid()[[1]]$value)
            lgb.save(lgbmd,fname)
            l <- l+1
        }
    }
}
scoredf <- data.table(scoredf)
fwrite(scoredf,paste(folder,"score",".csv",sep=""))
