item <- fread("inbr1009998.csv")
ttrain <- item[date<as.Date("2017-08-16"),]
lab <- ttrain[["unit_sales"]]
ttrain <- ttrain[,-c("date","unit_sales")]
ntrain <- nrow(ttrain)

valididx <- sample(1:ntrain,ceiling(.2*ntrain))
train <- ttrain[-valididx,]
valid <- ttrain[valididx,]
dtrain <- lgb.Dataset(data=as.matrix(train),label=lab[-valididx],categorical_feature=c("store_nbr"))
dvalid <- lgb.Dataset(data=as.matrix(valid),label=lab[valididx],categorical_feature=c("store_nbr"))
param <- list(learning_rate=.05,objective="regression",metric="l2_root")
lgbmd <- lgb.train(param,dtrain,nrounds=500,valids=list(test=dvalid),early_stoping_round=20)

