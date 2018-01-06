cat("Setting up Enviornment\n")
library(data.table)
library(lightgbm)
library(dplyr)

cat("Loading test and train")
train <- fread("../train.csv/train.csv", sep=",", na.strings="", skip=110000000,
               col.names=c("id","date","store_nbr","item_nbr","unit_sales","onpromotion"))
test <- fread("../test.csv/test.csv", sep=",", na.strings = "")


cat("Setting onpromotion to integer")
train[, ":=" ( onpromotion = as.integer(onpromotion)
               ,unit_sales = log1p(ifelse(unit_sales>0, unit_sales, 0)) )]
test[, onpromotion := as.integer(onpromotion)]

df_2017 <- train[date >= "2017-05-31"]

# Making wide table of onpromotion information
vars = c("item_nbr", "store_nbr", "date","onpromotion")
promo_2017_train <- dcast(df_2017[,vars,with=F], store_nbr+item_nbr~date, value.var="onpromotion")
promo_2017_test <- dcast(test[,vars,with=F], store_nbr+item_nbr~date, value.var="onpromotion")
promo_2017 <- left_join(promo_2017_train, promo_2017_test,by=c("store_nbr","item_nbr")) %>% as.data.table
promo_2017[is.na(promo_2017)] = 0
rm(promo_2017_train, promo_2017_test)

# Make wide table of unit_sales information
df_2017 <- dcast(df_2017, store_nbr+item_nbr ~ date, value.var = "unit_sales")
df_2017[is.na(df_2017)] = 0

# - Functions
prepare_dataset <- function(t2017, is_train = T) {
  m3_names = as.character(seq.Date(t2017-3, by = "days", length.out = 3))
  m7_names = as.character(seq.Date(t2017-7, by = "days", length.out = 7))
  m14_names = as.character(seq.Date(t2017-14, by = "days", length.out = 14))
  p14_names = as.character(seq.Date(t2017-14, by = "days", length.out = 14))

  X <- data.table( store_nbr = df_2017[, store_nbr]
                   ,item_nbr = df_2017[, item_nbr]
                   ,mean_3_2017 = df_2017[, rowMeans(.SD), .SDcols = m3_names]
                   ,mean_7_2017 = df_2017[,rowMeans(.SD), .SDcols = m7_names]
                   ,mean_14_2017 = df_2017[,rowMeans(.SD), .SDcols = m14_names]
                   ,promo_14_2017 = promo_2017[,rowSums(.SD), .SDcols = p14_names]
  )
  # Get onpromotion information for t2017 and the following 15 days.
  for(i in 0:15) {
    new_var = paste0("promo_",i)
    var = as.character(t2017+i)
    X[, (new_var) := promo_2017[, var, with=F]]
  }
  # Get unit_sales information for t2017 and following 15 days
  if(is_train) {
    y_dates = as.character(seq.Date(t2017, by="days", length.out = 16))
    y = df_2017[,c("store_nbr", "item_nbr", y_dates), with=F]
    colnames(y) = c("store_nbr", "item_nbr", "y_1","y_2","y_3","y_4","y_5","y_6","y_7","y_8",
                    "y_9","y_10","y_11","y_12","y_13","y_14","y_15","y_16")
    return(list(X = X, y = y))
  } else {
    return(list(X = X))
  }
}

cat("Making X_train and y_train\n")
t2017 = as.Date("2017-6-21")
# Four 'sets' of train data. With i=0, we have:
# t2017 = 2017-06-21
# X will have unit_sales info from 2017-06-07 --> 2017-06-20
# X will have onpromotion info through 2017-07-06
# Y will have unit_sales info from 2017-06-21 --> 2017-07-06
for(i in 0:3) {
  # Make X_tmp and Y_tmp
  delta <- 7*i
  results <- prepare_dataset(t2017+delta)
  X_tmp = results$X
  y_tmp = results$y
  # Concatenating X_l
  if(i == 0) {
    X_train <- X_tmp
    y_train <- y_tmp
  } else {
    X_train = rbindlist(list(X_train, X_tmp))
    y_train = rbindlist(list(y_train, y_tmp))
  }
}
val_date <- as.Date("2017-07-26")
results = prepare_dataset(val_date)
X_val = results$X
y_val = results$y
# Can't produce results$y since that is what we are predicting
results = prepare_dataset(as.Date("2017-08-16"), is_train = F)
X_test = results$X
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Preparing data for model
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat("Ordering all data.tables")
# Train set
X_train <- X_train[order(store_nbr, item_nbr)]
y_train <- y_train[order(store_nbr, item_nbr)]
# Val set
X_val <- X_val[order(store_nbr, item_nbr)]
y_val <- y_val[order(store_nbr, item_nbr)]
# Test set
X_test <- X_test[order(store_nbr,item_nbr)]

cat("Initializing prediction tables")
val_pred <- X_val[,.(store_nbr,item_nbr)]
test_pred <- X_test[,.(store_nbr,item_nbr)]

cat("Setting up X and Y matrices")
X_train_in = as.matrix(X_train[,-c("store_nbr", "item_nbr"), with=F])
X_val_in = as.matrix(X_val[,-c("store_nbr", "item_nbr"), with=F])
X_test_in = as.matrix(X_test[,-c("store_nbr","item_nbr"), with=F])

cat("Cleaning up before building model")
rm(train, df_2017, promo_2017); gc()
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Build Model and Getting Predictions
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat("Setting Parameters\n")
MAX_ROUNDS = 1000
params <- list(num_leaves = 2^5-2
               ,objective = "regression_l2"
               ,max_depth = 8
               ,min_data_in_leaf = 50 # Might be min_data
               ,learning_rate = .05
               ,feature_fraction = .75
               ,bagging_fraction = .75
               ,bagging_freq = 1
               ,metric = "l2"
               ,num_threads = 4
)

cat("Building Model")
for(i in 0:15) {
  cat("Round:",i,"of 15\n")
  y_col = paste0("y_",i+1)
  dtrain <- lgb.Dataset(X_train_in, label = as.matrix(y_train[, y_col, with=F]))
  dval <- lgb.Dataset(X_val_in, label = as.matrix(y_val[, y_col, with=F]), reference=dtrain)
  bst <- lgb.train(params, dtrain, nrounds = MAX_ROUNDS
                   ,valids = list(test=dval), early_stopping_rounds = 50
                   , verbose = 0
  )
  # Make predictions
  val_pred[, eval(y_col) := bst$predict(X_val_in, num_iteration = bst$best_iter)]
  test_pred[, eval(y_col) := bst$predict(X_test_in, num_iteration = bst$best_iter)]
}

# colnames are originally y_1, y_2, ... y_16. Switching to test dates for submission
colnames(test_pred) <- c("store_nbr", "item_nbr", sort(unique(test$date)))
# Make wide table long (columns are: store_nbr, item_nbr, date, unit_sales)
test_pred_long <- melt(test_pred, measure.vars = sort(unique(test$date)),
                       variable.name="date", value.name="unit_sales")
test_pred_long[, dates := as.Date(date)]
test[, unit_sales := 0]
# Merge with original test set
test[test_pred_long, unit_sales := i.unit_sales, on=c("store_nbr","item_nbr", "date")]
test[unit_sales<0, unit_sales:=0]
test[, unit_sales := expm1(unit_sales)]

cat("Making submission")
fwrite(test[,.(id, unit_sales)], "LGBM_sub.csv", sep=",", dec=".", quote=FALSE, row.names=FALSE)
