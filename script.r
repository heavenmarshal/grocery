library("data.table")
library("zoo")
library("RSQLite")
source("utility.r")
con <- dbConnect(SQLite(),dbname="sqldb/grocery.sqlite")
stores <- fread("stores.csv/storerv.csv")
items <- fread("items.csv/items.csv")
holcity <- fread("holidays_events.csv/holcity.csv")
trans <- fread("transactions.csv/transactions.csv")
## remove duplicated (caused by both nation and local)
redstore <- stores[,.(store_nbr,city_nbr,state_nbr,store_type_nbr,cluster)]
holcity <- holcity[!duplicated(holcity,by=c("date","city_nbr")),]

## global variables
lags <- c(1,2,3,7,14,21,28,35)
malags <- c(2,3,7,14,21,28,35,70,140)
pmlags <- c(14,21,28,35,70,140,210)
## process holiday
hdseq <- seq.Date(as.Date("2012-03-02"),as.Date("2017-12-26"),1)
hcityal <- alignmts(hdseq,holcity,id="city_nbr",field="isholiday")
hcityalt <- melt(hcityal,value.name="isholiday",id.vars="city_nbr",
                 variable.name="date",variable.factor=FALSE,value.factor=FALSE)
hcityalt <- dcast(hcityalt,date~city_nbr,value.var="isholiday")
itemno <- 103665

rs <- dbSendQuery(con,sqlquery)
item <- data.table(fetch(rs,n=Inf))
item[,onpromotion := procpromotion(onpromotion)]
dseq <- seq.Date(as.Date("2013-01-01"),as.Date("2017-08-15"),1)
transal <- alignmts(dseq,trans,"store_nbr","transactions")
transstore <- melt(transal,value.name="transactions",id.vars="store_nbr",
                   variable.name="date",variable.factor=FALSE,value.factor=FALSE)
transalt <- dcast(transstore,date~store_nbr)
itemalsale <- alignmts(dseq,item)
itemalprom <- alignmts(dseq,item,field="onpromotion")
itemprom <- melt(itemalprom,value.name="onpromotion",id.vars="store_nbr",
                    variable.name="date",variable.factor=FALSE,value.factor=FALSE)

res <- melt(itemalsale,value.name="unit_sales",id.vars="store_nbr",
            variable.name="date",variable.factor=FALSE,value.factor=FALSE)
itemalsalet <- dcast(res,date~store_nbr,value.var="unit_sales")

procitem(itemno,con,dseq,transalt,stores,hcityalt,lags,malags,pmlags,folder)
