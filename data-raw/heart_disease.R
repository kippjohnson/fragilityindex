# data downloaded from http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data

library(car)

mydata <- read.csv("data-raw/heartdisease.txt")
colnames(mydata) =c("age","sex","cp","trestbps","chol","fbs","restecg","thalach",
                    "exang","oldpeak","slope","ca","thal","num")
mydata$ca=as.numeric(mydata$ca)
mydata$thal=as.numeric(mydata$thal)
mydata$num <- recode(mydata$num, "1:4=1")
heart_disease <- mydata

devtools::use_data(heart_disease, overwrite = TRUE)
