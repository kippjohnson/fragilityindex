# data downloaded from http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data

mydata <- read.csv("data-raw/heartdisease.txt")
colnames(mydata) =c("age","sex","cp","trestbps","chol","fbs","restecg","thalach",
                    "exang","oldpeak","slope","ca","thal","num")
mydata$ca=as.numeric(mydata$ca)
mydata$thal=as.numeric(mydata$thal)
heart_disease <- mydata

devtools::use_data(heart_disease, overwrite = TRUE)
