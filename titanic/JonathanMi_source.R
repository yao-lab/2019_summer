#?read.csv
#Useless functions: sum()..True, ncol(), names()
setwd("~/Hong Kong Summer")
titanic.train<-read.csv("train.csv",stringsAsFactors = FALSE,header=TRUE)
titanic.test<-read.csv("test.csv",stringsAsFactors = FALSE,header=TRUE)


#Data Seperation
titanic.train$Test<-TRUE
titanic.test$Test<-FALSE
titanic.test$Survived<-NA
titanic.full<-rbind(titanic.train,titanic.test)

#Embark Cleaning
table(titanic.full$Embarked)
titanic.full[titanic.full$Embarked=='', "Embarked"]<-"S"


#Age Cleaning
#age.med<-median(titanic.full$Age, na.rm=TRUE)
#titanic.full[is.na(titanic.full$Age), "Age"]<-age.med

boxplot.stats(titanic.full$Age)
upper.whisker.age<-boxplot.stats(titanic.full$Age)$stats[5]
age.filter<-titanic.full$Age<upper.whisker.age

#age linear model
str(titanic.full)
age.equation<-"Age~Pclass+SibSp"
age.model<-lm(formula=age.equation,data=titanic.full[age.filter,])
age.row<-titanic.full[is.na(titanic.full$Age),c("Pclass","SibSp")]
age.predictions<-predict(age.model, newdata=age.row)
titanic.full[is.na(titanic.full$Age),"Age"]<-age.predictions

#Fare Cleaning
#fare.med<-median(titanic.full$Fare, na.rm=TRUE)
#titanic.full[is.na(titanic.full$Fare), "Fare"]<-fare.med
boxplot(titanic.full$Fare)
upper.whisker<-boxplot.stats(titanic.full$Fare)$stats[5]
fare.filter<-titanic.full$Fare<upper.whisker
titanic.full[fare.filter,]

#Fare linear model
str(titanic.full)
fare.equation<-"Fare~Pclass+Sex+Age+SibSp+Parch+Embarked"
fare.model<-lm(formula=fare.equation,data=titanic.full[fare.filter,])
fare.row<-titanic.full[is.na(titanic.full$Fare), c("Pclass","Sex","Age","SibSp","Parch","Embarked")]
fare.predictions<-predict(fare.model, newdata=fare.row)
titanic.full[is.na(titanic.full$Fare),"Fare"]<-fare.predictions

#Categorical Casting
titanic.full$Pclass<-as.factor(titanic.full$Pclass)
titanic.full$Sex<-as.factor(titanic.full$Sex)
titanic.full$Embarked<-as.factor(titanic.full$Embarked)

#Data Seperation
titanic.train<-titanic.full[titanic.full$Test==TRUE,]
titanic.test<-titanic.full[titanic.full$Test==FALSE,]
titanic.train$Survived<-as.factor(titanic.train$Survived)
str(titanic.train)

#randomForest modeling
survived.equation<-"Survived~Pclass+Sex+Age+SibSp+Parch+Fare+Embarked"
survived.formula<-as.formula(survived.equation)
install.packages("randomForest")
library(randomForest)
titanic.model<-randomForest(formula = survived.formula, data = titanic.train, ntree = 500, mtry = 3, nodesize = 0.01 * nrow(titanic.test))
#featured.equation<-"Pclass+Sex+Age+SibSp+Parch+Fare+Embarked"
Survived<-predict(titanic.model, newdata=titanic.test)

PassengerId<-titanic.test$PassengerId
titanic.data.frame<-as.data.frame(PassengerId)
titanic.data.frame$Survived<-Survived
write.csv(titanic.data.frame, file="kaggle titanic csv 3", row.names = FALSE)


importance(titanic.model)
print(titanic.model)

