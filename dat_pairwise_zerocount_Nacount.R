library(tidyr)
library(dplyr)
library(ggplot2)

#Read_data

dat<-read.csv(file = "GSE67758_Raw_gene_wise_quantifications_Human_and_Salmonella.csv", header = TRUE,sep = "\t")
dat=dat[grep("NC",dat$Sequence.name,fixed = TRUE),]
dat<-dat[,-c(1:9,11:28,47:61)]

Attributes <-dat[,c(1)]
Attributes<-str_split_fixed(Attributes,":", n=Inf)
Attributes<-Attributes[,c(2)]
Attributes<-str_split_fixed(Attributes,",", n=Inf)
Attributes<-Attributes[,-c(2)]
dat<-cbind(dat[,1],Attributes,dat[,2:19])
dat<-dat[,-c(1)]
dat=dat[grep("E1W",dat$Attributes,fixed = TRUE),]

dat<-aggregate(.~ Attributes,dat, FUN="mean")
rownames(dat) <- dat[,1]
dat <- dat[,-1]

dat=dat[rowSums(is.na(dat[,-1])) !=ncol(dat[,-1]),]


##code for zero count and NA count

dat?????

count0NA<-function(dat,col){
  x=dat[,col]
  x0=x==0
  xNA=is.na(x)
  n_NA=apply(xNA,1,sum)
  n_0=apply(x0,1,sum)
  out=data.frame(n_0=n_0,n_NA=n_NA)
  out
}

y=count0NA(dat, col = 2:7, )
dat=cbind(dat,zeros02=y[,1],NAs02=y[,2])

##pairwise zero count and NA count

#2h
y=count0NA(dat, col = 2:7)
dat=cbind(dat,zeros02=y[,1],NAs02=y[,2])

data_2h_0h=dat[,c(1:7,35:36)]
write.csv(data_2h_0h, file = "data_2h_0h.csv", row.names = FALSE)


#4h
y=count0NA(dat, c(2:4,8:10))
dat=cbind(dat,zeros04=y[,1],NAs04=y[,2])

data_4h_0h=dat[,c(1:7,37:38)]
write.csv(data_4h_0h, file = "data_4h_0h.csv", row.names = FALSE)


#8h
y=count0NA(dat, c(2:4,11:13))
dat=cbind(dat,zeros08=y[,1],NAs08=y[,2])

data_8h_0h=dat[,c(1:7,39:40)]
write.csv(data_8h_0h, file = "data_8h_0h.csv", row.names = FALSE)


#16h
y=count0NA(dat, c(2:4,14:16))
dat=cbind(dat,zeros016=y[,1],NAs016=y[,2])

data_16h_0h=dat[,c(1:7,41:42)]
write.csv(data_16h_0h, file = "data_16h_0h.csv", row.names = FALSE)


#24h
y=count0NA(dat, c(2:4,17:19))
dat=cbind(dat,zeros024=y[,1],NAs024=y[,2])

data_24h_0h=dat[,c(1:7,43:44)]
write.csv(data_24h_0h, file = "data_24h_0h.csv", row.names = FALSE)


#Write csv containing all columns
dat_pairwise_zerocount_Nacount= dat
write.csv(dat_pairwise_zerocount_Nacount, file = "dat_pairwise_zerocount_Nacount.csv")

dat_pairwise_zerocount_Nacount= read.csv(file = "dat_pairwise_zerocount_Nacount.csv", header = TRUE)
Dat_pairwise_zerocount_Nacount= dat_pairwise_zerocount_Nacount[ , -c(1)]

write.csv(Dat_pairwise_zerocount_Nacount, file = "Dat_pairwise_zerocount_Nacount.csv")

conditions<-gsub('HeLa.S3.','',colnames(dat_pairwise_zerocount_Nacount))
c<-gsub('.R[0-9]','',conditions)
