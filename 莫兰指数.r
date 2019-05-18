library(xlsx)
library(spdep)
theta=0.3
num=29  #省份
#long lat  经纬度
location<-read.xlsx("经纬度.xlsx",sheetName="2",header=T,encoding='UTF-8')
location=apply(location,2,as.numeric)
lat<-location[,2]
long<-location[,3]
##U数据矩阵
x<-read.xlsx("数据.xlsx",sheetName="unem",header=T,encoding='UTF-8')
x=apply(x,2,as.numeric)
#weight权重矩阵
weight <- matrix(0, nrow=num, ncol=num)
for (i in 1:num){
    for(j in 1:num){
        weight[i,j] = 1/sqrt(2*pi)*exp(-0.5*theta^2*((long[i]-long[j])^2+(lat[i]-lat[j])^2))
    }
}
I=matrix(0,15,2)
for(q in 1:15){
	a<-moran.test(x[,q+1],mat2listw(weight))
	I[q,1]=as.numeric(a$estimate[1])
	I[q,2]=as.numeric(a$p.value)
}
write.xlsx(I,"指数.xlsx")