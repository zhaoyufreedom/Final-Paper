library(MASS) ###求逆的时候用
m=8
n=m*m  ###地点
T=5   ###时间（5年，设置了5组变系数）
k=50 ##循环次数
l=4  ##beta的维度

num=n*T

###个体效应(保证可识别，求和为0，所以先生成n-1个，再用0减)
alpha<-matrix(0,num,1)
c=matrix(rnorm(n,0,1),nrow=n,ncol=1)
for (i in 1:n){
	for (t in 1:T){
		alpha[t+(i-1)*T]=c[i]
	}
}

###时间效应(保证可识别，求和为0，所以先生成T-1个，再用0减)
lamda<-matrix(0,num,1)
c<-matrix(rnorm(T,0,0.5),nrow=T,ncol=1)
for (i in 1:n){
	for (t in 1:T){
		lamda[t+(i-1)*T]=c[t]
	}
}

###第一种生成方式个体效应与时间效应的综合效应真值
zonghe1<-alpha+lamda

##格子点坐标
lat <- as.matrix(c(1:n-1)%/%m)
long<- as.matrix(c(1:n-1)%%m)

##把每个格子点重复T次，为变系数部分的估计做准备
varlat<-kronecker(lat,matrix(1,T,1))
varlong<-kronecker(long,matrix(1,T,1))


##T=5,5个变系数。变系数部分，同样用0-前四个函数，保证解唯一
lamda1<-as.matrix(1/6*(lat+long))
lamda2<-as.matrix(1/3*lat)
lamda3<-as.matrix(4*sin(1/12*pi*lat))
lamda4<-as.matrix(sin(pi*lat/(sqrt(n)-1)))
lamda5<-as.matrix(2/(sqrt(n)-1)*sqrt(((sqrt(n)-1)/2-abs((sqrt(n)-1)/2-lat))*((sqrt(n)-1)/2-abs((sqrt(n)-1)/2-long))))

com<-cbind(lamda1,lamda2,lamda3,lamda4,lamda5)


###第二种生成方式个体效应与时间效应的综合效应真值
zonghe2<-matrix(0,num,1)
for (i in 1:n){
	for (t in 1:T){
	zonghe2[(i-1)*T+t]=com[i,t]
}
}


####数据处理，加入虚拟变量，M1代表个体效应的系数,M2表示时间效应的系数（恰好与变系数函数的系数是一样的）
M1<-matrix(0,num,n)
M2<-matrix(0,num,T)
for (i in 1:n){
	for (j in 1:T){
	M1[(i-1)*T+j,i]=1
}
}
for (i in 1:num){
	if(i%%T >0){
	M2[i,i%%T]=1}
	if(i%%T ==0){
	M2[i,T]=1
}
}


##第一种数据生成方式的估计结果

beta11<-matrix(0,l,k)##原模型设定下beta的普通估计
comb11<-matrix(0,num,k)##原模型设定下的综合效应的普通估计

beta12<-matrix(0,l,k)##新模型设定下beta的普通估计
comb12<-matrix(0,num,k)##新模型设定下的综合效应的GWR估计

##第二种数据生成方式的估计结果

beta21<-matrix(0,l,k)##原模型设定下beta的普通估计
comb21<-matrix(0,num,k)##原模型设定下的综合效应的普通估计


beta22<-matrix(0,l,k)##新模型设定下beta的GWR估计
comb22<-matrix(0,num,k)##新模型设定下的综合效应的GWR估计





##系数设定为4维的，真值为（2,3,6,13），再加上常数项2.5###
beta<-matrix(c(2,3,6,13),nrow=l,ncol=1)

####循环k=1000次
for (jishu in 1:k)
{

###四个自变量，两个正态分布，两个个均匀分布,##
xx<-matrix(c(rnorm(num, 10, 1), runif(num, 2, 4),rnorm(num, 1, 1), runif(num, -1, 1)), nrow=num, ncol=4)

##将自变量与一列1合并（加上常数项）##
###xx<-cbind(matrix(1,num,1),xx)

###误差项设定为4种分布###
ee <- matrix(c(rnorm(num, 0, 0.5)), nrow=num, ncol=1)
## ee <- matrix(c(runif(num, -sqrt(3)/2,sqrt(3)/2)), nrow=num, ncol=1)
 #ee <- matrix(c(sqrt(3)/4*rt(num,8)), nrow=num, ncol=1)
 #ee <- matrix(c(1/8*rchisq(num,8)-1), nrow=num, ncol=1)

###alpha+lamda,原模型设定下生成的数据是yy1
yy1<-xx%*%beta+zonghe1+ee



####变系数设定,新模型设定下生成的数据是yy2
yy2<-xx%*%beta+zonghe2+ee


### 一、 普通估计
xishu1<-cbind(xx,M1,M2)

###第一种生成方式的普通估计##
guji1<-ginv(t(xishu1)%*%xishu1)%*%t(xishu1)%*%yy1

###加入可识别条件，即求和为0，在这简化为减去均值
##chat1<-guji1[1,]
bhat1<-guji1[1:l,]
singhat1<-guji1[(l+1):(n+l),]
timehat1<-guji1[(n+l+1):(n+T+l),]

##chat1<-chat1+mean(singhat1)+mean(timehat1)
##singhat1<-singhat1-mean(singhat1)
###timehat1<-timehat1-mean(timehat1)

##C和beta的估计值
##chat1  和  bhat1

###综合影响的估计值
zonghehat1<-matrix(0,num,1)
for (p in 1:n){
	for (q in 1:T){
		zonghehat1[(p-1)*T+q]=singhat1[p]+timehat1[q]
	}
}
	
##第二种生成方式的普通估计
guji2<-ginv(t(xishu1)%*%xishu1)%*%t(xishu1)%*%yy2

###加入可识别条件，即求和为0，在这简化为减去均值
##chat2<-guji2[1,]
bhat2<-guji2[1:l,]
singhat2<-guji2[(l+1):(n+l),]
timehat2<-guji2[(n+l+1):(n+T+l),]
##chat2<-chat2+mean(singhat2)+mean(timehat2)
##singhat2<-singhat2-mean(singhat2)
##timehat2<-timehat2-mean(timehat2)

##C和beta的估计值
##chat2  和  bhat2

###综合影响的估计值
zonghehat2<-matrix(0,num,1)
for (p in 1:n){
	for (q in 1:T){
		zonghehat2[(p-1)*T+q]=singhat2[p]+timehat2[q]
	}
}


####  二、GWR估计方法  （两组数据的theta值是不一样的）##
weight1<-matrix(0,num,num)
###选取第一种生成方式所得数据的theta
cv<-function(theta){
	for (i in 1:num){
		for (j in 1:num){
			weight1[i,j] = 1/sqrt(2*pi)*exp(-0.5*theta^2*((varlong[i]-varlong[j])^2+(varlat[i]-varlat[j])^2))
		}
	}
	sqr1<-0
	sqr2<-0
	for (i in 1:n){
		weight<-weight1
		for(j in 1:T){
			for(jj in 1:num){
				weight[(i-1)*T+j,jj]=0
				weight[jj,(j-1)*T+j]=0
			}
		}
		L<-ginv(t(M2)%*%diag(weight[1,])%*%M2)%*%t(M2)%*%diag(weight[1,])
		for (i in 2:n){
			L<-rbind(L,ginv(t(M2)%*%diag(weight[(i-1)*T+1,])%*%M2)%*%t(M2)%*%diag(weight[(i-1)*T+1,]))
		}
		betahat1<-ginv(t(xx)%*%t(diag(num)-L)%*%(diag(num)-L)%*%xx)%*%t(xx)%*%t(diag(num)-L)%*%(diag(num)-L)%*%yy1
		betahat2<-ginv(t(xx)%*%t(diag(num)-L)%*%(diag(num)-L)%*%xx)%*%t(xx)%*%t(diag(num)-L)%*%(diag(num)-L)%*%yy2
		for (j in 1:T){
			sqr1<-sqr1+(yy1[(i-1)*T+j]-xx[(i-1)*T+j,]%*%betahat1)^2
			sqr2<-sqr2+(yy2[(i-1)*T+j]-xx[(i-1)*T+j,]%*%betahat2)^2
		}	
	}
	value1<-sqr1
	value2<-sqr2
	list(value1=value1, value2=value2,theta=theta)
}

theta = seq(0,1,by=.01)
table<-matrix(0,101,3)
for (i in 1:100){
	table[i,1]<-cv(theta[i])$value1
	table[i,2]<-cv(theta[i])$value2
	table[i,3]<-cv(theta[i])$theta
}
theta1<-0
theta2<-2.5



weight11 <- matrix(0, nrow=num, ncol=num)
for (i in 1:num){
            for(j in 1:num){
                weight11[i,j] = 1/sqrt(2*pi)*exp(-0.5*theta1^2*((varlong[i]-varlong[j])^2+(varlat[i]-varlat[j])^2))
            }
        }
L11<-ginv(t(M2)%*%diag(weight11[1,])%*%M2)%*%t(M2)%*%diag(weight11[1,])
for (i in 2:n){
	L11<-rbind(L11,ginv(t(M2)%*%diag(weight11[(i-1)*T+1,])%*%M2)%*%t(M2)%*%diag(weight11[(i-1)*T+1,]))
	
}


###第一种生成方式生成的GWR估计，即xx,yy1
betahat1<-ginv(t(xx)%*%t(diag(num)-L11)%*%(diag(num)-L11)%*%xx)%*%t(xx)%*%t(diag(num)-L11)%*%(diag(num)-L11)%*%yy1
lamdahat1<-L11%*%(yy1-xx%*%betahat1)
##cc1<-betahat1[1,]
bb1<-betahat1[1:l,]
##cc1<-cc1+mean(lamdahat1)
##lamdahat1<-lamdahat1-mean(lamdahat1)

weight22 <- matrix(0, nrow=num, ncol=num)
for (i in 1:num){
            for(j in 1:num){
                weight22[i,j] = 1/sqrt(2*pi)*exp(-0.5*theta2^2*((varlong[i]-varlong[j])^2+(varlat[i]-varlat[j])^2))
            }
        }
L22<-ginv(t(M2)%*%diag(weight22[1,])%*%M2)%*%t(M2)%*%diag(weight22[1,])
for (i in 2:n){
	L22<-rbind(L22,ginv(t(M2)%*%diag(weight22[(i-1)*T+1,])%*%M2)%*%t(M2)%*%diag(weight22[(i-1)*T+1,]))
	
}

####第二种生成方式生成的GWR估计,即xx，yy2
betahat2<-ginv(t(xx)%*%t(diag(num)-L22)%*%(diag(num)-L22)%*%xx)%*%t(xx)%*%t(diag(num)-L22)%*%(diag(num)-L22)%*%yy2
lamdahat2<-L22%*%(yy2-xx%*%betahat2)
##cc2<-betahat2[1,]
bb2<-betahat2[1:l,]
##cc2<-cc2+mean(lamdahat2)
##lamdahat2<-lamdahat2-mean(lamdahat2)




###结果整理

##第一种数据生成方式，两种估计方法的结果

##beta11[1,jishu]<-chat1
for (i in 1:l){
	beta11[i,jishu]<-bhat1[i]
}
comb11[,jishu]=zonghehat1


##beta12[1,jishu]<-cc1
for (i in 1:l){
	beta12[i,jishu]<-bb1[i]
}
comb12[,jishu]=lamdahat1


###第二种数据生成方式，两种估计的结果
##beta21[1,jishu]<-chat2
for (i in 1:l){
	beta21[i,jishu]<-bhat2[i]
}
comb21[,jishu]=zonghehat2

##beta22[1,jishu]<-cc2
for (i in 1:l){
	beta22[i,jishu]<-bb2[i]
}
comb22[,jishu]=lamdahat2
}




####结果展示

##第一种数据生成方式
####比较均值，标准差，mse
mean11<-rbind(as.matrix(rowMeans(beta11)),as.matrix(rowMeans(comb11)))
sd11<-rbind(as.matrix(apply(beta11,1,sd)),as.matrix(apply(comb11,1,sd)))
mean12<-rbind(as.matrix(rowMeans(beta12)),as.matrix(rowMeans(comb12)))
sd12<-rbind(as.matrix(apply(beta12,1,sd)),as.matrix(apply(comb12,1,sd)))


###综合效应的mse
mse11<-matrix(0,num,1)
mse12<-matrix(0,num,1)

for (i in 1:num){
	sum1=0
	sum2=0
	for (j in (1:k)){
	sum1=sum1+(comb11[i,j]-zonghe1[i])*(comb11[i,j]-zonghe1[i])
	sum2=sum2+(comb12[i,j]-zonghe1[i])*(comb12[i,j]-zonghe1[i])
	}
	mse11[i]=sum1/k
	mse12[i]=sum2/k
}

####截距项与系数的mse

sum11<-matrix(0,l,k)
sum12<-matrix(0,l,k)
for (i in 1:k){
	sum11[,i]=(beta11[,i]-beta)*(beta11[,i]-beta)
	sum12[,i]=(beta12[,i]-beta)*(beta12[,i]-beta)
}

###把mse合并起来
mse11<-rbind(as.matrix(rowSums(sum11))/k,mse11)
mse12<-rbind(as.matrix(rowSums(sum12))/k,mse12)
result1<-cbind(mean11,sd11,mse11,mean12,sd12,mse12)
true1<-rbind(beta,zonghe1)





###第二种数据生成方式两种估计的结果
mean21<-rbind(as.matrix(rowMeans(beta21)),as.matrix(rowMeans(comb21)))
sd21<-rbind(as.matrix(apply(beta21,1,sd)),as.matrix(apply(comb21,1,sd)))
mean22<-rbind(as.matrix(rowMeans(beta22)),as.matrix(rowMeans(comb22)))
sd22<-rbind(as.matrix(apply(beta22,1,sd)),as.matrix(apply(comb22,1,sd)))


###综合效应的mse
mse21<-matrix(0,num,1)
mse22<-matrix(0,num,1)

for (i in 1:num){
	sum1=0
	sum2=0
	for (j in (1:k)){
	sum1=sum1+(comb21[i,j]-zonghe2[i])*(comb21[i,j]-zonghe2[i])
	sum2=sum2+(comb22[i,j]-zonghe2[i])*(comb22[i,j]-zonghe2[i])
	}
	mse21[i]=sum1/k
	mse22[i]=sum2/k
}

####截距项与系数的mse

sum21<-matrix(0,l,k)
sum22<-matrix(0,l,k)
for (i in 1:k){
	sum21[,i]=(beta21[,i]-beta)*(beta21[,i]-beta)
	sum22[,i]=(beta22[,i]-beta)*(beta22[,i]-beta)
}

###把mse合并起来
mse21<-rbind(as.matrix(rowSums(sum21))/k,mse21)
mse22<-rbind(as.matrix(rowSums(sum22))/k,mse22)
result2<-cbind(mean21,sd21,mse21,mean22,sd22,mse22)
true2<-rbind(beta,zonghe2)

write.csv(cbind(true1,result1),"1.csv")
write.csv(cbind(true2,result2),"2.csv")