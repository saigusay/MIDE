
#-----------------------------------------------------------------------------------------------------#
#   Simulation studies for contaminated PPP                                        		              #
#-----------------------------------------------------------------------------------------------------#

### Note:
  # The following represents a sample simulation script, which returns .csv files summarizing
  # the robust estimates of coefficient parameters under the Poisson point process model.
  # The cumulative distribution for the robust estimation is the Pareto distribution with the shape parametr of 1

library(MASS)
source("code/functions.R")



#-----------------------------------------------------------------------------------------------------#
#   (Sample) Simulation parameter values															  #
#-----------------------------------------------------------------------------------------------------#

lab = "sim"

# Candidate values of tuning parameters for the grid search
Eta = 1
Tau= c(Inf,20,10,5,1,0.1)
TPL = NULL #PT
for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1])
	TPL = c(TPL,paste("eta=",eta,"_tau=",tau,sep=""))

# Number of simulations and interations of estimation
n.sim = 10
iter.cal<-20; mit<-1

# Value of q of the root trimmed means square prediction error. The trimmed proportion is 1-q/100.
qq=c(80,90,99)

# Number of divisions for the penalty parameter of the L1 penalty term
nphi = 10

# True value of the coefficient parameters for the target and contamination distribution
be.true = matrix(c(-2,1,1,-1,-1),ncol=1)
be.ctm = matrix(c(-Inf,rep(0,length(be.true)-1)),ncol=1) #no-contamination
#be.ctm = -be.true; be.ctm[1] = -4.2 #light-contamination
#be.ctm = -be.true; be.ctm[1] = -3.4 #heavy-contamination
p <- nrow(be.true)

# True value of the coefficient parameters for the detection probability
p0.true = 0.5
al1.true = c(1,-1) #runif(1,0,2)
al.true = matrix(c(logit(p0.true),al1.true),ncol=1); int.al = T
al.true=al.true[-1]; int.al = F

# Number of divisions for the study area
nr <- 10; nc <- 20; nd <- nr*nc


#-----------------------------------------------------------------------------------------------------#
#   Simulation run																					  #
#-----------------------------------------------------------------------------------------------------#

issave=T
if(issave){ 
dir = paste("result/",format(Sys.time(),"%y%m%d%H%M%S"),"_",lab,sep="")
dir.create(dir)

fn = paste(dir,"/settings.csv",sep="")
write.table(cbind(paste("be",1:length(be.true),sep=""),c(be.true)),fn,sep=",",row.name=F,col.name=F,append=F)
write.table(cbind(paste("bectm",1:length(be.ctm),sep=""),c(be.ctm)),fn,sep=",",row.name=F,col.name=F,append=T)

write.table(cbind(paste("al",1:length(al.true),sep=""),c(al.true)),fn,sep=",",row.name=F,col.name=F,append=T)

A=""; for(i in Eta)A=paste(A,i)
write.table(cbind("Eta",A),fn,sep=",",row.name=F,col.name=F,append=T)
A=""; for(i in Tau)A=paste(A,i)
write.table(cbind("Tau",A),fn,sep=",",row.name=F,col.name=F,append=T)
write.table(t(c("n.cell",nd)),fn,sep=",",row.name=F,col.name=F,append=T)
write.table(t(c("n.sim",n.sim)),fn,sep=",",row.name=F,col.name=F,append=T)

write.table(be.true,paste(dir,"/be.true.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
write.table(t(TPL),paste(dir,"/summary_it2.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
}

for(i in 1:length(TPL))assign(paste("be",i,sep=""),matrix(nrow=0,ncol=length(be.true)))
for(it in 1:n.sim){
print(paste("   it=",it))
set.seed(it)
	x = rep(1,nd); for(i in 2:length(be.true))x=cbind(x,scale(rnorm(nd)))
	z = rep(1,nd); for(i in 2:length(al.true))z=cbind(z,scale(rnorm(nd)))
	if(!int.al){ z = scale(rnorm(nd)); if(length(al.true)>1)for(i in 2:length(al.true))z=cbind(z,scale(rnorm(nd))) }
	dimbe = ncol(x)

	lam.tar = exp(x%*%be.true)
	lam.ctm = exp(x%*%be.ctm)
	lam.gene = lam.tar + lam.ctm
	dp = expit(z%*%al.true)
	N = rpois(1,sum(lam.gene*dp))
	y = rmultinom(1,N,lam.gene*dp/sum(lam.gene*dp))
	
	n0 = floor(sum(y==0)/2); n1 = floor(sum(y>0)/2)
	for(i in 1:10){
		id1r = sample(which(y>0)); id0r = sample(which(y==0))
		assign(paste("si",i,sep=""),list(c(id0r[1:n0],id1r[1:n1]),c(setdiff(id0r,id0r[1:n0]),setdiff(id1r,id1r[1:n1]))))
	}
	idl.cv = list(si1,si2,si3,si4,si5,si6,si7,si8,si9,si10)

	for(q in qq)assign(paste("mspel",q,sep=""),NULL)
	it2l = NULL
	for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1]){
#eta=1;tau=0.1
		tpl = paste("eta=",eta,"_tau=",tau,sep="")
print(" ");print(tpl)
		out.xi = xiest(eta,tau,y,x,z);be0 = out.xi$be; al0 = out.xi$al
print("initial value");print(cbind(be.true,be0));print(cbind(al.true,al0))
if(issave){
			write.table(t(be0),paste(dir,"/",tpl,"_be0.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
			write.table(t(al0),paste(dir,"/",tpl,"_al0.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
}
		g0 = matrix(0,dimbe,1)
		for(id in 1:nd){
			yi = y[id]; wi = 1/max(c(1,yi))
			bx = x[id,]%*%be0; lam = c(exp(bx))
			az = z[id,]%*%al0; dp = c(1/(1+exp(-az)))
			if(tau==Inf){
				g0 = g0 + max(c(1,yi))*((yi>0) - wi*lam*dp)%x%matrix(x[id,],ncol=1)
				#g0 = g0 + max(c(1,yi))*((yi>0) - wi*lam)%x%matrix(x[id,],ncol=1)
			}else{
				g0 = g0 + max(c(1,yi))*c(1-(1+eta*tau*lam*dp)^(-1/eta))*((yi>0) - wi*lam*dp)%x%matrix(x[id,],ncol=1)
				#g0 = g0 + max(c(1,yi))*c(1-(1+eta*tau*lam)^(-1/eta))*((yi>0) - wi*lam)%x%matrix(x[id,],ncol=1)
			}
		}
		phi1 = max(abs(g0[-1]))
		if(is.nan(phi1)|phi1<=1)phi1=10^4; if(phi1>10^10)phi1=10^10
		pp = sort(c(exp(seq(log(phi1),0,length.out=nphi)),0))
		bep = be0 #=be.true
		for(q in qq)assign(paste("mm",q,sep=""),NULL)
		bemat = matrix(nrow=dimbe,ncol=0); almat = matrix(nrow=ncol(z),ncol=0)
		convl = NULL; wwmat = matrix(nrow=nd,ncol=0)
		for(ip in pp){
			out.xip = xitpest(eta,tau,y,x,z,bep,phi=ip)
			bep = out.xip$be; alp = out.xip$al; wwmat=cbind(wwmat,out.xip$weight)
			if(sum(is.nan(bep))>0)bep=be0
			for(q in qq)assign(paste("mmcv",q,sep=""),NULL)
			for(ib in 1:10){
				y1 = matrix(y[idl.cv[[ib]][[1]]],ncol=1); y2 = matrix(y[idl.cv[[ib]][[2]]],ncol=1)
				x1 = x[idl.cv[[ib]][[1]],]; x2 = x[idl.cv[[ib]][[2]],]
				z1 = z[idl.cv[[ib]][[1]],]; z2 = z[idl.cv[[ib]][[2]],]
				out.ib = xitpest(eta,tau,y1,x1,z1,bep,phi=ip)
				lam = exp(x2%*%out.ib$be); lam[lam>10^10] = 10^10
				dp = c(1/(1+exp(-z2%*%out.ib$al)))
				rsd = (y2-lam*dp)^2
				out.ib = xitpest(eta,tau,y2,x2,z2,bep,phi=ip)
				lam = exp(x1%*%out.ib$be); lam[lam>10^10] = 10^10
				dp = c(1/(1+exp(-z1%*%out.ib$al)))
				rsd = c(rsd,(y1-lam*dp)^2)
				for(q in qq){
					mi = mean(rsd[order(rsd)[1:ceiling(q/100*nd)]])
					assign(paste("mmcv",q,sep=""),c(get(paste("mmcv",q,sep="")),mi))
				}
			}
			for(q in qq)
				assign(paste("mm",q,sep=""),c(get(paste("mm",q,sep="")),mean(get(paste("mmcv",q,sep="")))))
			bemat = cbind(bemat,bep); almat=cbind(almat,alp)
			convl = c(convl,out.xip$conv)
		}
		for(q in qq){
			mm = get(paste("mm",q,sep=""))
			assign(paste("mspel",q,sep=""),c(get(paste("mspel",q,sep="")),min(mm)))
			is = which(mm==min(mm[convl]))[1]
			be.best = matrix(bemat[,is],ncol=1)
			al.best = matrix(almat[,is],ncol=1)
print(paste("q=",q));print(cbind(be.true,be.best));print(cbind(al.true,al.best));print(paste("RTMSPE",min(mm)))
if(issave){
			write.table(t(be.best),paste(dir,"/",tpl,"_q=",q,"_be.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
			write.table(t(al.best),paste(dir,"/",tpl,"_q=",q,"_al.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
			write.table(t(matrix(wwmat[,is],ncol=1)),paste(dir,"/",tpl,"_q=",q,"_weight.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
			write.table(sum(convl)>0,paste(dir,"/",tpl,"_q=",q,"_conv.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
}
		}
	}#tpl
if(issave){
	write.table(sum(lam.ctm)/sum(lam.gene),paste(dir,"/expected contamination proportion.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
	for(q in qq)
		write.table(t(get(paste("mspel",q,sep=""))),paste(dir,"/summary_mspel",q,".csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
}
}#it

for(q in qq){
	for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1]){
		tpl = paste("eta=",eta,"_tau=",tau,sep="")
		be = read.csv(paste(dir,"/",tpl,"_q=",q,"_be.csv",sep=""),header=F)
		assign(paste("be",which(tpl==TPL),sep=""),be)
	}
	labq = paste("mipde_q=",q,sep="")
	SM = read.csv(paste(dir,"/summary_mspel",q,".csv",sep=""),header=F)
	fs = function(v)order(v)[1]
	tplx = apply(SM,1,fs)
	bex = matrix(nrow=0,ncol=length(be.true))
n.sim=nrow(SM)
	for(i in 1:n.sim)
		bex = rbind(bex,get(paste("be",tplx[i],sep=""))[i,])
	write.table(Tau[tplx],paste(dir,"/",labq,"_tau.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
	write.table(bex,paste(dir,"/",labq,"_be.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
}


# Output a boxplot of the parameter estimates
ylim = c(-3,3)
f=function(ip){
	if(ip==1)return(expression(MLE))
	if(ip==2)return(expression(MIDE))
}

qqp = 90
TPLp = c(paste("eta=1_tau=Inf_q=",qqp,sep=""),paste("mipde_q=",qqp,sep=""))
par(mfrow=c(1,length(qqp)*2))
for(tpl in TPLp){
	be = read.csv(paste(dir,"/",tpl,"_be.csv",sep=""),header=F)
	boxplot(be[-1],ylim=ylim,xaxt="n",las=2,main=f(which(tpl==TPLp)),cex.main=2.5,cex.axis=2.5)
	axis(1,at=1:length(be[-1]),tick=F,line=1,labels=paste("beta",1:(ncol(be)-1)),cex.axis=2.5)
	points(1:length(be.true[-1]),be.true[-1],col="red",pch=20,cex=2,cex.axis=1.8)
}
dev.copy2eps(file=paste(dir,"/boxplot_beta.eps",sep=""),family="Times"); dev.off()


# Create a table of the percentages of the tuning parameter tau which was selected
tab = matrix(nrow=0,ncol=length(Tau))
for(q in qq){
	tt = read.csv(paste(dir,"/mipde_q=",q,"_tau.csv",sep=""),header=F)
	A = NULL; for(tau in sort(Tau))A = c(A,sum(tau==tt))
	A = as.character(round(A/sum(A),3))
	for(i in 1:length(A))while(nchar(A[i])<5)A[i] = paste(A[i],0,sep="")
	tab = rbind(tab,A)
}
tab = cbind(qq/100,tab)
colnames(tab) = c("q",sort(Tau))

tab

