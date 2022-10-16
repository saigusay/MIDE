
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
Tau= c(Inf,10,1,0.1)
TPL = NULL #PT

for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1])
	TPL = c(TPL,paste("eta=",eta,"_tau=",tau,sep=""))

# Number of simulations and interations of estimation
n.sim = 100
iter.cal<-10; mit<-1

# Values of q of the root trimmed mean square prediction error. The trimmed proportion is 1-q/100.
qq = c(70,80,90)

# True value of the coefficient parameters of the target and contamination distribution
be.true=matrix(c(-1.5,1,1.5),ncol=1)
be.ctm = matrix(c(-Inf,0,0),ncol=1) #no-contamination
#be.ctm = matrix(c(-2.5,-1,-0.5),ncol=1) #light-contamination
#be.true=matrix(c(-2,1,1.5),ncol=1); be.ctm =matrix(c(-2,-1,-0.5),ncol=1) #heavy-contamination

p <- nrow(be.true)



#-----------------------------------------------------------------------------------------------------#
#   Simulation run																					  #
#-----------------------------------------------------------------------------------------------------#

dir = paste("result/",format(Sys.time(),"%y%m%d%H%M%S"),"_",lab,sep="")
dir.create(dir)

fn = paste(dir,"/settings.csv",sep="")
write.table(cbind(paste("be",1:p,sep=""),c(be.true)),fn,sep=",",row.name=F,col.name=F,append=F)
write.table(cbind(paste("bectm",1:length(be.ctm),sep=""),c(be.ctm)),fn,sep=",",row.name=F,col.name=F,append=T)

A=""; for(i in Eta)A=paste(A,i)
write.table(cbind("Eta",A),fn,sep=",",row.name=F,col.name=F,append=T)
A=""; for(i in Tau)A=paste(A,i)
write.table(cbind("Tau",A),fn,sep=",",row.name=F,col.name=F,append=T)
write.table(t(c("n.sim",n.sim)),fn,sep=",",row.name=F,col.name=F,append=T)

write.table(be.true,paste(dir,"/be.true.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
write.table(t(TPL),paste(dir,"/summary_it2.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)

for(i in 1:length(TPL))assign(paste("be",i,sep=""),matrix(nrow=0,ncol=p))
for(it in 1:n.sim){
set.seed(it)
	x = cbind(1,rnorm(nd,0,5),rnorm(nd,0,5))
	x[,2:3] = scale(x[,2:3])

	x.ctm = x
	lam.tar = exp(x%*%be.true)
	lam.ctm = exp(x.ctm%*%be.ctm)
	lam.gene = lam.tar + lam.ctm
	N = rpois(1,sum(lam.gene))
	y = rmultinom(1,N,lam.gene/sum(lam.gene))

	for(q in qq)assign(paste("mspel",q,sep=""),NULL)
	msetl = it2l = NULL
	for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1]){
		tpl = paste("eta=",eta,"_tau=",tau,sep="")
		out.est = xiest(eta,tau,y,x)
		be = out.est$be
		lam = exp(x%*%be)
print(paste(eta,tau)); print(out.est$it2); print(cbind(be.true,be))
		rsd = (y-lam)^2
		for(q in qq){
			mi = mean(rsd[order(rsd)[1:ceiling(q/100*(nd+1))]])
			assign(paste("mspel",q,sep=""),c(get(paste("mspel",q,sep="")),mi))
		}
		msetl = c(msetl,mean((lam.tar-lam)^2))
		it2l = c(it2l,out.est$it2)
		write.table(t(be),paste(dir,"/",tpl,"_be.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
		write.table(out.est$conv,paste(dir,"/",tpl,"_conv.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
	}#tpl
	write.table(sum(lam.ctm)/sum(lam.gene),paste(dir,"/expected contamination proportion.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
	for(q in qq)
		write.table(t(get(paste("mspel",q,sep=""))),paste(dir,"/summary_mspel",q,".csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
	write.table(t(msetl),paste(dir,"/summary_msetl.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
	write.table(t(it2l),paste(dir,"/summary_it2.csv",sep=""),sep=",",row.name=F,col.name=F,append=T)
}#it

for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1]){
	tpl = paste("eta=",eta,"_tau=",tau,sep="")
	be = read.csv(paste(dir,"/",tpl,"_be.csv",sep=""),header=F)
	assign(paste("be",which(tpl==TPL),sep=""),be)
}

for(q in qq){
	labq = paste("mipde",q,sep="")
	SM = read.csv(paste(dir,"/summary_mspel",q,".csv",sep=""),header=F)
	fs = function(v)order(v)[1]
	tplx = apply(SM,1,fs)
	bex = matrix(nrow=0,ncol=p)
	for(i in 1:n.sim)
		bex = rbind(bex,get(paste("be",tplx[i],sep=""))[i,])
	write.table(Tau[tplx],paste(dir,"/",labq,"_tau.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
	write.table(bex,paste(dir,"/",labq,"_be.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
}


# Output a boxplot of the parameter estimates

ylim = c(-5,5)
TPLp = c("eta=1_tau=Inf",paste("mipde",qq,sep=""))
par(mfrow=c(1,1+length(qq)))
for(tpl in TPLp){
	be = read.csv(paste(dir,"/",tpl,"_be.csv",sep=""),header=F)
	labt = "MLE"
	if(which(tpl==TPLp)>1)labt = paste("MIDE_q=",qq[which(tpl==TPLp)-1]/100,sep="")
	boxplot(be,ylim=ylim,xaxt="n",las=2,main=labt,cex.main=2,cex.axis=1.8)
	axis(1,at=1:ncol(be),tick=F,line=1,labels=c(TeX(paste("$\\beta_{",0:(ncol(be)-1),"}$",sep=""))),cex.axis=1.5)
	points(1:length(be.true),be.true,col="red",pch=20,cex=2,cex.axis=1.8)
}
dev.copy2eps(file=paste(dir,"/betahat.eps",sep=""),family="Times"); dev.off()


# Create a table of the percentages of the tuning parameter tau which was selected

tab = matrix(nrow=0,ncol=length(Tau))
for(q in qq){
	tt = read.csv(paste(dir,"/mipde",q,"_tau.csv",sep=""),header=F)
	A = NULL; for(tau in sort(Tau))A = c(A,sum(tau==tt))
	A = as.character(A/sum(A))
	for(i in 1:length(A))while(nchar(A[i])<5)A[i] = paste(A[i],0,sep="")
	tab = rbind(tab,A)
}
tab = cbind(qq/100,tab)
colnames(tab) = c("q",sort(Tau))
tab

