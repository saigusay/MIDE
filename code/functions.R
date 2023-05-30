
logit=function(x)log(x/(1-x))

expit=function(x)1/(1+exp(-x))

xiest <- function(eta=1,tau,y,x,z,be=matrix(0,nrow=ncol(x),ncol=1),al=matrix(0,nrow=ncol(z),ncol=1)){
	ne = nrow(y)
	be.pre=be; dimbe = length(be)
	al.pre=al; dimal = length(al)
	fp = function(t)1-(1+t)^-1
	difmax=Inf; difth = 10^-2; it2=0; le.pre = Inf; ll = dd = Inf; ka0=NaN
	while((!is.nan(difmax) & difmax>difth) & it2<iter.cal){
		it2 = it2+1; upd = F
		gb = matrix(0,dimbe,1); hbb = matrix(0,dimbe,dimbe)
		ga = matrix(0,dimal,1); haa = matrix(0,dimal,dimal); hba = matrix(0,dimbe,dimal)
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			bw = x[id,]%*%be; lam = c(exp(bw))
			az = z[id,]%*%al; dp = c(1/(1+exp(-az))); pp = c(1/(1+exp(az)))
			if(tau==Inf){
				gb = gb + max(c(1,yi)) * ((yi>0) - wi*lam*dp)%x%matrix(x[id,],ncol=1)
				ga = ga + max(c(1,yi)) * pp*((yi>0) - wi*lam*dp)%x%matrix(z[id,],ncol=1)
				hbb = hbb - - max(c(1,yi)) * (wi*lam*dp)%x%(t(x[id,])%x%x[id,])
				hba = hba - - max(c(1,yi)) * (pp*wi*lam*dp)%x%((x[id,])%x%t(z[id,]))
				haa = haa - - max(c(1,yi)) * (dp*pp*(wi*lam*pp + (yi>0)-wi*lam*dp))%x%(t(z[id,])%x%z[id,])
			}else{
				gb = gb + max(c(1,yi)) * (fp(tau*lam*dp))*((yi>0) - wi*lam*dp)%x%matrix(x[id,],ncol=1)
				ga = ga + max(c(1,yi)) * (fp(tau*lam*dp))*pp*((yi>0) - wi*lam*dp)%x%matrix(z[id,],ncol=1)
				hbb = hbb - max(c(1,yi)) * (fp(tau*lam*dp)*(((yi>0)-wi*lam*dp)/(1+tau*lam*dp) - wi*lam*dp))%x%(t(x[id,])%x%x[id,])
				hba = hba - max(c(1,yi)) * (fp(tau*lam*dp)*pp*(((yi>0)-wi*lam*dp)/(1+tau*lam*dp)-wi*lam*dp))%x%((x[id,])%x%t(z[id,]))
				haa = haa - max(c(1,yi)) * (fp(tau*lam*dp)*pp*(-dp*((yi>0)-wi*lam*dp) + pp/(1+tau*lam*dp)*((yi>0)-wi*lam*dp) - pp*wi*lam*dp))%x%(t(z[id,])%x%z[id,])
			}
		}
		G = matrix(c(gb,ga),ncol=1)
		H = as.matrix(rbind(cbind(hbb,hba),cbind(t(hba),haa)))
		Hinv = try(ginv(H),silent=T)
		if(class(Hinv)[1]!="try-error"){
			HG = Hinv%*%G
			be = be + HG[1:dimbe,] * mit
			al = al + HG[(dimbe+1):(dimbe+dimal),] * mit
			upd = T
		}else{
			upd = F
		}
		le = 0
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			bw = x[id,]%*%be; lam = c(exp(bw))
			az = z[id,]%*%al; dp = c(1/(1+exp(-az))); pp = c(1/(1+exp(az)))
			if(tau==Inf){
				le = le - max(c(1,yi))*((yi>0)*log(lam*dp)-wi*lam*dp)
			}else{
				le = le - max(c(1,yi))*(((yi>0)+wi/tau)*log(1+tau*lam*dp)-wi*lam*dp)
			}
		}
		ll = c(ll,le)
		difmax = abs(ll[length(ll)-1]-le)
		#difmax = max(abs(c(be-be.pre,al-al.pre)))
		dd = c(dd,difmax)
		be.pre=be; al.pre=al
	} #while
	conv = upd & difmax<difth
	if(is.na(conv))conv=F

	return(list(be=be,al=al,conv=conv,difmax=difmax,it2=it2))
}

xitpest <- function(eta=1,tau,y,x,z,be=matrix(0,nrow=ncol(x),ncol=1),phi,al=matrix(0,nrow=ncol(z),ncol=1)){
	be[abs(be)<10^-10] = 0
	np = sum(y); ne = nrow(y)
	be.pre = be; dimbe = length(be)
	al.pre = al; dimal = length(al)
	fp = function(t)1-(1+t)^-1
	lxp = function(rho){
		l = 0; beopt = be + rho*gp
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			bx = x[id,]%*%beopt; lam = exp(bx)
			az = z[id,]%*%al; dp = c(1/(1+exp(-az)))
			if(tau==Inf){
				l = l - max(c(1,yi)) * ((yi>0)*log(lam*dp) - wi*lam*dp)
			}else{
				l = l - max(c(1,yi)) * (((yi>0)+wi/tau)*log(1+tau*lam*dp) - wi*lam*dp)
			}
		}
		l = l + phi*sum(abs(beopt[-1]))
		return(l)
	}
	difmax=-Inf; difth = 10^-2; it2=0; le.pre = Inf; dd = NULL; llp = Inf; ka0=NaN
	while((!is.nan(difmax) & -difmax>difth) & it2<iter.cal){
		it2 = it2+1; upd = T

		ga = matrix(0,dimal,1); haa = matrix(0,dimal,dimal)
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			bx = x[id,]%*%be; lam = c(exp(bx))
			az = z[id,]%*%al; dp = c(1/(1+exp(-az))); pp = c(1/(1+exp(az)))
			if(tau==Inf){
				ga = ga + max(c(1,yi)) * pp*((yi>0) - wi*lam*dp)%x%matrix(z[id,],ncol=1)
				haa = haa - - max(c(1,yi)) * (dp*pp*(wi*lam*pp + (yi>0)-wi*lam*dp))%x%(t(z[id,])%x%z[id,])
			}else{
				ga = ga + max(c(1,yi)) * fp(tau*lam*dp)*pp*((yi>0) - wi*lam*dp)%x%matrix(z[id,],ncol=1)
				haa = haa - max(c(1,yi)) * (fp(tau*lam*dp)*pp*(-dp*((yi>0)-wi*lam*dp) + pp/(1+tau*lam*dp)*((yi>0)-wi*lam*dp) - pp*wi*lam*dp))%x%(t(z[id,])%x%z[id,])
			}
		}
		A = try(ginv(haa),silent=T)
		if(class(A)[1]!="try-error"){
			al <- al + A%*%ga
		}
		gb = matrix(0,dimbe,1)
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			bx = x[id,]%*%be; lam = c(exp(bx))
			az = z[id,]%*%al; dp = c(1/(1+exp(-az)))
			if(tau==Inf){
				gb = gb + max(c(1,yi)) * ((yi>0) - wi*lam*dp)%x%matrix(x[id,],ncol=1)
			}else{
				gb = gb + max(c(1,yi)) * fp(tau*lam*dp)*((yi>0) - wi*lam*dp)%x%matrix(x[id,],ncol=1)
			}
		}
		if(sum(is.nan(gb))>0){ gb<-rep(0,length(be)); upd = F}
		gp <- NULL
		for(i in 1:dimbe){
			if(i==1){
				gp <- c(gp,gb[i])
			}else if(be[i]!=0){
				gp <- c(gp, gb[i]-phi*sign(be[i]))
			}else if(abs(gb[i])>phi){
				gp <- c(gp, gb[i]-phi*sign(gb[i]))
			}else{
				gp <- c(gp, 0)
			}
		}
		A = sign(be)==-sign(gp) & be!=0; A[1] = F
		rh = -be/gp; rh[abs(rh) < 10^-10] = 0
		rhedg = ifelse(sum(A)>0,min(rh[A]),10^4)
		A = try(optim(0,lxp,lower=0,upper=rhedg,method="L"),silent=T)
		rhopt = ifelse(class(A)[1]!="try-error",A$par,0)
		llp = c(llp,lxp(rhopt))
		difmax = lxp(rhopt) - lxp(0)
		be <- be + rhopt*gp
		be[abs(be)<10^-10] <- 0
		dd = c(dd,difmax)
		be.pre=be
	} #while
	conv = upd & -difmax<difth
	if(is.na(conv))conv=F
if(!conv)print(dd)

	weight = NULL
	for(id in 1:ne){
		bw = x[id,]%*%be; lam = c(exp(bw))
		az = z[id,]%*%al; dp = c(1/(1+exp(-az)))
		if(tau==Inf){
			weight = c(weight,1)
		}else{
			weight = c(weight,fp(tau*lam*dp))
		}
	}
	return(list(be=be,al=al,weight=weight,conv=conv,difmax=difmax,it2=it2))
}
