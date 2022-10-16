
xiest <- function(eta=1,tau,y,x,be=matrix(0,nrow=ncol(x),ncol=1)){
	np = sum(y); ne = nrow(y)
	be.pre=be; p = length(be)
	difmax=Inf; difth = 10^-2; it2=0; le.pre = Inf; ll = dd = NULL; ka0=NaN
	while((!is.nan(difmax) & difmax>difth) & it2<iter.cal){
		it2 = it2+1; upd = F
		g1 = g1n = g2n = h2n = matrix(0,p,1); g1d = g2d = 0
		h1 = h1n = h2n = matrix(0,p,p)
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			lp = x[id,]%*%be; lam = exp(lp)
			for(i in 1:max(c(1,yi))){
				if(tau==Inf){
					g1 = g1 + ((yi>0) - wi*lam)%x%matrix(x[id,],ncol=1)
					h1 = h1 - (-wi*lam)%x%(t(x[id,])%x%x[id,])
				}else{
					g1 = g1 + c(1-(1+eta*tau*lam)^(-1/eta))*((yi>0) - wi*lam)%x%matrix(x[id,],ncol=1)
					h1 = h1 - (-eta*tau*lam*(1+eta*tau*lam)^(-(1+eta)/eta)*((yi>0) - wi*lam) + (1-(1+eta*tau*lam)^(-1/eta))*(-wi*lam))%x%(t(x[id,])%x%x[id,])
				}
			}
		}
		A=try(ginv(h1),silent=T)
		if(class(A)[1]!="try-error"){
			be = be + A%*%g1 * mit; upd = T
		}else{
			upd = F
		}
		le = NA
		ll = c(ll,le)
		difmax = max(abs(c(be-be.pre)))
		dd = c(dd,difmax)
		be.pre=be
	} #while
	conv = upd & difmax<difth
	if(is.na(conv))conv=F

	#return(list(be=be,b=b,th=th,lgg=lgg,conv=conv,difmax=difmax,cll=cll,Lgg=Lgg,rhh=rhh))
	return(list(be=be,conv=conv,difmax=difmax,it2=it2))
}

xipest <- function(eta=1,tau,y,x,be0=matrix(0,nrow=ncol(x),ncol=1),phi){
	be0[abs(be0)<10^-10] = 0
	np = sum(y); ne = nrow(y)
	be.pre = be = be0; p = length(be0)
	lxp = function(rho){
		l = 0; beopt = be+rho*gp
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			lp = x[id,]%*%beopt; lam = exp(lp)
			for(i in 1:max(c(1,yi))){
				if(tau==Inf){
					l = l - (yi>0)*log(lam) + wi*lam
				}else{
					l = l - ((yi>0)+wi/tau)*log(1+tau*lam) + wi*lam
				}
			}
		}
		l = l + phi*sum(abs(beopt[-1]))
		return(l)
	}
	difmax=-Inf; difth = 10^-2; it2=0; le.pre = Inf; dd = NULL; llp = Inf; ka0=NaN
	while((!is.nan(difmax) & -difmax>difth) & it2<iter.cal){
		it2 = it2+1
		g1 = g1n = g2n = h2n = matrix(0,p,1); g1d = g2d = 0
		h1 = h1n = h2n = matrix(0,p,p)
		for(id in 1:ne){
			yi = y[id]; wi = 1/max(c(1,yi))
			lp = x[id,]%*%be; lam = exp(lp)
			for(i in 1:max(c(1,yi))){
				if(tau==Inf){
					g1 = g1 + ((yi>0) - wi*lam)%x%matrix(x[id,],ncol=1)
					#h1 = h1 - (-wi*lam)%x%(t(x[id,])%x%x[id,])
				}else{
					g1 = g1 + c(1-(1+eta*tau*lam)^(-1/eta))*((yi>0) - wi*lam)%x%matrix(x[id,],ncol=1)
					#h1 = h1 - (-eta*tau*lam*(1+eta*tau*lam)^(-(1+eta)/eta)*((yi>0) - wi*lam) + (1-(1+eta*tau*lam)^(-1/eta))*(-wi*lam))%x%(t(x[id,])%x%x[id,])
				}
			}
		}
		gp <- NULL
		for(i in 1:p){
			if(i==1){
				gp <- c(gp,g1[i])
			}else if(be[i]!=0){
				gp <- c(gp, g1[i]-phi*sign(be[i]))
			}else if(abs(g1[i])>phi){
				gp <- c(gp, g1[i]-phi*sign(g1[i]))
			}else{
				gp <- c(gp, 0)
			}
		}
		A = sign(be)==-sign(gp) & be!=0; A[1] = F
		rh = be/gp; rh[abs(rh) < 10^-10] = 0
		rhedg = ifelse(sum(A)>0,min(-rh[A]),10^4)
		A = try(optim(0,lxp,lower=0,upper=rhedg,method="L"))
		rhopt = ifelse(class(A)[1]!="try-error",A$par,0)

		llp = c(llp,lxp(rhopt))
		difmax = lxp(rhopt) - lxp(0)
		be <- be + rhopt*gp * mit
		be[abs(be)<10^-10] = 0
		dd = c(dd,difmax)
		be.pre=be
	} #while
	conv = -difmax<difth
	if(is.na(conv))conv=F
if(!conv)print(dd[(length(dd)-10):length(dd)])
	return(list(be=be,conv=conv,difmax=difmax,it2=it2))
}
