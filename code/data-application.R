
#-----------------------------------------------------------------------------------------------------#
#   Simulation studies for contaminated PPP                                        		              #
#-----------------------------------------------------------------------------------------------------#

library(MASS)
library(qPPP)
library(sf)
library(spatstat)
library(pROC)
library(ggplot2); library(reshape2); library(ggpubr); library(gridExtra)
source("code/functions.R")


#-----------------------------------------------------------------------------------------------------#
#   Pre-process of vascular plants data in qPPP package                           		              #
#-----------------------------------------------------------------------------------------------------#

lx = 0.125  
ly = 0.25/3 

xxg.po = seq(min(POenv[,"lon"]),max(POenv[,"lon"]),lx)
yyg.po = seq(min(POenv[,"lat"]),max(POenv[,"lat"]),ly)
ixenv.po = iyenv.po = NULL
for(i in 1:nrow(POenv)){
	dx = abs(xxg.po-POenv[i,"lon"]); ixenv.po = c(ixenv.po,order(dx)[1])
	dy = abs(yyg.po-POenv[i,"lat"]); iyenv.po = c(iyenv.po,order(dy)[1])
if(min(dx)>lx | min(dy)>ly)print(paste(i,min(dx),min(dy)))
}
res.pa = 4
lx = 0.125/res.pa  
ly = 0.25/3/res.pa 
xxg.pa = seq(min(PAenv[,"lon"]),max(PAenv[,"lon"]),lx)
yyg.pa = seq(min(PAenv[,"lat"]),max(PAenv[,"lat"]),ly)
ixenv.pa = iyenv.pa = ixenvpo.pa = iyenvpo.pa = NULL
for(i in 1:nrow(PAenv)){
	dx = abs(xxg.pa-PAenv[i,"lon"]); ixenv.pa <- c(ixenv.pa,order(dx)[1])
	dy = abs(yyg.pa-PAenv[i,"lat"]); iyenv.pa <- c(iyenv.pa,order(dy)[1])
if(min(dx)>lx | min(dy)>ly)print(paste(i,min(dx),min(dy)))
	dx = abs(xxg.po-PAenv[i,"lon"]); ixenvpo.pa <- c(ixenvpo.pa,order(dx)[1])
	dy = abs(yyg.po-PAenv[i,"lat"]); iyenvpo.pa <- c(iyenvpo.pa,order(dy)[1])
}

jpn = st_read("N03-20210101_GML/N03-21_210101.shp")
jpn_pref = jpn$N03_001
A = st_centroid(jpn$geometry)
jpn_cent = t(matrix(unlist(A), nrow=2))
AA = unique(jpn_pref)

area.po = matrix(nrow=0,ncol=10)
for(i in 1:nrow(POenv)){
	dd = abs(POenv$lon[i] - jpn_cent[,1]) + abs(POenv$lat[i] - jpn_cent[,2])
	pi = jpn_pref[order(dd)[1]]
	ai = c(is.element(pi,AA[1]),is.element(pi,AA[2:7]),is.element(pi,AA[8:14]),is.element(pi,AA[15:23]),is.element(pi,AA[24:30]),is.element(pi,AA[31:39]),is.element(pi,AA[40:47]))
	ai2 = c(sum(ai[1:4])>0,sum(ai[5:7])>0)
	area.po = rbind(area.po,c(T,ai,ai2))
}
area.pa = matrix(nrow=0,ncol=10)
for(i in 1:nrow(PAenv)){
	dd = abs(PAenv$lon[i] - jpn_cent[,1]) + abs(PAenv$lat[i] - jpn_cent[,2])
	pi = jpn_pref[order(dd)[1]]
	ai = c(is.element(pi,AA[1]),is.element(pi,AA[2:7]),is.element(pi,AA[8:14]),is.element(pi,AA[15:23]),is.element(pi,AA[24:30]),is.element(pi,AA[31:39]),is.element(pi,AA[40:47]))
	ai2 = c(sum(ai[1:4])>0,sum(ai[5:7])>0)
	area.pa = rbind(area.pa,c(T,ai,ai2))
}
areas = c("area.all","area.hkd","area.thk","area.knt","area.cb","area.knk","area.cs","area.ks","area.east","area.west")
colnames(area.po) = colnames(area.pa) = areas


#-----------------------------------------------------------------------------------------------------#
#   Parameter settings																				  #
#-----------------------------------------------------------------------------------------------------#

# Candidate values of tuning parameters for the grid search
Eta = 1
Tau= c(Inf,10,1,0.1)
TPL = NULL
for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1])TPL = c(TPL,paste("eta=",eta,"_tau=",tau,sep=""))

# Number of interations of estimation
iter.cal <- 50; mit <- 1

# Values of q of the root trimmed mean square prediction error. The trimmed proportion is 1-q/100.
qq = c(70,80,90,95,99)#quantiles for RTMSPE

# Number of divisions of candidate values for the tuning parameter phi
nphi = 100


#-------------------------------------------------------------------------------------------------------#
#   Data analysis run																					#
#-------------------------------------------------------------------------------------------------------#

vvori = colnames(PAenv)[4:ncol(PAenv)]

fnall = "result/result_all_mit=1.csv"
write.table(t(c("Area","q","tau","phi","RTMSPE_PO_MLE","AUC_PA_MLE","RTMSPE_PO_MIDE","AUC_PA_MIDE","Dif-AUC_PA")),fnall,sep=",",row.name=F,col.name=F,append=T)

for(iy in 1:length(presence_list40)){

xy.pre = PO40[names(presence_list40)[iy]==PO40[,"st.species"],c("lon","lat")]
ypo = matrix(0,nrow=nrow(POenv),ncol=1)
for(i in 1:nrow(xy.pre)){
	id = which(xy.pre[i,"lon"]==POenv[,"lon"] & xy.pre[i,"lat"]==POenv[,"lat"])
	ypo[id] = ypo[id] + 1
if(length(id)!=1)print(paste("unmatch",i))
}
rry = NULL; for(i in 1:length(ypo))rry = c(rry,rep(i,ypo[i]))

lab = paste(names(presence_list40)[iy],"_iy=",iy,sep="")
dir = paste("result/",format(Sys.time(),"%y%m%d%H%M%S"),"_",lab,sep="")
dir.create(dir)

write.table(vvori,paste(dir,"/vvori.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)

ypo.all = ypo
ypa.all = matrix(presence_list40[[iy]],ncol=1)

ss.po = NULL
for(i in 1:(ncol(area.po)-2)){
	if(sum(ypo[area.po[,i]])>= 50 & sum(ypa.all[area.pa[,i]])>0)
		ss.po = c(ss.po,i)
}

for(is in ss.po){
	idsub.po = area.po[,is]
	idsub.pa = area.pa[,is]
	ari = colnames(area.po)[is]
	fnalls = paste(dir,"/",ari,"_RTMSPE_AUC.csv",sep="")
	write.table(t(c("Area","tau","phi","q","RTMSPE_PO","AUC_PA")),fnalls,sep=",",row.name=F,col.name=F,append=F)

	ypo = matrix(ypo.all[idsub.po,],ncol=1)
	xpo = as.matrix(cbind(1,POenv[idsub.po,vvori]))
	xpa = as.matrix(cbind(1,PAenv[idsub.pa,vvori]))
	xpo[,vvori] = scale(xpo[,vvori])
	xpa[,vvori] = scale(xpa[,vvori])
	ne = nrow(xpo); p = ncol(xpo)
	ypa = matrix(presence_list40[[iy]][idsub.pa],ncol=1)

	if(is==1){
		out.xi = xiest(1,Inf,ypo,xpo); be0=out.xi$be
	}
	for(q in qq){
		assign(paste("pp.best",q,sep=""),NULL)
		assign(paste("mspel.po",q,sep=""),NULL)
		assign(paste("aucl.pa",q,sep=""),NULL)
	}
	for(eta in Eta)for(tau in Tau)if(tau!=Inf|eta==Eta[1]){
		tpl = paste("eta=",eta,"_tau=",tau,sep="")
		fnsp = paste(dir,"/",ari,"_",tpl,"_solution-path.csv",sep="")
		spmat = matrix(nrow=2+1+length(vvori),ncol=0)
		mspe.best = Inf
		for(q in qq)
			assign(paste("mspel",q,sep=""),NULL)
		bep = be0
		g0 = matrix(0,p,1)
		for(id in 1:ne){
			yi = ypo[id]; wi = 1/max(c(1,yi))
			lp = xpo[id,]%*%bep; lam = exp(lp)
			for(i in 1:max(c(1,yi))){
				if(tau==Inf){
					g0 = g0 + ((yi>0) - wi*lam)%x%matrix(xpo[id,],ncol=1)
				}else{
					g0 = g0 + c(1-(1+eta*tau*lam)^(-1/eta))*((yi>0) - wi*lam)%x%matrix(xpo[id,],ncol=1)
				}
			}
		}
		phi1 = max(abs(g0[-1])); pp = seq(phi1,0,length.out=nphi)
		for(ip in rev(pp)){
			out.xip = xipest(eta,tau,ypo,xpo,bep,phi=ip)
			bep = out.xip$be
			lam.po = exp(xpo%*%bep); lam.po[lam.po>10^10] = 10^10
			rsd.po = (ypo-lam.po)^2
			for(q in qq){
				mi = mean(rsd.po[order(rsd.po)[1:ceiling(q/100*length(rsd.po))]])
				assign(paste("mspel",q,sep=""),c(get(paste("mspel",q,sep="")),mi))
			}
			spmat = cbind(spmat,c(ip,out.xip$conv,out.xip$be))
		}
		for(q in qq){
			A = get(paste("mspel",q,sep=""))
			ip = which(A==min(A))[1]
			assign(paste("pp.best",q,sep=""),c(get(paste("pp.best",q,sep="")),pp[ip]))
			assign(paste("mspel.po",q,sep=""),c(get(paste("mspel.po",q,sep="")),min(A)))
			be.best = matrix(spmat[-c(1:2),ip],ncol=1)
			lam.pa = exp(xpa%*%be.best); lam.pa[lam.pa>10^10] = 10^10
			#rsd.pa = (ypa-lam.pa)^2
			auc = pROC::roc(c(ypa>0),c(lam.pa),direction="<")$auc
			assign(paste("aucl.pa",q,sep=""),c(get(paste("aucl.pa",q,sep="")),auc))
			write.table(cbind(c("Intercept",vvori),be.best),paste(dir,"/",ari,"_",tpl,"_q=",q,"_be.csv",sep=""),sep=",",row.name=F,col.name=F,append=F)
		}#q
		write.table(cbind(c("phi","Convergence","Intercept",vvori),spmat),fnsp,sep=",",row.name=F,col.name=F,append=F)
	}#tpl
	laba = paste(names(presence_list40)[iy],"_iy=",iy,"_",ari,sep="")
	for(q in qq){
		res = cbind(Tau,get(paste("pp.best",q,sep="")),rep(q,length(Tau)),get(paste("mspel.po",q,sep="")),get(paste("aucl.pa",q,sep="")))
		write.table(cbind(rep(ari,length(Tau)),res),fnalls,sep=",",row.name=F,col.name=F,append=T)
		write.table(t(c(laba,q,res[which(res[,4]==min(res[,4]))[1],1:2],res[which(TPL=="eta=1_tau=Inf"),4:5],res[which(res[,4]==min(res[,4]))[1],4:5],res[which(res[,4]==min(res[,4]))[1],5]-res[which(TPL=="eta=1_tau=Inf"),5])),fnall,sep=",",row.name=F,col.name=F,append=T)
	}
}#is

}#iy


#-------------------------------------------------------------------------------------------------------#
#   Output result summaries and figures in a subfolder													#
#-------------------------------------------------------------------------------------------------------#

# condition of the result to outout

# Number of plant to output
iy = 29

# Study area to output
ari = "area.ks"

# Values of tuning parameters to output
tpl0 = "eta=1_tau=Inf"
tau = 1
tplx = paste("eta=1_tau=",tau,sep="")
qstr = "q=90"
is = which(ari==colnames(area.po))
xxlim=NULL; yylim=NULL
if(ari=="area.knt"){ xxlim=NULL;yylim=c(100,160) }

A = list.files(dir)
dirs = paste(dir,"/",A[grep(names(presence_list40)[iy],A)][2],sep="") # subfolder to output


hplot = function(val,idyy,idxx,vlim=c(-Inf,Inf),xxlim=NULL,yylim=NULL,fn="test.png",lab=NULL,pal=NULL,islg=T){
	if(!is.null(xxlim)|!is.null(yylim)){
		#rxyi = ((xxlim[2]-xxlim[1])/lx) / ((yylim[2]-yylim[1])/ly)
		#if(rxyi<rxy)xxlim[1] = xxlim[1]-(rxy-rxyi)*(yylim[2]-yylim[1])
		if(is.null(xxlim))xxlim=c(min(idxx),max(idxx))
		if(is.null(yylim))yylim=c(min(idyy),max(idyy))
		idnew = which(xxlim[1]<=idxx & idxx<=xxlim[2] & yylim[1]<=idyy & idyy<=yylim[2])
		val = val[idnew]; idyy = idyy[idnew]; idxx = idxx[idnew]
	}
	if(vlim[2]==Inf){
		vlim = c(min(val),max(val))
	}else{
		val[val<vlim[1]] = vlim[1] 
		val[val>vlim[2]] = vlim[2] 
	}

	lat = sort(unique(yyg.po[idyy])); lon = sort(unique(xxg.po[idxx]))
	A = cbind(idyy,yyg.po[idyy])
	latint = latidx = lonint = lonidx = NULL
	for(i in 1:(length(lat)-1)){
		if(floor(lat[i])<floor(lat[i+1])){
			latint = c(latint,floor(lat[i+1]))
			latidx = c(latidx,(A[which(A[,2]==lat[i])[1],1]*abs(lat[i+1]-floor(lat[i+1]))+A[which(A[,2]==lat[i+1])[1],1]*abs(floor(lat[i+1])-lat[i]))/(abs(lat[i+1]-floor(lat[i+1]))+abs(floor(lat[i+1])-lat[i])) - min(idyy) + 1)
		}
	}
	A = cbind(idxx,xxg.po[idxx])
	for(i in 1:(length(lon)-1)){
		if(floor(lon[i])<floor(lon[i+1])){
			lonint = c(lonint,floor(lon[i+1]))
			lonidx = c(lonidx,(A[which(A[,2]==lon[i])[1],1]*abs(lon[i+1]-floor(lon[i+1]))+A[which(A[,2]==lon[i+1])[1],1]*abs(floor(lon[i+1])-lon[i]))/(abs(lon[i+1]-floor(lon[i+1]))+abs(floor(lon[i+1])-lon[i])) - min(idxx) + 1)
		}
	}
	idyy = idyy - min(idyy) + 1; idxx = idxx - min(idxx) + 1
	ylimits = quantile(idyy,c(0,1)); xlimits = quantile(idxx,c(0,1))
	r = (ylimits[2]-ylimits[1])/(xlimits[2]-xlimits[1])
	if(r > 1.5){
		cx = ceiling((((ylimits[2]-ylimits[1])/1.5)-xlimits[2]-xlimits[1])/2)
		xlimits = c(xlimits[1]-cx,xlimits[2]+cx)
	}else if(r < 1/1.5){
		cy = ceiling((((xlimits[2]-xlimits[1])/1.5)-ylimits[2]-ylimits[1])/2)
		ylimits = c(ylimits[1]-cx,ylimits[2]+cy)
	}

	xmat = matrix(NA,nrow=max(idyy),ncol=max(idxx))
	for(i in 1:length(val))
		xmat[idyy[i],idxx[i]] = val[i]
	xmatd = melt(xmat)
	if(is.null(pal))
		pal = c("#FFFFFF","#0000FF","#00FFFF","#00FF00","#FFFF00","#FF7700","#FF0000")
	xg = ggplot(xmatd,aes(Var2,Var1)) +
			geom_tile(aes(fill=value))+
			scale_fill_gradientn(colors=pal,limits=vlim) + 
			xlab("Longitude") + ylab("Latitude") +
			scale_y_continuous(breaks=latidx,labels=latint,expand=c(0,0),limits=ylimits) + 
			scale_x_continuous(breaks=lonidx,labels=lonint,expand=c(0,0),limits=xlimits) + 
			theme(legend.position=ifelse(islg,"right","none"))
	g = ggarrange(xg,labels=c(lab))
	ggsave(file=fn, plot=g, dpi=100, width=8, height=8)
}

# Draw maps of the PO and PA data

idsub.po = area.po[,is]
idsub.pa = area.pa[,is]
xy.pre = PO40[names(presence_list40)[iy]==PO40[,"st.species"],c("lon","lat")]
ypo = matrix(0,nrow=nrow(POenv),ncol=1)
for(i in 1:nrow(xy.pre)){
	id = which(xy.pre[i,"lon"]==POenv[,"lon"] & xy.pre[i,"lat"]==POenv[,"lat"])
	ypo[id] = ypo[id] + 1
if(length(id)<1)print(paste("unmatch",i))
}
ypo.all = ypo
ypo.sub = matrix(ypo.all[idsub.po,],ncol=1)
hplot(ypo.sub,iyenv.po[idsub.po],ixenv.po[idsub.po],xxlim=xxlim,yylim=yylim,fn=paste(dirs,"/",ari,"_ypo.eps",sep=""),pal=c("#B9FFF8","#cc2828"),islg=F)
#hplot(ypo.sub,iyenv.po[idsub.po],ixenv.po[idsub.po],xxlim=NULL,yylim=NULL,fn=paste(dirs,"/",ari,"_ypo_po_xylim=null.eps",sep=""),pal=c("#B9FFF8","#cc2828"),islg=F)

ypa.sub = matrix(presence_list40[[iy]][idsub.pa],ncol=1)
#hplot(ypa.sub,iyenv.pa[idsub.pa],ixenv.pa[idsub.pa],xxlim=xxlim,yylim=yylim,fn=paste(dirs,"/",ari,"_ypa_res.pa=4.eps",sep="")) #original resolution

iypa.sub = iyenvpo.pa[idsub.pa]; ixpa.sub = ixenvpo.pa[idsub.pa]
uniy = unique(iypa.sub); unix = unique(ixpa.sub)
iypa.new = ixpa.new = ypa.new = NULL
for(jy in uniy)for(jx in unix){
	id = jy==iypa.sub & jx==ixpa.sub
	if(sum(id)>0){
		iypa.new = c(iypa.new,jy); ixpa.new = c(ixpa.new,jx)
		ypa.new = c(ypa.new,as.numeric(sum(ypa.sub[id])>0))
	}
}
hplot(ypa.new,iypa.new,ixpa.new,xxlim=xxlim,yylim=yylim,fn=paste(dirs,"/",ari,"_ypa_res.pa=res.po.eps",sep=""),pal=c("#FFFFFF","#cc2828"),islg=F) #PO data resolution
#hplot(ypa.new,iypa.new,ixpa.new,xxlim=NULL,yylim=NULL,fn=paste(dirs,"/",ari,"_ypa_res.pa=res.po_xylim=null.eps",sep=""),pal=c("#FFFFFF","#cc2828"),islg=F)



# Draw heatmap of the predicted suitability

vv = unlist(read.csv(paste(dirs,"/vvori.csv",sep=""),header=F))
xpo = as.matrix(cbind(1,POenv[idsub.po,vv]))
xpa = as.matrix(cbind(1,PAenv[idsub.pa,vv]))
xpo[,vv] = scale(xpo[,vv])
xpa[,vv] = scale(xpa[,vv])

be0 = read.csv(paste(dirs,"/",ari,"_",tpl0,"_",qstr,"_be.csv",sep=""),header=F)
be0 = matrix(unlist(be0[,2]),ncol=1)
bex = read.csv(paste(dirs,"/",ari,"_",tplx,"_",qstr,"_be.csv",sep=""),header=F)
bex = matrix(unlist(bex[,2]),ncol=1)
lam0 = xpa%*%be0; lamx = xpa%*%bex
pROC::roc(ypa.sub>0,lam0,direction="<")$auc; pROC::roc(ypa.sub>0,lamx,direction="<")$auc

lam0.new = lamx.new = NULL
for(jxy in 1:length(iypa.new)){
	id = iypa.new[jxy]==iypa.sub & ixpa.new[jxy]==ixpa.sub
	lam0.new = c(lam0.new,mean(lam0[id]))
	lamx.new = c(lamx.new,mean(lamx[id]))
}
A = lam0.new-min(lam0.new); A = A/max(A)
hplot(A,iypa.new,ixpa.new,xxlim=xxlim,yylim=yylim,fn=paste(dirs,"/",ari,"_scaled-log-lam_",tpl0,"_",qstr,"_res.po.eps",sep=""))
A = lamx.new-min(lamx.new); A = A/max(A)
hplot(A,iypa.new,ixpa.new,xxlim=xxlim,yylim=yylim,fn=paste(dirs,"/",ari,"_scaled-log-lam_",tplx,"_",qstr,"_res.po.eps",sep=""))


# Draw a barplot of the parameter estimates

matbe = as.matrix(rbind(c(tpl0,be0),c(tplx,bex)))
#levels = c(tpl0,tplx)
levels = c("MLE","MIDE")
vvp = c("Intercept",vv)

colnames(matbe) = c("lab",vvp)
nt = nrow(matbe)
cc = heat.colors(nt)

dg = matrix(nrow=0,ncol=3)
for(l in levels){
	#dg = rbind(dg,cbind(paste(colnames(x),ifelse(be==0,"",i),sep=""),apply(d,2,median),be))
	dg = rbind(dg,cbind(paste(vvp,ifelse(which(l==levels)==1,"",paste("___",which(l==levels),sep="")),sep=""),unlist(matbe[l==levels,vvp]),rep(l,length(vvp))))
}
dg = as.data.frame(dg)
colnames(dg) = c("par","coef","Estimates")
llv = NULL
for(p in rev(vvp)){#rev(c("Intercept",vvnew))){
	for(it in c("",paste("___",2:length(levels),sep="")))
		llv = c(llv,paste(p,it,sep=""))
}
dg[,"par"] = factor(dg[,"par"],levels=llv)
dg[,"coef"] = as.numeric(dg[,"coef"])
dg[,"Estimates"] = factor(dg[,"Estimates"],levels=levels)
#cols = c("#0000FF","#00FFFF","#00FF00","#FFFF00","#FF0000","#000000")[1:nt]
cols = c("#F44336","#536DFE")
g = ggplot(dg,aes(par,coef,fill=Estimates))+geom_bar(stat="identity") + 
		scale_x_discrete(breaks=vvp) + 
		coord_flip() + 
		scale_fill_manual(values=cols)+
		theme(legend.position="top")
print(g)
ggsave(file=paste(dirs,"/",ari,"_coef_",tplx,"_",qstr,".eps",sep=""), plot=g, dpi=100, width=8, height=8); dev.off()


