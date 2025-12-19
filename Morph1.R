## analyses of genitalic and wing variation, testing for intermediacy of artifical hybrids
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


## male genitalia
hyb<-read.csv("artificial_hybrids/analyses/F1vsParent_ween_7vii13_labonly.csv",header=TRUE)
pcaout<-prcomp(hyb[,5:9],center=TRUE,scale=TRUE)

## male wings
hyb2<-read.csv("artificial_hybrids/analyses/F1vsParents_malespotsize_9vii13_standardized.csv",header=TRUE)
pcaout2<-prcomp(hyb2[,7:23],center=TRUE,scale=TRUE)

## female wings
fhyb2<-read.csv("artificial_hybrids/analyses/F1vsParents_femalespotsize_16vii13_standardized.csv",header=TRUE)
fpcaout2<-prcomp(fhyb2[,7:23],center=TRUE,scale=TRUE)

## PCA plots
cl<-1.75;cm<-1.65;ca<-1.35
pdf("fig_morphPCA.pdf",width=6,height=12)
layout(matrix(1:6,nrow=3,ncol=2,byrow=TRUE),widths=c(4,2),heights=c(4,4,4))

# m genitalia
## size is really length, not wdith, which loads negative, shape is mostly H and widths vs F and U
par(mar=c(4.5,5.5,2.5,1.5))
plot(pcaout$x[,1],pcaout$x[,2],type="n",xlab="PC1 (69.5%)",ylab="PC2 (16.9%)",cex.lab=cl,cex.axis=ca)
title(main="(a) Male genitalia",adj=0,cex.main=cm)
points(pcaout$x[hyb[,3] == 'yg',1],pcaout$x[hyb[,3] == 'yg',2],col="firebrick",pch=19)
points(pcaout$x[hyb[,3] == 'fcr',1],pcaout$x[hyb[,3] == 'fcr',2],col="darkorange",pch=19)
points(pcaout$x[hyb[,3] == 'sla',1],pcaout$x[hyb[,3] == 'sla',2],col="lightskyblue",pch=19)
points(pcaout$x[hyb[,3] == 'vcp',1],pcaout$x[hyb[,3] == 'vcp',2],col="cadetblue",pch=19)
points(pcaout$x[hyb[,3] == 'ygxvcp',1],pcaout$x[hyb[,3] == 'ygxvcp',2],col="gray30",pch=0)
points(pcaout$x[hyb[,3] == 'fcrxvcp',1],pcaout$x[hyb[,3] == 'fcrxvcp',2],col="gray30",pch=1)
points(pcaout$x[hyb[,3] == 'fcrxsla',1],pcaout$x[hyb[,3] == 'fcrxsla',2],col="gray30",pch=2)
points(pcaout$x[hyb[,3] == 'ygxsla',1],pcaout$x[hyb[,3] == 'ygxsla',2],col="gray30",pch=5)
points(pcaout$x[hyb[,3] == 'annxmel',1],pcaout$x[hyb[,3] == 'annxmel',2],pch=20)
points(pcaout$x[hyb[,3] == 'vcpxfcr',1],pcaout$x[hyb[,3] == 'vcpxfcr',2],pch=3)
points(pcaout$x[hyb[,3] == 'vcpxyg',1],pcaout$x[hyb[,3] == 'vcpxyg',2],pch=4)
points(pcaout$x[hyb[,3] == 'slaxfcr',1],pcaout$x[hyb[,3] == 'slaxfcr',2],pch=8)
points(pcaout$x[hyb[,3] == 'slaxyg',1],pcaout$x[hyb[,3] == 'slaxyg',2],pch=11)
points(pcaout$x[hyb[,3] == 'melxann',1],pcaout$x[hyb[,3] == 'melxann',2],pch=21)

par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),xlab="",ylab="",axes=FALSE,type='n')
legend(x=0,y=.9,c("L. anna (YG)","L. anna (FCR)","L. melissa (SLA)","L. melissa (VCP)","L. anna X L. melissa","L. melissa X L. anna","YG X VCP","FCR X VCP","FCR X SLA","YG X SLA","VCP X FCR","VCP X YG","SLA X FCR","SLA X YG"),col=c("firebrick","darkorange","lightskyblue","cadetblue","gray30","black","gray30","gray30","gray30","gray30","black","black","black","black"),pch=c(19,19,19,19,20,21,0,1,2,5,3,4,8,11),cex=1.35,bty="n")


# m wings, for wings (m and f) PC1 clearly element size and PC2 small aurorae vs big spots
par(mar=c(4.5,5.5,2.5,1.5))
plot(pcaout2$x[,1],pcaout2$x[,2],type="n",xlab="PC1 (70.0%)",ylab="PC2 (12.1%)",cex.lab=cl,cex.axis=ca)
title(main="(b) Male wing pattern",adj=0,cex.main=cm)
points(pcaout2$x[hyb2[,4] == 'yg',1],pcaout2$x[hyb2[,4] == 'yg',2],col="firebrick",pch=19)
points(pcaout2$x[hyb2[,4] == 'fcr',1],pcaout2$x[hyb2[,4] == 'fcr',2],col="darkorange",pch=19)
points(pcaout2$x[hyb2[,4] == 'sla',1],pcaout2$x[hyb2[,4] == 'sla',2],col="lightskyblue",pch=19)
points(pcaout2$x[hyb2[,4] == 'vcp',1],pcaout2$x[hyb2[,4] == 'vcp',2],col="cadetblue",pch=19)
points(pcaout2$x[hyb2[,4] == 'ygXvcp',1],pcaout2$x[hyb2[,4] == 'ygXvcp',2],col="gray30",pch=0)
points(pcaout2$x[hyb2[,4] == 'fcrXvcp',1],pcaout2$x[hyb2[,4] == 'fcrXvcp',2],col="gray30",pch=1)
points(pcaout2$x[hyb2[,4] == 'fcrXsla',1],pcaout2$x[hyb2[,4] == 'fcrXsla',2],col="gray30",pch=2)
points(pcaout2$x[hyb2[,4] == 'ygXsla',1],pcaout2$x[hyb2[,4] == 'ygXsla',2],col="gray30",pch=5)
points(pcaout2$x[hyb2[,4] == 'vcpXfcr',1],pcaout2$x[hyb2[,4] == 'vcpXfcr',2],pch=3)
points(pcaout2$x[hyb2[,4] == 'vcpXyg',1],pcaout2$x[hyb2[,4] == 'vcpXyg',2],pch=4)
points(pcaout2$x[hyb2[,4] == 'slaXfcr',1],pcaout2$x[hyb2[,4] == 'slaXfcr',2],pch=8)
points(pcaout2$x[hyb2[,4] == 'slaXyg',1],pcaout2$x[hyb2[,4] == 'slaXyg',2],pch=11)

par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),xlab="",ylab="",axes=FALSE,type='n')
legend(x=0,y=.9,c("L. anna (YG)","L. anna (FCR)","L. melissa (SLA)","L. melissa (VCP)","YG X VCP","FCR X VCP","FCR X SLA","YG X SLA","VCP X FCR","VCP X YG","SLA X FCR","SLA X YG"),col=c("firebrick","darkorange","lightskyblue","cadetblue","gray30","gray30","gray30","gray30","black","black","black","black"),pch=c(19,19,19,19,0,1,2,5,3,4,8,11),cex=1.35,bty="n")

## f wings
par(mar=c(4.5,5.5,2.5,1.5))
plot(fpcaout2$x[,1],fpcaout2$x[,2],type="n",xlab="PC1 (71.7%)",ylab="PC2 (11.7%)",cex.lab=cl,cex.axis=ca)
title(main="(c) Female wing pattern",adj=0,cex.main=cm)
points(fpcaout2$x[fhyb2[,4] == 'yg',1],fpcaout2$x[fhyb2[,4] == 'yg',2],col="firebrick",pch=19)
points(fpcaout2$x[fhyb2[,4] == 'fcr',1],fpcaout2$x[fhyb2[,4] == 'fcr',2],col="darkorange",pch=19)
points(fpcaout2$x[fhyb2[,4] == 'sla',1],fpcaout2$x[fhyb2[,4] == 'sla',2],col="lightskyblue",pch=19)
points(fpcaout2$x[fhyb2[,4] == 'vcp',1],fpcaout2$x[fhyb2[,4] == 'vcp',2],col="cadetblue",pch=19)
points(fpcaout2$x[fhyb2[,4] == 'ygXvcp',1],fpcaout2$x[fhyb2[,4] == 'ygXvcp',2],col="gray30",pch=0)
points(fpcaout2$x[fhyb2[,4] == 'fcrXvcp',1],fpcaout2$x[fhyb2[,4] == 'fcrXvcp',2],col="gray30",pch=1)
points(fpcaout2$x[fhyb2[,4] == 'ygXsla',1],fpcaout2$x[fhyb2[,4] == 'ygXsla',2],col="gray30",pch=5)
points(fpcaout2$x[fhyb2[,4] == 'vcpXfcr',1],fpcaout2$x[fhyb2[,4] == 'vcpXfcr',2],pch=3)
points(fpcaout2$x[fhyb2[,4] == 'vcpXyg',1],fpcaout2$x[fhyb2[,4] == 'vcpXyg',2],pch=4)

par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),xlab="",ylab="",axes=FALSE,type='n')
legend(x=0,y=.9,c("L. anna (YG)","L. anna (FCR)","L. melissa (SLA)","L. melissa (VCP)","YG X VCP","FCR X VCP","YG X SLA","VCP X FCR","VCP X YG"),col=c("firebrick","darkorange","lightskyblue","cadetblue","gray30","gray30","gray30","black","black"),pch=c(19,19,19,19,0,1,5,3,4),cex=1.35,bty="n")

dev.off()


# getting specific population combos
h<-grep(hyb2$pop,pattern="X")
ids<-unlist(strsplit(hyb2$pop[h],split="X"))
id1<-rep(NA,dim(hyb2)[1])
id2<-rep(NA,dim(hyb2)[1])

id1[h]<-ids[seq(1,length(ids),2)]
id2[h]<-ids[seq(2,length(ids),2)]
id1[-h]<-hyb2$pop[-h]
id2[-h]<-hyb2$pop[-h]
unique(id1)
#[1] "yg"  "fcr" "vcp" "sla"
unique(id2)
#[1] "yg"  "fcr" "vcp" "sla"

# id1 index/id: 1=fcr, 2=sla, 3=vcp, 4=yg
nid1<-as.numeric(as.factor(id1))
nid2<-as.numeric(as.factor(id2))
hid<-as.numeric(nid1!=nid2)
annaM<-as.numeric(nid1==1 | nid1==4)

## stan model
D<-list(y=pcaout2$x[,1],nid1=nid1,nid2=nid2,hid=hid,annaM=annaM,Npop=4,N=length(nid1))
fit2Pc1<-stan("wingF1mod.stan",data=D)
fit2Pc1
D<-list(y=pcaout2$x[,2],nid1=nid1,nid2=nid2,hid=hid,annaM=annaM,Npop=4,N=length(nid1))
fit2Pc2<-stan("wingF1mod.stan",data=D)
fit2Pc2
## works great, two anna pops similar, two melissa pops similar, credible effect of hybrid, shifted towards anna
## not really an effect of anna male
