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

#### male genitalia #############

# getting specific population combos
h<-grep(hyb$pop,pattern="x")
ids<-unlist(strsplit(hyb$pop[h],split="x"))
id1<-rep(NA,dim(hyb)[1])
id2<-rep(NA,dim(hyb)[1])

id1[h]<-ids[seq(1,length(ids),2)]
id2[h]<-ids[seq(2,length(ids),2)]
id1[-h]<-hyb$pop[-h]
id2[-h]<-hyb$pop[-h]
unique(id1)
#[1] "ann" "fcr" "mel" "sla" "vcp" "yg"
unique(id2)
#[1] "mel" "fcr" "sla" "ann" "vcp" "yg" 
## drop anna and mel ids
drop<-which(id1 %in% c("ann","mel"))

id1<-id1[-drop]
id2<-id2[-drop]
h<-h[-drop]

# id1 index/id: 1=fcr, 2=sla, 3=vcp, 4=yg
nid1<-as.numeric(as.factor(id1))
nid2<-as.numeric(as.factor(id2))
hid<-as.numeric(nid1!=nid2)
annaM<-as.numeric(nid1==1 | nid1==4)

## stan model
D<-list(y=pcaout$x[-drop,1],nid1=nid1,nid2=nid2,hid=hid,annaM=annaM,Npop=4,N=length(nid1))
fit1Pc1<-stan("wingF1mod.stan",data=D)
fit1Pc1
D<-list(y=pcaout$x[-drop,2],nid1=nid1,nid2=nid2,hid=hid,annaM=annaM,Npop=4,N=length(nid1))
fit1Pc2<-stan("wingF1mod.stan",data=D)
fit1Pc2

## PC1
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#
#          mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#alpha[1] -3.18    0.00 0.16 -3.50 -3.28 -3.18 -3.07 -2.86  3939    1
#alpha[2]  1.68    0.00 0.11  1.46  1.61  1.68  1.75  1.88  4834    1
#alpha[3]  1.70    0.00 0.11  1.49  1.62  1.69  1.77  1.91  3783    1
#alpha[4] -3.17    0.00 0.21 -3.57 -3.31 -3.17 -3.03 -2.75  3719    1
#betaH    -0.08    0.00 0.19 -0.44 -0.21 -0.08  0.04  0.28  2780    1
#betaAM   -0.23    0.00 0.23 -0.68 -0.39 -0.23 -0.07  0.23  2745    1
#sig       0.60    0.00 0.04  0.52  0.57  0.60  0.62  0.69  4172    1
#lp__      2.25    0.04 1.92 -2.22  1.19  2.60  3.66  4.99  1922    1

## PC2
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.

#           mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#alpha[1]  -0.32    0.00 0.25  -0.80  -0.48  -0.31  -0.16   0.17  3964    1
#alpha[2]   0.33    0.00 0.16   0.02   0.22   0.33   0.44   0.65  4233    1
#alpha[3]  -0.05    0.00 0.16  -0.38  -0.15  -0.04   0.06   0.28  4018    1
#alpha[4]   0.32    0.01 0.32  -0.33   0.11   0.32   0.53   0.94  4091    1
#betaH      0.12    0.01 0.28  -0.43  -0.07   0.12   0.31   0.66  3078    1
#betaAM    -0.23    0.01 0.35  -0.91  -0.47  -0.23   0.01   0.45  3141    1
#sig        0.91    0.00 0.06   0.80   0.87   0.91   0.96   1.04  4023    1
#lp__     -44.74    0.05 2.02 -49.58 -45.76 -44.39 -43.28 -41.92  1823    1





#### male wings #############

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

#PC 1, all rhat 1
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#
#            mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#alpha[1]   -5.07    0.01 0.46   -5.98   -5.38   -5.08   -4.76   -4.17  3276
#alpha[2]    3.23    0.00 0.31    2.64    3.03    3.23    3.42    3.83  4714
#alpha[3]    4.32    0.00 0.30    3.71    4.13    4.33    4.52    4.89  3790
#alpha[4]   -5.09    0.01 0.52   -6.13   -5.43   -5.08   -4.75   -4.07  3825
#betaH      -0.96    0.01 0.35   -1.65   -1.19   -0.95   -0.72   -0.28  2516
#betaAM      0.18    0.01 0.37   -0.54   -0.06    0.18    0.42    0.91  3800
#sig         1.70    0.00 0.10    1.52    1.63    1.70    1.76    1.90  3640
#lp__     -167.35    0.05 1.95 -172.05 -168.42 -167.01 -165.89 -164.56  1652

## PC2, all rhat 1
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#
#            mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#alpha[1]    1.65    0.01 0.31    1.01    1.45    1.65    1.86    2.26  3027
#alpha[2]   -0.54    0.00 0.21   -0.96   -0.68   -0.54   -0.39   -0.13  3768
#alpha[3]    0.37    0.00 0.21   -0.04    0.23    0.37    0.50    0.78  3790
#alpha[4]    1.78    0.01 0.37    1.05    1.54    1.78    2.01    2.49  3508
#betaH      -1.78    0.00 0.24   -2.25   -1.94   -1.79   -1.62   -1.31  2417
#betaAM      1.25    0.00 0.27    0.72    1.07    1.25    1.43    1.76  2966
#sig         1.19    0.00 0.07    1.07    1.15    1.19    1.24    1.33  4362
#lp__     -109.35    0.04 1.88 -113.90 -110.38 -109.02 -107.97 -106.65  1839

#### female wings #############

# getting specific population combos
h<-grep(fhyb2$pop,pattern="X")
ids<-unlist(strsplit(fhyb2$pop[h],split="X"))
id1<-rep(NA,dim(fhyb2)[1])
id2<-rep(NA,dim(fhyb2)[1])

id1[h]<-ids[seq(1,length(ids),2)]
id2[h]<-ids[seq(2,length(ids),2)]
id1[-h]<-fhyb2$pop[-h]
id2[-h]<-fhyb2$pop[-h]
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
D<-list(y=fpcaout2$x[,1],nid1=nid1,nid2=nid2,hid=hid,annaM=annaM,Npop=4,N=length(nid1))
fit3Pc1<-stan("wingF1mod.stan",data=D)
fit3Pc1
D<-list(y=fpcaout2$x[,2],nid1=nid1,nid2=nid2,hid=hid,annaM=annaM,Npop=4,N=length(nid1))
fit3Pc2<-stan("wingF1mod.stan",data=D)
fit3Pc2
## works great, two anna pops similar, two melissa pops similar, credible effect of hybrid, anna male has an effect of shifting towards melissa, PC2 more complicated more population variation

## PC1, all rhat = 1, but print wraps to new line
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#
#            mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#alpha[1]   -3.89    0.01 0.39   -4.65   -4.15   -3.89   -3.63   -3.14  4298
#alpha[2]    3.69    0.00 0.31    3.08    3.47    3.68    3.90    4.31  4685
#alpha[3]    4.61    0.01 0.32    3.99    4.38    4.60    4.83    5.24  3754
#alpha[4]   -4.44    0.01 0.32   -5.07   -4.65   -4.43   -4.22   -3.81  3960
#betaH      -1.17    0.01 0.35   -1.87   -1.41   -1.17   -0.93   -0.50  2491
#betaAM      1.32    0.01 0.40    0.55    1.04    1.32    1.59    2.13  3321
#sig         1.50    0.00 0.10    1.32    1.43    1.50    1.56    1.70  3350
#lp__     -119.91    0.05 1.94 -124.71 -120.91 -119.60 -118.49 -117.19  1650

## PC 2
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#
#           mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#alpha[1]   0.81    0.00 0.27   0.29   0.62   0.81   1.00   1.35  4557    1
#alpha[2]  -0.56    0.00 0.23  -1.02  -0.72  -0.55  -0.40  -0.10  4187    1
#alpha[3]   0.37    0.00 0.23  -0.09   0.21   0.37   0.53   0.83  3953    1
#alpha[4]   1.20    0.00 0.22   0.77   1.05   1.20   1.35   1.64  3802    1
#betaH     -2.22    0.00 0.24  -2.69  -2.38  -2.21  -2.05  -1.74  2754    1
#betaAM     2.08    0.00 0.28   1.54   1.90   2.08   2.27   2.64  3249    1
#sig        1.08    0.00 0.07   0.96   1.03   1.08   1.12   1.22  4204    1
#lp__     -75.83    0.05 1.91 -80.56 -76.81 -75.49 -74.46 -73.16  1772    1

save(list=ls(),file="morph1.rdat")

