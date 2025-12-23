## analyses of genitalic and wing variation, comparing means and SD for admixed populations
## vs lab hybrids
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


## male genitalia
weens2<-read.csv("artificial_hybrids/analyses/F1vsAlpine_ween_7vii13b.csv",header=T)
pcaout3<-prcomp(weens2[,4:8],center=T,scale.=T)
mns2<-tapply(X=pcaout3$x[,1],INDEX=weens2[,1],mean)


## male wings
hyb3<-read.csv("artificial_hybrids/analyses/F1vsAlpine_malespotsize_13vii13_standardized.csv",header=TRUE)
hyb3[,7:23]
pcaout4<-prcomp(hyb3[,7:23],center=T,scale.=T)
mns3<-tapply(X=pcaout4$x[,1],INDEX=hyb3[,5],mean)
mns4<-tapply(X=pcaout4$x[,2],INDEX=hyb3[,5],mean)


## female wings
hyb5<-read.csv("artificial_hybrids/analyses/F1vsAlpine_femalespotsize_16vii13_standardized.csv",header=TRUE)
hyb5[,7:23]
pcaout5<-prcomp(hyb5[,7:23],center=T,scale.=T)
mns5<-tapply(X=pcaout5$x[,1],INDEX=hyb5[,5],mean)
mns6<-tapply(X=pcaout5$x[,2],INDEX=hyb5[,5],mean)


## Lauren's original plot
pdf("PCs_F1_HHS_Fig5_31vii13.pdf",width=10,height=7)
par(mfrow=c(2,3),mar=c(5,4.5,1.5,1))

plot(mns2[order(mns2)],axes=F,xlab="",ylab="Genitalia shape (mean PC1, sd)",ylim=c(-3,3))
title(main="A.",adj=0)
pops<-sort(unique(as.character(weens2[,1])))
sds<-tapply(X=pcaout3$x[,1],INDEX=weens2[,1],sd)
segments(1:8,mns2[order(mns2)] + sds[order(mns2)],1:8,mns2[order(mns2)] - sds[order(mns2)])
axis(1,at=1:8,pops[order(mns2)],cex=0.7,las=3)
axis(2)
box()

plot(mns3[order(mns3)],axes=F,xlab="",ylab="Wing element size (mean PC1, sd)",ylim=c(-6,6))
title(main="B.",adj=0)
pops<-sort(unique(as.character(hyb3[,5])))
sds<-tapply(X=pcaout4$x[,1],INDEX=hyb3[,5],sd)
segments(1:8,mns3[order(mns3)] + sds[order(mns3)],1:8,mns3[order(mns3)] - sds[order(mns3)])
axis(1,at=1:8,pops[order(mns3)],cex=0.7,las=3)
axis(2)
box()


plot(mns4[order(mns4)],axes=F,xlab="",ylab="Wing aurorae vs. spot size (mean PC2, sd)",ylim=c(-2.2,2.1))
title(main="C.",adj=0)
pops<-sort(unique(as.character(hyb3[,5])))
sds<-tapply(X=pcaout4$x[,2],INDEX=hyb3[,5],sd)
segments(1:8,mns4[order(mns4)] + sds[order(mns4)],1:8,mns4[order(mns4)] - sds[order(mns4)])
axis(1,at=1:8,pops[order(mns4)],cex=0.7,las=3)
axis(2)
box()

plot(mns5[order(mns5)],axes=F,xlab="",ylab="Wing element size (mean PC1, sd)",ylim=c(-6,6))
title(main="D.",adj=0)
pops<-sort(unique(as.character(hyb5[,5])))
sds<-tapply(X=pcaout5$x[,1],INDEX=hyb5[,5],sd)
segments(1:8,mns5[order(mns5)] + sds[order(mns5)],1:8,mns5[order(mns5)] - sds[order(mns5)])
axis(1,at=1:8,pops[order(mns5)],cex=0.7,las=3)
axis(2)
box()

plot(mns6[order(mns6)],axes=F,xlab="",ylab="Wing aurorae vs. spot size (mean PC2, sd)",ylim=c(-3,2))
title(main="E.",adj=0)
pops<-sort(unique(as.character(hyb5[,5])))
sds<-tapply(X=pcaout5$x[,2],INDEX=hyb5[,5],sd)
segments(1:8,mns6[order(mns6)] + sds[order(mns6)],1:8,mns6[order(mns6)] - sds[order(mns6)])
axis(1,at=1:8,pops[order(mns6)],cex=0.7,las=3)
axis(2)
box()

dev.off()

### Model fit
#### male genitalia #############

table(weens2$pop)
#
#annaXmel       CP      FCR melXanna       MR      SLA      VCP       YG 
#      25       20       28       25       20       49       45       37 
id<-as.numeric(as.factor(weens2$pop))
## stan model
D<-list(y=pcaout3$x[,1],id=id,Npop=8,N=length(id))
fitPc1<-stan("morpHHmod.stan",data=D)
fitPc1

## PC1
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#          mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#alpha[1]  0.44    0.00 0.15  0.13  0.34  0.44  0.53  0.75  6646    1
#alpha[2] -0.29    0.00 0.12 -0.52 -0.36 -0.29 -0.21 -0.05  5957    1
#alpha[3]  2.49    0.00 0.09  2.31  2.43  2.49  2.55  2.66  8475    1
#alpha[4]  0.18    0.00 0.14 -0.10  0.09  0.18  0.27  0.47  7302    1
#alpha[5]  0.92    0.00 0.14  0.63  0.83  0.92  1.01  1.20  6863    1
#alpha[6] -2.12    0.00 0.08 -2.28 -2.18 -2.12 -2.07 -1.97  7534    1
#alpha[7] -1.88    0.00 0.10 -2.06 -1.95 -1.88 -1.82 -1.70  8371    1
#alpha[8]  2.46    0.00 0.08  2.30  2.41  2.46  2.51  2.61  6982    1
#sig[1]    0.74    0.00 0.12  0.55  0.66  0.73  0.81  1.01  4829    1
#sig[2]    0.54    0.00 0.09  0.39  0.47  0.52  0.59  0.75  5196    1
#sig[3]    0.46    0.00 0.07  0.35  0.41  0.45  0.50  0.63  6573    1
#sig[4]    0.70    0.00 0.11  0.53  0.63  0.69  0.77  0.95  5614    1
#sig[5]    0.61    0.00 0.11  0.45  0.54  0.60  0.67  0.87  4832    1
#sig[6]    0.54    0.00 0.06  0.44  0.50  0.53  0.57  0.67  6767    1
#sig[7]    0.64    0.00 0.07  0.51  0.59  0.63  0.68  0.79  6611    1
#sig[8]    0.47    0.00 0.06  0.38  0.44  0.47  0.51  0.60  7691    1
#lp__     15.63    0.07 2.92  9.01 13.89 15.96 17.71 20.44  1562    1


#### male wings #############

table(hyb3$pop2)
#
#annaXmel       CP      FCR melXanna       MR      SLA      VCP       YG 
#      41       20       21       47       20       40       44       35 
id<-as.numeric(as.factor(hyb3$pop2))
## stan model
D<-list(y=pcaout4$x[,1],id=id,Npop=8,N=length(id))
fitMwPc1<-stan("morpHHmod.stan",data=D)
fitMwPc1
D<-list(y=pcaout4$x[,2],id=id,Npop=8,N=length(id))
fitMwPc2<-stan("morpHHmod.stan",data=D)
fitMwPc2

## PC1 
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.

#            mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#alpha[1]   -1.22    0.00 0.21   -1.62   -1.35   -1.22   -1.08   -0.79  7575
#alpha[2]    2.36    0.01 0.43    1.52    2.08    2.35    2.64    3.22  6360
#alpha[3]   -5.14    0.00 0.21   -5.56   -5.28   -5.14   -5.01   -4.72  8126
#alpha[4]   -1.42    0.00 0.27   -1.95   -1.60   -1.42   -1.23   -0.88  7339
#alpha[5]    1.96    0.01 0.60    0.77    1.57    1.98    2.35    3.13  7886
#alpha[6]    3.34    0.00 0.35    2.63    3.11    3.34    3.57    4.01  7653
#alpha[7]    3.76    0.00 0.22    3.30    3.62    3.76    3.91    4.20  7683
#alpha[8]   -4.60    0.00 0.19   -4.99   -4.73   -4.60   -4.47   -4.23  7764
#sig[1]      1.31    0.00 0.15    1.05    1.21    1.29    1.40    1.66  8004
#sig[2]      1.83    0.00 0.33    1.33    1.60    1.78    2.01    2.63  6226
#sig[3]      0.99    0.00 0.16    0.72    0.87    0.97    1.08    1.37  6452
#sig[4]      1.84    0.00 0.20    1.50    1.70    1.81    1.96    2.28  7474
#sig[5]      2.60    0.01 0.47    1.86    2.28    2.54    2.86    3.70  6366
#sig[6]      2.17    0.00 0.26    1.73    1.99    2.15    2.33    2.74  6302
#sig[7]      1.45    0.00 0.16    1.18    1.33    1.43    1.55    1.81  6880
#sig[8]      1.16    0.00 0.14    0.92    1.06    1.15    1.25    1.49  5875
#lp__     -247.62    0.07 2.91 -254.17 -249.28 -247.27 -245.48 -242.86  1572

## PC2
#4 chains, each with iter=2000; warmup=1000; thin=1; 
#post-warmup draws per chain=1000, total post-warmup draws=4000.
#
#            mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#alpha[1]   -0.03    0.00 0.15   -0.32   -0.14   -0.04    0.06    0.26  6610
#alpha[2]    0.19    0.00 0.28   -0.37    0.00    0.19    0.37    0.76  5353
#alpha[3]    1.19    0.00 0.11    0.97    1.12    1.19    1.27    1.42  5927
#alpha[4]   -1.18    0.00 0.16   -1.49   -1.28   -1.18   -1.07   -0.88  6541
#alpha[5]    0.86    0.00 0.30    0.28    0.67    0.86    1.06    1.44  6020
#alpha[6]   -0.58    0.00 0.22   -1.00   -0.72   -0.58   -0.44   -0.16  6081
#alpha[7]    0.06    0.00 0.15   -0.23   -0.03    0.07    0.16    0.35  7271
#alpha[8]    0.89    0.00 0.15    0.59    0.79    0.89    0.99    1.18  6336
#sig[1]      0.96    0.00 0.12    0.76    0.88    0.95    1.03    1.21  6354
#sig[2]      1.24    0.00 0.22    0.90    1.08    1.20    1.36    1.77  4681
#sig[3]      0.52    0.00 0.09    0.38    0.46    0.51    0.57    0.72  5606
#sig[4]      1.08    0.00 0.12    0.88    0.99    1.07    1.15    1.34  5533
#sig[5]      1.33    0.00 0.24    0.95    1.16    1.30    1.46    1.88  5055
#sig[6]      1.38    0.00 0.17    1.10    1.27    1.36    1.48    1.75  6335
#sig[7]      0.98    0.00 0.11    0.79    0.90    0.97    1.04    1.22  5322
#sig[8]      0.89    0.00 0.11    0.71    0.81    0.88    0.96    1.15  6280
#lp__     -133.59    0.07 2.99 -140.43 -135.32 -133.22 -131.42 -128.88  1606


#### female wings #############

table(hyb5$pop2)
#
#annaXmel       CP      FCR melXanna       MR      SLA      VCP       YG 
#      30       20       37       30       20       48       49       50
id<-as.numeric(as.factor(hyb5$pop2))
## stan model
D<-list(y=pcaout5$x[,1],id=id,Npop=8,N=length(id))
fitFwPc1<-stan("morpHHmod.stan",data=D)
fitFwPc1
D<-list(y=pcaout5$x[,2],id=id,Npop=8,N=length(id))
fitFwPc2<-stan("morpHHmod.stan",data=D)
fitFwPc2

## PC1
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.

#            mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#alpha[1]   -0.14    0.00 0.30   -0.73   -0.33   -0.13    0.07    0.45  7567
#alpha[2]   -1.36    0.00 0.36   -2.07   -1.59   -1.36   -1.13   -0.63  6187
#alpha[3]    4.34    0.00 0.20    3.94    4.20    4.34    4.48    4.74  7247
#alpha[4]    0.92    0.00 0.28    0.38    0.74    0.92    1.10    1.47  7632
#alpha[5]   -1.49    0.00 0.36   -2.21   -1.73   -1.50   -1.25   -0.79  5999
#alpha[6]   -3.62    0.00 0.25   -4.11   -3.79   -3.63   -3.46   -3.13  7709
#alpha[7]   -3.58    0.00 0.26   -4.08   -3.75   -3.58   -3.41   -3.08  8164
#alpha[8]    4.44    0.00 0.16    4.13    4.33    4.44    4.54    4.75  7644
#sig[1]      1.63    0.00 0.22    1.25    1.47    1.61    1.77    2.12  7489
#sig[2]      1.63    0.00 0.29    1.19    1.43    1.59    1.79    2.31  6270
#sig[3]      1.23    0.00 0.15    0.98    1.12    1.21    1.32    1.56  6383
#sig[4]      1.53    0.00 0.21    1.18    1.38    1.50    1.65    2.02  5965
#sig[5]      1.55    0.00 0.28    1.12    1.35    1.51    1.71    2.17  5532
#sig[6]      1.75    0.00 0.19    1.43    1.62    1.74    1.87    2.16  7530
#sig[7]      1.81    0.00 0.19    1.48    1.68    1.79    1.92    2.23  6462
#sig[8]      1.11    0.00 0.12    0.91    1.03    1.10    1.18    1.37  6371
#lp__     -247.68    0.07 2.89 -254.31 -249.44 -247.32 -245.56 -243.15  1621


## PC2
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.

#            mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#alpha[1]   -0.05    0.00 0.20   -0.45   -0.18   -0.05    0.09    0.36  7723
#alpha[2]    0.20    0.00 0.20   -0.19    0.07    0.20    0.33    0.62  7056
##alpha[3]    0.49    0.00 0.16    0.18    0.39    0.49    0.59    0.79  8417
#alpha[4]   -1.74    0.00 0.18   -2.09   -1.85   -1.74   -1.63   -1.39  8914
#alpha[5]    0.90    0.00 0.24    0.43    0.74    0.89    1.05    1.38  6403
#alpha[6]   -0.44    0.00 0.22   -0.87   -0.57   -0.43   -0.30   -0.02  7657
#alpha[7]    0.02    0.00 0.14   -0.25   -0.07    0.02    0.11    0.30  7750
#alpha[8]    0.66    0.00 0.11    0.44    0.59    0.66    0.73    0.87  8013
#sig[1]      1.11    0.00 0.16    0.86    1.00    1.09    1.20    1.46  6501
#sig[2]      0.90    0.00 0.16    0.65    0.78    0.87    0.99    1.26  6094
#sig[3]      0.96    0.00 0.12    0.77    0.88    0.95    1.03    1.23  5639
#sig[4]      0.98    0.00 0.14    0.76    0.88    0.97    1.06    1.29  6188
#sig[5]      1.05    0.00 0.18    0.77    0.92    1.03    1.16    1.47  5022
#sig[6]      1.49    0.00 0.16    1.22    1.38    1.48    1.59    1.85  6202
#sig[7]      0.97    0.00 0.11    0.79    0.89    0.96    1.03    1.21  7540
#sig[8]      0.75    0.00 0.08    0.62    0.70    0.74    0.80    0.92  7061
#lp__     -138.86    0.08 2.92 -145.24 -140.65 -138.60 -136.70 -134.22  1509
save(list=ls(),file="morph2.rdat")


