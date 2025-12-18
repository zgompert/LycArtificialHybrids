## analyses of genitalic and wing variation, testing for intermediacy of artifical hybrids
library(rstan)
options(mc.cores = parallel::detectCores()).
rstan_options(auto_write = TRUE)


## need to figure out to what to do with male genitalia data, some are only anna x melissa
hyb<-read.csv("artificial_hybrids/analyses/F1vsParent_ween_7vii13_labonly.csv",header=TRUE)
pcaout<-prcomp(hyb[,5:9],center=TRUE,scale=TRUE)

## trying the Bayesian model with male wings, mostly a test to get it working
hyb2<-read.csv("artificial_hybrids/F1wings_malespotsize_17xii12_standardized.csv",header=TRUE)
pcaout2<-prcomp(hyb2[,6:22],center=TRUE,scale=TRUE)

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
fit2<-stan("wingF1mod.stan",data=D)
fit2
## works great, two anna pops similar, two melissa pops similar, credible effect of hybrid, shifted towards anna
## not really an effect of anna male
