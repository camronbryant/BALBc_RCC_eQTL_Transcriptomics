#----Load Packages----
library(XLConnect)
library(qtl2)
library(qtl)
library(yaml)
library(xlsx)
QT <- function(x) orderNorm(x)$x.t

names(BALB$pheno)
.
#------rQTL file input and formatting----
BALB <-read.cross(format=c("csv"), file="BALBc QTL Workbook ChiSq 5-4-2021.csv",na.strings=c("-"),
                  genotypes=c("J","H","B"),alleles=c("J","B"),estimate.map=TRUE,
                  convertXdata=TRUE, error.prob=0.0001)
BALB<-jittermap(BALB)


#----Final Genetic QC----
#DELETING MARKERS W CALL RATE < 95%
#total # individuals=378 5% of 378 = 18.9
na<-nmissing(BALB, what="mar")
na<-as.data.frame(na)
na<-data.frame(names = row.names(na), na)
row.names(na)=NULL
na$names<-as.character(na$names)
i = 1
for (i in 1:nrow(na)){
  if(na$na[i] > 18.9){
    BALB<-drop.markers(BALB, na$names[i])
  } else{NULL}
}

#Counting number of crossover events. identifying and removing outliers
nxo <- countXO(BALB)
plot(nxo, ylab="# crossovers")
nxo[nxo<10]
#nxo[nxo>30]

BALB<-BALB[ , -c(291:293)]
BALB<-BALB[ , -c(259:261)]
BALB<-BALB[ , -c(46)]
BALB<-BALB[ , -c(34:40)]

nxo <- countXO(BALB)
plot(nxo, ylab="# crossovers")

howmany<-markernames(BALB)
length(howmany)

#genetic map (cM)
map<-pull.map(BALB)
par(cex=2.2, font.lab=2,mar=c(5,4,4,5), lwd=2)
plot(map,alternate.chrid = TRUE)

#calculating probability of genotypes
BALB<-calc.genoprob(BALB, step=1)


#----QTL Model: Sex as an additive Covariate, Analysis for Whole Brain OMOR---------
Sex<-as.matrix(BALB$pheno$Sex)
colnames(Sex)<- 'Sex'
BALB.SexAdd<-scanone(BALB, pheno.col=c(11:54), addcovar= Sex, method="hk")
BALB.SexAdd.perm<-scanone(BALB, pheno.col=c(11:54), addcovar= Sex, n.perm=1000, method="hk")
#Summary Table of LOD and p values
BALB.SexAdd.summary05<-summary(BALB.SexAdd, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.SexAdd.perm)

#-----Sex Specific, No Covatiate models for Whole brain OMOR-----
BALB.F<-subset(BALB,ind=BALB$pheno$Sex=="0")
BALB.M<-subset(BALB,ind=BALB$pheno$Sex=="1")

BALB.OXY<-subset(BALB,ind=BALB$pheno$Tx=="1")
BALB.SAL<-subset(BALB,ind=BALB$pheno$Tx=="0")
names(BALB$pheno)

BALB.NoCov.F<-scanone(BALB.F, pheno.col=c(11:54), method="hk")
BALB.NoCov.F.perm<-scanone(BALB.F, pheno.col=c(11:54), n.perm=1000, method="hk")
#Summary Table of Female LOD scores and p values
BALB.NoCov.F.summary05<-summary(BALB.NoCov.F, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.NoCov.F.perm)

BALB.NoCov.M<-scanone(BALB.M, pheno.col=c(11:54), method="hk")
BALB.NoCov.M.perm<-scanone(BALB.M, pheno.col=c(11:54), n.perm=1000, method="hk")

#----JPET OMOR Figure 3----
#Figure 3A
par(cex=1.9,cex.axis=1.25, cex.lab= 1.25, ann=TRUE,tck=NA,xaxt="s",yaxt="s")
plot(BALB.SexAdd, lodcolumn=c(37, 38, 39), ylim = c(0,8), alternate.chrid=TRUE,
     col=c("slateblue", "maroon","black"),lwd=2,ylab="LOD score",  main="Whole Brain [OXY], [NOR], and [OMOR]")

add.threshold(BALB.SexAdd, perms=BALB.SexAdd.perm, alpha=0.05, lodcolumn = c(37,38,39), 
              col = c("slateblue", "maroon","black"))
legend("topright",
       c( "[OXY]", "[NOR]","[OMOR]"), lwd = 2,# inset=c(-0.21,0),
       col = c("slateblue", "maroon","black"),xpd=TRUE,bty="n")



#Figure 3B
#Females
par(cex=1.9,cex.axis=1.25, cex.lab= 1.25, ann=TRUE,tck=NA,xaxt="s",yaxt="s")
plot(BALB.NoCov.F, lodcolumn=c(39), ylim = c(0,8), alternate.chrid=TRUE,
     col=c("red3"),lwd=2,ylab="LOD score",  main="Whole Brain [OMOR]: Chr 15", chr=15)

add.threshold(BALB.NoCov.F, perms=BALB.NoCov.F.perm, alpha=0.05, lodcolumn = c(39), 
              col = c("red3"))
#Males
par(cex=1.9,cex.axis=1.25, cex.lab= 1.25, ann=TRUE,tck=NA,xaxt="s",yaxt="s")
plot(BALB.NoCov.M, lodcolumn=c(39), ylim = c(0,8), alternate.chrid=TRUE,
     col=c("blue3"),lwd=2,ylab="LOD score",  main="Whole Brain Oxymorphone", chr=15, add = TRUE)

add.threshold(BALB.NoCov.M, perms=BALB.NoCov.M.perm, alpha=0.05, lodcolumn = c(39), 
              col = c("blue3"))
#Collapsed
plot(BALB.SexAdd, lodcolumn=c(39), ylim = c(0,8), alternate.chrid=TRUE,
     col=c("black"),lwd=2,ylab="LOD score",  main="Whole Brain Oxymorphone", chr=15, add = TRUE)

add.threshold(BALB.SexAdd, perms=BALB.SexAdd.perm, alpha=0.05, lodcolumn = c(39), 
              col = c("black"))
par(cex=1.5)
legend("topleft",
       c( "All", "Female","Male"), lwd = 2,# inset=c(-0.21,0),
       col = c("black", "red", "blue"),xpd=TRUE,bty="n")


#Figure 3C
#effect plot by sex
find.marker(BALB, chr=15, pos=24)

par(cex=1.8,ann=TRUE,tck=NA,xaxt="s",yaxt="s")
effectplot(BALB, mname1="rs264203947", pheno.col=49, col="Black", ylim = c(2,8), ylab="Average Latency",xlab="Genotype",
           main= "Hotplate BL, Chr 13 @ 31.938cM",  geno1=c("cJ/cJ", "Het","ByJ/ByJ"))
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.F, mname1="rs264203947", pheno.col=49,  ylim = c(2,8),  col="Red", ylab="ng OMOR / g Brain",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.M, mname1="rs264203947", pheno.col=49, col="Blue", ylim = c(2,8), ylab="ng OMOR / g Brain",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
legend("topright",
       c("All", "Females", "Males"), lwd = 2,# inset=c(-0.21,0),
       col = c("Black", "red","blue"),xpd=TRUE,bty="n")

#----JPET OMOR Table 1----
find.marker(BALB, chr=15, pos=24)

#Summary Table of Sex Collapsed LOD scores and p values
BALB.SexAdd.summary05<-summary(BALB.SexAdd, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.SexAdd.perm)
#Summary Table of Female LOD scores and p values
BALB.NoCov.F.summary05<-summary(BALB.NoCov.F, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.NoCov.F.perm)

#95% Bayesian Interval and 1.5 LOD drop interval, Sex Collapsed
bayesint(BALB.SexAdd, chr=15, 0.95, lodcolumn=39, expandtomarkers=FALSE)
lodint(BALB.SexAdd, chr=15, drop = 1.5, lodcolumn=39, expandtomarkers=FALSE)
#95% Bayesian Interval and 1.5 LOD drop interval, Female Only
bayesint(BALB.NoCov.F, chr=15, 0.95, lodcolumn=39, expandtomarkers=FALSE)
lodint(BALB.NoCov.F, chr=15, drop = 1.5, lodcolumn=39, expandtomarkers=FALSE)

#Sex Collapsed Variance
all_15<-makeqtl(BALB, 15, 29.95, what='prob')
summary(fitqtl(BALB, pheno.col=49, all_15, covar=Sex, method='hk', formula=y~Q1+Sex,
               get.ests=TRUE, run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE))
#Female Specific Variance
female_15<-makeqtl(BALB.F, 15, 29.95, what='prob')
summary(fitqtl(BALB.F, pheno.col=49, all_15, method='hk', formula=y~Q1,
               get.ests=TRUE, run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE))




#-----QTL Model: Sex and Tx as Additive Covatiate, Analysis for Hot Plate----------------
Sex<-as.numeric(BALB$pheno$Sex)
Tx<-BALB$pheno$Tx
SexTx<-cbind(Sex, Tx)
colnames(SexTx)<- c('Sex', 'Tx')

BALB.SexTxAdd<-scanone(BALB, pheno.col=c(11:56), addcovar= SexTx, method="hk")
BALB.SexTxAdd.perm<-scanone(BALB, pheno.col=c(11:56), addcovar= SexTx, n.perm=1000, method="hk")

BALB.SexTxAdd.summary05<-summary(BALB.SexTxAdd, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.SexTxAdd.perm)

#-----QTL Model: Sex, Tx, Age as Additive Covatiates, Analysis for Brain Weight----------------
Sex<-as.numeric(BALB$pheno$Sex)
Age<-as.numeric(BALB$pheno$D1.Age)
Tx<-as.numeric(BALB$pheno$Tx)
SexAgeTx<-cbind(Sex,Age, Tx)
colnames(SexAgeTx)<- c('Sex', 'Age','Tx')

names(BALB$pheno)

BALB.SexAgeTxAdd<-scanone(BALB, pheno.col=c(11:56), addcovar= SexAgeTx, method="hk")
BALB.SexAgeTxAdd.perm<-scanone(BALB, pheno.col=c(11:56), addcovar= SexAgeTx, n.perm=1000, method="hk")

BALB.SexAgeTxAdd.summary05<-summary(BALB.SexAgeTxAdd, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.SexTxAdd.perm)


#----Molecular Pain Figure 4----
#Figure 4A
par(cex=2, cex.main= 1.5, cex.lab=1.3, cex.axis=1.3, ann=TRUE,tck=NA,xaxt="s",yaxt="s")
plot(BALB.SexTxAdd, lodcolumn=c(46), ylim = c(0,12), alternate.chrid=TRUE,
     col=c("black"),lwd=2,ylab="LOD score",  main="Hot Plate BL Latency")
add.threshold(BALB.SexTxAdd, perms=BALB.SexTxAdd.perm, alpha=0.05, lodcolumn = c(46), 
              col = c("black"))


#Figure 4B
par(cex=2, cex.main= 1.5, cex.lab=1.3, cex.axis=1.3, ann=TRUE,tck=NA,xaxt="s",yaxt="s")
plot(BALB.SexTxAdd, lodcolumn=c(46), ylim = c(0,12), alternate.chrid=TRUE, chr = 13,
     col=c("black"),lwd=2,ylab="LOD score",  main="Hot Plate BL Latency")
add.threshold(BALB.SexTxAdd, perms=BALB.SexTxAdd.perm, alpha=0.05, lodcolumn = c(46), 
              col = c("black"))

#Figure 4C
find.marker(BALB, chr=13, pos=34.098)

par(cex=1.8,ann=TRUE,tck=NA,xaxt="s",yaxt="s")
effectplot(BALB, mname1="SBB132408413", pheno.col=55, col="Black", ylim = c(15,40), ylab="Average Latency",xlab="Genotype",
           main= "Hotplate BL, Chr 13 @ 31.938cM",  geno1=c("cJ/cJ", "Het","ByJ/ByJ"))
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.F, mname1="SBB132408413", pheno.col=55,  ylim = c(15,40),  col="Red", ylab="ng OMOR / g Brain",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.M, mname1="SBB132408413", pheno.col=55, col="Blue", ylim = c(15,40), ylab="ng OMOR / g Brain",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
legend("topright",
       c("All", "Females", "Males"), lwd = 2,# inset=c(-0.21,0),
       col = c("Black", "red","blue"),xpd=TRUE,bty="n")

#Figure 4D
par(cex=1.8,ann=TRUE,tck=NA,xaxt="s",yaxt="s")
effectplot(BALB, mname1="SBB132408413", pheno.col=56, col="Black", ylim = c(-1,1), ylab="Normalized Latency",xlab="Genotype",
           main= "Hotplate BL, Chr 13 @ 31.938cM",  geno1=c("cJ/cJ", "Het","ByJ/ByJ"))
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.F, mname1="SBB132408413", pheno.col=56,  ylim = c(-1,1),  col="Red", ylab="ng OMOR / g Brain",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.M, mname1="SBB132408413", pheno.col=56, col="Blue", ylim = c(-1,1), ylab="ng OMOR / g Brain",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
legend("topright",
       c("All", "Females", "Males"), lwd = 2,# inset=c(-0.21,0),
       col = c("Black", "red","blue"),xpd=TRUE,bty="n")

#----Molecular Pain Figure 5----
#Figure 5A
par(cex=2, cex.main= 1.5, cex.lab=1.3, cex.axis=1.3, ann=TRUE,tck=NA,xaxt="s",yaxt="s")
plot(BALB.SexAgeTxAdd, lodcolumn=c(36), ylim = c(0,10), alternate.chrid=TRUE,
     col=c("black"),lwd=2,ylab="LOD score",  main="Brain Weight")
add.threshold(BALB.SexAgeTxAdd, perms=BALB.SexAgeTxAdd.perm, alpha=0.05, lodcolumn = c(36), 
              col = c("black"))


#Figure 5B
par(cex=2, cex.main= 1.5, cex.lab=1.3, cex.axis=1.3, ann=TRUE,tck=NA,xaxt="s",yaxt="s")
plot(BALB.SexAgeTxAdd, lodcolumn=c(35), ylim = c(0,10), alternate.chrid=TRUE, ch = 5,
     col=c("black"),lwd=2,ylab="LOD score",  main="Brain Weight: Chr 5")
add.threshold(BALB.SexAgeTxAdd, perms=BALB.SexAgeTxAdd.perm, alpha=0.05, lodcolumn = c(36), 
              col = c("black"))


#Figure 5C
find.marker(BALB, chr=5, pos=58.6)

par(cex=1.7,ann=TRUE,tck=NA,xaxt="s",yaxt="s")
effectplot(BALB, mname1="SBJ054695332", pheno.col=45, col="Black", ylim = c(0.4,0.48), ylab="Grams",xlab="Genotype",
           main= "Brain Weight, Chr 5 @ 56.95cM",  geno1=c("cJ/cJ", "Het","ByJ/ByJ"))
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.F, mname1="SBJ054695332", pheno.col=45,  ylim = c(0.4,0.48),  col="Red", ylab="",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.M, mname1="SBJ054695332", pheno.col=45, col="Blue", ylim = c(0.4,0.48), ylab="",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
legend("topright",
       c("All", "Females", "Males"), lwd = 2,# inset=c(-0.21,0),
       col = c("Black", "red","blue"),xpd=TRUE,bty="n")


#Figure 5D
par(cex=1.7,ann=TRUE,tck=NA,xaxt="s",yaxt="s")
effectplot(BALB, mname1="SBJ054695332", pheno.col=46, col="Black", ylim = c(-1.5, 1.5), ylab="Normalized Brian Weight",xlab="Genotype",
           main= "Brain Weight, Chr 5 @ 56.95cM",  geno1=c("cJ/cJ", "Het","ByJ/ByJ"))
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.F, mname1="SBJ054695332", pheno.col=46,  ylim = c(-1.5, 1.5),  col="Red", ylab="",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(BALB.M, mname1="SBJ054695332", pheno.col=46, col="Blue", ylim = c(-1.5, 1.5), ylab="",xlab="Genotype",
           main= "Hotplate Baseline Average, Chr 13 @ 34.098cM")
legend("topright",
       c("All", "Females", "Males"), lwd = 2,# inset=c(-0.21,0),
       col = c("Black", "red","blue"),xpd=TRUE,bty="n")
.

#----Molecular Pain Table 1 - QTL metrics----
#HotPlate
#LOD Scores and Pvals
BALB.SexTxAdd.summary05<-summary(BALB.SexTxAdd, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.SexTxAdd.perm)
#Confidence Intervals
bayesint(BALB.SexTxAdd, chr=13, 0.95, lodcolumn=46, expandtomarkers=FALSE)
lodint(BALB.SexTxAdd, chr=13, drop = 1.5, lodcolumn=46, expandtomarkers=FALSE)
#%Variance Explained
all_13<-makeqtl(BALB, 13, 31.1, what='prob')
summary(fitqtl(BALB, pheno.col=56, all_13, covar=SexTx, method='hk', formula=y~Q1+Sex+Tx,
               get.ests=TRUE, run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE))


#Brain Weight
#LOD Scores and Pvals
BALB.SexAgeTxAdd.summary05<-summary(BALB.SexAgeTxAdd, alpha=0.05, format='allpheno', pvalues=TRUE, perms=BALB.SexTxAdd.perm)
#Confidence Intervals
bayesint(BALB.SexAgeTxAdd, chr=5, 0.95, lodcolumn=36, expandtomarkers=FALSE)
lodint(BALB.SexAgeTxAdd, chr=5, drop = 1.5, lodcolumn=36, expandtomarkers=FALSE)
#%Variance Explained
all_5<-makeqtl(BALB, 5, 58.6, what='prob')
summary(fitqtl(BALB, pheno.col=46, all_5, covar=SexAgeTx, method='hk', formula=y~Q1+Sex+Tx+Age,
               get.ests=TRUE, run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE))