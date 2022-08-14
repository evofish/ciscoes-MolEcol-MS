##script for estimating the molecular diversity indeces

library(vcfR)
library(adegenet)
library(pcadapt)
library(ggplot2)
library(pegas)
library (ape)
library(hierfstat)


setwd('~/Desktop/coregonids-transcriptomes/cisco-snps/')
cisco<-read.vcfR('cisco-snps-filter3-thin-NOCZ8.recode.vcf') 

cisco.gen<-vcfR2genind(cisco)
cisco.gen

#including populations in genind
pops<-c(rep("Zen",7), rep("Kiy",8), rep("Hoy",8), rep("Art",8))
cisco.gen@pop<-as.factor(rep(pops))
cisco.gen@pop

cisco.gen

sum<-summary(cisco.gen) #very large file, confirms pops are there 
head(sum, n=5)

base<-basic.stats(cisco.gen) #gives indexes for whole dataset
#Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.2695  0.2245  0.2274  0.0030  0.2284  0.0039  0.0130  0.0173 -0.2006  0.0051 

boot.vc(cisco.h[,1],cisco.h[,-c(1:2)]) #gives actual confidence intervals for Fst and Fis, report the average of all measures
#H-Total F-Pop/Total F-Ind/Total H-Pop F-Ind/Pop   Hobs
#2.5%   0.2661      0.0139     -0.2527 0.262   -0.2720 0.3283
#50%    0.2702      0.0155     -0.2399 0.266   -0.2594 0.3352
#97.5%  0.2743      0.0174     -0.2285 0.270   -0.2481 0.3419

base2<-as.matrix(base.pop) #allows exporting as table
write.table(base2, file="diversity-indexes-cisco-May19.txt", quote=FALSE, sep='\t')

##analyses with HierFstat in R

diff<-fstat(cisco.gen)
diff
#pop        Ind
#Total 0.01753703 -0.1795115
#pop   0.00000000 -0.2005659

gtest<- gstat.randtest(cisco.gen,nsim=99) #should be increased to 500 - This G-statistic Monte Carlo procedure tests for population structuring at different levels.
gtest
#Monte-Carlo test
#Call: gstat.randtest(x = cisco.gen, nsim = 99)
#Observation: 24222.84 f
#Based on 99 replicates
#Simulated p-value: 0.01 
#Alternative hypothesis: greater 
#Std.Obs Expectation    Variance 
#13.45813 20141.39456 91972.42984 

cisco.h<-genind2hierfstat(cisco.gen) #converts from geneind to hierfstat

#nei's distance by sample 
genet.dist(cisco.gen)

freqs<-pop.freq(cisco.h) #population allele frequences by each locus

#this is the conversion to genpop file
cisco.pop<- genind2genpop(cisco.gen, pop = c(rep("Zen",7), rep("Kiy",8), rep("Hoy",8), rep("Art",8)))
popNames(cisco.pop)

Hs(cisco.pop) #heterozygosity per population 
#Zen       Kiy       Hoy       Art 
#0.2351241 0.2181278 0.2325045 0.2165995 
Ht(cisco.pop)

#allele depth graph
depth<-extract.gt(cisco, element='AD', as.numeric = TRUE) #extracting the depth of each marker
par(mar=c(1,1,1,1))
boxplot(depth, las=3, col=c("#C0C0C0", "#808080"), ylab="Allele Depth", las=2)

#Qual histogram 
plot(cisco.sub)

##PCadapt
    #PCA
cisco<-read.pcadapt('cisco-snps-overlapping-orthos_NoCZ08-may17-2019.vcf', type='vcf')
x <- pcadapt(cisco, K = 10) 
plot(x, option='screeplot', K=10) #suggests best k=3
pop.cisco<-c(rep("Zen",7), rep("Kiy",8), rep("Hoy",8), rep("Art",8))
print(pop.cisco)
plot(x, option = "scores", pop = pop.cisco, label=TRUE, label.size = 4) #PC1=21.02% PC2=19.98
x #type x to get the values of the first 2 PC scores. Copy and paste and export as a table - draw PCA along 2 other datasets

###for finding outliers 
x <- pcadapt(cisco, K =3) #based on the result above
summary(x)
plot(x , option = "manhattan", k=3, threshold=0.10)
plot(x, option = "qqplot")

y <-outlier(x,K=3,threshold=0.1)

##Adegenet DAPC and Fstat
cisco<-read.vcfR('cisco-snps-overlapping-orthos_NoCZ08-may17-2019.vcf')

cisco.gen<-vcfR2genlight(cisco) #converts the matrix to Genlight format for adegenet
pop(cisco.gen)<- c(rep("Zen",7), rep("Kiy",8), rep("Hoy",8), rep("Art",8)) #stores population information in Adegenet
x.dist<-dist(cisco.gen)

#dapc 
grp <- find.clusters(cisco.gen, max.n.clust=10) #retained 10 PCs and K kept at 4
dapc1<-dapc(cisco.gen, grp$grp) # retained 10PCs , Chose 3DF 
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:3), col=c("red","blue", "gold", "green")) #graph test1 dapc
compoplot(dapc1, col=c("red","blue", "gold", "green"),lab="", txt.leg=paste("group", 1:4), ncol=4) #test of composition plot similar to structure 

