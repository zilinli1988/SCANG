## ------------------------------------------------------------------------
library("Matrix")
library(SCANG)
genotype <- example$genotype
maf <- example$maf
snploc <- example$snploc
sum((maf>0)*(maf<0.05))

## ------------------------------------------------------------------------
## kinship
library(kinship2)
library(MASS)
grid <- 1
Npercell <- 10000
ndiv <- 1
vfam <- 0.5
N <- round(grid*grid*Npercell/ndiv)

unitmat <- matrix(0.5, 4, 4)
diag(unitmat) <- 1
unitmat[1,2] <- unitmat[2,1] <- 0
ped <- data.frame(famid = rep(as.integer(1:2500), each=4), id = as.integer(1:10000L), fa = rep(0, 10000), mo = rep(0, 10000))
for(i in 1:2500) {
        ped$fa[4*i-(0:1)] <- ped$id[4*i-3]
	ped$mo[4*i-(0:1)] <- ped$id[4*i-2]
}
kins <- makekinship(ped$famid, ped$id, ped$fa, ped$mo)

## ------------------------------------------------------------------------
set.seed(168)
samplesize <- dim(genotype)[1]
X1 <- rnorm(samplesize,0,1)
X2 <- rbinom(samplesize,1,1/2)
epi <- rnorm(samplesize,0,1)

## ------------------------------------------------------------------------
## number of signal regions
n0 <- 2

## generate signal region location
pp <- floor(dim(genotype)[2]/n0)

## location od signal region
# In snese of variants order
sigloc <- matrix(rep(0,3*n0),ncol=3)
# In sense of location
sigloc_bp <- matrix(rep(0,3*n0),ncol=3)
for(i in 1:n0)
{
  begnum <- (i-1)*pp + 1
  endnum <- i*pp - 1000
  sigloc[i,1] <- sample(begnum:endnum,1)
  sigloc_bp[i,1] <- snploc$CHROM_POS[sigloc[i,1]]
  length_r <- runif(1,0,1)
  sigloc_bp[i,3] <- (length_r<=1/4)*3000 + (length_r>1/4)*(length_r<=2/4)*4000  + (length_r>2/4)*(length_r<=3/4)*5000 + (length_r>3/4)*6000
  sigloc_bp[i,2] <- sigloc_bp[i,1] + sigloc_bp[i,3] - 1
  sigloc[i,2] <- which.max(snploc$CHROM_POS>sigloc_bp[i,2]) - 1
  if(sigloc[i,2] == 0)
  {
    sigloc[i,2] <- dim(genotype)[2]
  }

  sigloc[i,3] <- sigloc[i,2] - sigloc[i,1] + 1
}

## ------------------------------------------------------------------------
sigloc

## ------------------------------------------------------------------------
sigloc_bp

## ------------------------------------------------------------------------
percen <-  0.10

sigploc <- matrix(rep(0,28*n0),ncol=28)
sigloctemp <- c()

for(ii in 1:n0)
{
  maftemp <- maf[(sigloc[ii,1]+1):(sigloc[ii,2]-1)]
  mafid <- (sigloc[ii,1]+1):(sigloc[ii,2]-1)
  mafid <- mafid[mafid>0]
  p0 <- floor(sigloc[ii,3]*percen)
  sigploc[ii,1:p0] <- c(sigloc[ii,1],sort(sample(mafid,p0-2,replace=FALSE)),sigloc[ii,2])
}

sigloc <- cbind(sigloc,sigploc)

## ------------------------------------------------------------------------
protect <- 0
c0 <- 0.18

## generate phenotype
phenotype <- 0.5*X1 + 0.5*X2
for(ii in 1:n0)
{
  sigloctemp <- sigloc[ii,4:(dim(sigloc)[2]-1)]
  sigloctemp <- sigloctemp[sigloctemp>0]
  # beta
  beta <- c0*log(maf[sigloctemp])/log(10)
  betadir <- rbinom(length(beta),1,protect)*2-1
  beta <- beta*betadir
  phenotype <- phenotype + genotype[,sigloctemp]%*%beta
  phenotype <- as.vector(phenotype)
}

phenotype <- phenotype + epi
X <- cbind(X1,X2)
phenotypedata <- data.frame(phenotype=phenotype,X=X)

## ------------------------------------------------------------------------
Lmax <- 250         
Lmin <- 70         
family <- "gaussian"
filter <- 2e-5
f <- 0

## ------------------------------------------------------------------------
# obj_nullmodel <- fit_null_glm_SCANG(phenotype~-1+X,data=phenotypedata,family=gaussian(link = "identity"))
# res_lm <- SCANG(genotype,obj_nullmodel,Lmin,Lmax,filter=filter,f=f)

## ------------------------------------------------------------------------
res_lm$SCANG_O_res

## ------------------------------------------------------------------------
sigloc[,1:20]

## ------------------------------------------------------------------------
res_lm$SCANG_S_res

## ------------------------------------------------------------------------
res_lm$SCANG_B_res

## ------------------------------------------------------------------------
protect <- 0
c0 <- 0.21
## generate phenotype
phenotype <- 0.5*X1 + 0.5*X2
for(ii in 1:n0)
{
  sigloctemp <- sigloc[ii,4:(dim(sigloc)[2]-1)]
  sigloctemp <- sigloctemp[sigloctemp>0]
  # beta
  beta <- c0*log(maf[sigloctemp])/log(10)
  betadir <- rbinom(length(beta),1,protect)*2-1
  beta <- beta*betadir
  phenotype <- phenotype + genotype[,sigloctemp]%*%beta
  phenotype <- as.vector(phenotype)
}

## random effect
randomfam <- mvrnorm(N/4, rep(0, 4), vfam * unitmat)
randomfam <- as.vector(t(randomfam))
id0 <- 1:N


phenotype <- phenotype + randomfam + epi
X <- cbind(X1,X2)

phenotypedata_relatedness <- data.frame(phenotype=phenotype,X=X,id=id0)

## ------------------------------------------------------------------------
# obj_nullmodel <- fit_null_glmmkin_SCANG(phenotype ~ -1+X, data=phenotypedata_relatedness, kins=kins, id="id",family=gaussian(link="identity"))

# res_lmm <- SCANG(genotype=genotype,obj_nullmodel=obj_nullmodel,Lmin=Lmin,Lmax=Lmax,filter=filter,f=f)

## ------------------------------------------------------------------------
res_lmm$SCANG_O_res

## ------------------------------------------------------------------------
res_lmm$SCANG_S_res

## ------------------------------------------------------------------------
res_lmm$SCANG_B_res

## ------------------------------------------------------------------------
kins <- as.matrix(kins)

## ------------------------------------------------------------------------
# obj_nullmodel <- fit_null_glmmkin_SCANG(phenotype ~ -1+X, data=phenotypedata_relatedness, kins=kins, id="id",family=gaussian(link="identity"),use_sparse = FALSE)
# res_lmm_dense <- SCANG(genotype,obj_nullmodel,Lmin,Lmax,filter=filter,f=f)

## ------------------------------------------------------------------------
res_lmm_dense$SCANG_O_res

## ------------------------------------------------------------------------
res_lmm_dense$SCANG_S_res

## ------------------------------------------------------------------------
res_lmm_dense$SCANG_B_res

