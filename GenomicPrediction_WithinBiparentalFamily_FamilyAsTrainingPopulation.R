#This code takes phenotype + genotype files from biparental families
#This code asseses genomic prediciton accuracy within biparental families, one focal family at a time
#This code estimates marker effects using all families except focal family to predict within focal family

## Part 1.Load libraries
########### 
library(BGLR)
library(rrBLUP)
############


## Part 2: Load input files
#################

### Read in phenotypes + genotypes from families 

#family_phenotypes <- read.csv("./PathToPhenotypeFile.csv")
#the first column of this file is "Line" which is the unique genotype identifier
#the second column of this file is "Family" which groups genotypes into families 
# *family name **MUST** be included in genotype identifer and can be retreived through e.g., stringr or grepl; grepl used in this code
#the third through nth column are names of phenotypes to predict

#family_genotypes <-  readRDS("./PathToGenotypeFile.RDS")
#rownames are Lines (above); columns are marker names

### Ensure that same individuals in phenotypes + genotypes from families 

#family_phenotypes<- family_phenotypes[which(family_phenotypes$Line %in% rownames(family_genotypes)),]
#family_genotypes<- family_genotypes[which(rownames(family_genotypes) %in% family_phenotypes$Line),]
#cat("\n Family phenotype data loaded and formatted with dimensions", dim(family_phenotypes) )
#cat("\n Family genotype data loaded and formatted with dimensions", dim(family_genotypes) )

### get vector of phenotype and family names

#phenos <- colnames(family_phenotypes)[2:ncol(family_phenotypes)]
#cat("\n There are", length(phenos), "phenotypes")

#fams  <- unique(family_phenotypes$Family)
#cat("\n There are", length(fams), "families")

###################


## Part 3: Make marker and relationship matrices for estimating marker effects + validation
# following https://github.com/MarcooLopez/Genomic-Selection/blob/master/single_environment.md
#################
validFam_X <- list()
TPFam_X <- list()
TPFam_M <- list()

for (i in 1:length(fams)) { #cycle through families
  # get M, G matrices for training populations to estimated snp effect for rest of families
  fam_tp_temp <- family_genotypes[-c(which(grepl(pattern = as.character(fams[i]), x = rownames(family_genotypes)))),]
  fam_tp_temp <- fam_tp_temp[,-c(which(apply(fam_tp_temp,FUN =  var, 2) == 0))] #remove zero var so can scale
  fam_tp_temp_M <- scale(fam_tp_temp) #when column has zero variance, this fails
  
  TPFam_X[[as.character(fams[i])]] <- fam_tp_temp
  TPFam_M[[as.character(fams[i])]] <- fam_tp_temp_M
  
  #marker matrix for single family to use in validation
  focal_fam_temp <- family_genotypes[which(grepl(pattern = as.character(fams[i]), x = rownames(family_genotypes))),
                                which(colnames(family_genotypes) %in% colnames(fam_tp_temp_M))]
  validFam_X[[as.character(fams[i])]] <- focal_fam_temp
  
  rm(fam_tp_temp, fam_tp_temp_M, focal_fam_temp)
  
  cat("\n Done with X M matrices for ", as.character(fams[i]))
}

cat("\n X, M matrices made for families")

FamMrkMatrices <- list(validFam_X = validFam_X,
                       TPFam_X = TPFam_X,
                       TPFam_M = TPFam_M)

#saveRDS(FamMrkMatrices, "./PathToSaveMarkerMatrices.RDS")
rm(validFam_X,TPFam_X, TPFam_M)

cat("\n X, M  matrices saved for families")


#########################


## Part 4: Estimate marker effects
# following https://github.com/MarcooLopez/Genomic-Selection/blob/master/single_environment.md
# This code uses GBLUP + Bayes B
#################
cat("\n Starting mrkr effect estimates from families")

mkrMeans <- list(GBLUP = list(), BB = list())
mkrEffs <- list(GBLUP = list(), BB = list())

nIter <- 20000
burnIn <-5000

## estimate marker effects without effect of family
for (j in 3:ncol(family_phenotypes)) { #phenos
  for(i in 1:length(FamMrkMatrices$TPFam_X)) { #for each training population
    Z = FamMrkMatrices$TPFam_X[[i]]
    X = FamMrkMatrices$TPFam_M[[i]]
    
    y=family_phenotypes[which(family_phenotypes$Line %in% rownames(Z)),]
    y<- y %>% arrange(factor(Line, levels=rownames(Z)))
    y <- y[,j]
    
    # G-BLUP model using rrBLUP
    fm_rrblup  <- mixed.solve(y=y , Z=Z, SE=FALSE, return.Hinv = F)
    mkrMeans$GBLUP[[as.character(names(FamMrkMatrices$TPFam_X)[i])]][[as.character(colnames(family_phenotypes)[j])]] <- fm_rrblup$beta
    mkrEffs$GBLUP[[as.character(names(FamMrkMatrices$TPFam_X)[i])]][[as.character(colnames(family_phenotypes)[j])]] <- as.matrix(fm_rrblup$u)
    
    # Bayes B model using 'BGLR' package
    fm_bb <- BGLR(y = y,ETA=list(list(X=X,model="BayesB")),nIter=nIter,burnIn=burnIn, verbose = F)
    mkrMeans$BB[[as.character(names(FamMrkMatrices$TPFam_X)[i])]][[as.character(colnames(family_phenotypes)[j])]] <- fm_bb$mu
    mkrEffs$BB[[as.character(names(FamMrkMatrices$TPFam_X)[i])]][[as.character(colnames(family_phenotypes)[j])]] <- as.matrix(fm_bb$ETA[[1]]$b)
    
    rm(fm_rrblup,fm_bb)
    
    cat("\n Done with marker eff est for fam:", as.character(names(FamMrkMatrices$TPFam_X)[i]), "\n y:", colnames(family_phenotypes)[j])
    
    Mrk_eff_for_GP = list(mkrMeans = mkrMeans, mkrEffs = mkrEffs)
    #saveRDS(Mrk_eff_for_GP, "./PathToSaveMarkerEffects.RDS")
    
  }
  cat("\n Done with marker effect est for all families for y:", colnames(family_phenotypes)[j])
}



#########################


## Part 5: Correlation between predicted and observed
#################

## predict family
corRes <- as.data.frame(matrix(ncol=5)) ; colnames(corRes) <- c("FocalFam", "Model", "TrainingPop", "pheno", "r")


### cycle through model, training population, focal families, phenotypes
for (i in 1:length(FamMrkMatrices$validFam_X)) { #i = focal fam to test
  X_FocalFam <- FamMrkMatrices$validFam_X[[i]]
  dim(X_FocalFam)
  
  for (j in 1:length( Mrk_eff_for_GP$mkrEffs)) { #j = model
    
    for (k in 1:length( Mrk_eff_for_GP$mkrEffs[[j]])) { #k=panel
      
      for (m in 1:length(phenos)) { #m = phenos
        
        TmpB <- as.matrix(Mrk_eff_for_GP$mkrEffs[[j]][[k]][[m]])
        
        X_FocalFam <- X_FocalFam[,which(colnames(X_FocalFam) %in% rownames(TmpB))]
        TmpB <- as.matrix(TmpB[which(rownames(TmpB) %in% colnames(X_FocalFam)),])
        gHat<-  X_FocalFam  %*% TmpB
        
        y=family_phenotypes[which(family_phenotypes$Line %in% rownames(X_FocalFam)),]
        y<- y %>% arrange(factor(y$Line, levels=rownames(gHat)))
        cor.test(gHat, y[,which(colnames(y) %in% phenos[m])])
        
        corResTemp <- as.data.frame(matrix(ncol=5)) ; colnames(corResTemp) <- c("FocalFam", "Model", "TrainingPop", "pheno", "r")
        
        corResTemp[1,] <- c(as.character(names(FamMrkMatrices$validFam_X)[i]),
                            as.character(names( Mrk_eff_for_GP$mkrEffs)[j]),
                            as.character(names( Mrk_eff_for_GP$mkrEffs[[j]])[k]),
                            as.character(names( Mrk_eff_for_GP$mkrEffs[[j]][[k]])[m]),
                            cor.test(gHat, y[,which(colnames(y) %in% phenos[m])])$estimate)
        
        corRes <- rbind(corRes, corResTemp)
      } 
    } 
  }
  corRes$r <- as.numeric(corRes$r)
} 


corRes$r <- as.numeric(corRes$r)
# ** important ** 
corRes <- corRes[corRes$TrainingPop == corRes$FocalFam,] ## training pop for specific family has same name as the focal fam
#write.csv(corRes, "./PathToPredictionAccResults.csv", row.names = F)


