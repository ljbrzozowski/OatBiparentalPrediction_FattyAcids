#This code takes phenotype + genotype files from biparental families + panels
#This code asseses genomic prediciton accuracy within biparental families, one focal family at a time
#This code estimates marker effects using the germplasm panels to predict within focal family

## Part 1.Load libraries
########### 
library(BGLR)
library(rrBLUP)
############


## Part 2: Load input files
# all individuals must have phenotypes + genotypes
# phenotype names must match between files
# both panels and families must have same marker set
# order of markers and phenotypes must match between panels and families
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


### Read in phenotypes + genotypes from panels
#panel_phenotypes <- read.csv("./PathToPhenotypeFile.RDS")
#list of phenotype files where named by panel
#rownames are Lines (above); columns are phenotypes

#panel_genotypes <-  readRDS("./PathToGenotypeFile.RDS")
#list of genotype files where named by panel
#rownames are Lines (above); columns are marker names

###################


## Part 3: Make marker and relationship matrices for estimating marker effects + validation
# following https://github.com/MarcooLopez/Genomic-Selection/blob/master/single_environment.md
## get all marker and relationship matrices for set of markers in all panels
## get marker matrix for each family for validation
#########################

#panels
X_Panel1 <- panel_genotypes$Panel1[,which(colnames(panel_genotypes$Panel1) %in% allMrkrs)]
M_Panel1 <- scale(X_Panel1)
X_Panel2 <- panel_genotypes$Panel2[,which(colnames(panel_genotypes$Panel2) %in% allMrkrs)]
M_Panel2 <- scale(X_Panel2)

PanelMrkrMatrices <- list(X = list(DP = X_Panel1, EP=X_Panel2),
                          M = list(DP = M_Panel1, EP=M_Panel2))

rm(X_Panel1,M_Panel1, X_Panel2,M_Panel2)
#saveRDS(PanelMrkrMatrices, "./PathToSavePanelMarkerMatrices.RDS")

cat("\n saved panel marker matrices")


#families
all_fam_gbs <- all_fam_gbs[,which(colnames(all_fam_gbs) %in% allMrkrs)]

validFam_X <- list()

for (i in 1:length(fams)) {
  #marker matrix for single family to use in validation
  ## this is focal family
  focal_fam_temp <- family_genotypes[which(grepl(pattern = as.character(fams[i]), x = rownames(family_genotypes))),
                                     which(colnames(family_genotypes) %in% colnames(fam_tp_temp_M))]
  validFam_X[[as.character(fams[i])]] <- focal_fam_temp

  rm(focal_fam_temp)
  cat("\n Done with X matrices for ", as.character(fams[i]))
}

cat("\n X matrices made for families")

FamMrkMatrices <- list(validFam_X = validFam_X )

rm(validFam_X)
#saveRDS(FamMrkMatrices, "./PathToSaveMarkerMatrices.RDS")

cat("\n X matrices saved for families")


################


## Part 4: Estimate marker effects
# following https://github.com/MarcooLopez/Genomic-Selection/blob/master/single_environment.md
# This code uses GBLUP + Bayes B
#################
cat("\n Starting mrkr effect estimates from panels")

mkrMeans <- list(GBLUP = list(), BB = list())
mkrEffs <- list(GBLUP = list(), BB = list())

nIter <- 20000
burnIn <-5000

#panels
for (j in 1:length(phenos)) { #for each penotype
  for(i in 1:length(PanelMrkrMatrices$X)) { #for each germplasm panel
    
    y=panel_phenotypes[[i]][,which(colnames(panel_phenotypes[[i]]) %in% phenos[j])]
    Z = PanelMrkrMatrices$X[[i]]
    X = PanelMrkrMatrices$M[[i]]
    
    # G-BLUP model using rrBLUP
    fm_rrblup  <- mixed.solve(y=y , Z=Z, SE=FALSE, return.Hinv = F)
    mkrMeans$GBLUP[[as.character(names(PanelMrkrMatrices$X)[i])]][[as.character(phenos[j])]] <- fm_rrblup$beta
    mkrEffs$GBLUP[[as.character(names(PanelMrkrMatrices$X)[i])]][[as.character(phenos[j])]] <- as.matrix(fm_rrblup$u)
    
    # Bayes B model using 'BGLR' package
    fm_bb <- BGLR(y = y,ETA=list(list(X=X,model="BayesB")),nIter=nIter,burnIn=burnIn)
    mkrMeans$BB[[as.character(names(PanelMrkrMatrices$X)[i])]][[as.character(phenos[j])]] <- fm_bb$mu
    mkrEffs$BB[[as.character(names(PanelMrkrMatrices$X)[i])]][[as.character(phenos[j])]] <- as.matrix(fm_bb$ETA[[1]]$b)
    
    rm(fm_rrblup,fm_bb)
    
    Mrk_eff_for_GP_panel = list(mkrMeans = mkrMeans, mkrEffs = mkrEffs)
    #saveRDS(Mrk_eff_for_GP_panel, "./PathToSavePanelMrkEffects.RDS")
    
  }
  
  cat("\n Done with marker effect est for", as.character(phenos[j]), "in panels")
  
}

###################


## Part 5: Correlation between predicted and observed
#################

## predict family
corRes <- as.data.frame(matrix(ncol=5)) ; colnames(corRes) <- c("FocalFam", "Model", "TrainingPop", "pheno", "r")


for (i in 1:length(FamMrkMatrices$validFam_X)) { #i = focal fam to test
  X_FocalFam <- FamMrkMatrices$validFam_X[[i]]
  dim(X_FocalFam)
  
  for (j in 1:length( Mrk_eff_for_GP_panel$mkrEffs)) { #j = model
    
    for (k in 1:length( Mrk_eff_for_GP_panel$mkrEffs[[j]])) { #k=panel
      
      for (m in 1:length(phenos)) { #m = phenos
        
      TmpB <- as.matrix(Mrk_eff_for_GP_panel$mkrEffs[[j]][[k]][[m]])
      
      X_FocalFam <- X_FocalFam[,which(colnames(X_FocalFam) %in% rownames(TmpB))]
      TmpB <- as.matrix(TmpB[which(rownames(TmpB) %in% colnames(X_FocalFam)),])
      gHat<-  X_FocalFam  %*% TmpB
      
      y=bipTotFAMEs[which(bipTotFAMEs$newLine %in% rownames(X_FocalFam)),]
      y<- y %>% arrange(factor(y$newLine, levels=rownames(gHat)))
      cor.test(gHat, y[,which(colnames(y) %in% phenos[m])])
      
      corResTemp <- as.data.frame(matrix(ncol=5)) ; colnames(corResTemp) <- c("FocalFam", "Model", "TrainingPop", "pheno", "r")
      
      corResTemp[1,] <- c(as.character(names(FamMrkMatrices$validFam_X)[i]),
                          as.character(names( Mrk_eff_for_GP_panel$mkrEffs)[j]),
                          as.character(names( Mrk_eff_for_GP_panel$mkrEffs[[j]])[k]),
                          as.character(names( Mrk_eff_for_GP_panel$mkrEffs[[j]][[k]])[m]),
                          cor.test(gHat, y[,which(colnames(y) %in% phenos[m])])$estimate)
      
      corRes <- rbind(corRes, corResTemp)
      } 
    } 
  }
} 

corRes
corRes$r <- as.numeric(corRes$r)
#write.csv(corRes, "../PathToPredAccResults.csv, row.names=F")

#######