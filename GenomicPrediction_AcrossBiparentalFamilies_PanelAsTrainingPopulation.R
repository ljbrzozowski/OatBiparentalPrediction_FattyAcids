#This code takes phenotype + genotype files from  panels, pedgiree file, and phenotype file of family mean
#This code asseses genomic prediciton accuracy across biparental families (e.g., of family mean)
#This code estimates marker effects using the germplasm panels to predict across families

## Part 1.Load libraries
########### 
library(BGLR)
library(rrBLUP)
############


## Part 2: Load input files
# all individuals must have phenotypes + genotypes
# phenotype names must match between files
# individual names must match between files
#################

### Read in phenotypes + genotypes from panels
#panel_phenotypes <- read.csv("./PathToPhenotypeFile.RDS")
#list of phenotype files where named by panel
#rownames are Lines (above); columns are phenotypes

#panel_genotypes <-  readRDS("./PathToGenotypeFile.RDS")
#list of genotype files where named by panel
#rownames are Lines (above); columns are marker names



### Read in phenotypes from families -- this is either family mean / family BLUE
#family_phenotypes <- read.csv("./PathToPhenotypeFile.csv")
#the first column of this file is "Family" which is the unique family identifier
# family identifier matches pedigree file family name
#the second through nth column are names of phenotypes to predict

#phenos <- colnames(family_phenotypes)[2:ncol(family_phenotypes)]
#cat("\n There are", length(phenos), "phenotypes")

#fams  <- unique(family_phenotypes$Family)
#cat("\n There are", length(fams), "families")


### Read in parent pedigree file
#parent_df <- read.csv("./PathToParentPedigreeFile.csv")
#the first column of this file is "Family" which is the unique family identifier
#the second column of this file is "SeedParent" which is the seed parent genotype identifer that is found in the panel genotype file
#the third column of this file is "PollenParent" which is the pollen parent genotype identifer that is found in the panel genotype file

###################


## Part 3: Make marker and relationship matrices for estimating marker effects + validation
# following https://github.com/MarcooLopez/Genomic-Selection/blob/master/single_environment.md
## get marker matrix for each family mean for validation
#########################

## this panel has all parent genotpyes 
X_Panel1 <- panel_genotypes$Panel1

# then make marker matrix of mean parent value
X_parents <- matrix(nrow = length(unique(parent_df$family)), ncol=ncol(X_Panel1))
dim(X_parents )
rownames(X_parents) <- parent_df$family ; colnames(X_parents) <- colnames(X_Panel1)
X_parents[,1:10]

for (i in 1:nrow(X_parents)) {
  
  X_Panel1_tmp <- X_Panel1[c(which(rownames(X_Panel1) %in% c(as.character(parent_df[i,2]), as.character(parent_df[i,3])))),]
  
  if (nrow(X_Panel1_tmp) == 2) {
    X_parents[i,] <- colMeans(X_Panel1_tmp)
  }
  
}
#X_parents is marker matrix

##################


## Part 4: Estimate marker effects
#################
## these were computed is in GenomicPrediction_WithinBiparentalFamily_PanelAsTrainingPopulation.R
#Mrk_eff_for_GP_panel <- readRDS("./PathToSavePanelMrkEffects.RDS")
cat("\n Marker effects loaded")
################


## Part 5: Correlation between predicted and observed
#################

## predict across family
corRes <- as.data.frame(matrix(ncol=4)) ; colnames(corRes) <- c("Model", "TrainingPop", "pheno", "r")

  for (j in 1:length( Mrk_eff_for_GP_panel$mkrEffs)) { #j = model
    for (k in 1:length( Mrk_eff_for_GP_panel$mkrEffs[[j]])) { #k=panel
      for (m in 1:length(phenos)) { #m = phenos
        
        TmpB <- as.matrix(Mrk_eff_for_GP_panel$mkrEffs[[j]][[k]][[m]])
        
        X_parents_tmp <- X_parents[,which(colnames(X_parents) %in% rownames(TmpB))]
        TmpB <- as.matrix(TmpB[which(rownames(TmpB) %in% colnames( X_parents_tmp)),])
        gHat<-  X_parents_tmp  %*% TmpB
        
        y<- fam_phenos %>% arrange(factor(fam_phenos$Family, levels=rownames(gHat)))
        cor.test(gHat, y[,which(colnames(y) %in% phenos[m])])
        
        corResTemp <- as.data.frame(matrix(ncol=4)) ; colnames(corResTemp) <- c("Model", "TrainingPop", "pheno", "r")
        
        corResTemp[1,] <- c( as.character(names( Mrk_eff_for_GP_panel$mkrEffs)[j]),
                            as.character(names( Mrk_eff_for_GP_panel$mkrEffs[[j]])[k]),
                            as.character(names( Mrk_eff_for_GP_panel$mkrEffs[[j]][[k]])[m]),
                            cor.test(gHat, y[,which(colnames(y) %in% phenos[m])])$estimate)
        
        corRes <- rbind(corRes, corResTemp)
      } 
    } 
  }

corRes
corRes$r <- as.numeric(corRes$r)

#write.csv(corRes, "./PathToPredictionAccuracyResults.csv", row.names=F)
###################
