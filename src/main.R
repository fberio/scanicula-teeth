# Clean and clear everything in the working env
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
if(!is.null(dev.list())) dev.off() # clean all previous figures
# unload and detach all packages
lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))

install.packages("devtools")
require(devtools)

install_version("geomorph", version = "3.3.1", upgrade="never")
install_version("ade4", version = "1.7-15", upgrade="never")
install_version("ggplot2", version = "3.3.2", upgrade="never")
install_version("MASS", version = "7.3-53", upgrade="never")
install_version("mclust", version = "5.4.6", upgrade="never")
install_version("plotly", version = "4.9.2.1", upgrade="never")
install_version("gridExtra", version = "2.3", upgrade="never")
install_version("tidyr", version = "1.1.2", upgrade="never")
install_version("plyr", version = "1.8.6", upgrade="never")
install_version("abind", version = "1.4-5", upgrade="never")
install_version("randomForest", version = "4.6-14", upgrade="never")
install_version("reshape2", version = "1.4.4", upgrade="never")
install_version("dplyr", version = "1.0.2", upgrade="never")
install_version("MKinfer", version = "0.6", upgrade="never")
install_version("REdaS", version = "0.9.3", upgrade="never")
#install_version("RVAideMemoire", version = "0.9-80")
install.packages("RVAideMemoire")
devtools::install_github("vbonhomme/mevolCVP@fe43754fa42ee764d561f89a6e80054756ebec61", upgrade="never")

require(geomorph)
require(ade4)
require(ggplot2)
require(MASS)
require(mclust)
require(plotly)
require(gridExtra)
require(tidyr)
require(plyr)
require(abind)
require(randomForest)
require(reshape2)
require(dplyr)
require(mevolCVP)
require(RVAideMemoire)
require(REdaS)
getCurrentFileLocation <-  function() {
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
currentWorkingDir = getCurrentFileLocation()
setwd(currentWorkingDir)
source(paste(currentWorkingDir, "/utils.R", sep=""))

set.seed(42)

################################## LOADING DATASET #################################
filename <- paste(currentWorkingDir, "/../data/canicula_dataset.txt", sep="")
canicula <- read.table(filename, header=T, dec=".", sep=",")
caniculaRDataFilename <- paste(currentWorkingDir, "/../data/canicula.RData", sep="")
save(canicula, file = caniculaRDataFilename)
load(caniculaRDataFilename)

semilands <- read.table(paste(currentWorkingDir, "/../data/semilandmarks.txt", sep=""), header=F, sep="")

################################ PROCRUSTES SUPERIMPOSITION #################################
rawcoord <- mat2array(canicula[,2:115], dimx = 3)
sup <- gpagen(rawcoord, curves=semilands)

############################# PCA ###########################
mypca <- Rpca(array2mat(sup$coords))
mypca$percent[1:20]
mypca$cpercent[1:20]

################################# DEFINING FACTORS ################################

Indiv <- canicula[,116] #collection reference of specimens 
Pop <- canicula[,118] #populations
TL <- canicula[,119] #total length of specimens
Sex <- canicula[,120] #sex of specimens
Stage <- canicula[,121] #ontogenetic stage of specimens
Jaw <- canicula[,122] #palatoquadrate (upper) or meckelian (lower) teeth
t_pos <- canicula[,123] #tooth position along the jaw (e.g., 1 symphyseal, 20 distal)
Size <- sup$Csize #tooth centroid size

# Creating the general dataset
full_dataset <- cbind(array2mat(sup$coords), canicula[,116:123], sup$Csize)

#Visualisations
## PCA
mypca <- Rpca(full_dataset[,1:114])
mypca$percent[1:2]

graph.pca <- as.data.frame(cbind(mypca$scores[,1:4], full_dataset[,115:123]))
colnames(graph.pca)[1:4] <- c("pc1", "pc2", "pc3", "pc4")

ggplot(data = graph.pca, aes(x=pc1, y=pc2, color=pop)) +
  geom_point() +
  scale_colour_manual(values = c("#808080", "black"),
                      breaks=c("atl", "med"), 
                      labels=c("Atlantic", "Mediterranean")) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = NULL,
                                        linetype = "solid",
                                        color = NULL,
                                        inherit.blank = FALSE))+
  theme(panel.grid.major = element_line(size = 0.5,
                                        linetype = "solid",
                                        colour = "white"))+
  theme(panel.grid.minor = element_line(size = 0.25,
                                        linetype = "solid",
                                        colour = ("white")))+
  theme(strip.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = NULL, 
                                        linetype = "solid",
                                        color = NULL, 
                                        inherit.blank = FALSE))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.text.x = element_text(size =10), 
        axis.text.y = element_text(size =10)) +
  theme(legend.title = element_text(face = "bold"))+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Principal Component 1 (40.50%)", 
       y = "Principal Component 2 (20.66%)")

# Extreme shapes
sup.pca <- gm.prcomp(sup$coords)
canicula.pca.plot <- plot(sup.pca, axis1=1, axis2=2)
picknplot.shape(canicula.pca.plot)

# Centroid sizes
# General data frame with centroid sizes
graph.pca <- as.data.frame(cbind(mypca$scores[,1:2], full_dataset[,115:123]))
# Subsampling by sex and jaw
graph.csize.sex <- filter(graph.pca, sex=="male", jaw=="sup")
graph.csize.sex <- as.data.frame(cbind(as.character(graph.csize.sex[,5]), #population
                                       as.character(graph.csize.sex[,8]), #stage
                                       substr(graph.csize.sex[,10], 6, 7), #tooth position
                                       graph.csize.sex[,11])) #centroid size
colnames(graph.csize.sex) <- c("pop", "stage", "tooth_pos", "csize")
graph.csize.sex$csize <- as.numeric(graph.csize.sex$csize)
#Calculation of mean and sd for each tooth position, per population and stage
csize.means <- aggregate(graph.csize.sex$csize~as.character(graph.csize.sex$pop)*as.character(graph.csize.sex$stage)*graph.csize.sex$tooth_pos, 
                         data = graph.csize.sex, function(x) c(mean=mean(x),sd=sd(x)))
csize.means <- as.data.frame(cbind(csize.means[,1:3], csize.means$`graph.csize.sex$csize`[,1],
                                   csize.means$`graph.csize.sex$csize`[,2]))
colnames(csize.means) <- c("pop", "stage", "tooth_pos", "mean", "sd")

#Generation of the plot
ggplot(data = csize.means, aes(x = as.numeric(tooth_pos), y = mean, color = pop)) +
  geom_point(aes(group = paste(pop, stage), shape = stage)) +
  geom_line(aes(group = paste(pop, stage))) +
  scale_linetype_discrete(breaks=c("atl", "med"), 
                          labels = c("Atlantic", "Mediterranean")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  scale_y_continuous(limits = c(0, 5)) +
  scale_x_continuous(breaks = c(1,5,10,15,20,25)) +
  scale_colour_manual(values=c("#808080", "black"),
                      breaks=c("atl", "med"),
                      labels=c("Atlantic", "Mediterranean")) +
  labs(x = "Tooth proximo-distal position", y = "Centroid size") +
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid")) +
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10)) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = NULL,
                                        color = NULL,
                                        inherit.blank = FALSE))+
  theme(panel.grid.major = element_line(size = 0.5,
                                        linetype = "solid",
                                        colour = "white"))+
  theme(panel.grid.minor = element_line(size = 0.25,
                                        linetype = "solid",
                                        colour = ("white"))) +
  theme(strip.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = NULL,
                                        linetype = "solid",
                                        color = NULL,
                                        inherit.blank = FALSE)) +
  theme(legend.title = element_blank()) +
  theme(legend.key.size = unit(6, "mm"), 
        legend.key = element_rect(fill = "white"),
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13))

## Creating dataset for RF process: aligned coords + pop (RF on SHAPE DATA, i.e. /!\ no CSize)
dataset_forest <- as.data.frame(cbind(full_dataset[,1:114], full_dataset[,117], full_dataset[,119:120]))
colnames(dataset_forest)[115:117] <- c("pop", "sex", "stage")

########### TUNING PARAMETERS #############
## Formatting train & test sets
train <- dataset_forest %>% sample_frac(0.8) 
test <- anti_join(dataset_forest, train)

# Formating factors for the randomForest function
dataset_forest[,115] <- as.factor(dataset_forest[,115])
dataset_forest[,116] <- as.factor(dataset_forest[,116])
dataset_forest[,117] <- as.factor(dataset_forest[,117])

## Step 1: Train
nbFolds <- 5
nbReplication <- 1
#getAccuracy
accuracies_model <- vector(mode="numeric", length=nbFolds * nbReplication)
accuracyCpt_model <- 1
accuracies_pred <- vector(mode="numeric", length=nbFolds * nbReplication)
accuracyCpt_pred <- 1
#getRecall
recall_model <- matrix(0, nrow=nbReplication,ncol=2)
colnames(recall_model) <- c("atl", "med")
recall_pred <- matrix(0, nrow=nbReplication,ncol=2)
colnames(recall_pred) <- c("atl", "med")
#getPrecision
prec_model <- matrix(0, nrow=nbReplication,ncol=2)
colnames(prec_model) <- c("atl", "med")
prec_pred <- matrix(0, nrow=nbReplication,ncol=2)
colnames(prec_pred) <- c("atl", "med")

list.preds <- vector(mode = "list", length = nbReplication)

# The new randomForest package does not handle column name that start with a numeric value.
# Add a letter to the colnames starting with a numeric value.
for (i in 1:114) {
  colnames(dataset_forest)[i] <- paste("l", i, sep="")
  colnames(test)[i] <- paste("l", i, sep="")
}

for (i in 1:nbReplication) {
  print(i)
  metrics.preds <- matrix(0, nrow=12, ncol=3)
  if(exists("tmp")) {
    rm(tmp)
  }
  for (test in getTestSets(dataset_forest, nbFolds)) {
    model.rf <- randomForest(pop ~ ., data = suppressMessages(anti_join(dataset_forest, test)[,1:115]), 
                             ntree = 500, importance = TRUE, mtry = 114, 
                             nodesize = 1)
    accuracies_model[accuracyCpt_model] <- computeAccuracy(model.rf$confusion[1:2,1:2])
    accuracyCpt_model <- accuracyCpt_model + 1
    recall_model[i,1] <- computeRecall(model.rf$confusion, "atl")
    recall_model[i,2] <- computeRecall(model.rf$confusion, "med")
    prec_model[i,1] <- computePrecision(model.rf$confusion, "atl")
    prec_model[i,2] <- computePrecision(model.rf$confusion, "med")
    
    model.pred <- predict(model.rf, test[,1:115])
    confMat.pred <- table(observed=test$pop, predicted=model.pred)
    accuracies_pred[accuracyCpt_pred] <- computeAccuracy(confMat.pred)
    accuracyCpt_pred <- accuracyCpt_pred + 1
    recall_pred[i,1] <- computeRecall(confMat.pred, "atl")
    recall_pred[i,2] <- computeRecall(confMat.pred, "med")
    prec_pred[i,1] <- computePrecision(confMat.pred, "atl")
    prec_pred[i,2] <- computePrecision(confMat.pred, "med")
    if(exists("tmp")) {
      tmp <- table(actual = paste(as.character(test[,115]), 
                                  as.character(test[,116]),
                                  as.character(test[,117])), 
                   predicted = paste(as.character(model.pred), 
                                     as.character(test[,116]), 
                                     as.character(test[,117]))) 
      + tmp 
    }
    else{
      tmp <- table(actual = paste(as.character(test[,115]), 
                                  as.character(test[,116]),
                                  as.character(test[,117])), 
                   predicted = paste(as.character(model.pred), 
                                     as.character(test[,116]), 
                                     as.character(test[,117])))
    }
  }
  for (popsexstage in 1:length(colnames(tmp))) {
    metrics.preds[popsexstage,1] <- colnames(tmp)[popsexstage]
    metrics.preds[popsexstage,2] <- computePrecision(tmp, popsexstage)
    metrics.preds[popsexstage,3] <- computeRecall(tmp, popsexstage)
  }
  list.preds[[i]] <- metrics.preds
  results <- list(accuracies_model, accuracies_pred, prec_model, prec_pred,
                  recall_model, recall_pred)
}

# results for Random Forests on shape data
list.preds # 12 values for precision and recall
results # (accuracy/precision/recall) * (atl/med) 

#Visualise accuracies, recalls, and precisions train/test for both populations
boxplot(results[[1]], results[[2]], names = c("model", "pred"), main = "Accuracies")
par(mfrow=c(1,2))
boxplot(recall_model, main = "Recalls model", border = "brown", col = "orange")
boxplot(recall_pred, main = "Recalls prediction")
boxplot(prec_model, main = "Precisions model")
boxplot(prec_pred, main = "Precisions prediction")
par(mfrow=c(1,1))

#Evaluate feature importance
importance.feats <- as.data.frame(cbind(model.rf$importance, seq(1,114,1))) #114 for shape data, 115 for form data
colnames(importance.feats)[5] <- "lmk_semilmk" #x, y or z coordinate of landmarks and semilandmarks
acc <- importance.feats[order(importance.feats[,3], decreasing = TRUE),]
what.to.plot <- acc
column.to.plot <- 3
par(mfrow=c(1,1))
plot(what.to.plot[,column.to.plot], type = "l", xlab = "Feature",
     ylab = "Feature importance")
text(what.to.plot[,column.to.plot], labels = what.to.plot[,5], cex=0.7)

#### RF on Form (shape + centroid size data) ####
set.seed(42)

dataset_forest.csize <- as.data.frame(cbind(full_dataset[,1:114], full_dataset[,123], full_dataset[,117], full_dataset[,119:120]))
colnames(dataset_forest.csize)[115:118] <- c("csize", "pop", "sex", "stage")

nbFolds <- 5
nbReplication <- 1
#getAccuracy
accuracies_model <- vector(mode="numeric", length=nbFolds * nbReplication)
accuracyCpt_model <- 1
accuracies_pred <- vector(mode="numeric", length=nbFolds * nbReplication)
accuracyCpt_pred <- 1
#getRecall
recall_model <- matrix(0, nrow=nbReplication,ncol=2)
colnames(recall_model) <- c("atl", "med")
recall_pred <- matrix(0, nrow=nbReplication,ncol=2)
colnames(recall_pred) <- c("atl", "med")
#getPrecision
prec_model <- matrix(0, nrow=nbReplication,ncol=2)
colnames(prec_model) <- c("atl", "med")
prec_pred <- matrix(0, nrow=nbReplication,ncol=2)
colnames(prec_pred) <- c("atl", "med")

list.preds <- vector(mode = "list", length = nbReplication)

for (i in 1:114) {
  colnames(dataset_forest.csize)[i] <- paste("l", i, sep="")
  colnames(test)[i] <- paste("l", i, sep="")
}

dataset_forest.csize[,116] <- as.factor(dataset_forest.csize[,116])

for (i in 1:nbReplication) {
  print(i)
  metrics.preds <- matrix(0, nrow=12, ncol=3)
  if(exists("tmp")) {
    rm(tmp)
  }
  for (test in getTestSets(dataset_forest.csize, nbFolds)) {
    model.rf <- randomForest(pop ~ ., data = suppressMessages(anti_join(dataset_forest.csize, test)[,1:116]), 
                             ntree = 500, importance = TRUE, mtry = 115,
                             nodesize = 1)
    accuracies_model[accuracyCpt_model] <- computeAccuracy(model.rf$confusion[1:2,1:2])
    accuracyCpt_model <- accuracyCpt_model + 1
    recall_model[i,1] <- computeRecall(model.rf$confusion, "atl")
    recall_model[i,2] <- computeRecall(model.rf$confusion, "med")
    prec_model[i,1] <- computePrecision(model.rf$confusion, "atl")
    prec_model[i,2] <- computePrecision(model.rf$confusion, "med")
    
    model.pred <- predict(model.rf, test[,1:116])
    confMat.pred <- table(observed=test$pop, predicted=model.pred)
    accuracies_pred[accuracyCpt_pred] <- computeAccuracy(confMat.pred)
    accuracyCpt_pred <- accuracyCpt_pred + 1
    recall_pred[i,1] <- computeRecall(confMat.pred, "atl")
    recall_pred[i,2] <- computeRecall(confMat.pred, "med")
    prec_pred[i,1] <- computePrecision(confMat.pred, "atl")
    prec_pred[i,2] <- computePrecision(confMat.pred, "med")
    if(exists("tmp")) {
      tmp <- table(actual = paste(as.character(test[,116]), 
                                  as.character(test[,117]),
                                  as.character(test[,118])), 
                   predicted = paste(as.character(model.pred), 
                                     as.character(test[,117]), 
                                     as.character(test[,118]))) 
      + tmp 
    }
    else{
      tmp <- table(actual = paste(as.character(test[,116]), 
                                  as.character(test[,117]),
                                  as.character(test[,118])), 
                   predicted = paste(as.character(model.pred), 
                                     as.character(test[,117]), 
                                     as.character(test[,118])))
    }
  }
  for (popsexstage in 1:length(colnames(tmp))) {
    metrics.preds[popsexstage,1] <- colnames(tmp)[popsexstage]
    metrics.preds[popsexstage,2] <- computePrecision(tmp, popsexstage)
    metrics.preds[popsexstage,3] <- computeRecall(tmp, popsexstage)
  }
  list.preds[[i]] <- metrics.preds
  results <- list(accuracies_model, accuracies_pred, prec_model, prec_pred,
                  recall_model, recall_pred)
}

# results for Random Forests on form data
list.preds # 12 values for precision and recall
results # (accuracy/precision/recall) * (atl/med) 

# LDA on shape and form data
performLDA <-  function(shapeOrForm, useRawCoord) {
  shapeMat <- Rpca(array2mat(sup$coords))
  if (useRawCoord) {
    myData <- array2mat(sup$coords)
  } else {
    myData <- shapeMat$scores[,1:12]
  }
  if (shapeOrForm == "shape") {
    dataset <- myData #shape data
  } else {
    dataset <- cbind(myData, sup$Csize) #form data
  }
  population <- droplevels(as.factor(Pop))
  gpe <- paste(Sex, Pop)
  shapeMatDetails <- as.data.frame(cbind(shapeMat$scores[,1:ncol(myData)], 
                                         sup$Csize, TL, Pop, Sex, Stage,
                                         gpe=paste(Sex,Pop)))
  
  set.seed(42)
  # Formatting data
  lda.data <- as.matrix(shapeMatDetails[,c(1:ncol(myData))])
  lda_raw <- cbind(lda.data, as.character(shapeMatDetails$gpe))
  colnames(lda_raw)[ncol(myData)+1] <- "group"
  
  nbFolds <- 5
  nbReplication <- 1
  #getAccuracy
  accuracies_model <- vector(mode="numeric", length=nbFolds * nbReplication)
  accuracyCpt_model <- 1
  accuracies_pred <- vector(mode="numeric", length=nbFolds * nbReplication)
  accuracyCpt_pred <- 1
  #getRecall
  recall_model <- matrix(0, nrow=nbReplication,ncol=2)
  colnames(recall_model) <- c("atl", "med")
  recall_pred <- matrix(0, nrow=nbReplication,ncol=2)
  colnames(recall_pred) <- c("atl", "med")
  #getPrecision
  prec_model <- matrix(0, nrow=nbReplication,ncol=2)
  colnames(prec_model) <- c("atl", "med")
  prec_pred <- matrix(0, nrow=nbReplication,ncol=2)
  colnames(prec_pred) <- c("atl", "med")
  
  list.preds <- vector(mode = "list", length = nbReplication)
  
  for (i in 1:nbReplication) {
    metrics.preds <- matrix(0, nrow=12, ncol=3)
    for (testIndexes in getTestSets(dataset, nbFolds)) {
      model.lda <- lda(x = dataset[-testIndexes,],
                       grouping=population[-testIndexes])
      # lda predictions on test set
      predict.lda <- predict(model.lda, newdata=dataset[testIndexes,])
      confMat.pred <- table(observed=population[testIndexes], predicted=predict.lda$class)
      accuracies_pred[accuracyCpt_pred] <- computeAccuracy(confMat.pred)
      accuracyCpt_pred <- accuracyCpt_pred + 1
      recall_pred[i,1] <- computeRecall(confMat.pred, "atl")
      recall_pred[i,2] <- computeRecall(confMat.pred, "med")
      prec_pred[i,1] <- computePrecision(confMat.pred, "atl")
      prec_pred[i,2] <- computePrecision(confMat.pred, "med")
      if(exists("tmp")) {
        tmp <- table(actual = paste(as.character(population[testIndexes]),
                                    as.character(Sex[testIndexes]),
                                    as.character(Stage[testIndexes])),
                     predicted = paste(as.character(predict.lda$class),
                                       as.character(Sex[testIndexes]),
                                       as.character(Stage[testIndexes])))
        + tmp
      } else {
        tmp <- table(actual = paste(as.character(population[testIndexes]),
                                    as.character(Sex[testIndexes]),
                                    as.character(Stage[testIndexes])),
                     predicted = paste(as.character(predict.lda$class),
                                       as.character(Sex[testIndexes]),
                                       as.character(Stage[testIndexes])))
      }
    }
    for (popsexstage in 1:length(colnames(tmp))) {
      metrics.preds[popsexstage,1] <- colnames(tmp)[popsexstage]
      metrics.preds[popsexstage,2] <- computePrecision(tmp, popsexstage)
      metrics.preds[popsexstage,3] <- computeRecall(tmp, popsexstage)
    }
    list.preds[[i]] <- metrics.preds
    results <- list(accuracies_model, accuracies_pred, prec_model, prec_pred,
                    recall_model, recall_pred)
  }
  
  write.table(metrics.preds, paste(currentWorkingDir, "/../data/AccPrecRec_AtlMed_SexStage_", shapeOrForm, ".txt", sep=""),
              sep = ",", dec = ".")
  print(mean(accuracies_pred))
  print(list.preds)
  print(results)
  shapeMat <- Rpca(array2mat(sup$coords))
}

# LDA results 
performLDA("shape", useRawCoord =  FALSE)
performLDA("form", useRawCoord = FALSE)

#Tooth file number
data_nb_files <- read.table(paste(currentWorkingDir, "/../data/nb_files.txt", sep=""), header = T, sep = ",") 
data_pergroup <- filter(data_nb_files, jaw == "pq", 
                        sex == "f",
                        stage == "mat",
                        pop == "med")
mean(data_pergroup$nb_max)
sd(data_pergroup$nb_max)

#Visualization of tooth file number per sex, stage, and population
data_pergroup <- filter(data_nb_files, sex == "f") # f for females, m for males
ggplot(data_pergroup, aes(x = jaw, y = nb_max, fill = paste(stage))) +
  geom_boxplot(aes(linetype = pop)) +
scale_fill_manual(values=c("firebrick1", "dodgerblue2", "purple1"),
                  breaks=c("hat", "juv", "mat"),
                  labels=c("Hatchling", "Juvenile", "Mature")) +
  scale_y_continuous(name = "Number of tooth files",
                     breaks = seq(15, 30, 5), # (15, 30, 5) for meckelian; (17, 28, 5) for palatoquadrate
                     limits=c(15, 30)) + # (15, 30) for meckelian; (17, 28) for palatoquadrate
  scale_x_discrete(name = "",
                   breaks = c("pq", "mc"),
                   labels = c("Palatoquadrate", "Meckelian")) +
  theme(axis.line = element_line(colour = "black",
                                 size = 1, linetype = "solid")) +
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10)) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = NULL,
                                        color = NULL,
                                        inherit.blank = FALSE))+
  theme(panel.grid.major = element_line(size = 0.5,
                                        linetype = "solid",
                                        colour = "white"))+
  theme(panel.grid.minor = element_line(size = 0.25,
                                        linetype = "solid",
                                        colour = ("white"))) +
  theme(strip.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = NULL,
                                        linetype = "solid",
                                        color = NULL,
                                        inherit.blank = FALSE)) +
  theme(legend.title = element_blank()) +
  #theme_bw() +
  theme(legend.key.size = unit(6, "mm"),
        legend.key = element_rect(fill = "white"),
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13))

############ GRAPHS ACCURACY, PRECISION, AND RECALL
model <- c("LDA shape", "LDA form", "RF shape", "RF form")
accuracies <- c(64.5, 74.5, 81.8, 87.4)
sd_accuracies <- c(0.7, 1.2, 1.4, 1.7)

population <- c("Atlantic", "Mediterranean", "Atlantic", "Mediterranean", 
                "Atlantic", "Mediterranean","Atlantic", "Mediterranean",
                "Atlantic", "Mediterranean","Atlantic", "Mediterranean",
                "Atlantic", "Mediterranean","Atlantic", "Mediterranean")
data <- c("shape", "shape", "shape", "shape", "shape", "shape", "shape", "shape",
          "form", "form", "form", "form", "form", "form", "form", "form")
metric <- c("precision", "precision", "precision", "precision", "recall", "recall", "recall",
            "recall", "precision", "precision", "precision", "precision", "recall", "recall",
            "recall", "recall")


mydf <- data.frame(xp=c("1LDA shapeM", "2LDA shapeA", "3break", "4LDAformM", "5LDAformA", 
                        "6break", "7RF shapeM", "8RF shapeA", "91break", "92RF formM", 
                        "93RF formA"), 
                   population=c("Mediterranean", "Atlantic", "", "Mediterranean", "Atlantic", "",
                                "Mediterranean", "Atlantic", "", "Mediterranean", "Atlantic"),
                   precisions=c(65.7, 63.5, NA, 74.6, 74.5, NA, 83.5, 80.3, NA, 85.1, 86.7),
                   recalls=c(53.6, 74.3, NA, 71.0, 77.8, NA, 75.6, 86.9, NA, 82.8, 88.6))
                  

ggplot(data=mydf, aes(x=xp, y=recalls, fill=population)) +
  geom_bar(stat="identity", position = "dodge") + 
  geom_text(size=3, aes(label=recalls), position=position_dodge(width=0.9), vjust=-0.4) +
  theme_classic() + 
  labs(y="Recall", x="Condition") + 
  scale_x_discrete(name = "", breaks = c("1LDA shapeM", "2LDA shapeA", "3break", "4LDAformM", "5LDAformA", 
                                         "6break", "7RF shapeM", "8RF shapeA", "91break", "92RF formM", 
                                         "93RF formA"), 
                                         labels = c("LDA \n shape", "", "", 
                                                    "LDA \n form", "", "",
                                                    "RF \n shape", "", "",
                                                    "RF \n form", "")) +
     scale_y_continuous(name = "Recall (%)",
                     breaks = seq(0, 100, 20),
                     limits=c(0, 100)) +
     scale_fill_manual(values=c("white", "#808080", "black"))

# Differences in tooth file number between populations?
set.seed(42)
atl <- as.data.frame(filter(data_nb_files, sex == "f", jaw == "mc",
                            pop == "atl", stage == "juv")[,5:6])
med <- as.data.frame(filter(data_nb_files, sex == "f", jaw == "mc",
                            pop == "med", stage == "juv")[,5:6])
dat.t.test <- rbind(atl, med)
perm.t.test(dat.t.test$nb_max~dat.t.test$pop, nperm = 100)
