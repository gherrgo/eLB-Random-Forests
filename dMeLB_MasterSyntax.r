#######################################################################################
## Diagnostic-Meningioma epigenetic Liquid Biopsy (d-MeLB) Construction & Validation ##
#######################################################################################
# All IDATs may be downloaded using :: DOI: 10.17632/zrc982rvjm.1
# All excel files referenced are located within the Supplementary Data

#Necessary package loading & function definition
lapply(c('readxl', 'caret', 'parallel', 'doMC', 'pROC', 'ggrepel', 'impute',
         'DescTools', 'magrittr', 'dplyr', 'minfi', 'tibble', 'data.table'), library, character.only = TRUE)
'%notin%' <- Negate("%in%")
'%notlike%' <- Negate('%like%')
row_to_column <- function(df){
  var = colnames(df)[1]
  df <- df %>% remove_rownames() %>% column_to_rownames(var = var)
}

###########################################################################
### Step #1: Definition of parent signature set: Tissue Supervised Analysis.
###########################################################################

### Preprocessing of internal meningioma ______
RGset <- read.metharray.exp(base = 'Tissue IDATs/', force = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
Mset <- preprocessIllumina(RGset)
pan.mm <- data.frame(getBeta(Mset, type = 'EPIC'))
colnames(pan.mm) <- sub('X', '', colnames(pan.mm))

### Preprocessing of publicly available nontumor controls ACCESSION #: GSE111165 ______
RGset <- read.metharray.exp(base = 'GSE111165/', force = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
Mset <- preprocessIllumina(RGset)
epilep <- data.frame(getBeta(Mset, type = 'EPIC'))
colnames(epilep) <- sub('X', '', colnames(epilep))

total <- row_to_column(merge(pan.mm, epilep, by = 0))

### Unmasking procedures and NA-value removal of entire genome _____
## Function built with consideration of EPIC hg38 Manifest
unmask <- function(df){
  #Loading EPIC manifest :: available for download at https://zwdzwd.github.io/InfiniumAnnotation 
  probe.anno <- read.table("EPIC.hg38.manifest.tsv.gz", sep = "\t", header = T)
  #Removal of masked probes
  probe.anno <- probe.anno[!probe.anno$MASK_general & probe.anno$CpG_chrm %in% paste0("chr", 1:22),]
  # Subsetting of probes found after unmasking
  df <- df[rownames(df) %in% probe.anno$probeID,]
  #Removal of any NA values
  df <- df[rowSums(is.na(df)) == 0,]
  
}
total <- unmask(total)

### Supervised Analysis: meningioma and non-tumor controls
## Base Wilcoxon-test function: returns p-values across whole genome
my.wilcox.test.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error"))
    return(NA)
  else
    return(obj$p.value)
}
## Group definition:
group1 <- colnames(pan.mm)
group2 <- colnames(epilep)

## Function running and p-value recording
p <- apply(total, 1, function(x) {
  zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
  return(zz)
})
total$p.value.raw <- p
total$p.value.adj <- p.adjust(p, method = "fdr")

## Need to save CpG names and p-values for later incorporation into random forest
cpgs <- data.frame('CpG' = c(rownames(total)),
           'p.value.adj' = c(total$p.value.adj))

###########################################################################
### Step #2: Machine-driven randomization of samples into Discovery (80%)## 
### & Independent (20%) Sets.
###########################################################################

# Loading Matching Tissue & Serum (MTS) cases clinical data
mts.cases <- read_excel('SupplementalData1.xlsx', sheet = 'MTS Case IDs')
# Subsetting only MTS cases across the tissue-meningioma methylation matrix and formatting column names appropriately
pan.mm <- unmask(pan.mm[,mts.cases$Tissue.Barcode])
colnames(pan.mm) <- mts.cases$`Tissue ID`

# Loading serum-based clinical data
pd <- read_excel('SupplementalData1.xlsx')
pd2 <- read_excel('SupplementalData1.xlsx', sheet = 'Additional Samples')
pd2 <- pd2[pd2$`Sample Type` == "Non-neoplastic disease",]

pd <- data.frame('Sample.ID' = c(pd$`Sample ID`, pd2$`Sample ID`),
           'Sentrix.ID' = c(pd$`Sentrix ID`, pd2$`Sentrix ID`),
           'Sample Type' = c(pd$`Sample Type`, pd2$`Sample Type`),
           'Histotype' = c(pd$Histotype, pd2$Histotype),
           'Time of Sample Collection' = c(pd$`Time of Sample Collection`, pd2$`Time of Sample Collection`))
names(pd) <- c("Sample ID", "Sentrix ID", "Sample Type", 'Histotype', "Time of Sample Collection")

# Appropriate exclusions
excluded.samples <- pd[pd$`Sample Type` == "Glioma" & pd$`Time of Sample Collection` == "Recurrence" |
                         pd$Histotype %in% c('Infection', 'Necrosis', 'Abscess'),]
pd <- pd[pd$`Sample ID` %notin% excluded.samples$`Sample ID`,]

#### Labeling MNG by treatment
pd3 <- read_excel('SupplementalData1.xlsx', sheet = 'MNG Specific Features')
treatment <- pd3$`Sample ID`[pd3$`Pre-Surgical Radiotherapy (Y/N)` == "Yes" | 
                               pd3$`Pre-Surgical Embolization (Y/N)` == "Yes"]
pd$Treatment <- ifelse(pd$`Sample ID` %in% treatment, "Treated", "Untreated")
pd$`Sample Type` <- ifelse(pd$`Sample Type` == "Meningioma" & pd$`Time of Sample Collection` == "Primary" & pd$Treatment == "Treated", "Meningioma (IT)", pd$`Sample Type`)
pd$`Sample Type` <- ifelse(pd$`Sample Type` == "Meningioma" & pd$`Time of Sample Collection` == "Recurrence" & pd$Treatment == "Treated", "Meningioma (RT)", pd$`Sample Type`)
pd$`Sample Type` <- ifelse(pd$`Sample Type` == "Meningioma" & pd$`Time of Sample Collection` == "Primary" & pd$Treatment == "Untreated", "Meningioma (IU)", pd$`Sample Type`)
pd$`Sample Type` <- ifelse(pd$`Sample Type` == "Meningioma" & pd$`Time of Sample Collection` == "Recurrence" & pd$Treatment == "Untreated", "Meningioma (RU)", pd$`Sample Type`)

###### Discovery & Validation derivation 
discovery <- pd %>%
  group_by(`Sample Type`) %>% sample_frac(0.8) %>%
  ungroup()
validation <- pd[pd$`Sample ID` %notin% discovery$`Sample ID`,]

# Conglomeration of total validation (n=122 samples)
# Additional MNG: formatting
pd4 <- read_excel('SupplementalData1.xlsx', sheet = 'Additional Samples')
pd4 <- pd4[pd4$`Sample Type` == "Meningioma",]
pd4$`Sample Type` <- ifelse(pd4$`Time of Sample Collection` == 'Initial' & pd4$`Post-Surgical Radiotherapy (Y/N)` == 'Yes', 'Meningioma (IT)', pd4$`Sample Type`)
pd4$`Sample Type` <- ifelse(pd4$`Time of Sample Collection` == 'Recurrent' & pd4$`Post-Surgical Radiotherapy (Y/N)` == 'Yes', 'Meningioma (RT)', pd4$`Sample Type`)
pd4$`Sample Type` <- ifelse(pd4$`Time of Sample Collection` == 'Initial' & pd4$`Post-Surgical Radiotherapy (Y/N)` == 'No', 'Meningioma (IU)', pd4$`Sample Type`)
pd4$`Sample Type` <- ifelse(pd4$`Time of Sample Collection` == 'Recurrent' & pd4$`Post-Surgical Radiotherapy (Y/N)` == 'No', 'Meningioma (RU)', pd4$`Sample Type`)

validation <- data.frame(`Sample ID` = c(validation$`Sample ID`, excluded.samples$`Sample ID`, pd4$`Sample ID`),
                         `Sentrix ID` = c(validation$`Sentrix ID`, excluded.samples$`Sentrix ID`, pd4$`Sentrix ID`),
                         `Sample Type` = c(validation$`Sample Type`, excluded.samples$`Sample Type`, pd4$`Sample Type`),
                         `Histotype` = c(validation$Histotype, excluded.samples$Histotype, pd4$Histotype),
                         `Time of Sample Collection` = c(validation$`Time of Sample Collection`, excluded.samples$`Time of Sample Collection`, pd4$`Time of Sample Collection`))
names(validation) <- names(discovery)[-6]

##### Preprocessing of Quality Checked IDATs into methylation matrices (discovery & validation separate)
RGset <- read.metharray.exp(base = 'LB IDATs/', force = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
Mset <- preprocessIllumina(RGset)
beta <- data.frame(getBeta(Mset, type = 'EPIC'))
colnames(beta) <- sub('X', '', colnames(beta))

discovery.beta <- beta[,discovery$`Sentrix ID`] # 865k CpGs x 117 samples
colnames(discovery.beta) <- discovery$`Sample ID`

validation.beta <- beta[,validation$`Sentrix ID`] # 865k CpGs x 30 samples
colnames(validation.beta) <- validation$`Sample ID`

### Optional Garbage collection: free up unused memory
# gc()

#### Setting up discovery matrix and including binary classification of samples 
discovery.data <- data.frame(t(discovery.beta[,discovery$`Sample ID`]))
discovery.data$class <- ifelse(discovery$`Sample Type` %like% "Meningioma", "MNG", "Non")

###########################################################################
### Steps #3-9: Defining the master function ###
###########################################################################
# Miscellaneous Hyperparameters -- tissue similarity min and maximum; minimum # of signatures desired
range_min <- 0.0005
range_max <- 0.002
min_probes <- 20


menin_RF <- function(discovery.data, seed1, ntree = 1000, verbose=FALSE, Ncores=10){
  
  # Set your seed so your work is repeatable
  set.seed(seed1)
  ########################################################################################
  ### Step #3: Create a subset of your data for Training and the other for Model Selection
  ########################################################################################
  
  inTraining <- createDataPartition(discovery.data$class, p=0.8, list=FALSE, times=1) 
  # Training Set : 80% of the total discovery set
  myTrain <- discovery.data[inTraining, ]
  pdata.train <- discovery[discovery$`Sample ID` %in% rownames(myTrain),]
  # Model Selection Set : 20% of the total discovery set
  myMS <- discovery.data[-inTraining, ]
  pdata.MS <- discovery[discovery$`Sample ID` %in% rownames(myMS),]
  
  ### Misc: Tuning random parameters
  # Setting d-MeLB signature set size (random for each separate iteration)
  range <- c(20:30)
  N_probe <- sample(range, 1)
  
  # Setting MNG & Nontumor comparison p-value cutoff (significance; random for each separate iteration)
  range <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
  coeff <- cpgs$CpG[cpgs$p.value.adj <= sample(range, 1)]
  
  # Setting matching tumor tissue and liquid biopsy similarity cutoff (aka diff.mean < cutoff)
  similarity <- round(runif(1, min = range_min, max = range_max), 4)
  
  ###############################################################################################
  ### Step #4: Dimensionality reduction using matching tumor tissue & LB Cases w/in training set.
  ###############################################################################################
  
  mts.cases.tmp <- mts.cases[mts.cases$`Serum ID` %in% rownames(myTrain),]
  dimension.reduc <- function(x){ 
    ### IF statement: if there are more than 2 pairs of matching tissue & serum (MTS) in the randomized training set
    if(length(mts.cases.tmp$`Serum ID`) >=2){
      # Subsetting serum samples found in the MTS cases
      beta.tmp <- discovery.beta[x, mts.cases.tmp$`Serum ID`]
      # Subsetting tissue samples found in the MTS cases
      pan.tmp <- pan.mm[x, mts.cases.tmp$`Tissue ID`]
      # Merging serum and tissue
      beta.tmp <- row_to_column(merge(beta.tmp, pan.tmp, by = 0))
      # Calculating the differences in average serum and tissue methylation across the methylome 
      beta.tmp$diff.mean <- abs(rowMeans(beta.tmp[,mts.cases.tmp$`Serum ID`])-rowMeans(beta.tmp[,mts.cases.tmp$`Tissue ID`]))
      # Defining the dimension-reduced probes based on the similarity cutoff that is prior randomized
      sigProbes <- rownames(beta.tmp)[beta.tmp$diff.mean<similarity]
      return(sigProbes)
    }
    ### IF statement: if there is only 1 pair of matching tissue and serum in the randomized training set
    if(length(mts.cases.tmp$`Serum ID`)==1){
      beta.tmp <- row_to_column(merge(beta, pan.mm, by = 0))
      sigProbes <- rownames(beta.tmp)[abs(beta.tmp[x,mts.cases.tmp$`Serum ID`]-beta.tmp[x,mts.cases.tmp$`Tissue ID`])<similarity]
      return(sigProbes)
    }
    ### IF statement: if there are no pairs of matching tissue and serum in the randomized training set
    if(length(mts.cases.tmp$`Serum ID`)<1){
      sigProbes <- rownames(beta)
      return(sigProbes)
    }
    else(return(sigProbes))
  }
  
  ### Dimension Reduction: Conducted across Nontumor probes (labeled as: coeff)
  sigProbes <- dimension.reduc(coeff)
  
  ###############################################################################################
  ### Step #5: Supervised analysis series to define signature set (DMPs); Wilcoxon rank-sum test.
  ###############################################################################################
  
  #Supervised Analyses: Group 1 - Untreated MNG; Group 2- Non-MNG (Excluding Glioblastoma)
  group1 <- pdata.train$`Sample ID`[pdata.train$`Sample Type` == "Meningioma (IU)" | pdata.train$`Sample Type` == 'Meningioma (RU)']
  group2 <- pdata.train$`Sample ID`[pdata.train$`Sample Type` %notlike% "Meningioma" & pdata.train$Histotype != "Glioblastoma"]
  #Sampling of non-MNG in effort to potentially exclude unknown outliers
  group2 <- sample(group2, length(group2)*0.8)
  
  # Creation of temporary beta-matrix for conducting statistical tests
  beta.tmp <- discovery.beta[sigProbes, c(group1, group2)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(beta.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
    return(zz)
  })
  
  # P-value FDR adjustment
  p.adj <- p.adjust(p, method = "fdr")
  
  # Assignment of results to methylation matrix
  beta.tmp$p.value.raw <- p
  beta.tmp$p.value.adj <- p.adj
  # Mean methylation difference calculations 
  beta.tmp$mean.group1 <- apply(beta.tmp[, group1], 1, mean, na.rm = T)
  beta.tmp$mean.group2 <- apply(beta.tmp[, group2], 1, mean, na.rm = T)
  beta.tmp$diff.mean <- beta.tmp$mean.group1 - beta.tmp$mean.group2
  
  # Sampling of probes according to p-value cutoff and highest diff.mean
  sigProbes2 <- beta.tmp %>%
    filter(p.value.adj <= 0.01) %>% rownames_to_column(var = 'rowname') %>% arrange(desc(abs(diff.mean))) %>%
    slice_head(n = N_probe) %>% pull(rowname)
  
  ### Reassurance for number of probes - in case we do not have sufficient CpG sets sizes
  # Less stringent significances: only run in case probe set size is lower than 'min_probes' required for RF building
  less.sig = function(x){
    beta.tmp %>%
      dplyr::filter(p.value.adj <= 0.05) %>% rownames_to_column(var = 'rowname') %>% dplyr::arrange(desc(abs(diff.mean))) %>%
      slice_head(n = N_probe) %>% pull(rowname)
  }
  reassur.probes <- less.sig(beta.tmp)
  if(length(sigProbes2)<min_probes){
    sigProbes2 <- reassur.probes
  }
  
  # Seed definition
  set.seed(seed1)
  
  ###############################################################################################
  ### Step #6: Generation of predictive Random Forest classifier.
  ###############################################################################################
  
  # Subsetting significant probes across the training data set
  myTrain1 <- discovery.data[inTraining, c(sigProbes2, "class")]
  
  # Setting up k-fold cross validation
  fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv", number = 10,
    ## repeated five times
    repeats = 5)
  
  # Mtry value setting -- defining minimum number of signatures required for decision trees
  median.mtry <- floor(sqrt(ncol(myTrain1)-1))
  mtryVals <- seq(median.mtry-3, median.mtry+4, by=1)
  mtryGrid <- data.frame(.mtry=mtryVals)
  
  # Confirm seed again
  set.seed(seed1)
  
  # Set number of computer cores - this can be altered
  registerDoMC(cores = 12)
  
  # Construction of actual random forest object 
  RF.obj <- train(class ~ ., # variable to be trained on
                  data = myTrain1, # Data we are using
                  method = "rf", # Method we are using
                  trControl = fitControl, # How we validate
                  ntree = ntree, # number of trees is dependent on training data size
                  importance = TRUE, # calculate variable importance
                  tuneGrid = mtryGrid # set mtrys
  )
  
  ###############################################################################################
  ### Step #7: Using signature set derived in #5 to subset Model Selection set.
  ###############################################################################################
  
  myMS <- myMS[,sigProbes2]
  
  ###############################################################################################
  ### Step #8: Application of predictive classifier to MS set.
  ###############################################################################################
  
  RF.pred1 <- predict(RF.obj, myMS, type="prob")
  
  ###############################################################################################
  ### Step #9.1: Storage of dMeLB scores & Model selection set classifications
  ###############################################################################################
  
  ## Ensuring class and diagnosis are recorded in the model selection set results
  class.new = pdata.MS$`Sample Type`[match(row.names(myMS), pdata.MS$`Sample ID`)]
  diag.new = pdata.MS$Histotype[match(row.names(myMS), pdata.MS$`Sample ID`)]
  
  ## Recorded result :: Model Selection classifications
  pred.rslt <- data.frame(sample = row.names(myMS), class = class.new, diagnosis = diag.new, 
                          prob_menin = RF.pred1$MNG)
  
  ## Returning of necessary results: Random Forest Object (RF.obj), MS predictions (pred.rslt),
  ## tissue similarity level (similarity), coeff (Nontumor DMPs), sigProbes (dimension reduced CpGs), sigProbes2 (serum-supervised DMPs)
  if(verbose) return(list(RF.obj, pred.rslt,  similarity, coeff, sigProbes, sigProbes2))
  else return(pred.rslt)
}

###############################################################################################
### Step #9.2: Repetition of Stages #2-8 for 1,000 iterations
###############################################################################################

# Defining unique seed #'s for each iteration -- allows work to be repeatable
seed_1000 = sample(1:999999, 1000)
seeds <- seed_1000[1:1000]

## Start of iterative process 
res.list1 = list() 
start_time <- Sys.time()
for (i in 1:length(seeds)) {
  cat(i, "\n")
  seed.tmp = seeds[i]
  res.tmp = menin_RF(discovery.data = discovery.data, seed1 = seed.tmp, ntree = 1000)
  res.list1[[i]] = res.tmp
}
# Calculations for time-elapsed and average iteration run time :: estimated runtime per 1,000 iterations ~8 hours
end_time <- Sys.time()
total.time <- end_time - start_time
avg.time <- total.time/length(seeds)

###############################################################################################
### Step #9.3: Generating ROC curves for 1,000 iterations
###############################################################################################

# Creating a pseudo-results list so we preserve initial results
res.list <- c(res.list1)
# ROC generation and function definition
rocs <- data.frame(cutoff = rep(NA, length(res.list)),
                   accuracy = rep(NA, length(res.list)),
                   sensitivity = rep(NA, length(res.list)),
                   specificity = rep(NA, length(res.list)))
for (i in seq_along(res.list)) {
  
  ## Modifying Meningioma & Non
  res.list[[i]] <- res.list[[i]] %>%
    mutate(class = if_else(class %like% "Meningioma", "Meningioma", "Control"))
  
  ## Calculations of optimal cutoffs
  roc_data <- roc(res.list[[i]]$class, res.list[[i]]$prob_menin)
  optimal_cutoff <- coords(roc_data, 'best', ret = 'threshold')
  if(length(optimal_cutoff$threshold)>1){
    optimal_cutoff$threshold <- optimal_cutoff[1,]
  }
  ## Sensitivity and specificity calculation
  sensitivity <- roc_data$sensitivities[which.min(abs(roc_data$specificities - 0.95))]
  specificity <- roc_data$specificities[which.min(abs(roc_data$sensitivities - sensitivity))]
  
  ## Accuracy calculation
  res.list[[i]]$cutoff <- optimal_cutoff$threshold[1]
  classified <- ifelse(res.list[[i]]$prob_menin>=res.list[[i]]$cutoff, "Meningioma", "Control")
  accuracy <- sum(classified == res.list[[i]]$class)/length(classified)
  
  ## Storage of the necessary parameters
  rocs$cutoff[i] <- optimal_cutoff$threshold
  rocs$accuracy[i] <- accuracy
  rocs$sensitivity[i] <- sensitivity
  rocs$specificity[i] <- specificity
}

## Model slimming:: isolating models with acceptable performances 
rocs <- rownames_to_column(rocs, var = 'model')
rocs <- rocs[rocs$sensitivity > 0 & rocs$specificity > 0,]

## Optimal models
high_param_models <- rocs %>%
  filter(accuracy >= 0.9 & sensitivity >= 0.8 & specificity >= 0.8)

###############################################################################################
### Step #10: Selection of the best performing model and validation across Independent Samples (n=122)
###############################################################################################

## Best models
get_optimal_results <- function(rocs){
  rocs <- rocs %>% rownames_to_column('Model')
  best_accuracy <- rocs[which.max(rocs$accuracy),]
  best_sensitivity <- rocs[which.max(rocs$sensitivity),]
  best_specificity <- rocs[which.max(rocs$specificity),]
  best_results <- rbind(best_accuracy, best_sensitivity, best_specificity)
  row.names(best_results) <- c('ACC', 'SE', 'SP')
  return(best_results)
}
best_results <- get_optimal_results(rocs)
model_number <- as.numeric(best_results$Model[1]) # Default selection: highest accuracy; this can be changed

## Retrieval of selected model
RF.obj <- menin_RF(discovery.data = discovery.data, seed1 = seeds[model_number], ntree = 1000, verbose = TRUE)

#Specify model parameter - cutoff
cutoff <- rocs$cutoff[rocs$model == model_number]

########################################################################################################################
### Step #10: Selection of best performing model and validation of chosen model across independent samples (n=122).
########################################################################################################################

#### Prediction using validation sources generated in step #1
valid.preds <- predict(RF.obj[[1]], t(validation.beta), type = 'prob')

## Marking of correctly predicted validation samples using custom function 
valid_correct <- function(dt.res){
  # Default value: Incorrect
  dt.res$correct <- 0
  # If sample is a Meningioma and score is >= the chosen cutoff, the sample is marked as correct
  dt.res$correct[dt.res$`Sample Type` %like% "Meningioma" & dt.res$MNG >= cutoff] <- 1
  # If sample is a Non-Meningioma and score is < the chosen cutoff, the sample is marked as correct
  dt.res$correct[dt.res$`Sample Type` %notlike% "Meningioma" & dt.res$MNG < cutoff] <- 1
  # Return data
  return(dt.res)
}
# Final liquid biopsy based validation - 122 samples; Column 'correct' = binary Yes (1) /No (0)
valid <- valid_correct(merge(valid.preds, validation, by.x = 0, by.y = 1))