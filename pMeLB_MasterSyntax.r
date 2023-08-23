#######################################################################################
## Prognostic-Meningioma epigenetic Liquid Biopsy (p-MeLB) Construction & Validation ##
#######################################################################################
# All IDATs may be downloaded using :: DOI: 10.17632/zrc982rvjm.2
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
## Base Wilcoxon-test function: returns p-values across whole genome
my.wilcox.test.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error"))
    return(NA)
  else
    return(obj$p.value)
}
## Marking of correctly predicted validation samples using custom function 
add_correct <- function(dt.res){
  dt.res$correct <- 0
  dt.res$correct[dt.res$class == "CR" & dt.res$prob_CR >= 0.5] <- 1
  dt.res$correct[dt.res$class != "CR" & dt.res$prob_CR < 0.5] <- 1
  
  return(dt.res)
}

### Preprocessing of internal meningioma tissue samples ______
RGset <- read.metharray.exp(base = 'Tissue IDATs/', force = TRUE, recursive = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
Mset <- preprocessIllumina(RGset)
pan.mm <- data.frame(getBeta(Mset, type = 'EPIC'))
colnames(pan.mm) <- sub('X', '', colnames(pan.mm))

### Unmasking procedures and NA-value removal of entire genome _____
## Function built with consideration of EPIC hg38 Manifest
unmask <- function(df){
  #Loading EPIC manifest :: available for download at https://zwdzwd.github.io/InfiniumAnnotation 
  probe.anno <- read.table("EPIC.hg38.manifest.tsv.gz", sep = "\t", header = T)
  #Removal of masked probes
  probe.anno <- probe.anno[!probe.anno$MASK_general & probe.anno$CpG_chrm %in% paste0("chr", 1:22),]
  # Subsetting of probes found after unmasking
  df <- df[rownames(df) %in% probe.anno$Probe_ID,]
  #Removal of any NA values
  df <- df[rowSums(is.na(df)) == 0,]
  
}
pan.mm <- unmask(pan.mm) #741k CpGs x 103 samples
# Retrieval of tissue clinical data
pd <- read_excel('SupplementalData1.xlsx', sheet = 'Tissue Information')
# Creation of binary classification
pd$Class_Bi <- ifelse(pd$Category == 'Confirmed Recurrence', 'CR', 'Others')
# Reformatting column names of tissue methylation matrix
pan.mm <- pan.mm[,pd$`Sentrix ID`]
colnames(pan.mm) <- pd$`Sample ID`

### Optional Garbage collection: free up unused memory
# gc()

##### Preprocessing of Quality Checked IDATs into methylation matrix.
lb.pd <- read_excel('Revisions/Revisions3/Clin_info.xlsx', sheet = 'MNG Specific Features')
sentrix <- read_excel('Revisions/Revisions3/Clin_info.xlsx', sheet = 'Whole Cohort Clinical Features')
lb.pd <- merge(sentrix[,c(1, 2)], lb.pd, by = 1)

RGset <- read.metharray.exp(base = 'LB IDATs/', force = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
Mset <- preprocessIllumina(RGset)
lb.beta <- data.frame(getBeta(Mset, type = 'EPIC'))

# Reformatting column names of the methylation matrix to sample ID's
colnames(lb.beta) <- sub('X', '', colnames(lb.beta))
lb.beta <- unmask(lb.beta[,lb.pd$`Sentrix ID`]) #739k CpGs x 63 samples
colnames(lb.beta) <- lb.pd$`Sample ID`

# Construction of iterative function
rr_RF <- function(training.data, seed1, ntree = 1000, verbose=FALSE, Ncores=12){
  
  ###########################################################################
  ### Step #1: Machine-driven randomization of tissue samples into Training (66.6%) & Independent Validation (33.4%).
  ###########################################################################
  
  # Randomization of clinical data into separate cohorts
  training <- pd %>%
    group_by(Category) %>% sample_frac(0.7) %>% ungroup()
  validation <- pd[pd$`Sample ID` %notin% training$`Sample ID`,]
  
  # Construction of training data methylation matrix & training data
  train.beta <- pan.mm[,training$`Sample ID`] # 741k CpGs x  34 samples
  train.data <- data.frame(t(train.beta))
  train.data$class <- training$Class_Bi
  
  # Construction of validation data methylation matrix
  valid.beta <- pan.mm[,validation$`Sample ID`] # 741k CpGs x 16 samples
  
  ###########################################################################
  ### Step #2: Series of machine-driven genome-wide supervised analyses w/ resampling of CNR
  ### training samples to define signature set DMPs.
  ###########################################################################
  
  # Set your seed so your work is repeatable
  set.seed(seed1)
  
  ###########################################################################
  ### Step #2a: Disimilarity: Tissue groups - supervised analysis
  ###########################################################################
  
  #Supervised Analyses: Group 1 - Confirmed Recurrence; Group 2- Confirmed Non-Recurrence
  group1 <- training$`Sample ID`[training$Class_Bi == "CR"]
  group2 <- training$`Sample ID`[training$Class_Bi == "Others"]

  # Creation of temporary beta-matrix for conducting statistical tests
  beta.tmp <- train.beta
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(beta.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
    return(zz)
  })
  # P-value FDR adjustment
  p.adj <- p.adjust(p, method = "fdr")
  
  # Assignment of results to methylation matrix
  beta.tmp$tiss.p <- p.adj
  # Subsetting of significant probes according to pre-assigned N_probe value
  sigProbes = names(sort(p.adj))[1:N_probe]
  
  # Recording of tissue probe values -- mean methylation differences and recording of methylation changes (hyper- or hypo-methylated).
  tiss.probes <- beta.tmp[sigProbes,]
  # Mean methylation difference calculations 
  tiss.probes$mean.group1 <- apply(tiss.probes[, group1], 1, mean, na.rm = T)
  tiss.probes$mean.group2 <- apply(tiss.probes[, group2], 1, mean, na.rm = T)
  tiss.probes$diff.mean <- (tiss.probes$mean.group1 - tiss.probes$mean.group2)
  # Labeling of probe behavior 
  tiss.probes$category <- "hyper"
  tiss.probes[tiss.probes$diff.mean<0, "category"] <- "hypo"
  # Creation of probe-level description data frame
  tiss.probes <- data.frame('probe' = rownames(tiss.probes), 'category' = tiss.probes$category,
                            'tiss.sig' = tiss.probes$tiss.p)
  
  ###########################################################################
  ### Step #2b & 3: Disimilarity: Serum k-mean clusters - supervised analysis & repetition
  ### for each of 4 clusters resulting in four unique DMP sets.
  ###########################################################################
  
  #_____________________ Round 1 _______________________#
  #Supervised Analyses: Group 1 - k1; Group 2- k2, k3, k4
  group1 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k1"]
  group2 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` != "k1"]
  # Creation of temporary beta-matrix for conducting statistical tests
  lb.tmp <- lb.beta[sigProbes, c(group1, group2)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(lb.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are significant at 95% confidence
  k1.probes = names(p[p<0.05])
  # Creation of probe-level description data frame for k1-specific probes
  k1.probes <- unmask(lb.tmp[k1.probes,])
  # Mean methylation difference calculations 
  k1.probes$mean.group1 <- apply(k1.probes[, group1], 1, mean, na.rm = T)
  k1.probes$mean.group2 <- apply(k1.probes[, group2], 1, mean, na.rm = T)
  k1.probes$diff.mean <- (k1.probes$mean.group1 - k1.probes$mean.group2)
  # Labeling of probe behavior 
  k1.probes$category <- "hyper"
  k1.probes[k1.probes$diff.mean<0, "category"] <- "hypo"
  # Creation of probe-level description data frame
  k1.probes <- data.frame('probe' = rownames(k1.probes), 'category' = k1.probes$category,
                          'cluster' = "k1")
  
  #_____________________ Round #2 _______________________#
  #Supervised Analyses: Group 1 - k2; Group 2- k1, k3, k4
  group1 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k2"]
  group2 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` != "k2"]
  # Creation of temporary beta-matrix for conducting statistical tests
  lb.tmp <- lb.beta[sigProbes, c(group1, group2)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(lb.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are significant at 95% confidence
  k2.probes = names(p[p<0.05])
  # Creation of probe-level description data frame for k2-specific probes
  k2.probes <- unmask(lb.tmp[k2.probes,])
  # Mean methylation difference calculations 
  k2.probes$mean.group1 <- apply(k2.probes[, group1], 1, mean, na.rm = T)
  k2.probes$mean.group2 <- apply(k2.probes[, group2], 1, mean, na.rm = T)
  k2.probes$diff.mean <- (k2.probes$mean.group1 - k2.probes$mean.group2)
  # Labeling of probe behavior 
  k2.probes$category <- "hyper"
  k2.probes[k2.probes$diff.mean<0, "category"] <- "hypo"
  # Creation of probe-level description data frame
  k2.probes <- data.frame('probe' = rownames(k2.probes), 'category' = k2.probes$category,
                          'cluster' = "k2")
  
  
  #_____________________ Round #3 _______________________#
  #Supervised Analyses: Group 1 - k3; Group 2- k1, k2, k4
  group1 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k3"]
  group2 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` != "k3"]
  # Creation of temporary beta-matrix for conducting statistical tests
  lb.tmp <- lb.beta[sigProbes, c(group1, group2)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(lb.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are significant at 95% confidence
  k3.probes = names(p[p<0.05])
  # Creation of probe-level description data frame for k3-specific probes
  k3.probes <- unmask(lb.tmp[k3.probes,])
  # Mean methylation difference calculations 
  k3.probes$mean.group1 <- apply(k3.probes[, group1], 1, mean, na.rm = T)
  k3.probes$mean.group2 <- apply(k3.probes[, group2], 1, mean, na.rm = T)
  k3.probes$diff.mean <- (k3.probes$mean.group1 - k3.probes$mean.group2)
  # Labeling of probe behavior 
  k3.probes$category <- "hyper"
  k3.probes[k3.probes$diff.mean<0, "category"] <- "hypo"
  # Creation of probe-level description data frame
  k3.probes <- data.frame('probe' = rownames(k3.probes), 'category' = k3.probes$category,
                          'cluster' = "k3")
  
  #_____________________ Round #4 _______________________#
  #Supervised Analyses: Group 1 - k4; Group 2- k1, k2, k3
  group1 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k4"]
  group2 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` != "k4"]
  # Creation of temporary beta-matrix for conducting statistical tests
  lb.tmp <- lb.beta[sigProbes, c(group1, group2)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(lb.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[group1]), as.matrix(x[group2]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are significant at 95% confidence
  k4.probes = names(p[p<0.05])
  # Creation of probe-level description data frame for k4-specific probes
  k4.probes <- unmask(lb.tmp[k4.probes,])
  # Mean methylation difference calculations 
  k4.probes$mean.group1 <- apply(k4.probes[, group1], 1, mean, na.rm = T)
  k4.probes$mean.group2 <- apply(k4.probes[, group2], 1, mean, na.rm = T)
  k4.probes$diff.mean <- (k4.probes$mean.group1 - k4.probes$mean.group2)
  # Labeling of probe behavior 
  k4.probes$category <- "hyper"
  k4.probes[k4.probes$diff.mean<0, "category"] <- "hypo"
  # Creation of probe-level description data frame
  k4.probes <- data.frame('probe' = rownames(k4.probes), 'category' = k4.probes$category,
                          'cluster' = "k4")
  
  ### Creation of serum-specific probe-level descriptions: all 4 rounds of analyses
  serum.probes <- rbind(k1.probes, k2.probes, k3.probes, k4.probes)
  ### Merging of serum probe and tissue probe information data.frames
  serum.probes<- merge(serum.probes, tiss.probes, by = 1)
  ### Labeling of methylation changes across sample mediums: tissue and serum
  serum.probes$methyl <- "dissimilar"
  serum.probes[serum.probes$category.x ==  serum.probes$category.y, 'methyl'] <- "similar"
  
  ###########################################################################
  ### Step #2c & 3:Translation: tissue & serum & repetition
  ### for each of 4 clusters resulting in four unique DMP sets.
  ###########################################################################
  
  ### Creating sample ID vectors based on k-means cluster
  k1 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k1"]
  k2 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k2"]
  k3 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k3"]
  k4 <- lb.pd$`Sample ID`[lb.pd$`K-means Cluster` == "k4"]
  
  ### Creating sample ID vectors based on confirmed follow-up information, with active resampling
  CR <- training$`Sample ID`[training$Class_Bi == "CR"]
  CR <- sample(CR, 30, replace = TRUE)
  CNR <- training$`Sample ID`[training$Class_Bi == "Others"]
  CNR <- sample(CNR, 20, replace = TRUE) 
  
  # Creating liquid biopsy and tumor tissue combined methylation matrix
  total.beta <- row_to_column(merge(lb.beta, train.beta, by = 0))
  
  #_____________________ Round #1 _______________________#
  
  # Creation of temporary matrix with groups of interest across k1-specific similarly methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k1" & serum.probes$methyl == "similar"], c(k1, CR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k1]), as.matrix(x[CR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k1.sim = names(p[p>0.5])
  # Assurance function - if there are no similar probes, grab the top 10 most similar between tissue and serum
  missing.probes <- function(df){
    if (is_empty(df)) { 
      df <- names(sort(-p))[1:10]
    }else{ 
      return(df)
    }}
  k1.sim <-missing.probes(k1.sim)
  
  # Creation of temporary matrix with groups of interest across k1-specific differentially methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k1" & serum.probes$methyl == "dissimilar"], c(k1, CNR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k1]), as.matrix(x[CNR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k1.dis = names(p[p>0.5])
  
  # Assurance function - if there are no significant probes, grab the top 10 most similar
  k1.dis <-missing.probes(k1.dis)
  k1.total <- c(k1.sim, k1.dis)

  #_____________________ Round #2 _______________________#
  
  # Creation of temporary matrix with groups of interest across k2-specific similarly methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k2" & serum.probes$methyl == "similar"], c(k2, CR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k2]), as.matrix(x[CR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k2.sim = names(p[p>0.5])
  # Assurance function - if there are no similar probes, grab the top 10 most similar between tissue and serum
  missing.probes <- function(df){
    if (is_empty(df)) { 
      df <- names(sort(-p))[1:10]
    }else{ 
      return(df)
    }}
  k2.sim <-missing.probes(k2.sim)
  
  # Creation of temporary matrix with groups of interest across k2-specific differentially methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k2" & serum.probes$methyl == "dissimilar"], c(k2, CNR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k2]), as.matrix(x[CNR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k2.dis = names(p[p>0.5])
  
  # Assurance function - if there are no significant probes, grab the top 10 most similar
  k2.dis <-missing.probes(k2.dis)
  k2.total <- c(k2.sim, k2.dis)
  
  #_____________________ Round #3 _______________________#
  
  # Creation of temporary matrix with groups of interest across k3-specific similarly methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k3" & serum.probes$methyl == "similar"], c(k3, CR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k3]), as.matrix(x[CR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k3.sim = names(p[p>0.5])
  # Assurance function - if there are no similar probes, grab the top 10 most similar between tissue and serum
  missing.probes <- function(df){
    if (is_empty(df)) { 
      df <- names(sort(-p))[1:10]
    }else{ 
      return(df)
    }}
  k3.sim <-missing.probes(k3.sim)
  
  # Creation of temporary matrix with groups of interest across k3-specific differentially methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k3" & serum.probes$methyl == "dissimilar"], c(k3, CNR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k3]), as.matrix(x[CNR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k3.dis = names(p[p>0.5])
  
  # Assurance function - if there are no significant probes, grab the top 10 most similar
  k3.dis <-missing.probes(k3.dis)
  k3.total <- c(k3.sim, k3.dis)
  
  #_____________________ Round #4 _______________________#
  
  # Creation of temporary matrix with groups of interest across k4-specific similarly methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k4" & serum.probes$methyl == "similar"], c(k4, CR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k4]), as.matrix(x[CR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k4.sim = names(p[p>0.5])
  # Assurance function - if there are no similar probes, grab the top 10 most similar between tissue and serum
  missing.probes <- function(df){
    if (is_empty(df)) { 
      df <- names(sort(-p))[1:10]
    }else{ 
      return(df)
    }}
  k4.sim <-missing.probes(k4.sim)
  
  # Creation of temporary matrix with groups of interest across k4-specific differentially methylated probes
  total.tmp <- total.beta[serum.probes$probe[serum.probes$cluster == "k4" & serum.probes$methyl == "dissimilar"], c(k4, CNR)]
  # Wilcoxon rank sum: P-value calculation across dimension reduced cpgs
  p <- apply(total.tmp, 1, function(x) {
    zz <- my.wilcox.test.p.value(as.matrix(x[k4]), as.matrix(x[CNR]), na.action = na.omit)
    return(zz)
  })
  # Subsetting of p-values which are not significant at 95% confidence
  k4.dis = names(p[p>0.5])
  
  # Assurance function - if there are no significant probes, grab the top 10 most similar
  k4.dis <-missing.probes(k4.dis)
  k4.total <- c(k4.sim, k4.dis)
  
  ###########################################################################
  ### Step #4:Selection of the top 25% DMPs w/ highest dissimilarity between
  ### CR & CNR tumor tissue groups.
  ###########################################################################
  
  ### Combining results :: similar and dissimilar methylated probes - similar across serum and tissue
  total.sim <- c(k1.sim, k2.sim, k3.sim, k4.sim)
  total.dis <- c(k1.dis, k2.dis, k3.dis, k4.dis)
  total.probes <- data.frame("Probe.ID" = c(k1.total, k2.total, k3.total, k4.total), 
                             "Group" = c(rep("k1", length(k1.total)), rep("k2", length(k2.total)),
                                         rep("k3", length(k3.total)), rep("k4", length(k4.total))))
  ### Labeling probes
  total.probes[total.probes$Probe.ID %in% total.sim, "category"] <- "sim"
  total.probes[total.probes$Probe.ID %in% total.dis, "category"] <- "dis"
  ### Removing any duplicates
  total.probes <- total.probes[!duplicated(total.probes$Probe.ID),]
  ### Retrieving tissue-based information (results of Step #2a).
  total.probes <- merge(total.probes, tiss.probes[,c(1,3)], by = 1)
  ### Sampling probes based on tissue-oriented group significances (CR vs CNR)
  sigProbes <- total.probes$Probe.ID[total.probes$tiss.sig<quantile(total.probes$tiss.sig)[2]]
  
  ### Assurance function :: if # of signatures is not sufficient for RF construction, run random sampling
  min.set <- function(df) {
    # If statement: if there are less than 16 CpGs
    if (length(df)<16) { 
      # Sampling variable: combination of Group and Category
      total.probes$sampling <-paste0(total.probes$Group, total.probes$category)
      # Sampling of probes
      total.probes <- data.table(total.probes)[,.SD[sample(.N, min(3,.N))],by=sampling]
      # Retrieval of randomized signatures
      df <- total.probes$Probe.ID
    } else{ 
      return(df)
    }}
  # Run assurance function
  sigProbes <- min.set(sigProbes)
  
  ###############################################################################################
  ### Step #5: Generation of predictive Random Forest classifier.
  ###############################################################################################
  
  # Subsetting significant probes across the training data set
  myTrain <- train.data[, c(sigProbes, "class")]
  
  # Setting up k-fold cross validation
  fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv", number = 10,
    ## repeated five times
    repeats = 5)
  
  # Mtry value setting -- defining minimum number of signatures required for decision trees
  median.mtry <- floor(sqrt(ncol(myTrain)-1))
  mtryVals <- seq(median.mtry-3, median.mtry+3, by=1)
  mtryGrid <- data.frame(.mtry=mtryVals)
  
  # Confirm seed again
  set.seed(seed1)
  
  # Set number of computer cores - this can be altered
  registerDoMC(cores = Ncores)
  
  # Construction of actual random forest object 
  RF.obj <- train(class ~ ., # variable to be trained on
                  data = myTrain, # Data we are using
                  method = "rf", # Method we are using
                  trControl = trainControl(method = "oob"), # How we validate
                  ntree = ntree, # number of trees is dependent on training data size
                  importance = TRUE, # calculate variable importance
                  tuneGrid = mtryGrid # set mtrys
  )
  
  ###############################################################################################
  ### Step #6.1: Storage of pMeLB scores across training set & OOB scores.
  ###############################################################################################
  
  train.results <- data.frame(RF.obj$finalModel$votes)
  train.results <- merge(train.results, pd, by.x = 0, by.y = 1)
  oob <-RF.obj$finalModel$err.rate[1000,1]
  pred.rslt <- data.frame(sample = train.results$Row.names, class = train.results$Class_Bi,  prob_CR = train.results$CR, 
                          oob = oob, row.names = NULL)
  
  ## Returning of necessary results: Random Forest Object (RF.obj), training predictions (pred.rslt),
  ## probe information from sequential supervised analyses (total.probes), validation clinical (validation) and methylation matrix (valid.beta).
  if(verbose) return(list(RF.obj, pred.rslt, total.probes, validation, valid.beta))
  else return(pred.rslt)
  
}

###############################################################################################
### Step #9.2: Repetition of Stages #1-5 for 1,000 iterations
###############################################################################################

# Defining hyperparamters and seeds
N_probe <- 500
seed_1000 = sample(1:999999, 1000)
seeds <- seed_1000[1:1000]

## Start of iterative process 
res.list1 = list() 
start_time <- Sys.time()
for (i in 1:length(seeds)) {
  cat(i, "\n")
  seed.tmp = seeds[i]
  res.tmp = rr_RF(seed1 = seed.tmp, ntree = 1000)
  res.list1[[i]] = res.tmp
}

# Calculations for time-elapsed and average iteration run time :: estimated runtime per 1,000 iterations ~ 64.5 hours
end_time <- Sys.time()
total.time <- end_time - start_time
avg.time <- total.time/length(seeds)

###############################################################################################
### Step #7: Selection of pMeLB score & classifier with the smalles Out of Bag (OOB) Error.
###############################################################################################

# Retrieval of results
res.list = list()
for (i in 1:1000){
  res.list[[i]] <- res.list1[[i]]
}

#######################################################################################################################
### Step #7: Selection of the pMeLB score & classifier with the smallest OOB Error across training set CR & CNR cases.
#######################################################################################################################

### Checking oob
df.oob <- (lapply(res.list1,function(i) i[,4]) %>% do.call(rbind,.) %>% data.frame)[,1:2]
head(df.oob[order(df.oob$X1),], 15) # Visualizing the out of bag errors for the first 15 models
# Retrieving the model with the minimal OOB Error
RF.obj <- rr_RF(seed1 = seed_1000[as.numeric(rownames(df.oob)[order(df.oob$X1)][[1]])], ntree = 1000, verbose = TRUE)
# Retrieval of other necessary items from the chosen Random Forest
RF.obj1 <- RF.obj[[1]]
train.results <- add_correct(RF.obj[[2]])
validation <- RF.obj[[4]]
valid.beta <- RF.obj[[5]]

#######################################################################################################################
### Step #7.5: Cohort Addition: retrospectively collected Independent Validation cases.
#######################################################################################################################
### Preprocessing of internal meningioma tissue and serum samples ______
RGset <- read.metharray.exp(base = 'Tissue IDATs/', force = TRUE, recursive = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
Mset <- preprocessIllumina(RGset)
add.mm <- data.frame(getBeta(Mset, type = 'EPIC'))
colnames(add.mm) <- sub('X', '', colnames(add.mm))

RGset <- read.metharray.exp(base = 'LB IDATs/', force = TRUE, recursive = TRUE)
RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
Mset <- preprocessIllumina(RGset)
add.lb <- data.frame(getBeta(Mset, type = 'EPIC'))
colnames(add.lb) <- sub('X', '', colnames(add.lb))

### Retrieval of additional sample information
additional.pd <- read_excel('SupplementalData1.xlsx', sheet = 'pMeLB Validation')
### Subsetting and merging of additional validation samples and original validation
add.mm <- add.mm[RF.obj1$coefnames , colnames(add.mm) %in% additional.pd$`Sentrix ID`]
add.lb <- add.lb[RF.obj1$coefnames , colnames(add.lb) %in% additional.pd$`Sentrix ID`]
add.mm <- row_to_column(merge(add.mm, add.lb, by = 0))
add.mm <- add.mm[,additional.pd$`Sentrix ID`]
colnames(add.mm) <- additional.pd$`Sample ID`

total.valid <- row_to_column(merge(valid.beta[RF.obj1$coefnames, ], add.mm , by= 0))
### Construction of final validation clinical data
valid.pd <- data.frame(`Sample ID` = c(validation$`Sample ID`, additional.pd$`Sample ID`),
                       `Category` = c(validation$Category, additional.pd$Category))

### Note: Bayley et al., 2022 samples may be added here according to clinical notes provided within their publication.
### Publication DOI: 10.1126/sciadv.abm6247; GEO Accession #: GSE189673

#######################################################################################################################
### Step #8: Application of the chosen model to the total Independent Validation Set.
#######################################################################################################################
### Function for labeling sample predictions according to pre-determined cutoff (here, 0.5)
valid_correct <- function(dt.res){
  dt.res$correct <- 0
  dt.res$correct[dt.res$Category == "Confirmed Recurrence" & dt.res$CR >= 0.5] <- 1
  dt.res$correct[dt.res$Category != "Confirmed Recurrence" & dt.res$CR < 0.5] <- 1
  
  return(dt.res)
}

valid.preds <- predict(RF.obj[[1]], t(total.valid), type = 'prob')
### Final performance across total validation 
valid.pd <- valid_correct(merge(valid.pd, valid.preds, by.x = 1, by.y = 0))