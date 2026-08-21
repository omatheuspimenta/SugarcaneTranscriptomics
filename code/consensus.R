## ---------------------------
##
## Script name: consensus.R
##
## Purpose of script: Get the consensus prediction of lncRNA and mRNA from 
##                    the five tools: CPC2, LncADeep2, RNAplonc, RNAsamba,
##                    and FEELnc.
##
## Author: Pimenta-Zanon, M. H.
##
## Date Created: 2026-07-27
##
## Copyright (c) Pimenta-Zanon, M. H., 2026
## Email: omatheuspimenta@outlook.com

## set working directory

setwd("~/Documents/ESALQ/SugarcaneTranscriptomics/code/lncRNA")

## ---------------------------

## load up the packages we will need:

library(readr)
library(dplyr)
library(purrr)
library(pROC)

## ---------------------------

## load up our functions into memory
# Standardize prediction labels to "mRNA" and "lncRNA"
standardize_pred <- function(x) {
  x <- as.character(x)
  case_when(
    # FEELnc
    x == "0" ~ "lncRNA",
    x == "1" ~ "mRNA",
    
    # RNAplonc
    grepl("^1:", x) ~ "lncRNA",
    grepl("^2:", x) ~ "mRNA",
    
    # Other tools
    grepl("^non[- ]?coding$", x, ignore.case = TRUE) ~ "lncRNA",
    grepl("^coding$", x, ignore.case = TRUE) ~ "mRNA",
    grepl("^lncrna$", x, ignore.case = TRUE) ~ "lncRNA",
    grepl("^mrna$", x, ignore.case = TRUE) ~ "mRNA",
    
    TRUE ~ NA_character_
  )
}

# Calculate entropy for a given probability
entropy <- function(x){
  entrop <- 0
  if (x > 0)
    entrop <- (-(x * log2(x)))
  if (is.nan(x) == TRUE)
    entrop <- 0
  return(entrop)
}

# Calculate maximum entropy and the corresponding threshold from a histogram of scores
maxentropy <- function(histogram) {
  totalofpixels <- sum(histogram)
  maximum_entropy <- 0
  threshold <- NULL
  descendinghistogram <- sort(histogram, decreasing = TRUE)
  curveofentropy <- NULL
  for (t in seq_len(length(descendinghistogram) - 1)) {
    P0 <- 0
    P1 <- 0
    for (i in seq_len(t)) {
      P0 <- P0 + descendinghistogram[i] / totalofpixels
    }
    for (i in (t + 1):length(descendinghistogram)) {
      P1 <- P1 + descendinghistogram[i] / totalofpixels
    }
    H0 <- 0
    H1 <- 0
    HT <- 0
    H0 <- H0 + entropy(P0)
    H1 <- H1 + entropy(P1)
    if (is.nan(H1) == TRUE)
      H1 <- 0
    if (is.nan(H0) == TRUE)
      H0 <- 0
    
    HT <- H0 + H1
    curveofentropy <- c(curveofentropy, HT)
    if (HT > maximum_entropy) {
      maximum_entropy <- HT
      threshold <- t
    }
    frequency <- descendinghistogram[threshold]
  }
  list <- list(maximum_entropy, threshold, frequency, curveofentropy)
  names(list) <- c("MaxEntropy",
                   "Threshold",
                   "Frequency",
                   "EntropyCurve")
  return(list)
}

# Create a plot of the entropy curve and mark the threshold
curveofentropy <- function(H, threshold) {
  H <- unlist(H)
  threshold <- unlist(threshold)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(
    bg = "white",
    mar = c(5, 5, 4, 2) + 0.1,
    las = 1,
    cex.lab = 1.2,
    cex.main = 1.3
  )
  
  plot(
    H,
    type = "l",
    lwd = 2.5,
    col = "steelblue4",
    xlab = "Consensus score Distribution",
    ylab = "Sum of entropies",
    main = "Entropy Curve",
    bty = "l"
  )
  
  grid(nx = NULL, ny = NULL, col = "grey90", lty = "dotted")
  
  abline(v = threshold, col = "#D55E00", lwd = 2, lty = 2)
  
  text(
    x = threshold,
    y = max(H, na.rm = TRUE) * 0.95,
    labels = paste0("Threshold = ", threshold),
    pos = 4,
    col = "#D55E00",
    cex = 0.9
  )
  
  x_coord <- grconvertX(1, "npc", "user")
  y_coord <- grconvertY(1, "npc", "user")
  
  par(xpd = NA)
  
  x_coord <- grconvertX(1, "nfc", "user")
  y_coord <- grconvertY(0.8, "nfc", "user")
  
  legend(
    x = x_coord,
    y = y_coord,
    legend = "Selected threshold",
    col = "#D55E00",
    lty = 2,
    lwd = 2,
    bty = "n",
    cex = 0.9,
    xjust = 1,
    yjust = 1
  )
}

## Load lncRNA predictors data
cpc2 <- read.table("../../data/lncRNA/CPC2/cpc2output.txt",
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE)

lncadeepv2 <- read_csv("../../data/lncRNA/LncADeep2/lncadeep2.csv")

rnaplonc <- read.table("../../data/lncRNA/RNAplonc/final_result.txt",
                       header = TRUE,
                       stringsAsFactors = FALSE)

rnasamba <- read.table("../..//data/lncRNA/RNAsamba/rnasamba.tsv",
                          header = TRUE,
                          stringsAsFactors = FALSE)

feelnc <- read.table("../../data/lncRNA/FEELnc/lncRNA_candidates.gtf_RF.txt",
                      header = TRUE,
                      stringsAsFactors = FALSE)

# Standardize column names for merging
cpc2 <- cpc2 %>%
  transmute(
    transcript = ID,
    cpc2_score = coding_probability,
    cpc2_pred = label
  )

lncadeep <- lncadeepv2 %>%
  transmute(
    transcript = SeqID,
    lncadeep_score = Coding_probability,
    lncadeep_pred = Prediction
  )

rnaplonc <- rnaplonc %>%
  transmute(
    transcript = Seq,
    rnaplonc_score = prediction,
    rnaplonc_pred = predicted
  )

rnasamba <- rnasamba %>%
  transmute(
    transcript = sequence_name,
    rnasamba_score = coding_score,
    rnasamba_pred = classification
  )

feelnc <- feelnc %>%
  transmute(
    transcript = name,
    feelnc_score = coding_potential,
    feelnc_pred = label
  )

## Merge all predictors data into a single data frame
merged <- list(
  cpc2,
  lncadeep,
  rnaplonc,
  rnasamba,
  feelnc
) %>%
  reduce(full_join, by = "transcript")

# Standardize prediction labels to "mRNA" and "lncRNA"
merged <- merged %>%
  mutate(across(ends_with("_pred"), standardize_pred))

## Consensus score
score_cols <- grep("_score$", names(merged), value = TRUE)

merged <- merged %>%
  mutate(
    consensus_score = rowMeans(across(all_of(score_cols)), na.rm = TRUE)
  )

## Determine consensus prediction based on majority vote
pred_cols <- grep("_pred$", names(merged), value = TRUE)

merged <- merged %>%
  rowwise() %>%
  mutate(
    n_lncRNA = sum(c_across(all_of(pred_cols)) == "lncRNA", na.rm = TRUE),
    n_mRNA   = sum(c_across(all_of(pred_cols)) == "mRNA", na.rm = TRUE)
  ) %>%
  ungroup()

## We don't use RNAplonc in the high confidence prediction because 
## it has a different scoring system and may not be directly comparable 
## to the other tools. 
## Select the samples where the scores are lower than 0.2 and 
## all votes are in lncRNA
high_conf_lncRNA <- merged %>%
  filter(
    cpc2_pred == "lncRNA",
    lncadeep_pred == "lncRNA",
    rnaplonc_pred == "lncRNA",
    rnasamba_pred == "lncRNA",
    feelnc_pred == "lncRNA",
    cpc2_score < 0.2,
    lncadeep_score < 0.2,
    rnasamba_score < 0.2,
    feelnc_score < 0.2,
  )

# same for mRNA
high_conf_mRNA <- merged %>%
  filter(
    cpc2_pred == "mRNA",
    lncadeep_pred == "mRNA",
    rnaplonc_pred == "mRNA",
    rnasamba_pred == "mRNA",
    feelnc_pred == "mRNA",
    cpc2_score > 0.8,
    lncadeep_score > 0.8,
    rnasamba_score > 0.8,
    feelnc_score > 0.8,
  )

merged <- merged %>%
  mutate(
    high_confidence = case_when(
      cpc2_pred == "lncRNA" & lncadeep_pred == "lncRNA" &
        rnasamba_pred == "lncRNA" & feelnc_pred == "lncRNA" &
        cpc2_score < 0.2 & lncadeep_score < 0.2 & rnasamba_score < 0.2 &
        feelnc_score < 0.2~ TRUE,
      
      cpc2_pred == "mRNA" & lncadeep_pred == "mRNA" &
        rnasamba_pred == "mRNA" & feelnc_pred == "mRNA" &
        cpc2_score > 0.8 & lncadeep_score > 0.8 & rnasamba_score > 0.8 &
        feelnc_score > 0.8 ~ TRUE,
      
      TRUE ~ FALSE
    )
  )

# Include a column for high confidence prediction based on the criteria above
merged <- merged %>%
  mutate(
    high_confidence_pred = case_when(
      cpc2_pred == "lncRNA" & lncadeep_pred == "lncRNA" &
        rnasamba_pred == "lncRNA" & feelnc_pred == "lncRNA" &
        cpc2_score < 0.2 & lncadeep_score < 0.2 & rnasamba_score < 0.2 &
        feelnc_score < 0.2 ~ "lncRNA",
      
      cpc2_pred == "mRNA" & lncadeep_pred == "mRNA" & # rnaplonc_pred == "mRNA" &
        rnasamba_pred == "mRNA" & feelnc_pred == "mRNA" &
        cpc2_score > 0.8 & lncadeep_score > 0.8 & rnasamba_score > 0.8 &
        feelnc_score > 0.8 ~ "mRNA",
    )
  )

## Maximum Entropy threshold approach
maxentropy_result <- maxentropy(merged$consensus_score)
best_thresh <- maxentropy_result$Frequency
curveofentropy(maxentropy_result$EntropyCurve, maxentropy_result$Threshold)

## Determine the final consensus prediction based on threshold on 'consensus_score'
threshold <- best_thresh
merged <- merged %>%
  mutate(
    consensus_pred = case_when(
      consensus_score >= threshold ~ "mRNA",
      consensus_score < threshold ~ "lncRNA",
      TRUE ~ NA_character_
    )
  )

## Create a "vote_pred" based in the majority vote of the predictions from the five tools
merged <- merged %>%
  mutate(
    vote_pred = case_when(
      n_lncRNA > n_mRNA ~ "lncRNA",
      n_mRNA > n_lncRNA ~ "mRNA",
      n_lncRNA == n_mRNA ~ consensus_pred,
      TRUE ~ NA_character_
    )
  )

## Save the final merged data frame with consensus predictions
write_csv(merged, "../../data/lncRNA/consensus_predictions.csv")

## Get the IDs of the high confidence lncRNA and mRNA predictions
# Column to extract the transcript IDs of high confidence lncRNA predictions
col2extract <- "high_confidence_pred" # "consensus_pred"
lncRNA_final_ids <- merged %>%
  filter(.data[[col2extract]] == "lncRNA") %>%
  pull(transcript)

## Save the high confidence lncRNA IDs to a text file
writeLines(lncRNA_final_ids, "../../data/lncRNA_final_ids.txt")
