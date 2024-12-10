#!/usr/bin/env Rscript

### Polygenic Score Catalog - Benchmarking analysis (v1.2)
### Written by Joel T. Gibson (2024)
### Email: jtg29@cam.ac.uk
### Cardiovascular Epidemiology Unit
### Department of Public Health and Primary Care
### Victor Phillip Dahdaleh Heart and Lung Research Institute
### University of Cambridge

start.time <- round(Sys.time(), "secs")


############################## Required Packages ###############################

# This script requires the "data.table", "R.utils", "pROC", "survival", "boot",
# "parallel", "dplyr", "ggplot2" and "r2redux" packages

# Specify directory for installing/loading packages
#.libPaths("<path to user-specified library>")

# Set repository for downloading packages
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install required packages
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")

if (!requireNamespace("R.utils", quietly = TRUE))
  install.packages("R.utils")

if (!requireNamespace("pROC", quietly = TRUE))
  install.packages("pROC")

if (!requireNamespace("survival", quietly = TRUE))
  install.packages("survival")

if (!requireNamespace("boot", quietly = TRUE))
  install.packages("boot")

if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel")

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!requireNamespace("r2redux", quietly = TRUE))
  install.packages("r2redux")

# Load required packages
suppressMessages(library("data.table"))
suppressMessages(library("R.utils"))
suppressMessages(library("pROC"))
suppressMessages(library("boot"))
suppressMessages(library("parallel"))
suppressMessages(library("dplyr"))
suppressMessages(library("r2redux"))
suppressMessages(library("ggplot2"))
suppressMessages(library("survival"))


################### Data, Parameters and Running Requirements ##################

# These options should normally be specified when calling the script from the
# command line but can be specified here if running R interactively. See the
# user manual for a description of each option. Use quotation marks and/or
# vector notation where appropriate

### Data description ###

SAMPLESET <- NULL
PHENOTYPE <- NULL

### Input and output ###

DEM.FILE <- NULL
PGS.FILE <- NULL
POP.FILE <- NULL
OUTDIR <- NA

### Required columns ###

ID <- "IID"
AGE <- "age"
SEX <- "sex"

### Binary phenotype data ###

ANY.CASE <- NA
PREV.CASE <- NA
INC.CASE <- NA
FOLLOW.TIME <- NA

### Quantitative phenotype data ###

TRAIT.VAL <- NA
UNITS <- NA

### Additional data ###

COVS.CAT <- NA
COVS.QUANT <- NA

### Ancestry stratification ###

ANCESTRY <- "all"
STRATIFIED <- FALSE

### Run specifications ###

SCORE.TYPE <- "all"
MIN.CASES <- 50
SKIP <- 0
MAX.PGS <- NA
BOOTSTRAPS <- 200
N.CORES <- 1
RUN.DEMS <- TRUE
RUN.METRICS <- TRUE
RUN.PLOTS <- TRUE


######################### Read Command Line Arguments ##########################

set.arg <- function(arg) {
  # Accept an option-argument character string in the form "--option=arg" and
  # creates a variable 'OPTION' with the value 'arg'
  
  # Extract variable name and value from argument
  arg <- strsplit(arg, "=")[[1]]
  var.name <- toupper(gsub("-", ".", substring(arg[1], first = 3)))
  value <- format.val(arg[2])
  
  do.call("<<-", list(var.name, value))
}

format.val <- function(val) {
  # Coerces a value to its intended data type
  
  if (grepl(",", val, fixed = TRUE)) {
    val <- strsplit(val, ",")[[1]]
  } else if (!is.na(as.logical(val))) {
    val <- as.logical(val)
  } else if (val == "NA") {
    val <- NA
  } else if (suppressWarnings(!is.na(as.numeric(val)))) {
    val <- as.numeric(val)
  }
  
  return(val)
}

# Read and set command line arguments
args <- commandArgs(trailingOnly = TRUE)
args <- lapply(args, set.arg)


############################ Directory Preparation #############################

# Set output directory
if (is.na(OUTDIR)) {
  OUTDIR <- "~"
}

### Create directories for storing output files ###

# Directory names
benchmark.dir <- "pgsc_benchmark_results"
study.dir <- gsub(" ", "_", SAMPLESET)
pheno.dir <- gsub(" ", "_", PHENOTYPE)

results.dir <- gsub(":", "-", start.time)
results.dir <- gsub(" ", "_", results.dir)

# Full paths to output directories
results.path <- file.path(OUTDIR, benchmark.dir, study.dir, pheno.dir,
                          results.dir)

# Create all directories
if (RUN.PLOTS) {
  plots.path <- file.path(results.path, "Plots")
  dir.create(plots.path, recursive = TRUE, showWarnings = FALSE)
} else {
  dir.create(results.path, recursive = TRUE, showWarnings = FALSE)
}

### Create file paths for writing output data ###

# Demographic data
dem.outfile <- "demographics.csv"
dem.outpath <- file.path(results.path, dem.outfile)

# Performance metrics
met.outfile <- "metrics.csv"
met.outpath <- file.path(results.path, met.outfile)

# Run specifications
run.outfile <- "run_specifications.txt"
run.outpath <- file.path(results.path, run.outfile)


########################## Data Preparation Functions ##########################

dem.cols <- function() {
  # Returns a vector specifying which columns to import from the file
  # containing demographic/phenotype data (needed to speed up import)
  
  # Extract column names
  one.row <- data.table::fread(DEM.FILE, nrows = 1, showProgress = FALSE)
  all.cols <- colnames(one.row)
  
  # Specify required columns
  col.classes <- rep("NULL", length(all.cols))
  
  col.classes[match(ID, all.cols)] <- "character"
  col.classes[match(c(SEX, COVS.CAT), all.cols)] <- "factor"
  col.classes[match(c(AGE, ANY.CASE, PREV.CASE, INC.CASE, TRAIT.VAL,
                      FOLLOW.TIME, COVS.QUANT), all.cols)] <- "numeric"
  
  return(col.classes)
}

pgs.cols <- function() {
  # Returns a vector specifying which columns to import from the file
  # containing PGS data (needed to speed up import)
  
  # Extract column names
  one.row <- data.table::fread(PGS.FILE, nrows = 1, showProgress = FALSE)
  all.cols <- colnames(one.row)
  num.cols <- ncol(one.row)
  
  # Specify required columns
  if (num.cols == 8) {
    col.classes <- c(rep("character", 2), "factor", rep("NULL", 5))
  } else {
    col.classes <- c(rep("character", 3), "factor", rep("NULL", 5))
  }
  col.classes[match(score.types, all.cols)] <- "numeric"
  
  return(col.classes)
}

pop.cols <- function() {
  # Returns a vector specifying which columns to import from the file
  # containing population similarity data (needed to speed up import)
  
  # find number of columns in population similarity file
  one.row <- data.table::fread(POP.FILE, nrows = 1, showProgress = FALSE)
  num.cols <- ncol(one.row)
   
  if (num.cols == 23) {
    col.classes <- c(rep("character", 2), rep("numeric", 10), rep("NULL", 7),
                     "character", rep("NULL", 3))
  } else {
    col.classes <- c(rep("character", 3), rep("numeric", 10), rep("NULL", 7),
                     "character", rep("NULL", 3))
  }
  
  return(col.classes)
}

load.dem.data <- function(filename) {
  # Loads and prepares demographic data
  
  cat("Loading demographic data\n")
  
  dem.data <- data.table::fread(
    filename, colClasses = dem.cols(), showProgress = FALSE,
    na.strings = c("", "NA")
  )
  
  # Rename columns
  names(dem.data)[names(dem.data) == ID] <- ID.COL
  names(dem.data)[names(dem.data) == AGE] <- AGE.COL
  names(dem.data)[names(dem.data) == SEX] <- SEX.COL
  
  if (!is.na(ANY.CASE)) {
    names(dem.data)[names(dem.data) == ANY.CASE] <- ANY.CASE.COL
  }
  
  if (!is.na(PREV.CASE)) {
    names(dem.data)[names(dem.data) == PREV.CASE] <- PREVALENT.COL
  }
  
  if (!is.na(INC.CASE)) {
    names(dem.data)[names(dem.data) == INC.CASE] <- INCIDENT.COL
  }
  
  if (!is.na(FOLLOW.TIME)) {
    names(dem.data)[names(dem.data) == FOLLOW.TIME] <- FOLLOW.UP.COL
  }
  
  if (!is.na(TRAIT.VAL)) {
    names(dem.data)[names(dem.data) == TRAIT.VAL] <- TRAIT.VALUE.COL
    
    # Exclude rows with missing quantitative phenotype data
    dem.data <- dem.data[!is.na(dem.data[[TRAIT.VALUE.COL]]),]
  }
  
  # Create a column for 'any.cases' using prevalent and incident data (if there
  # isn't already a column for this)
  if (is.na(ANY.CASE) && !is.na(PREV.CASE) && !is.na(INC.CASE)) {
    dem.data[[ANY.CASE.COL]] <- dem.data[[PREVALENT.COL]]
    dem.data[[ANY.CASE.COL]][!is.na(dem.data[[INCIDENT.COL]]) &
                            dem.data[[INCIDENT.COL]] == 1] <- 1
  }
  
  # Exclude rows with incomplete data in required columns
  req.cols <- match(c(AGE.COL, SEX.COL, COVS.CAT, COVS.QUANT),
                    colnames(dem.data))
  req.cols <- req.cols[!is.na(req.cols)]
  dem.data <- dem.data[complete.cases(dem.data[, ..req.cols]),]
  
  return(dem.data)
}

load.pop.data <- function(filename) {
  # Loads and prepares population similarity data
  
  cat("Loading population similarity data\n")
  
  pop.data <- data.table::fread(
    filename, sep = "\t", colClasses = pop.cols(), showProgress = FALSE,
    na.strings = c("", "NA")
  )
  
  pop.data <- pop.data[pop.data$sampleset != "reference",]
  pop.data <- subset(pop.data, select = -c(sampleset))
  
  # Remove participants with incomplete PC data
  req.cols <- match(PCs, colnames(pop.data))
  pop.data <- pop.data[complete.cases(pop.data[, ..req.cols]),]
  
  # Subset to specified ancestries
  if (ANCESTRY[1] != "all") {
    pop.data <- pop.data[pop.data$MostSimilarPop %in% ANCESTRY,]
  }
  
  return(pop.data)
}

load.pgs.data <- function(filename) {
  # Loads and prepares PGS data
  
  cat("Loading PGS data\n")
  
  pgs.data <- data.table::fread(
    filename, sep = "\t", colClasses = pgs.cols(), showProgress = FALSE,
    na.strings = c("", "NA")
  )
  
  # Exclude rows with no PGS (just check first score type) and rows containing
  # reference data
  pgs.data <- pgs.data[!is.na(pgs.data[[score.types[1]]]) &
                         pgs.data$sampleset != "reference",]
  pgs.data <- subset(pgs.data, select = -c(sampleset))
  
  return(pgs.data)
}

index.pgs.data <- function() {
  # Creates an index data frame containing the first row number of each PGS in
  # the 'pgs.data' data frame
  
  cat("Indexing PGS data\n")
  
  pgs.ids <- levels(pgs.data$PGS)
  starts <- match(pgs.ids, pgs.data$PGS)
  
  # Find the total number of observations for each PGS
  if (length(pgs.ids) > 1) {
    ends <- c(starts[2:length(pgs.ids)], nrow(pgs.data) + 1)
  } else {
    ends <- nrow(pgs.data) + 1
  }
  n.samples <- ends - starts
  
  # Extract PGS ID from full score string
  pgs.ids <- substr(pgs.ids, start = 1, stop = 9)
  
  index <- data.frame(
    pgs.id = pgs.ids,
    start = starts,
    rows = n.samples,
    row.names = NULL
  )
  
  # Add the multi-pgs to the index data frame
  if ("PGS003443" %in% index$pgs.id &&
      "PGS003444" %in% index$pgs.id &&
      "PGS003445" %in% index$pgs.id) {
    multi.pgs <- data.frame(pgs.id = "PGS00344X", start = NA, rows = NA)
    index <- rbind(index, multi.pgs)
  }
  
  # Clean up data
  pgs.data <<- subset(pgs.data, select = -c(PGS))
  
  return(index)
}

get.case.types <- function() {
  # Returns a vector of all case types present in the input data (any cases,
  # prevalent cases, incident cases, or quantitative trait)
  
  input <- c(ANY.CASE, PREV.CASE, INC.CASE, TRAIT.VAL)
  all.types <- c(ANY.CASE.COL, PREVALENT.COL, INCIDENT.COL, TRAIT.VALUE.COL)
  
  case.types <- all.types[!is.na(input)]
  
  # Add 'any cases' if it can be created from the other data
  if (is.na(ANY.CASE) && !is.na(PREV.CASE) && !is.na(INC.CASE)) {
    case.types <- c(ANY.CASE.COL, case.types)
  }
  
  return(case.types)
}

get.dataset <- function(pgs.id, ancestry) {
  # Creates a merged data set for the specified PGS
  
  # Extract scores for the current PGS ID
  start <- index[index$pgs.id == pgs.id, "start"]
  rows <- index[index$pgs.id == pgs.id, "rows"]
  curr.pgs <- pgs.data[start:(start + rows - 1),]
  
  # Merge demographic, PGS and population similarity data
  merged.data <- merge(x = curr.pgs, y = pop.data,
                       by.x = ID.COL, by.y = ID.COL)
  merged.data <- merge(x = dem.data, y = merged.data,
                       by.x = ID.COL, by.y = ID.COL)
  
  # Subset to specified ancestry
  if (STRATIFIED) {
    merged.data <- merged.data[merged.data$MostSimilarPop == ancestry,]
  }
  
  # Scale and center SUM
  if ("SUM" %in% score.types) {
    merged.data[["SUM"]] <- as.numeric(scale(merged.data[["SUM"]]))
  }
  
  # Scale and center quantitative phenotype values
  if (!is.na(TRAIT.VAL)) {
    merged.data[[TRAIT.SCALED.COL]] <- as.numeric(
      scale(merged.data[[TRAIT.VALUE.COL]])
    )
  }
  
  return(merged.data)
}

select.score <- function(data, score.type) {
  # Returns the input data set containing only the specified score type
  
  excluded.scores <- score.types[score.types != score.type]
  excluded.cols <- match(excluded.scores, colnames(data))
  names(data)[names(data) == score.type] <- SCORE.COL
  
  return(data[, -c(..excluded.cols)])
}


##################### Plotting PGS Distributions Functions #####################

plot.pgs.binary <- function(data, case.type, pgs.id, score.type, ancestry) {
  # Plots PGS distributions for the specified PGS and cohort
  
  # Define groups
  data[[case.type]] <- as.factor(data[[case.type]])
  cases <-  status.name(case.type, 1)
  controls <-  status.name(case.type, 0)
  
  # Create headings for plot
  heading <- paste(pgs.id, "distributions:", PHENOTYPE)
  subheading <- SAMPLESET
  if (STRATIFIED || ANCESTRY[1] != "all") {
    subheading <- paste(subheading, ancestry.list(ancestry, par = TRUE))
  }
  subheading <- paste0(subheading, ", ", cohort.name(case.type), " cohort")
  
  # Find x-axis bounds (mean +/- 4 sd)
  score.mean <- mean(data[[score.type]])
  score.sd <- sd(data[[score.type]])
  x.lower <- score.mean - 4 * score.sd
  x.upper <- score.mean + 4 * score.sd
  
  # Create the density plot
  figure <- ggplot(data, aes(x = get(score.type), colour = data[[case.type]])) +
    geom_density() +
    theme_bw() +
    labs(title = heading,
         subtitle = subheading,
         x = paste0("Score (", score.type, ")"),
         y = "Density") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom",
          panel.grid.minor.y = element_blank()) +
    coord_cartesian(xlim = c(x.lower, x.upper)) +
    scale_colour_manual(values = c("blue", "red"),
                        labels = c(controls, cases),
                        guide = guide_legend(reverse = TRUE))
  
  print(figure)
}

plot.pgs.quant <- function(data, pgs.id, score.type, ancestry) {
  # Plots mean trait values stratified by PGS deciles for the specified PGS and
  # cohort (for a quantitative phenotype)
  
  # Assign data to PGS percentile bins
  n.bins <- 10
  data$bin <- as.factor(dplyr::ntile(data[[score.type]], n.bins))
  
  # Create headings for plot
  heading <- paste("Mean", PHENOTYPE, "by", pgs.id, "decile")
  subheading <- SAMPLESET
  if (STRATIFIED || ANCESTRY[1] != "all") {
    subheading <- paste(subheading, ancestry.list(ancestry, par = TRUE))
  }
  
  # Calculate summary data to be plotted for each bin
  mean.pgs <- tapply(data[[score.type]], data$bin, mean)
  mean.trait <- tapply(data[[TRAIT.VALUE.COL]], data$bin, mean)
  ci <- tapply(data[[TRAIT.VALUE.COL]], data$bin, get.ci)
  lower95 <- sapply(ci, function(ci) ci[1])
  upper95 <- sapply(ci, function(ci) ci[2])
  plot.values <- data.frame(mean.pgs, mean.trait, lower95, upper95)
  
  # Find equation of regression line
  model <- lm(paste(TRAIT.VALUE.COL, "~", score.type), data)
  slope <- coef(summary(model))[score.type, "Estimate"]
  x.mid <- mean(data[[score.type]], na.rm = TRUE)
  y.mid <- mean(data[[TRAIT.VALUE.COL]], na.rm = TRUE)
  int <- y.mid - slope * x.mid
  
  # Create y-axis label
  y.lab <- PHENOTYPE
  if (!is.na(UNITS)) {
    y.lab <- paste0(y.lab, " (", UNITS, ")")
  }
  
  # Create the plot
  figure <- ggplot(plot.values, aes(x = mean.pgs, y = mean.trait)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower95, ymax = upper95)) +
    theme_bw() +
    labs(title = heading,
         subtitle = subheading,
         x = paste0("Score (", score.type, ")"),
         y = y.lab) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom",
          panel.grid.minor.y = element_blank()) +
    geom_abline(intercept = int,
                slope = slope,
                linetype = 2)
  
  print(figure)
}

get.ci <- function(values) {
  # Calculates the confidence interval from a vector of numeric values
  
  ci <- t.test(values)$conf.int
  c(ci[1], ci[2])
}

plot.pgs <- function(case.type, pgs.id, data, ancestry) {
  # Plots the score distributions for each PGS ID in the specified cohort
  
  # Only plot scores if the number of cases exceeds the required number
  if (!enough.cases(data, case.type)) {
    return(NULL)
  }
  
  # Exclude rows with missing phenotype data
  data <- data[!is.na(data[[case.type]]),]
  
  # Exclude rows without follow-up time (to be consistent with cohort used to
  # calculate performance metrics)
  if (case.type == INCIDENT.COL && !is.na(FOLLOW.TIME)) {
    data <- data[!is.na(data[[FOLLOW.UP.COL]]),]
  }
  
  for (score.type in score.types) {
    if (case.type == TRAIT.VALUE.COL) {
      # Plot average trait values for quantitative phenotype
      plot.pgs.quant(data, pgs.id, score.type, ancestry)
    } else {
      # Plot density curves for binary phenotype
      plot.pgs.binary(data, case.type, pgs.id, score.type, ancestry)
    }
    
    # Create png file name
    plot.outfile <- paste(pgs.id, score.name(score.type), sep = "_")
    if (STRATIFIED) {
      plot.outfile <- paste(plot.outfile, ancestry, sep = "_")
    }
    if (case.type != TRAIT.VALUE.COL) {
      plot.outfile <- paste(plot.outfile, cohort.name(case.type), sep = "_")
    }
    plot.outfile <- paste0(plot.outfile, ".png")
    
    # Save plot as png image
    plot.outpath <- file.path(plots.path, plot.outfile)
    ggsave(plot.outpath, width = 14, height = 10, units = "cm")
  }
}


####################### Demographic Statistics Functions #######################

get.base.dems <- function(ancestry) {
  # Calculates the baseline demographic statistics (not including PGS means)
  
  cat(" - Calculating demographic statistics\n")
  
  # Calculate demographics for a single cohort (using the data set for the
  # first PGS ID)
  first.dataset <- get.dataset(index$pgs.id[1], ancestry)
  split.dems <- lapply(
    case.types, stratify.cohort, first.dataset, index$pgs.id[1],
    "base.demographics", ancestry
  )
  demographics <- do.call(rbind, split.dems)
}

get.pgs.means <- function(pgs.id, data, ancestry) {
  # Calculates the mean and standard deviation of the specified PGS across all
  # case types and disease statuses
  
  split.means <- lapply(
    case.types, stratify.cohort, data, pgs.id, "pgs.means", ancestry
  )
  means <- do.call(rbind, split.means)
}

stratify.cohort <- function(case.type, data, pgs.id, analysis, ancestry) {
  # Applies the specified analysis to the given case type, stratified by
  # disease status
  
  # Exclude rows with missing data
  data <- data[!is.na(data[[case.type]]),]
  if (case.type == INCIDENT.COL && !is.na(FOLLOW.TIME)) {
    data <- data[!is.na(data[[FOLLOW.UP.COL]]),]
  }
  
  # Find all possible disease statuses for the case type
  if (case.type == TRAIT.VALUE.COL) {
    disease.status <- NA
  } else {
    disease.status <- c("Any", "1", "0")
  }
  
  if (analysis == "base.demographics") {
    split.output <- lapply(
      disease.status, dem.row, case.type, data, ancestry
    )
    output <- do.call(rbind, split.output)
  } else if (analysis == "pgs.means") {
    split.output <- lapply(
      disease.status, calc.pgs.mean, case.type, pgs.id, data
    )
    output <- do.call(rbind, split.output)
  }
  
  return(output)
}

calc.pgs.mean <- function(disease.status, case.type, pgs.id, data) {
  # Calculates mean and standard deviation for a single PGS in the specified
  # cohort and disease status
  
  # Only calculate statistics if the number of cases exceeds the required number
  if (!enough.cases(data, case.type)) {
    skip.stats <- TRUE
  } else {
    skip.stats <- FALSE
  }
  
  # Extract data for specified disease status
  if (!is.na(disease.status) && disease.status != "Any") {
    data <- data[data[[case.type]] == as.numeric(disease.status),]
  }
  
  mean.data <- list()
  
  # Find mean and SD for each score type
  for (score.type in score.types) {
    
    # New column names
    mean.col <- paste(pgs.id, score.name(score.type), "mean", sep = ".")
    sd.col <- paste(pgs.id, score.name(score.type), "sd", sep = ".")
    
    # Calculate the mean and standard deviation of the PGS
    if (skip.stats) {
      mean.data[[mean.col]] <- NA
      mean.data[[sd.col]] <- NA
    } else {
      mean.data[[mean.col]] <- mean(data[[score.type]], na.rm = TRUE)
      mean.data[[sd.col]] <- sd(data[[score.type]], na.rm = TRUE)
    }
  }
  
  as.data.frame(mean.data)
}

dem.row <- function(disease.status, case.type, data, ancestry) {
  # Creates a single row of demographics data using the specified input

  # Only calculate statistics if the number of cases exceeds the required number
  if (!enough.cases(data, case.type)) {
    skip.stats <- TRUE
  } else {
    skip.stats <- FALSE
  }
  
  # Extract status name and data for specified disease status
  if (!is.na(disease.status) && disease.status != "Any") {
    data <- data[data[[case.type]] == as.numeric(disease.status),]
    disease.status <- status.name(case.type, as.numeric(disease.status))
  }

  dem.data <- data.frame(
    sampleset = SAMPLESET,
    phenotype = PHENOTYPE,
    ancestry = ancestry.list(ancestry),
    cohort = cohort.name(case.type),
    status = disease.status,
    n.total = nrow(data),
    n.males = NA,
    age.recruit.mean = NA,
    age.recruit.sd = NA,
    age.event.mean = NA,
    age.event.sd = NA,
    trait.value.mean = NA,
    trait.value.sd = NA,
    trait.value.units = NA,
    row.names = NULL
  )
  
  if (skip.stats) {
    return(dem.data)
  }
  
  dem.data[["n.males"]] <- sum(is.male(data[[SEX.COL]]), na.rm = TRUE)
  dem.data[["age.recruit.mean"]] <- mean(data[[AGE.COL]], na.rm = TRUE)
  dem.data[["age.recruit.sd"]] <- sd(data[[AGE.COL]], na.rm = TRUE)

  # Calculate age at event/censoring for incident cohort
  if (case.type == INCIDENT.COL && !is.na(FOLLOW.TIME)) {
    dem.data[["age.event.mean"]] <- mean(
      data[[AGE.COL]] + data[[FOLLOW.UP.COL]], na.rm = TRUE
    )
    dem.data[["age.event.sd"]] <- sd(
      data[[AGE.COL]] + data[[FOLLOW.UP.COL]], na.rm = TRUE
    )
  }

  # Calculate mean value of quantitative phenotype
  if (case.type == TRAIT.VALUE.COL) {
    dem.data[["trait.value.mean"]] <- mean(
      data[[TRAIT.VALUE.COL]], na.rm = TRUE
    )
    dem.data[["trait.value.sd"]] <- sd(
      data[[TRAIT.VALUE.COL]], na.rm = TRUE
    )
    dem.data[["trait.value.units"]] <- UNITS
  }

  return(dem.data)
}

is.male <- function(sex) {
  # Returns TRUE if the input value represents a male (e.g. "Male", "M", 1)
  
  toupper(substr(sex, start = 1, stop = 1)) == "M" | sex == 1
}


##################### Logistic Regression Model Functions ######################

log.model <- function(case.type, data, baseline = FALSE) {
  # Creates a logistic regression model
  
  # Create formula
  form <- paste(case.type, "~ .")
  if (baseline) {
    form <- paste(form, "-", SCORE.COL)
  }
  
  glm(form, "binomial", data)
}

get.or <- function(model) {
  # Extracts the PGS odds ratio from a logistic regression model
  
  prs.data <- coef(summary(model))[SCORE.COL,]
  beta <- prs.data["Estimate"]
  se <- prs.data["Std. Error"]
  lower95 <- beta - qnorm(0.975) * se
  upper95 <- beta + qnorm(0.975) * se
  p.val <- prs.data["Pr(>|z|)"]
  
  c(exp(beta), exp(lower95), exp(upper95), p.val)
}

bstrap.r2 <- function(data, ind, case.type, sample.prev) {
  # Function used in bootstrapping to calculate R2 on the liability scale
  
  data <- data[ind, ]
  calc.r2.liab(data, case.type, sample.prev)
}

calc.r2.liab <- function(data, case.type, sample.prev) {
  # Calculates R2 on the liability scale for each score type for a single
  # bootstrap replicate
  
  # Calculate R2 values for each score type
  all.r2.liab <- NULL
  for (score.type in score.types) {
    
    # Create formula for regression
    excluded.scores <- score.types[score.types != score.type]
    form <- paste(case.type, "~ .")
    if (length(excluded.scores) > 0) {
      form <- paste(form, "-", paste(excluded.scores, collapse = " - "))
    }
    
    # Calculate R2 on the observed scale
    model <- lm(form, data)
    beta <- coef(summary(model))[score.type, "Estimate"]
    r2.obs <- beta ^ 2
    
    # Calculate R2 on the liability scale
    r2.liab <- transform.r2(sample.prev, r2.obs)
    all.r2.liab <- c(all.r2.liab, r2.liab)
  }
  
  return(all.r2.liab)
}

transform.r2 <- function(sample.prev, r2.obs) {
  # Transforms R2 from the observed scale to the liability scale. Adapted from
  # function by Xilin Jiang
  
  # Find height on normal density curve at point that truncates ratio of cases
  # in upper tail
  thd <- qnorm(1 - sample.prev)
  z <- dnorm(thd)
  
  # Transform R2 to liability scale using eq. (8) from Lee et al. 2012
  r2.obs * sample.prev * (1 - sample.prev) / (z ^ 2)
}

get.r2.liab <- function(data, case.type, boot.data, sample.prev, score.type) {
  # Returns R2 and 95% CI on the liability scale
  
  if (BOOTSTRAPS > 1) {
    # Extract values from bootstrap data
    pos <- match(score.type, score.types)
    r2.liab <- boot.data$t0[pos]
    ci <- boot::boot.ci(boot.data, type = "norm", index = pos)
    lower95 <- ci$normal[2]
    upper95 <- ci$normal[3]
  } else {
    # Calculate R2 on the observed scale
    model <- lm(paste(case.type, "~ ."), data)
    beta <- coef(summary(model))[SCORE.COL, "Estimate"]
    r2.obs <- beta ^ 2
    
    # Calculate R2 on the liability scale
    r2.liab <- transform.r2(sample.prev, r2.obs)
    lower95 <- NA
    upper95 <- NA
  }
  
  p.val <- NA
  
  c(r2.liab, lower95, upper95, p.val)
}

get.auc <- function(roc.object) {
  # Extracts the AUC from an ROC curve
  
  auc <- roc.object$auc
  lower95 <- roc.object$ci[1]
  upper95 <- roc.object$ci[3]
  p.val <- NA
  
  c(auc, lower95, upper95, p.val)
}

get.delta.auc <- function(model.full, model.covs, roc.comparison) {
  # Extracts the change in AUC from a roc.test object and the p-value from the
  # likelihood ratio test between the full and baseline logistic regression
  # models
  
  delta.auc <- roc.comparison$estimate[1] - roc.comparison$estimate[2]
  lower95 <- roc.comparison$conf.int[1]
  upper95 <- roc.comparison$conf.int[2]
  
  # Compare model fit using likelihood ratio Chi-squared test
  p.val <- anova(model.covs, model.full, test = "LRT")[2, "Pr(>Chi)"]
  
  c(delta.auc, lower95, upper95, p.val)
}

log.metrics <- function(case.type, pgs.id, data, model.covs, ancestry) {
  # Calculates performance metrics for the specified PGS ID in the given cohort
  # using logistic regression analysis (requires binary phenotype)
  
  # Only calculate metrics if the number of cases exceeds the required number
  if (!enough.cases(data, case.type)) {
    return(metric.row(case.type, NA, pgs.id, NA, NA, NA, ancestry))
  }
  
  # Select data for modelling
  data <- data[!is.na(data[[case.type]]),]
  data <- subset(data, select = c(case.type, all.covs, score.types))
  
  # Generate bootstrap data for each score type
  # Find ratio of cases in sample
  sample.prev <- mean(data[[case.type]], na.rm = TRUE)
  
  # Bootstrap data to find 95% CI for R2 on liability scale
  if (BOOTSTRAPS > 1) {
    
    # Expand factors and scale columns (for R2 liability calculation)
    data.scaled <- as.data.frame(model.matrix( ~ . - 1, data = data))
    data.scaled <- as.data.frame(scale(data.scaled))
    
    boot.data <- boot::boot(
      data.scaled, bstrap.r2, R = BOOTSTRAPS, case.type = case.type,
      sample.prev = sample.prev
    )
  } else {
    boot.data <- NA
  }
  
  # Calculate performance metrics for each score type
  combined.metrics <- list()
  for (score.type in score.types) {
    
    data.final <- select.score(data, score.type)
    
    if (BOOTSTRAPS <= 1) {
      # Expand factors and scale columns (for R2 liability calculation)
      data.scaled <- as.data.frame(model.matrix( ~ . - 1, data = data.final))
      data.scaled <- as.data.frame(scale(data.scaled))
    }
    
    # Create full logistic regression model using PGS + covariates
    model.full <- log.model(case.type, data.final)
    
    # Extract PGS odds ratio from full model
    or <- get.or(model.full)
    or <- metric.row(case.type, "Full", pgs.id, score.type, "OR", or, ancestry)
    
    # Extract PGS R2 on liability scale from full model
    r2.liab <- get.r2.liab(data.scaled, case.type, boot.data,
                           sample.prev, score.type)
    r2.liab <- metric.row(case.type, "Full", pgs.id, score.type,
                          "R2 liability", r2.liab, ancestry)
    
    # Calculate predicted responses
    pred.full <- predict(model.full, data.final, type = "response")
    pred.covs <- predict(model.covs, data.final, type = "response")
    
    roc.data <- data.frame(
      observed = data.final[[case.type]],
      predicted.full = pred.full,
      predicted.covs = pred.covs
    )
    
    # Create ROC curve for each model
    suppressMessages(
      roc.full <- pROC::roc(
        data = roc.data,
        response = observed,
        predictor = predicted.full,
        auc = TRUE, ci = TRUE
      )
    )
    
    suppressMessages(
      roc.covs <- pROC::roc(
        data = roc.data,
        response = observed,
        predictor = predicted.covs,
        auc = TRUE, ci = TRUE
      )
    )
    
    # Extract AUC for full and baseline (covariates-only) models
    auc.full <- get.auc(roc.full)
    auc.full <- metric.row(case.type, "Full", pgs.id, score.type, "AUC",
                           auc.full, ancestry)
    
    auc.covs <- get.auc(roc.covs)
    auc.covs <- metric.row(case.type, "Baseline", pgs.id, score.type, "AUC",
                           auc.covs, ancestry)
    
    # Calculate difference between AUCs
    roc.comparison <- pROC::roc.test(
      roc1 = roc.full,
      roc2 = roc.covs,
      response = observed,
      predictor1 = predicted.full,
      predictor2 = predicted.covs,
      data = roc.data,
      method = "delong"
    )
    
    delta.auc <- get.delta.auc(model.full, model.covs, roc.comparison)
    delta.auc <- metric.row(case.type, "Full", pgs.id, score.type, "Delta AUC",
                            delta.auc, ancestry)
    
    # Combine all data
    met.data <- list(or, r2.liab, auc.full, auc.covs, delta.auc)
    met.data <- do.call(rbind, met.data)
    combined.metrics[[score.type]] <- met.data
  }
  
  combined.metrics <- do.call(rbind, combined.metrics)
  rownames(combined.metrics) <- NULL
  return(combined.metrics)
}


################### Cox Proportional Hazards Model Functions ###################

cox.model <- function(case.type, data, baseline = FALSE) {
  # Creates a Cox proportional hazards regression model
  
  surv.object <- Surv(time = data[[FOLLOW.UP.COL]], event = data[[case.type]])
  
  # Create formula
  form <- paste("surv.object ~ .", case.type, FOLLOW.UP.COL, sep = "-")
  if (baseline) {
    form <- paste(form, SCORE.COL, sep = "-")
  }
  
  coxph(as.formula(form), data)
}

get.hr <- function(model) {
  # Extracts the PGS hazard ratio from a Cox proportional hazards regression
  # model
  
  prs.data <- summary(model)$conf.int[SCORE.COL,]
  hr <- prs.data["exp(coef)"]
  lower95 <- prs.data["lower .95"]
  upper95 <- prs.data["upper .95"]
  p.val <- coef(summary(model))[SCORE.COL, "Pr(>|z|)"]
  
  c(hr, lower95, upper95, p.val)
}

get.c <- function(model) {
  # Extracts the C-index from a Cox proportional hazards model
  
  c.data <- summary(model)$concordance
  c.index <- c.data["C"]
  se <- c.data["se(C)"]
  lower95 <- c.index - qnorm(0.975) * se
  upper95 <- c.index + qnorm(0.975) * se
  p.val <- NA
  
  c(c.index, lower95, upper95, p.val)
}

get.delta.c <- function(model.full, model.covs, data, case.type) {
  # Estimates the change in C-index between two Cox proportional hazards
  # regression models and the p-value from the likelihood ratio test.
  # Adapted from code by Sam Lambert
  
  # Calculate predicted responses
  pred.full <- -predict(model.full)
  pred.covs <- -predict(model.covs)
  
  # Contrast vector for full-baseline
  contrast.full <- c(1, -1)
  
  # Calculate concordance
  surv.object <- Surv(time = data[[FOLLOW.UP.COL]], event = data[[case.type]])
  concord <- concordance(surv.object ~ pred.full + pred.covs, data)
  
  # Calculate Delta C-index and 95% CI
  delta.c <- as.numeric(contrast.full %*% coef(concord))
  sd <- as.numeric(sqrt(contrast.full %*% vcov(concord) %*% contrast.full))
  lower95 <- delta.c - qnorm(0.975) * sd
  upper95 <- delta.c + qnorm(0.975) * sd
  
  # Compare model fit using likelihood ratio Chi-squared test
  p.val <- anova(model.covs, model.full, test = "LRT")[2, "Pr(>|Chi|)"]
  
  c(delta.c, lower95, upper95, p.val)
}

cox.metrics <- function(case.type, pgs.id, data, model.covs, ancestry) {
  # Calculates performance metrics for the specified PGS ID in the given cohort
  # using Cox proportional hazards regression analysis (requires binary
  # phenotype with time-to-event data)
  
  # Only calculate metrics if the number of cases exceeds the required number
  if (!enough.cases(data, case.type)) {
    return(metric.row(case.type, NA, pgs.id, NA, NA, NA, ancestry))
  }
  
  # Select data for modelling
  data <- data[!is.na(data[[case.type]]) & !is.na(data[[FOLLOW.UP.COL]]),]
  data <- subset(data, select = c(case.type, FOLLOW.UP.COL, all.covs,
                                  score.types))
  
  # Calculate performance metrics for each score type
  combined.metrics <- list()
  for (score.type in score.types) {
    
    data.final <- select.score(data, score.type)
    
    # Create full Cox proportional hazards model using PGS + covariates
    model.full <- cox.model(case.type, data.final)
    
    # Extract PGS hazard ratio from full model
    hr <- get.hr(model.full)
    hr <- metric.row(case.type, "Full", pgs.id, score.type, "HR", hr, ancestry)
    
    # Extract C-index for full and baseline (covariates-only) models
    c.full <- get.c(model.full)
    c.full <- metric.row(case.type, "Full", pgs.id, score.type,
                         "C-index", c.full, ancestry)
    
    c.covs <- get.c(model.covs)
    c.covs <- metric.row(case.type, "Baseline", pgs.id, score.type,
                         "C-index", c.covs, ancestry)
    
    # Calculate difference between C-indices
    delta.c <- get.delta.c(model.full, model.covs, data.final, case.type)
    delta.c <- metric.row(case.type, "Full", pgs.id, score.type,
                          "Delta C-index", delta.c, ancestry)
    
    # Combine all data
    met.data <- list(hr, c.full, c.covs, delta.c)
    met.data <- do.call(rbind, met.data)
    combined.metrics[[score.type]] <- met.data
  }
  
  combined.metrics <- do.call(rbind, combined.metrics)
  rownames(combined.metrics) <- NULL
  return(combined.metrics)
}


###################### Linear Regression Model Functions #######################

linear.model <- function(data, baseline = FALSE) {
  # Creates a linear regression model
  
  # Create formula
  form <- paste(TRAIT.SCALED.COL, "~ .")
  if (baseline) {
    form <- paste(form, "-", SCORE.COL)
  }
  
  lm(form, data)
}

get.beta <- function(model) {
  # Extracts the PGS effect size from a linear regression model
  
  prs.data <- coef(summary(model))[SCORE.COL,]
  beta <- prs.data["Estimate"]
  ci <-  confint(model, SCORE.COL)
  lower95 <- ci[1]
  upper95 <- ci[2]
  p.val <- prs.data["Pr(>|t|)"]
  
  c(beta, lower95, upper95, p.val)
}

get.r2 <- function(model, trait.vals) {
  # Extracts the R2 value from a linear regression model and estimates the
  # 95% confidence intervals
  
  # Calculate R2 statistics
  fitted.vals <- model$fitted.values
  data <- data.frame(trait.vals, fitted.vals)
  suppressWarnings(
    r2.data <- r2redux::r2_var(data, 1, nrow(data))
  )
  
  r2 <- r2.data$rsq
  lower95 <- r2.data$lower_r2
  upper95 <- r2.data$upper_r2
  p.val <- NA
  
  c(r2, lower95, upper95, p.val)
}

get.delta.r2 <- function(model.full, model.covs, trait.vals) {
  # Extracts the change in R2 between two linear regression models and the
  # p-value from the likelihood ratio test
  
  # Calculate delta R2 statistics
  fit.vals.full <- model.full$fitted.values
  fit.vals.covs <- model.covs$fitted.values
  data <- data.frame(trait.vals, fit.vals.full, fit.vals.covs)
  r2.data <- r2redux::r2_diff(data, 1, 2, nrow(data))
  
  delta.r2 <- r2.data$mean_diff
  lower95 <- r2.data$lower_diff
  upper95 <- r2.data$upper_diff
  
  # Compare model fit using likelihood ratio Chi-squared test
  p.val <- anova(model.covs, model.full, test = "LRT")[2, "Pr(>Chi)"]
  
  c(delta.r2, lower95, upper95, p.val)
}

linear.metrics <- function(case.type, pgs.id, data, model.covs, ancestry) {
  # Calculates performance metrics for the specified PGS ID in the given cohort
  # using linear regression analysis (requires quantitative phenotype)
  
  # Only calculate metrics if the number of cases exceeds the required number
  if (!enough.cases(data, case.type)) {
    return(metric.row(case.type, NA, pgs.id, NA, NA, NA, ancestry))
  }
  
  # Select data for modelling
  data <- data[!is.na(data[[case.type]]),]
  data <- subset(data, select = c(TRAIT.SCALED.COL, all.covs, score.types))
  
  # Calculate performance metrics for each score type
  combined.metrics <- list()
  for (score.type in score.types) {
    
    data.final <- select.score(data, score.type)
    
    # Create full linear regression model using PGS + covariates
    model.full <- linear.model(data.final)
    
    # Extract PGS beta from full model
    beta <- get.beta(model.full)
    beta <- metric.row(case.type, "Full", pgs.id, score.type, "Beta", beta,
                       ancestry)
    
    # Extract R2 for full and baseline (covariates-only) models
    r2.full <- get.r2(model.full, data.final[[TRAIT.SCALED.COL]])
    r2.full <- metric.row(case.type, "Full", pgs.id, score.type, "R2", r2.full,
                          ancestry)
    
    r2.covs <- get.r2(model.covs, data.final[[TRAIT.SCALED.COL]])
    r2.covs <- metric.row(case.type, "Baseline", pgs.id, score.type, "R2",
                          r2.covs, ancestry)
    
    # Calculate difference between both R2 values
    delta.r2 <- get.delta.r2(model.full, model.covs,
                             data.final[[TRAIT.SCALED.COL]])
    delta.r2 <- metric.row(case.type, "Full", pgs.id, score.type, "Delta R2",
                           delta.r2, ancestry)
    
    # Combine all data
    met.data <- list(beta, r2.full, r2.covs, delta.r2)
    met.data <- do.call(rbind, met.data)
    combined.metrics[[score.type]] <- met.data
  }
  
  combined.metrics <- do.call(rbind, combined.metrics)
  rownames(combined.metrics) <- NULL
  return(combined.metrics)
}

############################## General Functions ###############################

status.name <- function(case.type, status) {
  # Dictionary for disease status
  
  dict <- list(
    any.cases = c("Non-cases", "Cases"),
    prevalent.cases = c("Controls", "Cases"),
    incident.cases = c("Non-events", "Events")
  )
  
  dict[[case.type]][status + 1]
}

cohort.name <- function(case.type) {
  # Dictionary for cohort name
  
  dict <- list(
    any.cases = "Full",
    prevalent.cases = "Prevalent",
    incident.cases = "Incident",
    trait.value = "Full"
  )
  
  dict[[case.type]]
}

score.name <- function(score.type) {
  # Dictionary for score type abbreviations
  
  dict <- list(
    SUM = "SUM",
    Z_MostSimilarPop = "ZMSP",
    Z_norm1 = "ZN1",
    Z_norm2 = "ZN2"
  )
  
  dict[[score.type]]
}

get.ancestries <- function() {
  # Returns a character vector containing the required the ancestry codes
  
  return(levels(as.factor(pop.data$MostSimilarPop)))
}

ancestry.list <- function(ancestry, par = FALSE) {
  # Returns a string containing the specified ancestry codes (to display in
  # output tables and plots). Optionally surrounded by parentheses (`par`)
  
  if (ANCESTRY[1] == "all" && !STRATIFIED) {
    ancestry <- "All"
  } else {
    ancestry <- paste(ancestry, collapse = ", ")
  }
  
  # Add parentheses
  if (par) {
    ancestry <- paste0("(", ancestry, ")")
  }
  
  return(ancestry)
}

get.covariates <- function() {
  # Returns a character vector containing the column names of all covariates
  # to be included in the regression models
  
  covs <- c(AGE.COL, SEX.COL, PCs, COVS.CAT, COVS.QUANT)
  
  # Exclude sex as a covariate if only one sex is present in input data
  if (length(levels(dem.data[[SEX.COL]])) == 1) {
    covs <- covs[covs != SEX.COL]
  }
  
  return(covs[!is.na(covs)])
}

get.score.types <- function() {
  # Returns a character vector containing the score types to use in analysis
  
  all.types <- c("SUM", "Z_MostSimilarPop", "Z_norm1", "Z_norm2")
  abbrevs <- c("SUM", "ZMSP", "ZN1", "ZN2")
  
  if (is.na(SCORE.TYPE[1]) | SCORE.TYPE[1] == "all") {
    return(all.types)
  }
  
  included <- NULL
  for (score in SCORE.TYPE) {
    if (score %in% all.types) {
      included <- c(included, score)
    } else if (score %in% abbrevs) {
      included <- c(included, all.types[match(score, abbrevs)])
    }
  }
  
  return(included)
}

get.base.models <- function(ancestry) {
  # Calculates the required set of baseline models (covariates without PGS)
  
  if (!RUN.METRICS) {
    return(NA)
  }
  
  cat(" - Calculating baseline models\n")
  
  data <- get.dataset(index$pgs.id[1], ancestry)
  data <- select.score(data, score.types[1])
  
  models <- parallel::mclapply(
    case.types, calc.baseline, data, mc.cores = N.CORES
  )
  
  base.models <- list()
  for (model in models) {
    base.models[[model$case.type]] <- model$model
  }
  
  return(base.models)
}

calc.baseline <- function(case.type, data) {
  # Calculates the baseline model for the specified case type
  
  # Only calculate baseline model if the number of cases exceeds the required
  # number
  if (!enough.cases(data, case.type)) {
    return(list(case.type = case.type, model = NA))
  }
  
  data <- data[!is.na(data[[case.type]]),]
  
  # Calculate baseline model for a binary phenotype
  if (case.type == ANY.CASE.COL || case.type == PREVALENT.COL ||
      (case.type == INCIDENT.COL && is.na(FOLLOW.TIME))) {
    data <- subset(data, select = c(case.type, all.covs, SCORE.COL))
    model <- log.model(case.type, data, baseline = TRUE)
  }
  
  # Calculate baseline model for a binary phenotype with time-to-event data
  else if (case.type == INCIDENT.COL && !is.na(FOLLOW.TIME)) {
    data <- data[!is.na(data[[FOLLOW.UP.COL]]),]
    data <- subset(data, select = c(case.type, FOLLOW.UP.COL, all.covs,
                                    SCORE.COL))
    model <- cox.model(case.type, data, baseline = TRUE)
  }
  
  # Calculate baseline model for a quantitative phenotype
  else if (case.type == TRAIT.VALUE.COL) {
    data <- subset(data, select = c(TRAIT.SCALED.COL, all.covs, SCORE.COL))
    model <- linear.model(data, baseline = TRUE)
  }
  
  list(case.type = case.type, model = model)
}

enough.cases <- function(data, case.type) {
  # Checks if the number of cases (for a binary trait) or the total number of
  # individuals (for a quantitative trait) is greater than the minimum required
  
  if (case.type == TRAIT.VALUE.COL) {
    return(sum(!is.na(data[[case.type]]), na.rm = TRUE) >= MIN.CASES)
  } else {
    return(sum(data[[case.type]], na.rm = TRUE) >= MIN.CASES)
  }
}

metric.row <- function(case.type, model.type, pgs.id, score.type,
                       metric.type, estimates, ancestry) {
  # Creates a single row of metric data using the specified input
  
  # Abbreviate name of score type
  if (!is.na(score.type)) {
    score.type <- score.name(score.type)
  }
  
  data.frame(
    sampleset = SAMPLESET,
    phenotype = PHENOTYPE,
    ancestry = ancestry.list(ancestry),
    pgs.id = pgs.id,
    score.type = score.type,
    cohort = cohort.name(case.type),
    model = model.type,
    metric = metric.type,
    estimate = estimates[1],
    lower.95 = estimates[2],
    upper.95 = estimates[3],
    p.val = estimates[4],
    row.names = NULL
  )
}

analyse.ancestries <- function() {
  # Executes the main analysis pipeline for each ancestry
  
  if (STRATIFIED) {
    demographics <- list()
    metrics <- list()
    
    for (ancestry in all.ancestries) {
      cat("Stratified analysis: ", ancestry, "\n", sep = "")
      
      # Calculate baseline models (covariates only, no PGS)
      base.models <- get.base.models(ancestry)
      
      # Calculate demographics and performance metrics
      results <- analyse.data(SKIP, MAX.PGS, base.models, ancestry)
      demographics[[ancestry]] <- results$demographics
      metrics[[ancestry]] <- results$metrics
    }
    
    # Combine data from each ancestry
    demographics <- do.call(rbind, demographics)
    rownames(demographics) <- NULL
    
    metrics <- do.call(rbind, metrics)
    rownames(metrics) <- NULL
    
  } else {
    if (ANCESTRY[1] == "all") {
      cat("Combined analysis: all ancestries\n")
    } else {
      cat("Combined analysis:", ancestry.list(all.ancestries), "\n")
    }
    
    # Calculate baseline models (covariates only, no PGS)
    base.models <- get.base.models(all.ancestries)
    
    # Calculate demographics and performance metrics
    results <- analyse.data(SKIP, MAX.PGS, base.models, all.ancestries)
    demographics <- results$demographics
    metrics <- results$metrics
  }
  
  list(demographics = demographics, metrics = metrics)
}

analyse.data <- function(skip, max.pgs, base.models, ancestry) {
  # Executes the analysis pipeline using the specified number of PGS IDs
  
  n.pgs <- nrow(index)
  
  # Find the first and last PGS ID indices to analyse
  if (is.na(skip)) {
    start <- 1
  } else {
    start <- skip + 1
  }
  
  if (is.na(max.pgs)) {
    end <- n.pgs
  } else {
    end <- start + max.pgs - 1
    if (end > n.pgs) {
      end <- n.pgs
    }
  }
  
  # Calculate baseline demographic data
  if (RUN.DEMS) {
    base.demographics <- get.base.dems(ancestry)
  } else {
    base.demograhpics <- NA
  }
  
  cat(" - Analysing PGS data\n")
  
  # Calculate PGS means and performance metrics for each PGS ID
  if (start <= end) {
    pgs.ids <- index$pgs.id[start:end]
    raw.output <- parallel::mclapply(
      pgs.ids, analyse.pgs, base.models, ancestry, mc.cores = N.CORES
    )
    full.output <- format.output(raw.output, base.demographics)
  } else {
    full.output <- list(demographics = base.demographics, metrics = NA)
  }
  
  return(full.output)
}

analyse.pgs <- function(next.pgs, base.models, ancestry) {
  # Executes the demographics/performance metrics calculations and plotting
  # for the specified PGS ID. This is the main analysis loop
  
  if (next.pgs == "PGS00344X") {
    dataset <- multi.pgs.dataset(ancestry)
  } else {
    dataset <- get.dataset(next.pgs, ancestry)
  }
  
  # Calculate demographic data (PGS means)
  if (RUN.DEMS) {
    pgs.means <- get.pgs.means(next.pgs, dataset, ancestry)
  } else {
    pgs.means <- NA
  }
  
  # Plot PGS distributions
  if (RUN.PLOTS) {
    lapply(case.types, plot.pgs, next.pgs, dataset, ancestry)
  }
  
  # Calculate performance metrics
  if (RUN.METRICS) {
    split.metrics <- lapply(
      case.types, run.metrics, next.pgs, dataset, base.models, ancestry
    )
    metrics <- do.call(rbind, split.metrics)
  } else {
    metrics <- NA
  }
  
  list(metrics = metrics, means = pgs.means)
}

run.metrics <- function(case.type, next.pgs, data, base.models, ancestry) {
  # Executes the relevant analyses based on the specified case type
  
  # Calculate performance metrics for binary phenotype
  if (case.type == ANY.CASE.COL || case.type == PREVALENT.COL ||
      (case.type == INCIDENT.COL && is.na(FOLLOW.TIME))) {
    met.data <- log.metrics(
      case.type, next.pgs, data, base.models[[case.type]], ancestry
    )
  }
  
  # Calculate performance metrics for binary phenotype with time-to-event data
  else if (case.type == INCIDENT.COL && !is.na(FOLLOW.TIME)) {
    met.data <- cox.metrics(
      case.type, next.pgs, data, base.models[[case.type]], ancestry
    )
  }
  
  # Calculate performance metrics for quantitative phenotype (full cohort)
  else if (case.type == TRAIT.VALUE.COL) {
    met.data <- linear.metrics(
      case.type, next.pgs, data, base.models[[case.type]], ancestry
    )
  }
  
  return(met.data)
}

format.output <- function(raw.output, base.dems) {
  # Formats the raw output returned by the PGS analysis
  
  # Extract and combine the demographics data
  if (RUN.DEMS) {
    split.means <- lapply(raw.output, function(x) x$means)
    means <- do.call(cbind, split.means)
    full.dems <- cbind(base.dems, means)
  } else {
    full.dems <- NA
  }
  
  # Extract and combine the performance metrics data
  if (RUN.METRICS) {
    split.metrics <- lapply(raw.output, function(x) x$metrics)
    metrics <- do.call(rbind, split.metrics)
  } else {
    metrics <- NA
  }
  
  list(demographics = full.dems, metrics = metrics)
}


########################### Exporting Data Functions ###########################

get.run.spec <- function() {
  # Creates a run log containing the input parameters of the run
  
  # Report ancestry
  if (ANCESTRY[1] == "all" && !STRATIFIED) {
    ancestry <- "all"
  } else if (ANCESTRY[1] == "all" && STRATIFIED) {
    ancestry <- paste("all", ancestry.list(all.ancestries, par = TRUE))
  } else {
    ancestry <- ancestry.list(all.ancestries)
  }
  
  # Report bootstraps
  if (ANY.CASE.COL %in% case.types || PREVALENT.COL %in% case.types ||
      (INCIDENT.COL %in% case.types && is.na(FOLLOW.TIME))) {
    n.bootstraps <- BOOTSTRAPS
  } else {
    n.bootstraps <- NA
  }
  
  # Calculate run time
  end.time <- Sys.time()
  run.time <- difftime(end.time, start.time)
  
  run.spec <- c(
    "RUN",
    paste("    sampleset:", SAMPLESET),
    paste("    phenotype:", PHENOTYPE),
    paste("    start.time:", start.time),
    paste("    run.time:", round(run.time, 4), units(run.time)),
    "\nDATA",
    paste("    any.cases:", !is.na(ANY.CASE)),
    paste("    prevalent.cases:", !is.na(PREV.CASE)),
    paste("    incident.cases:", !is.na(INC.CASE)),
    paste("    incident.follow.up.time:", !is.na(FOLLOW.TIME)),
    paste("    quantitative.trait:", !is.na(TRAIT.VAL)),
    "\nPGS",
    paste("    pgs.detected:", nrow(index)),
    "\nPARAMETERS",
    paste("    regression.covariates:", report.covs(all.covs)),
    paste("    ancestry:", ancestry),
    paste("    stratified:", STRATIFIED),
    paste("    score.types:", paste(score.types, collapse = ", ")),
    paste("    min.cases:", MIN.CASES),
    paste("    skip:", SKIP),
    paste("    max.pgs:", MAX.PGS),
    paste("    bootstraps:", n.bootstraps),
    paste("    n.cores:", N.CORES),
    "\nRUN SPECIFICATIONS",
    paste("    run.demographics:", RUN.DEMS),
    paste("    run.metrics:", RUN.METRICS),
    paste("    run.plots:", RUN.PLOTS)
  )
}

report.covs <- function(covs) {
  # Returns a character string of covariates used in the regression analyses
  # to include in the run report
  
  covs <- covs[!(covs %in% PCs)]
  covs <- c(covs, "10 PCs")
  paste0(covs, collapse = ", ")
}

export.data <- function(demographic.data, metric.data) {
  # Writes output data to files
  
  cat("Exporting data\n")
  
  # Export demographic data as a CSV file
  if (RUN.DEMS) {
    write.csv(demographic.data, dem.outpath, row.names = FALSE)
  }
  
  # Export performance metrics data as a CSV file
  if (RUN.METRICS) {
    write.csv(metric.data, met.outpath, row.names = FALSE)
  }
  
  # Export running specifications as a text file
  run.spec <- get.run.spec()
  writeLines(run.spec, run.outpath)
}


############################# Multi-PGS Functions ##############################

multi.pgs.dataset <- function(ancestry) {
  # Creates and returns the data set for the multi-PGS (Huerta-Chagoya 2023)
  
  # Component scores
  score1 <- "PGS003443"
  score2 <- "PGS003444"
  score3 <- "PGS003445"
  
  # Extract data for each component score
  data1 <- get.dataset(score1, ancestry)
  data2 <- subset(get.dataset(score2, ancestry),
                  select = c(ID.COL, score.types))
  data3 <- subset(get.dataset(score3, ancestry),
                  select = c(ID.COL, score.types))
  
  for (score.type in score.types) {
    names(data1)[names(data1) == score.type] <- paste0(score.type, "_1")
    names(data2)[names(data2) == score.type] <- paste0(score.type, "_2")
    names(data3)[names(data3) == score.type] <- paste0(score.type, "_3")
  }
  
  # Merged the component scores into a single data frame
  data <- merge(x = data1, y = data2, by.x = ID.COL, by.y = ID.COL)
  data <- merge(x = data, y = data3, by.x = ID.COL, by.y = ID.COL)
  
  # Calculate the combined PGS for each ancestral group separately
  ancestries <- levels(as.factor(data$MostSimilarPop))
  data <- lapply(ancestries, calc.multi.pgs, data)
  data <- do.call(rbind, data)
  
  # Calculate the Z-score for the combined PGS
  for (score.type in score.types) {
    data[[score.type]] <- as.numeric(scale(data[[score.type]]))
  }
  
  return(data)
}

calc.multi.pgs <- function(ancestry, data) {
  # Calculates the multi-PGS score for the specified ancestral group
  
  # Linear PGS weights for each component score
  weight1 <- 0.531117
  weight2 <- 0.5690198
  weight3 <- 0.1465538
  
  data <- data[data$MostSimilarPop == ancestry,]
  
  # Calculate the Z-score for each component PGS and sum by their weights
  for (score.type in score.types) {
    score1 <- paste0(score.type, "_1")
    score2 <- paste0(score.type, "_2")
    score3 <- paste0(score.type, "_3")
    
    data[[score1]] <- as.numeric(scale(data[[score1]]))
    data[[score2]] <- as.numeric(scale(data[[score2]]))
    data[[score3]] <- as.numeric(scale(data[[score3]]))
    
    data[[score.type]] <- weight1 * data[[score1]] + weight2 * data[[score2]] +
      weight3 * data[[score3]]
  }
  
  return(data)
}


##################################### Main #####################################

# Desired column names
ID.COL <- "IID"
AGE.COL <- "age"
SEX.COL <- "sex"
ANY.CASE.COL <- "any.cases"
PREVALENT.COL <- "prevalent.cases"
INCIDENT.COL <- "incident.cases"
FOLLOW.UP.COL <- "follow.up"
TRAIT.VALUE.COL <- "trait.value"
TRAIT.SCALED.COL <- "trait.scaled"
SCORE.COL <- "score"
PCs <- paste0("PC", 1:10)

# Find which score types to use in analysis
score.types <- get.score.types()

# Load and prepare data
dem.data <- load.dem.data(DEM.FILE)
pop.data <- load.pop.data(POP.FILE)
pgs.data <- load.pgs.data(PGS.FILE)
index <- index.pgs.data()

# Find which case types are present in the input data (e.g. any, prevalent,
# incident, quantitative trait)
case.types <- get.case.types()

# Find which ancestries are present in input data
all.ancestries <- get.ancestries()

# Find which regression covariates should be used in the analysis
all.covs <- get.covariates()

# Analyse and export the data
all.output <- analyse.ancestries()
demographics <- all.output$demographics
metrics <- all.output$metrics

export.data(demographics, metrics)
