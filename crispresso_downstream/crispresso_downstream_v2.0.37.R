#crispresso_downstream_v2.0.37.R
# conda environment: crispresso_downstream_env
# last modified: 2020_07_21 Anne Shen

# for use with CRISPResso version 2.0.37

#load necessary R packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyselect))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(effsize))
suppressPackageStartupMessages(library(cowplot))

# "grid" is not found in anaconda
if (!require("grid", character.only = TRUE)) {
  install.packages("grid", dependencies = TRUE)
}
suppressPackageStartupMessages(library(grid))

#turn off scientific notation
options(scipen=999)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

#define user-input options (called in main())
get_option_list <- function(){
  
  #define command line user-input options
  option_list <- list(
    
    #mandatory inputs
    make_option(c("-c", "--crisp_out_dir"), action="store", type = "character",
                help="path to folder storing CRISPRessoPooled output directories"),
    make_option(c("-d", "--dataID"), action="store", type = "character", default = NULL,
                help="data identifier common to all names of CRISPRessoPooled input amplicons (regular expressions accepted)"),
    
    
    #optional inputs (with default settings)
    make_option(c("-m", "--mode"), action="store", type = "character", default = "collapse_only",
                help="analysis modes [default collapse_only]:\n
                \t collapse_only = generate collapsed allele table only \n
                \t collapse_BE = generate collapsed allele table & base editing summary (requires BE parameters) \n
                \t collapse_OT = generate collapsed allele table & off-target summary (requires OT parameters) \n
                \t collapse_BE_OT = generate collapsed allele table, base editing summary, and off-target summary \n
                \t \t \t (requires BE and OT parameters) \n
                \t BE_OT = generates base editing and off-target summary \n
                \t \t \t (requires BE parameters, OT parameters, and collapsed allele tables) \n
                \t BE_only = generates base editing summary (requires BE parameters and collapsed allele tables) \n
                \t OT_only = generates off-target summary (requires OT parameters and collapsed allele tables)"),
    make_option(c("-p", "--percent_freq_cutoff"), action="store", type = "double", default=0,
                help="minimum frequency an allele must appear in any sample/amplicon to be
                included in the collapsed output table [default %default]"),
    
    #collapse-mode specific inputs
    make_option(c("--CRISPRessoBatch"), action="store_true", type = "logical", default = FALSE,
                help="for any collapse mode: Use CRISPRessoBatch files as inputs"),
    make_option(c("--CRISPRessoPooled"), action="store_true", type = "logical", default = FALSE,
                help="for any collapse mode: Use CRISPRessoPooled files as inputs"),
    make_option(c("-n", "--noSub"), action="store_true", type = "logical", default = FALSE,
                help="do not include substitutions as edits in output table [default FALSE]"),
    
    #BE-mode specific inputs
    make_option(c("-f", "--conversion_nuc_from"), action="store", type = "character", default = "C",
                help="for any BE mode: the nucleotide targeted by the base editor [default C]"),
    make_option(c("-t", "--conversion_nuc_to"), action="store", type = "character", default = "T",
                help="for any BE mode: the nucleotide(s) produced by the base editor. If multiple nucleotides, \n
                \t \t \t enter all letters without separators (ex. ATCG) [default T]"),
    
    #OT-mode specific inputs
    make_option(c("-s", "--ot_sample_csv"), action="store", type = "character", default = "",
                help="for any OT mode: path to sample csv file (see manual for input files)"),
    make_option(c("-r", "--ref_seq_csv"), action="store", type = "character", default = "",
                help="for any OT mode: path to reference and guide sequences csv file (see manual for input files)"),
    make_option(c("-b", "--be_summary_exists"), action="store_true", type = "logical", default = FALSE,
                help="for OT_only mode: add this flag if BE summaries already exist AND set [-f, -t] if not default"),
    
    make_option(c("-e", "--scale_size_by_editing_freq"), action="store_true", type = "logical", default = FALSE,
                help="for OT_only mode: add this flag to separate the points in OT % editing scatterplot by size 
                according to sample read coverage [default FALSE]"),
    make_option(c("-l", "--low_coverage"), action="store", type = "double", default=1000,
                help="for OT_only mode with --editing_freq_scale flag: the upper read count cutoff for 
                \"low-coverage\" amplicons/samples  [default 1000]"),
    make_option(c("-u", "--high_coverage"), action="store", type = "double", default=10000,
                help="for OT_only mode with --editing_freq_scale flag: the lower read count cutoff for 
                \"high-coverage\" amplicons/samples [default 10000]")
  )
  
  return(option_list)
}


#check validity of user-input options and generate approriate error message (called in main())
#if options$crisp_out_dir exists, set working directory to options$crisp_out_dir
check_options <- function(options){
  
  valid_modes <- c("collapse_only", "collapse_BE", "collapse_OT", "collapse_BE_OT", 
                   "BE_OT", "BE_only", "OT_only")
  
  #stop if mode is not a valid option
  if (! options$mode %in% valid_modes) { 
    stop("Invalid analysis mode. See help [-h] for valid options.", call.=FALSE)
  }
  
  #stop if crisp_out_dir directory does not exist
  if(! dir.exists(options$crisp_out_dir)){
    stop("crisp_out_dir [-c] does not exist.", call.=FALSE)
  }else{
    setwd(options$crisp_out_dir)
  }
  
  #stop if in collapse mode & CRISPRessoBatch/CRISPRessoPooled option not input
  if(grepl("collapse", options$mode)){
    if(!(options$CRISPRessoBatch | options$CRISPRessoPooled)){
      stop("Must indicate either CRISPRessoBatch or CRISPRessoPooled files as inputs.", call.=FALSE)
    }
  }
  
  #if using OT mode:
  if(grepl("OT", options$mode)){
    
    #stop if running OT analysis and the ot_sample_csv does not exist
    if(!file.exists(options$ot_sample_csv)){
      stop("ot_sample_csv [-s] does not exist.", call.=FALSE)
    }
    
    #stop if running OT analysis and the ref_seq_csv does not exist
    if(!file.exists(options$ref_seq_csv)){
      stop("ref_seq_csv [-r] does not exist.", call.=FALSE)
    }
  }
  
  #stop if dataID is not provided
  if(is.null(options$dataID)){
    stop("dataID missing.", call.=FALSE)
  }
}



main <- function() {
  
  #source analysis functions and move to crisp_out_dir directory
  source("CRISPRessoPooled_alleles_multiguide_functs.R")
  source("CRISPResso_BE_analysis_functs.R")
  source("Summarize_off-target_editing_functs.R")
  source("CRISPRessoBatch_alleles_multiguide_multiallele_functs.R")
  #load extrafonts
  loadfonts()
  
  #define command line inputs/options
  option_list <- get_option_list()
    
  #parse command line options/arguments
  parser <- OptionParser(usage="%prog [options/input]", option_list=option_list)
  args <- parse_args(parser, 
                     args = commandArgs(trailingOnly = TRUE),
                     print_help_and_exit = TRUE,
                     positional_arguments = 0, 
                     convert_hyphens_to_underscores = TRUE)
  options <- args$options
  
  #check that inputs are valid AND setwd to options$crisp_out_dir if it is valid
  check_options(options)
  
  #store user command line inputs as variables 
  mode <- options$mode
  crisp_out_dir <- options$crisp_out_dir
  dataID <- options$dataID
  percent_freq_cutoff <- options$percent_freq_cutoff
  noSub <- options$noSub
  crispressoBatch <- options$CRISPRessoBatch
  crispressoPooled <- options$CRISPRessoPooled
  #BE-specific inputs
  conversion_nuc_from <- options$conversion_nuc_from
  conversion_nuc_to <- options$conversion_nuc_to
  #OT-specific inputs
  ot_sample_csv <- options$ot_sample_csv
  ref_seq_csv <- options$ref_seq_csv
  be_summary_exists <- options$be_summary_exists
  scale_size_by_editing_freq <- options$scale_size_by_editing_freq
  low_coverage <- options$low_coverage
  high_coverage <- options$high_coverage
  
  #if be_summary_exists and the user wishes to generate OT summaries only, append "BE" to mode
  if(be_summary_exists){
    mode <- paste("BE", mode, sep = "_")
  }
  
  #begin analysis run log
  sink("CRISPResso_downstream_analysis_log.txt", append=FALSE, split=TRUE)
  cat("Begin CRISPResso downstream analysis...\n\n")
  
  #if mode includes generating collapsed allele tables
  if(grepl("collapse", mode)){ 
    
    if(crispressoPooled){
      #get list of CRISPRessoPooled output directories
      crispresso2_pooled_files <- grep("CRISPRessoPooled", list.dirs(recursive = FALSE), value = TRUE)

      #generate CRISPResso2 summary allele tables
      get_CRISPRessoPooled_allele_tbs(crispresso2_pooled_files,
                                      dataID,
                                      percent_freq_cutoff,
                                      noSub)
    }else if(crispressoBatch){

      #get list of CRISPRessoBatch output directories
      crispresso2_batch_files <- grep("CRISPRessoBatch", list.dirs(recursive = FALSE), value = TRUE)
      
      #generate CRISPResso2 summary allele tables
      get_CRISPRessoBatch_allele_tbs(crispresso2_batch_files,
                                     dataID,
                                     percent_freq_cutoff,
                                     noSub)
    }
  }
  
  #if mode includes generating BE summary tables
  if(grepl("BE", mode)  & ! grepl("BE_OT_only", mode) ){
    
    #get list of collapsed allele tables
    summary_files <- grep(paste("collapsed_", percent_freq_cutoff, ".csv", sep = ""),
                          list.files(recursive = FALSE), value = TRUE)

    #generate CRISPResso2 BE summary allele tables
    get_BE_summary_tb(summary_files, dataID, percent_freq_cutoff,
                      ref_nucleotide = conversion_nuc_from, 
                      target_nucleotides = conversion_nuc_to)
  }
  
  #if mode includes generating OT summary tables
  if(grepl("OT", mode)){
    
    #generate off-target editing summary, statistics (t.test) table, and plots
    summarize_off_targets(mode, ref_seq_csv, ot_sample_csv, percent_freq_cutoff,
                          conversion_nuc_from, conversion_nuc_to,
                          scale_size_by_editing_freq, low_coverage, high_coverage)
  }

  #end run log
  sink()
  
  #remove empty Rplot.pdf if it was generated
  system("rm Rplots.pdf")
}

#call main() and run analysis
main()
