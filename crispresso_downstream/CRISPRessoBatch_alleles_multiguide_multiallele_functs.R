#CRISPRessoBatch_alleles_multiguide_multiallele_functs.R
# conda_environment: crispresso_downstream_env
# last modified: 2020_07_21 Anne Shen
# For use with CRISPResso version 2.0.37
#
# Modified to take multiple guides and alleles.
#
### Dependencies:
# library(tidyselect)
# library(tidyverse)
# library(RColorBrewer)
# library(grid)
# library(pheatmap)
# library(gtools)
#
### Options:
# options(scipen=999) #turn off scientific notation
#
#
# The functions are listed in the following order.
### Functions:
#     get_CRISPRessoBatch_allele_tbs() 
#           add_digit_names()
#       get_batch_allele_frequency_table() 
#           get_batch_run_info()
#           get_batch_run_table()
#           reverse_complement() 
#           get_batch_plot_allele_table_from_file()
#           select_alleles_for_plotting()
#               get_indel_columns()
#                 check_editing_pattern()
#                 get_deletion_summary()
#                 get_insertion_summary()
#                 get_substitution_summary()
#           collapse_duplicate_indels()
#                 get_duplicated_row_index()
#           order_indels()
#                 numeric_order()
#                 non_numeric_order()
#       adjust_percent_frequency()




###### get_CRISPRessoBatch_allele_tbs() #############################################################################
# This function takes a list of CRISPRessoBatch file names within the current working
# directory, extracts summary allele tables from each file, and saves the tables in the
# current directory.
#
# CALLS HELPERS: get_batch_allele_frequency_table(), adjust_percent_frequency()
# CALLED IN: scripts for CRISPResso2 analysis
# ARGUMENTS:  crispresso2_batch_files = list of sample directories (ex. list("./CRISPRessoBatch_on",...))
#             dataID = string idenitfier in the sample names that identifies a column as a column of 
#                      allele frequencies
#             percent_freq_cutoff = minimum percent frequency an allele must appear (in ANY sample within
#                                   the CRISPRessoBatch run) to be included in the collapsed allele table 
#             noSub = boolean indicating whether to include substitutions
# OUTPUT: none (directly saves allele tables to .csv files in current directory)
get_CRISPRessoBatch_allele_tbs <- function(crispresso2_batch_files,
                                           dataID,
                                           percent_freq_cutoff,
                                           noSub){
  
  #generate analysis log
  #sink("CRISPRessoBatch_collapsing_log.txt", append=FALSE, split=TRUE)
  
  cat("CRISPRessoBatch_allele_collapsing_log\n",
      paste(Sys.time(), "\n", sep = ""),
      paste(getwd(), "\n", sep = ""),
      "dataID: ", dataID, "\n",
      "percent_freq_cutoff: ", percent_freq_cutoff, "\n",
      "noSub: ", noSub, "\n",
      "\n")
  
  for(n in seq(1,length(crispresso2_batch_files))){

    #set working directory to correct batch directory)
    cat("***************************************************************************\n",
        "CRISPRessoBatch run: ", crispresso2_batch_files[n], "\n")
    setwd(crispresso2_batch_files[n])
    
    #get csv_prefix for saving file
    csv_prefix <- sub("./CRISPRessBatch_on_", "", crispresso2_batch_files[n])
    
    #select the list of CRISPresso2 output directories (may be 2-3 directories deep) with the sample_name of 
    # interest
    sample_dirs_list <- grep("/CRISPResso_on", list.dirs(), value = TRUE)
    
    #compile and format all allele frequency data from sample directories in sample_dirs_list into one 
    #  data table
    final_allele_table <- get_batch_allele_frequency_table(crisp_out_dir, sample_dirs_list, 
                                                     percent_frequency = percent_freq_cutoff, 
                                                     dataID , noSub)

    final_allele_table <- adjust_percent_frequency(final_allele_table, dataID)
    
    #remove Aligned_Sequence & Reference_Sequence columns
    final_allele_table <- select(final_allele_table, -c("Aligned_Sequence", "Reference_Sequence"))
    
    setwd("../") #save summary tables in outer directory
    
    #save allele tables
    save_to_csv(final_allele_table, ".", paste(csv_prefix,"collapsed", percent_freq_cutoff, sep = "_"))
  }
  
  #end analysis run log
  cat("\n")
  #sink()
}


###### get_batch_allele_frequency_table() #############################################################################
# This function loops through sample_dirs_list and compiles all the allele frequency data from samples in
# sample_dirs_list into one data table and formats the data. The function:
#     1. Loops through all directories in sample_dirs_list and 
#       2. Extracts CRISPResso2 sample/run information from RUNNING_LOG.txt
#       3. Extracts information on desired_allele from CRISPREsso Alleles_frequency_table.txt
#       4. Generates data frame that only contains alleles with valid editing patterns and includes their 
#          indel summaries in the "indel" column
#       5. Collapses ranked_indel_table duplicate indels so that there are no duplicate indels and filters
#          out alleles with no frequency > 0.1 percent in any sample
#     6. Collapse duplicate indels within the final compiled table
#     7. Order the alleles/rows by indel
#
# CALLS HELPERS: get_batch_plot_allele_table_from_file(), select_alleles_for_plotting(),
#                collapse_duplicate_indels(), order_indels()
# CALLED IN: CRISPResso2_analysis.R
# ARGUMENTS:  sample_dirs_list = list of sample directories (ex. list("./CRISPResso_on_BM_1_1",...))
#             percent_frequency = minimum allele frequency percent required for an allele to be included
#                                 in the final table
#             dataID = string idenitfier in the sample names that identifies a column as a column of 
#                      allele frequencies
# OUTPUT: ordered_alleles_table = data frame containing extracted and organized allele frequencies for all
#                                  samples listed in sample_dirs_list
get_batch_allele_frequency_table <- function(crisp_out_dir, sample_dirs_list, percent_frequency, dataID, noSub){
  #initalize ranking_table data frame
  ranking_table <- data.frame(Aligned_Sequence = as.character(c()), 
                              Reference_Sequence = as.character(c()), 
                              Unedited = as.logical(c()),
                              n_deleted = as.numeric(c()), 
                              n_inserted = as.numeric(c()), 
                              n_mutated= as.numeric(c()),
                              indel = as.character(c()))
  
  ### 1. Loops through directories in sample_dirs_list and compiles allele frequency data from all #########
  ### CRISPResso2 samples into one data frame for organization. ############################################
  for(sample_dir in sample_dirs_list){ #loop through directories in sample_dirs_list
    
    setwd(sample_dir) #set working directory to sample_dir to access the reference allele .txt file
    cat(sample_dir, "\n")
    
    ### 2. read CRISPResso2 run info (guide_seq, ref_seq, q_window, q_w_center, etc.) from running log .txt 
    ### file &organize in a data table (crispresso_run_info)
    run_log_string <- read.delim("CRISPResso_RUNNING_LOG.txt", sep = "\n", as.is = TRUE)
    crispresso_run_info <- get_batch_run_info(run_log_string[1])
    
    #if CRISPResso run had an error (and the run is not complete), skip the CRISPResso run and go on to the next 
    # directory
    if(grepl("ERROR", run_log_string)){
      cat("CRISPRessoRun error, no results generated.\n")
      
    }else if(typeof(crispresso_run_info) == "character"){
      #If the desired allele is not found among the reference sequences, print "XXX allele not found."
      cat(crispresso_run_info, "\n")
      
    }else{ #If the desired allele is found among the reference sequences, extract the allele frequency data.
      
      #this should be the table containing an entire run (including information from 2+ guides)
      indel_table_unique <- get_batch_run_table(sample_dir, crispresso_run_info, dataID, noSub)
      
      #compile data into ranking_table by merging allele_table and ranking_table
      ranking_table <- full_join(ranking_table, indel_table_unique,
                                 by = c("Aligned_Sequence", "Reference_Sequence", "Unedited",
                                        "n_deleted", "n_inserted", "n_mutated", "indel"))
      
    } 
    setwd("../") #set working directory back so to access other CRISPResso2 output directories
  }
  
  ### 6. Collapse duplicate indels within the final compiled table ##########################################
  ranking_table <- collapse_duplicate_indels(ranking_table, dataID)[1:ncol(ranking_table)]
  
  #keep only alleles that meet the percent_frequency threshold in at least one sample
  ranking_table_filtered <- filter_at(ranking_table, 
                             grep("reads", 
                                  grep(dataID, names(ranking_table), value = TRUE),
                                  invert = TRUE,
                                  value = TRUE), 
                             any_vars(. >= percent_frequency))
  
  ### 7. Order the alleles/rows by indel ####################################################################
  
  ordered_alleles_table <-order_indels(ranking_table_filtered)
  ordered_alleles_table$indel <- trimws(ordered_alleles_table$indel)
  
  return(ordered_alleles_table)
}

###### get_batch_run_info() #######################################################################################
# This function is used to extract the run data (guide_seq, ref_seq, q_window, qw_center, plot_window_size, 
# etc.) from the tab-delimited CRISPResso_RUNNING_LOG.txt file containing CRISPResso2 run parameters.
# The extracted data is later used to identify the mutations in the modified reads. 
#
# CALLS HELPERS: NA
# CALLED IN: get_batch_allele_frequency_table()
# ARGUMENTS:  run_log_string = string extracted from CRISPResso2 RUNNING_LOG.txt containing run parameters
# OUTPUT: allele_table = data frame containing extracted CRISPResso parameters for each allele
get_batch_run_info <- function(run_log_string){
  
  #extract the guide sequence
  guide_seqs <- str_extract(run_log_string[1],"--guide_seq [ACGT[:punct:]]{6,} --") %>% 
    str_remove( "--guide_seq ") %>% str_remove(" --") %>% strsplit(split = ",") %>%
    unlist()
  
  # extract the reference amplicon sequences
  ref_seqs <- str_extract(run_log_string[1],"-amplicon_seq [ACGTatcg[:punct:]]{1,} --") %>%
    str_remove("-amplicon_seq ") %>% str_remove(" --") %>% strsplit(split = ",") %>%
    unlist()
  
  # extract the amplicon names (should be in same order as ref_seqs)
  if(length(ref_seqs) > 1){
    ref_names <- str_extract(run_log_string[1],"--amplicon_name [[:alnum:][:punct:]]{1,} --") %>%
      str_remove("--amplicon_name ") %>% str_remove(" --") %>% strsplit(split = ",") %>%
      unlist()
  }else{
    ref_names <- ""
  }
  
  #extract the quantification window size (#bp around the predicted cut site)
  q_window <- str_extract(run_log_string[1],"--quantification_window_size [[:digit:]]{1,} --") %>% 
    str_remove("--quantification_window_size ") %>% str_remove(" --") %>% as.numeric()
  #extract the quantification window center/predicted cut site (relative to the right side of the guide 
  # sequence)
  qw_center <- str_extract(run_log_string[1],"--quantification_window_center [[:punct:][:digit:]]{1,} --") %>%
    str_remove("--quantification_window_center ") %>% str_remove(" --") %>% as.numeric()
  
  #extract plot window size (length of sequence shown in allele table)
  plot_window_size <- str_extract(run_log_string[1],"--plot_window_size [[:digit:]]{1,} --") %>%
    str_remove("--plot_window_size ") %>% str_remove(" --") %>% as.numeric()
  
  #get number of guides:
  n_guides <- length(guide_seqs)
  n_refs <- length(ref_seqs)
  n_guide_x_refs <- n_guides * n_refs
  
  run_info <- data.frame(guide_seq = rep(guide_seqs, each = n_refs), 
                         ref_seq = rep(ref_seqs, times = n_guides), 
                         ref_name = rep(ref_names,  times = n_guides),
                         q_window = rep(q_window, n_guide_x_refs), 
                         qw_center = rep(qw_center,  n_guide_x_refs), 
                         plot_window_size = rep( plot_window_size * 2, n_guide_x_refs))
  
  #return CRISPResso2 run/sample information as a data frame
  return(run_info)
}



###### get_batch_run_table() #######################################################################################
# This function is used to extract, summarize, and collapse all the alleles for a CRISPResso sample. It
# returns an allele frequency table for the CRISPResso input run (sample_dir).
#
# CALLS HELPERS: NA
# CALLED IN: get_batch_allele_frequency_table()
# ARGUMENTS:  sample_dir = directory name of CRISPResso sample run
#             crispresso_run_info = crispresso run parameters, amplicon, and guide table 
#                                   (output of get_batch_run_info())
#             dataID = CRISPResso run identifier
#             noSub = logical indicating whether substitutions are counted as real mutations
#
# OUTPUT: run_table = data frame containing summarized and collapsed indel frequency table for the
#                     CRISPResso run
get_batch_run_table <- function(sample_dir, crispresso_run_info, dataID, noSub){
  
  #initiate allele table structure
  run_table <- data.frame(Aligned_Sequence = as.character(c()), 
                          Reference_Sequence = as.character(c()), 
                          Unedited = as.logical(c()),
                          n_deleted = as.numeric(c()), 
                          n_inserted = as.numeric(c()), 
                          n_mutated= as.numeric(c()),
                          indel = as.character(c()),
                          substitutions = as.character(c()))
  
  for(n in seq(1, nrow(crispresso_run_info))){
    
    ### 1. extract info from crispresso_run_info data table ####################################################
    cat(crispresso_run_info$guide_seq[n])
    guide_seq <- crispresso_run_info$guide_seq[n]
    
    #get full reference sequence to align with guide_sequence to determine the cut site position
    ref_seq <-  toupper(crispresso_run_info$ref_seq[n])
    
    #get reference amplicon name of interest
    ref_name <-  crispresso_run_info$ref_name[n]
    
    q_window <- crispresso_run_info$q_window[n]
    qw_center <- crispresso_run_info$qw_center[n]
    plot_window_size <- crispresso_run_info$plot_window_size[n]
    
    ### 2. find bp position immediately anterior (-1) to the cut site (within the plot window (centered #######
    ### on cut site)) to check for indels #####################################################################
    cut_minus1 <- plot_window_size/2 
    
    ### 3. Order guide sequences and their bp location on the full amplicon reference sequence in a data ######
    ### frame. ################################################################################################
    if(regexpr(guide_seq, ref_seq)[1] < 0){
      #if the guide must be reverse-complemented, the qw_center is counted from the other end of the guide
      reverse_adjustment <- (-qw_center)      
      cut_site <- regexpr(reverse_complement(guide_seq), ref_seq) - qw_center
    }else{
      cut_site <- regexpr(guide_seq, ref_seq) + (str_length(guide_seq) + qw_center)
    }
    
    cut_index_adjustment <- cut_minus1 - (str_length(guide_seq) + qw_center) 
    
    # The indel indices are centered around the sgRNA. When counting from the PAM distal end in Cas9 editing,
    # the cut site is between bp 17-18, where PAM is 21-23. This cut_index_adjustment is the 
    # number to subtract from all bp indices taken from the CRISPResso2 plot window sequences to get the 
    # correct indel indices. (ex. cut_minus1 = 20, cut_index_adjustment = 3)
    # Using the first guide as a reference, the distance between the guides are subtracted from the 2nd guide onward. This allows all
    # the indels to be indexed the same way (calculates all edits as they are on a single read).
    
    ### 4. generate one allele table for each guide ###########################################################
    ### A. Extract information from CRISPREsso2 Alleles_frequency_table.txt ################
    allele_table <- get_batch_plot_allele_table_from_file(sample_dir, guide_seq, ref_name)
    
    #print information to check everything is correct
    cat(paste("guide_seq:", guide_seq, sep = " "), "\n")
    cat(paste("ref_name:", ref_name, sep = " "), "\n")
    cat(paste("ref_seq:", ref_seq, sep = " "), "\n")
    cat(paste("cut_minus1:", as.character(cut_minus1), sep = " "), "\n")
    cat(paste("cut_index_adjustment:", as.character(cut_index_adjustment), sep = " "), "\n")
    cat(paste("quantification_window_size:", as.character(q_window), sep = " "), "\n")
    cat(paste("quantification_window_center:", as.character(qw_center), sep = " "), "\n")
    
    ### 5. Generate data frame that only contains alleles with valid editing patterns and includes their ###
    ### indel summaries in the "indel" column ##############################################################
    ranked_indel_table <- select_alleles_for_plotting(allele_table, cut_minus1, cut_index_adjustment, ref_seq, 
                                                      q_window, noSub)
    
    #rename columns (remove "_percent_reads", leaving just sample names)
    colnames(ranked_indel_table) <- sub("_percent_reads", "", colnames(ranked_indel_table))
    
    cat("*************\n")
    cat(sum(filter(allele_table, Unedited == TRUE)[,ncol(allele_table)], na.rm = TRUE), "\n")
    cat(sum(filter(ranked_indel_table, Unedited == TRUE)[,ncol(ranked_indel_table)-2], na.rm = TRUE), "\n")
    cat("*************\n")
    
    #join tables for each allele together so that there is only one table per run
    run_table <- full_join(run_table, ranked_indel_table, 
                            by = c("Aligned_Sequence", "Reference_Sequence", "Unedited",
                                   "n_deleted", "n_inserted", "n_mutated", "indel", "substitutions"))
    
  }
  
  #clean up the run_table mutations:
  #   1. replace all NAs in columns containing "substitutions" with "" so that the display will be neater
  #   2. unite columns containing "substitutions" into one (subs on the same amplicon may not be picked up by both guide windows)
  #   3. unite mutations columns (indels + substitutions)
  run_table <-run_table %>% mutate_at(vars(contains("substitutions")), funs(replace(., is.na(.), ""))) %>%
      unite("substitutions", grep("^substitutions", names(run_table), value = TRUE), sep = " ", remove = TRUE) %>%
      unite("indel", c("substitutions", "indel"), sep = ", ", remove = TRUE)
  
  #remove ", " at the start of run_table$indel elements (residues from "unite" function above)
  run_table$indel <- gsub("^, ", "", run_table$indel)
  run_table$indel <- gsub(", $", "", run_table$indel)
  
  ### 6. Collapse ranked_indel_table duplicate indels so that there are no duplicate indels ######################################
  run_table <- collapse_duplicate_indels(run_table, dataID)[1:ncol(run_table)]

  return(run_table)
}


###### reverse_complement() #######################################################################################
# This function takes a DNA sequence (comprised of A, T, C, and G only) and returns the reverse complement
# sequence.
#
# CALLS HELPERS: NA
# CALLED IN: get_batch_allele_frequency_table()
# ARGUMENTS:  sequence = sequence containing upper-case ATCG
# OUTPUT: rev_comp = reverse complement of input sequence
reverse_complement<- function(sequence){
  
  bps <- str_split(sequence, pattern = "") #split sequence into individual base pairs
  rev_comp <- ""  #initalize empty string to hold reverse-complement sequence
  
  #loop through the bps in the sequence in reverse order
  for(n in seq(length(bps[[1]]), 1)){
    
    #get the complement of the base pair
    if(bps[[1]][n] == "A"){
      complement_bp <- "T"
    }else if(bps[[1]][n] == "T"){
      complement_bp <- "A"
    }else if(bps[[1]][n] == "G"){
      complement_bp <- "C"
    }else{
      complement_bp <- "G"
    }
    #join the complement base pairs in reverse order
    rev_comp <- paste(rev_comp, complement_bp, sep = "")
  }
  return(rev_comp)
}


###### get_batch_plot_allele_table_from_file() #############################################################################
# This function is used to extract the data from a tab-delimited .txt file containing allele frequencies 
# of the desired allele from from CRISPResso2 (Alleles_frequency_table.txt) into a data frame in R. 
# The extracted data is later compiled into a master data frame containg allele frequencies for all
# CRISPResso2 samples, which will be ordered and sorted for further analysis. 
#
# CALLS HELPERS: NA
# CALLED IN: get_batch_allele_frequency_table()
# ARGUMENTS:  sample_dir = name of CRISPResso2 run output directory, contains .txt file containing the
#                         reference allele frequency table to be extracted (ex. ./CRISPResso_on_DB5533_1)
#             guide_seq = sequence of guide RNA used for editing
# OUTPUT: allele_table = data frame containing extracted CRISPResso read data for each allele
get_batch_plot_allele_table_from_file <- function(sample_dir, guide_seq, ref_name){
  
  #generate the read_name (ex. DB5533_1_percent_reads) with which to label the frequency column,
  # which is necessary to keep track of samples once all frequencies are compiled into one data frame
  sample_handle <- str_split(sample_dir, "/") #file handle (includes directory above sample file)
  read_name_percent <-paste(sub("CRISPResso_on_", "", sample_handle[[1]][length((sample_handle[[1]]))]),
                            "percent_reads", paste0("_", ref_name), guide_seq, sep = "_")
  read_name_reads <-paste(sub("CRISPResso_on_", "", sample_handle[[1]][length((sample_handle[[1]]))]),
                          "reads", paste0("_", ref_name), guide_seq, sep = "_")
  
  ### generate alleles plot .txt file name based on desired allele & guide sequence
  #if ref_name = "", do not include "." before "Alleles_frequency_table_around_sgRNA_"
  if(ref_name == ""){
    filename <- paste("Alleles_frequency_table_around_sgRNA_",
                      guide_seq, ".txt", sep = "")
  }else{
    filename <- paste(ref_name, ".", "Alleles_frequency_table_around_sgRNA_",
                      guide_seq, ".txt", sep = "")
  }

  #check if Alleles frequency table exists (and has been unzipped)
  if(length(list.files(pattern = filename)) == 0){
    system(command = paste("gunzip ", filename, sep = ""), ignore.stdout = TRUE,
           wait = TRUE, timeout = 5)
  }
  
  #read the frequency file (should be the only file in alleles_list) into a data frame and format it
  allele_table <- read.delim(file = filename, header = TRUE) %>%
    select(Aligned_Sequence, Reference_Sequence, Unedited, n_deleted, n_inserted, n_mutated, X.Reads, X.Reads.1) %>%
    mutate(Aligned_Sequence = as.character(Aligned_Sequence),
           Reference_Sequence = as.character(Reference_Sequence),
           Unedited = as.logical(Unedited),
           n_deleted = as.numeric(n_deleted),
           n_inserted = as.numeric(n_inserted),
           n_mutated = as.numeric(n_mutated),
           X.Reads = as.numeric(X.Reads),
           X.Reads.1 = as.numeric(X.Reads.1)) %>%
    rename_at("X.Reads",~read_name_reads) %>%
    rename_at("X.Reads.1",~read_name_percent) #rename X.Reads.1 to match sample name
  
  return(allele_table)
}


###### select_alleles_for_plotting() #####################################################################
# This function is used to generate a data table (ranked_indel_table) that contains only alleles that:
#   1. have at least (percent_frequency) in any sample
#   2. display indels overlapping with the cut site (true edited alleles)
# Only selected alleles will be shown in the final data frame. 
#
# CALLS HELPERS: get_indel_columns() 
# CALLED IN: get_organized_allele_table()
# ARGUMENTS:  ranking_table = table containing alleles from all CRISPResso2 output directories
#             cut_minus1 = the bp position immediately downstream of cut site
#             ref_seq = reference sequence
#             guide_seq = guide sequence
# OUTPUT: ranked_indel_table = data table including only alleles that will be displayed in heatmap
select_alleles_for_plotting <- function(ranking_table, cut_minus1, cut_index_adjustment, ref_seq, q_window, noSub){
  #1. get_indel_columns() gets indel summaries for valid mutations and adds them as a column to ranking_table
  #2. only keep the alleles for which there is at least one valid mutation
  ranked_indel_table <- get_indel_columns(ranking_table, cut_minus1, cut_index_adjustment, ref_seq, q_window, noSub) %>%
    filter(grepl("[1-9|U]{1}", indel) | grepl("[1-9|U]{1}", substitutions)) 
  #all valid ranked_indel_table$indel values have numbers (ex. T17C, +1, -2) or "U" for "Unedited"
  
  #return the data table with the selected alleles
  return(ranked_indel_table)
}


###### get_indel_columns() #############################################################################
# This function takes a ranking_table containing allele frequency data and adds the indel column, which 
# contains indel summary information for each allele, ex. "T17C, +1 (T, 21)". This is done in 3 steps:
#   1. A mutation data frame (subsetted from ranking_table's 3 mutation columns) is generated.
#   2. The function loops through all alleles and checks each of the three mutation columns for n>0. If 
#      n>0, check_editing_pattern() is called. If editing_pattern==TRUE, the function to get the 
#      mutation summary is called, and n is replaced by the summary. If editing==FALSE, the mutation
#      summary is "invalid."
#   3. The 3 columns containing mutation summaries are collapsed into a single string for each allele, and
#      the vector of strings is added as the indel column to ranking_table.
#
# CALLS HELPERS: check_editing_pattern(), get_deletion_summary(), get_insertion_summary(), get_substitution_summary()
# CALLED IN: select_alleles_for_plotting()
# ARGUMENTS: ranking_table = data table containing all allele frequency data compiled from all
#                            CRISPResso2 samples
#            cut_minus1 = bp position immediately downstream (-1) of the cut site
#            ref_seq = full reference amplicon sequence
#            cut_index_adjustment = number to subtract from all bp indices taken from the CRISPResso2 plot 
#                                   window sequences to get correct bp indices relative to PAM-distal end
# OUTPUT: (ranking_table) data table with indel column added
get_indel_columns <- function(ranking_table, cut_minus1, cut_index_adjustment, ref_seq, q_window, noSub){

  #if don't want to count substitutions as edits:
  if(noSub){
    ranking_table$n_mutated <-rep(0, nrow(ranking_table)) #set all n_mutated to 0
  }
  
  ### 1. Get just the mutation columns from ranking_table (rows same order as ranking_table) as a new ######
  ### mutations data frame. The numeric entries will be replaced by the mutation summary, then collapsed ### 
  ### into strings contatining all mutation summaries for each allele. #####################################
  mutations <- select(ranking_table, c(n_deleted, n_inserted, n_mutated))
  
  #calculate quantification window bp boundaries (based on plot allele sequence lengths)
  q_window_start <- cut_minus1 + 1 - q_window
  q_window_end <- cut_minus1 + q_window 
  
  cat(paste("q_window_start: ", q_window_start, sep = ""), "\n")
  cat(paste("q_window_end: ", q_window_end, sep = ""), "\n")
  
  ### 2. Loop through all the alleles and replace the numbers in mutations (df) with mutation summaries. ###
  for (n in seq(1,nrow(ranking_table))){ #loop through all rows in table
    
    row_indel <- "" #initialize string to accumulate indel summaries
    row_editing_pattern <- TRUE
    
    if(! ranking_table$Unedited[n]){#if the amplicon sequence contains mutations
      
      if (ranking_table$n_mutated[n] > 0){#if the amplicon sequence contains a substitution
        #check the validity (editing pattern == TRUE) of the substitution
        substitution_valid <- check_editing_pattern(ranking_table[n,], 
                                                    cut_minus1, 
                                                    edited = "substitution", 
                                                    q_window,
                                                    q_window_start,
                                                    q_window_end)
        if(substitution_valid){#if the substitution matches the editing pattern
          #add the substitution summary to mutations$n_mutated at the correct row
          mutations$n_mutated[n] <- get_substitution_summary(ranking_table, n, 
                                                             cut_minus1, 
                                                             cut_index_adjustment, 
                                                             ref_seq, 
                                                             q_window,
                                                             q_window_start,
                                                             q_window_end)
        }else{#if the substitution does not match the editing pattern
          mutations$n_mutated[n] <- "invalid"
        }
      }
      
      if (ranking_table$n_inserted[n] > 0){#if the amplicon sequence contains an insertion
        #check the validity (editing pattern == TRUE) of the insertion
        insertion_valid <- check_editing_pattern(ranking_table[n,], 
                                                 cut_minus1, 
                                                 edited = "insertion", 
                                                 q_window,
                                                 q_window_start,
                                                 q_window_end)
        if(insertion_valid){#if the substitution matches the editing pattern
          #add the insertion (ex: "+2 (CT, 21-22)") to mutations$n_inserted at the correct row
          mutations$n_inserted[n] <- get_insertion_summary(ranking_table, n, 
                                                           cut_minus1, 
                                                           cut_index_adjustment, 
                                                           q_window,
                                                           q_window_start,
                                                           q_window_end)
        }else{#if the substitution does not match the editing pattern
          mutations$n_inserted[n] <- "invalid"
        }
      }
      
      if(ranking_table$n_deleted[n] > 0){#if the amplicon sequence contains a deletion
        #check the validity (editing pattern == TRUE) of the deletion
        deletion_valid <- check_editing_pattern(ranking_table[n,], 
                                                cut_minus1, 
                                                edited = "deletion", 
                                                q_window,
                                                q_window_start,
                                                q_window_end)
        if(deletion_valid){#if the substitution matches the editing pattern
          #add the deletion (ex: "-2 (21-22)") to mutations$n_deleted at the correct row
          mutations$n_deleted[n] <- get_deletion_summary(ranking_table, n, 
                                                         cut_minus1, 
                                                         cut_index_adjustment, 
                                                         q_window,
                                                         q_window_start,
                                                         q_window_end)
        }else{#if the substitution does not match the editing pattern
          mutations$n_deleted[n] <- "invalid"
        }
      }
      
    }else{ #if the amplicon sequence does not contain mutations (Unedited == TRUE)
      mutations$n_mutated[n] <- "Unedited"
    }
  }
  mutations[mutations == 0] <- NA #revert all the mutations == 0 to NA
  
  ### 3. Collapse the columns in mutations into a string, with each different mutation separated by ", " ###
  ### ex. T17C, +2 (GC 18-19). #############################################################################
  mutations$indel_list <- apply(mutations[ , c(1,2) ] , 1 , function(x) paste(x[!is.na(x)], collapse = ", "))
  
  #substitute an empty string "" for "invalid" within any of the mutations strings & remove beginning spaces
  mutations$indel_list <- sub("invalid, ", "", mutations$indel_list)
  mutations$indel_list <- sub(", invalid", "", mutations$indel_list)
  mutations$indel_list <- sub("invalid", "", mutations$indel_list)
  mutations$indel_list <- sub("^ ", "", mutations$indel_list)
  mutations$indel_list <- sub("^, ", "", mutations$indel_list)
  mutations$indel_list <- sub("$,", "", mutations$indel_list)
  mutations$n_mutated <- sub("^, ", "", mutations$n_mutated)
  mutations$n_mutated <- sub(", $", "", mutations$n_mutated)
  
  #add indel summary and match_editing_pattern as columns in ranking_table and return the data table
    return(ranking_table %>% mutate(indel = mutations$indel_list, substitutions = mutations[,3]))
}


###### check_editing_pattern() #############################################################################
# This function checks to see if the indels within an amplicon sequence overlap with the cut site 
# if indels exist. If not
# It returns TRUE if there is a bulge ("-") detected at the +1/-1 bp positions around the given
# cut site (meaning that the indel is likely due to nuclease cutting and not sequencing errors). 
# The function takes Aligned_sequence to check for deletions and Reference_sequence for insertions.
#
# CALLS HELPERS: NA
# CALLED IN: get_indel_colums() 
# ARGUMENTS: allele_row = the allele_row containing the Aligned_sequence and Reference_sequence
#            cut_minus1 = bp index to immediate left of cut site
#            edited = string indicating the type of mutation (or "false" for no mutation) in the aligned sequence
# OUTPUT: boolean (TRUE/FALSE) indicating whether or not the sequence matches a true cut pattern
check_editing_pattern <- function(allele_row, cut_minus1, edited, q_window, q_window_start, q_window_end){
  
  if(edited == "deletion"){ #if the indel is a deletion
    # 1. Take sequence within the quantification window (on Aligned_sequence).
    # 2. If "-" is detected as one of the characters in the substring, return TRUE. Else return FALSE.
    return(str_detect(substring(allele_row[1], q_window_start, q_window_end),"-"))
    
  } else if (edited == "insertion"){#if the indel is an insertion
    # 1. Take sequence within the quantification window (on Reference_sequence).
    # 2. If "-" is detected as one of the characters in the substring, return TRUE. Else return FALSE.
    return(str_detect( substring(allele_row[2], q_window_start, q_window_end),"-"))
    
  } else if (edited == "substitution" ) {#if there is a substitution, check that the bps in the 
    # sequence within the quantification window  are not the same
    return(!substring(allele_row[1], q_window_start, q_window_end) == substring(allele_row[2], q_window_start, q_window_end))
  }
}


###### get_deletion_summary() #############################################################################
# This function gets the deletion summary of row n in ranking_table. Whether or not the row contains a 
# deletion is determined in get_indel_columns(). The function adjusts the bp indices by subtracting
# the cut_index_adjustment to center the indel label cut site at bp 17-18 (PAM @ 21-23).
#
# CALLS HELPERS: NA
# CALLED IN: get_indel_columns()
# ARGUMENTS: ranking_table = data table containing all allele frequency data compiled from all
#                            CRISPResso2 samples
#            n = ranking_table row number (allele row) generated in get_indel_columns() for loop
#            cut_index_adjustment = the number to subtract from the deletion bp indices to 
#                                        get the cut index centered at sgRNA, PAM-distal 17-18
#                                        for heatmap labels
# OUTPUT: string of deletion summary ex. "-1 (T, 21)"
get_deletion_summary <- function(ranking_table, n, cut_minus1, cut_index_adjustment, q_window, q_window_start, q_window_end){
  
  #get the deletion sites & lengths within the plot allele sequence
  del_start_sites <- gregexpr('[-]{1,}', ranking_table$Aligned_Sequence[n])
  del_lengths <- attr(del_start_sites[[1]], "match.length")
  del_starts <-  `attributes<-`(del_start_sites[[1]],NULL) 
  
  for(i in seq(1,length(del_starts))){
    #checks to see which deletion is a true edit (overlaps with the cut site)
    if( (del_starts[i] + del_lengths[i]) >= q_window_start && del_starts[i] <=  q_window_end ){
      
      if(del_starts[i] == 1){ #if the deletion started before the beginning of the plotting window
        
        deletion_end <- del_starts[i] + del_lengths[i] - 1
        deletion_start <- deletion_end - ranking_table$n_deleted[n] + 1 #adjust for the bps not included in plotting window
        
        #if the deletion extends beyond the plotting window
      }else if(del_starts[i] + del_lengths[i] >= str_length(ranking_table$Aligned_Sequence[n])){
        
        deletion_start <- del_starts[i]
        deletion_end <- del_starts[i] + ranking_table$n_deleted[n] - 1 #adjust for the bps not included in plotting window
        
      }else{ #if the entire deletion falls within the plotting window
        deletion_start <- del_starts[i]
        deletion_end <- deletion_start + del_lengths[i] - 1
      }
      
      #calculate deletion length (need to count)
      deletion_length <- (deletion_end + 1) - deletion_start 
      
      #### generate the string containing the deletion summary
      if (deletion_length > 1){#if the deletion is longer than -1
        #use the adjusted bp indices (subtract cut_index_adjustment) for this label string
        deletion_base_positions <- paste(as.character(deletion_start - cut_index_adjustment), 
                                         as.character(deletion_end - cut_index_adjustment),
                                         sep = "-")
      }else{ #if the deletion is -1
        #use the adjusted bp indices (subtract cut_index_adjustment) for this label string
        deletion_base_positions <-as.character(deletion_start - cut_index_adjustment)
      }
      
      #if the correct deletion is found, break out of the loop by returning the deletion summary
      return(paste("-", deletion_length, 
                   " (", deletion_base_positions, ")", sep = "" ))
    }
    
  }
  return("invalid") #return "invalid" if no legitimate deletions found
  
}


###### get_insertion_summary() #############################################################################
# This function gets the insertion summary of row n in ranking_table. Whether or not the row contains an 
# insertion is determined in get_indel_columns(). The function adjusts the bp indices by subtracting
# the cut_index_adjustment to center the indel label cut site at bp 17-18 (PAM @ 21-23).
#
# CALLS HELPERS: NA
# CALLED IN: get_indel_columns()
# ARGUMENTS: ranking_table = data table containing all allele frequency data compiled from all
#                            CRISPResso2 samples
#            n = ranking_table row number (allele row) generated in get_indel_columns() for loop
#            cut_index_adjustment = the number to subtract from the insertion bp indices to 
#                                        get the cut index centered at sgRNA, PAM-distal 17-18
#                                        for heatmap labels
# OUTPUT: string of insertion summary ex. "+1 (T, 20)"
get_insertion_summary <- function(ranking_table, n, cut_minus1, cut_index_adjustment, q_window, q_window_start, q_window_end){
  
  #get the deletion site within the plot allele sequence
  insertion_start <- regexpr('-', substring(ranking_table$Reference_Sequence[n], q_window_start, q_window_end))[1] + q_window_start -1
  
  if (ranking_table$n_inserted[n] > 1){
    insertion_end <- insertion_start + ranking_table$n_inserted[n]-1
    #use the adjusted bp indices (subtract cut_index_adjustment) for this label string
    insertion_base_positions <- paste(as.character(insertion_start - cut_index_adjustment), 
                                      as.character(insertion_end - cut_index_adjustment),
                                      sep = "-")
  }else{ #if the insertion is +1
    insertion_end <- insertion_start
    #use the adjusted bp indices (subtract cut_index_adjustment) this label string
    insertion_base_positions <-as.character(insertion_start - cut_index_adjustment)
  }
  insertion_bases <- substr(ranking_table$Aligned_Sequence[n], insertion_start,
                            insertion_end)
  return(paste("+", as.character(ranking_table$n_inserted[n]), " (", insertion_bases, 
               " ", insertion_base_positions, ")", sep = "" ))
}


###### get_substitution_summary() #############################################################################
# This function gets the substitution summary of row n in ranking_table. Whether or not the row contains
# a substitution is determined in get_indel_columns(). The function adjusts the bp indices by subtracting
# the cut_index_adjustment to center the indel label cut site at bp 17-18 (PAM @ 21-23).
#
# CALLS HELPERS: NA
# CALLED IN: get_indel_columns()
# ARGUMENTS: ranking_table = data table containing all allele frequency data compiled from all
#                            CRISPResso2 samples
#            n = ranking_table row number (allele row) generated in get_indel_columns() for loop
#            cut_index_adjustment = the number to subtract from the insertion bp indices to 
#                                        get the cut index centered at sgRNA, PAM-distal 17-18
#                                        for heatmap labels
# OUTPUT: string of substitution summary ex. "T21C"
get_substitution_summary <- function(ranking_table, n, cut_minus1, cut_index_adjustment, ref_seq, q_window, q_window_start, q_window_end){
  
  ######check for mismatches in which the reference amplicon is not "-"
  
  #get both Reference and Aligned quantification window sequences & split into individual bps for comparison
  ref_qw <-str_split(substring(ranking_table$Reference_Sequence[n], q_window_start, q_window_end), "")[[1]]
  align_qw <-str_split(substring(ranking_table$Aligned_Sequence[n], q_window_start, q_window_end), "")[[1]]
  
  substitution <- "" #initialize empty substitution string
  
  #compare the ref_qw and align_qw bps at the same index to find substitutions
  for(i in seq(1, length(ref_qw))){
    
    if(ref_qw[i] != "-" && align_qw[i] != "-" && ref_qw[i] != align_qw[i]){
      substitution <- paste(substitution, 
                            paste(ref_qw[i], 
                                  #use the adjusted bp indices for label 
                                  as.character(i + q_window_start -1 - cut_index_adjustment), 
                                  align_qw[i], sep = ""),
                            sep = " ") #HERE
    }
  }
  return(substitution)
}


###### collapse_duplicate_indels() #############################################################################
# This function takes a data frame containing unsorted alleles and collapses all rows of duplicate indels into
# a single row containing the sum of frequencies for all the mis-sequenced allele columns. The function:
#     1. Roughly sorts alleles so that duplicate indels are grouped together in the data frame
#     2. Loops through the alleles (rows) in the data frame 
#     3. Calls get_duplicated_row_index() every time it encounters a duplicated indel to find the row index
#        of the last duplicate allele
#     4. Collapses all duplicates of the allele in question by replacing the frequency of the first allele 
#        instance with the sum of all the duplicate allele frequencies
#
# CALLS HELPERS: get_duplicated_row_index()
# CALLED IN: get_organized_allele_table()
# ARGUMENTS: ranked_indel_table = table containing allele that are true edits, have a minimum percent_frequency 
#                                 in at least one sample, and is ordered by indel (for aesthetics)
#            dataID = some string within each of the sample column headings that identify columns 
#                             containing allele frequency percentages
# OUTPUT: data table containing allele frequency data in which duplicated indels have been collapsed
#         into one row for each indel
collapse_duplicate_indels <- function(ranked_indel_table, dataID){
  
  ### 1. order the ranked_indel_table (will be re-ordered later through order_indels(), this is just to ###
  ### place duplicate indels immediately adjacent to one another for easy collapsing. #####################
  ranked_indel_table  <- ranked_indel_table[order(ranked_indel_table$indel, decreasing = TRUE),]
  
  #generate list of duplicate indels in ranked_indel_table
  duplicate_indels_list <- unique(ranked_indel_table[duplicated(ranked_indel_table$indel),]$indel)
  #add column indicating whether the indel is duplicated (FALSE also for 1st instance of indel)
  ranked_indel_table$duplicate <- duplicated(ranked_indel_table$indel)
  
  #get the first & last column indices of allele frequency data
  allele_freq_cols <- grep(dataID, names(ranked_indel_table), value = FALSE)
  allele_freq_col_first <- allele_freq_cols[1]
  
  num_row = 1 #initialize num_row interator for while loop
  
  ### 2. Loop through rows in ranked_indel_table (skip duplicate indel rows as necessary) ################
  while (num_row < nrow(ranked_indel_table)){
    dup_index <- 0 #re-set dup_index
    
    #if the row indel has duplicates (first iteration of indel)
    if((ranked_indel_table$indel[num_row] %in% duplicate_indels_list)){
      
      ### 3. get the row index of the last duplicate of this duplicated indel ############################
      dup_index <- get_duplicated_row_index(num_row, ranked_indel_table)
      
      ### 4. collapse duplicates: replace the percentage frequencies of this indel row with the sum of ###
      ### frequencies across all duplicate indel columns #################################################
      if (length(allele_freq_cols) == 2){ #if there are 2 allele frequency+reads  columns 
        # (called during table compilation when single tables are joined)
        
        ranked_indel_table[num_row,allele_freq_col_first] <-
          sum(ranked_indel_table[num_row:(num_row + dup_index), allele_freq_col_first], na.rm = TRUE)
        
      }else{ #if there are > 2 allele frequency+reads columns/guides (called after final table is compiled)
        
        allele_freq_col_last <- allele_freq_cols[length(allele_freq_cols)]

        ranked_indel_table[num_row,allele_freq_col_first:allele_freq_col_last] <-
          colSums(ranked_indel_table[num_row:(num_row + dup_index), allele_freq_col_first:allele_freq_col_last], na.rm = TRUE)
      }
    }
    #increase num_row iterator
    num_row = num_row + dup_index + 1
  }
  
  ranked_indel_table[ranked_indel_table == 0] <- NA #replace all 0 values with NA (colSUM of NAs = 0)
  
  return(ranked_indel_table[is.na(ranked_indel_table$duplicate),]) #return all non-duplicate & summed rows
}


###### get_duplicated_row_index() #############################################################################
# This function finds and returns the row index of the last duplicate indel row following a particular 
# indel in the plotting table.
#
# CALLS HELPERS: NA
# CALLED IN: collapse_duplicate_indels()
# ARGUMENTS: plotting_table = table containing allele that are true edits, have a minimum percent_frequency 
#                             in at least one sample, and is ordered by indel (for aesthetics)
#            num_row = the counter forr row number generated in collapse_duplicate_indels()
# OUTPUT: a number indicating the index of the last duplicate indel row; 
#         collapse_duplicate_indels() will use it to sum duplicate indel columns and to skip to the next
#         relevant row/allele
get_duplicated_row_index <- function(num_row, plotting_table){
  next_row = 1 #initiate next_row iterator for while loop
  
  # if the next row exists  && is a duplicated indel
  while((num_row + next_row <= nrow(plotting_table) && plotting_table$duplicate[num_row + next_row]) ){
    #accumulate the number of duplicated rows for colSums in collapse_indel_duplicates()
    next_row = next_row + 1 
  }
  
  #subtract 1 from next_row (to get last number that met the loop requirements)
  return(next_row - 1) 
}


###### order_indels() #################################################################################
# This function is used to order the alleles found in all CRISPResso2 samples in order so that the indels
# go from (undedited, substitutitons, +1/-1, +2/2...etc.). Function steps:
#     1. Generate a new sorting_col in indel_table_unique for sorting indels in desired manner. The 
#        function pulls out the number of deletions/insertions from the indels for numeric as opposed
#        to character sorting (1,2,3 vs. 1,10,11,2) into a new sorting_col column to use in the 
#        order() function.
#     2. Generate two data frames by character vs. numeric indel for ordering separately.
#     3. rbind the two data frames into one data frame.
#     4. Return final ordered data frame
#
# This function replaced rank_indels() in get_organized_allele_table().
# CALLS HELPERS: NA
# CALLED IN: get_organized_allele_table()
# ARGUMENTS:  ranking_table = table containing alleles that are true edits (indels match cut 
#                                  site) and have a minimum percent_frequency in at least one sample;
#                                  the indels have no duplicates (already collapsed)
# OUTPUT: ordered_table =  indel_table_unique sorted so that the alleles (rows) are in order of indels
#         (unmodified, substitutions, then +1, then –1, then +2, then –2, then +3, -3 etc.)
order_indels <- function(ranking_table){
  
  ### 1. Generate two new columns (mutation1, mutation2) in ranking_table for sorting mutations in ##########
  ### desired manner ########################################################################################
  mutation1 <- c() #initialize vector to make into mutation1 sorting column in ranking_table
  mutation2 <- c() #initialize vector to make into mutation2 sorting column in ranking_table
  
  for (indel in ranking_table$indel){ #loop through indels & save whether they are numeric (indel) or not in
    #mutation1 & mutation2 columns
    
    mutations_sort <- str_split_fixed(indel, ",", n=2)
    
    #if the indel is an insertion/deletion
    if (substring(indel, 1, 1) == "+" | substring(indel, 1, 1) == "-" ){
      #extract the indel number and add it to sorting_indel_list for numeric ordering
      mutation1 <- c(mutation1, regmatches(mutations_sort[,1], regexpr("[0-9]{1,}", mutations_sort[,1])))
      mutation2 <- c(mutation2, "")
      
    }else{#if the indel is a substitution or unmodified
      #extract the substitution as mutation1 and any additional indels as mutation2
      mutation1 <- c(mutation1, mutations_sort[,1])
      
      if(mutations_sort[,2] == ""){
        mutation2 <- c(mutation2, mutations_sort[,2])
      }else{
        mutation2 <- c(mutation2, regmatches(mutations_sort[,2], regexpr("[0-9]{1,}", mutations_sort[,2])))
      }
    }
  }
  
  #add mutation1 and mutation2 ordering columns to ranking_table
  ranking_table <- mutate(ranking_table, mutation1 = mutation1, mutation2 = mutation2)
  
  ### 2. Generate two data frames by character vs. numeric indel for ordering separately #################
  #filter out rows with character indels (ex. "T17C" or "Unmodified") into nonumeric_rows data frame
  #non_numeric_rows <- filter(ranking_table, grepl("^[UATCG]", ranking_table$mutation1) == TRUE)
  non_numeric_rows <- filter(ranking_table, grepl("[UATCG]{1}", mutation1))
  
  #order rows according to decreasing (Z --> A) indels in sorting_col
  non_numeric_rows <-non_numeric_order(non_numeric_rows, non_numeric_rows$mutation1)
  
  #get list of substitutions ("first mutation")
  unique_substitutions <- unique(non_numeric_rows$mutation1)
  
  ordered_non_numeric_rows <- data.frame()
  #loop through all substitutions and order them by indels (2nd mutations)
  for(subs in unique_substitutions){
    sub_table <- filter(non_numeric_rows, mutation1 == subs)
    sub_non_numeric_rows <- filter(sub_table, sub_table$mutation2 == "")
    sub_numeric_rows <- filter(sub_table, grepl("^[1-9]{1,}", sub_table$mutation2))
    sub_numeric_rows <- numeric_order(sub_numeric_rows, sub_numeric_rows$mutation2) 
    
    #bind the ordered data frames together
    ordered_non_numeric_rows <- rbind(ordered_non_numeric_rows,sub_non_numeric_rows,sub_numeric_rows)
  }
  
  #move the last non-numeric row (Unedited) to the top
  if(ordered_non_numeric_rows$mutation1[nrow(ordered_non_numeric_rows)] == "Unedited"){
    ordered_non_numeric_rows <- rbind(ordered_non_numeric_rows[nrow(ordered_non_numeric_rows),], 
                                      ordered_non_numeric_rows[1:nrow(ordered_non_numeric_rows)-1,])
  }
  
  #filter out rows with numeric indels (ex. "2" or "10") into numeric_rows data frame
  numeric_rows <- filter(ranking_table, grepl("^[1-9]{1,}", mutation1))
  
  #order rows according to increasing (1 --> 9) indels in sorting_col
  order_numeric_rows <- numeric_order(numeric_rows, numeric_rows$mutation1)
  
  # ### 3. Bind the two data frames into one data frame ####################################################
  # #rbind the two data frames together, with the character indels on top, and remove sorting_col
  # #ORDER: unmodified, substitutions, then +1, then –1, then +2, then –2, then +3, -3 etc.
  ordered_table <- rbind(ordered_non_numeric_rows, order_numeric_rows )[, 1:(ncol(order_numeric_rows)-2)]
  
  ### 4. Return final ordered data frame #################################################################
  return(ordered_table)
}


###### numeric_order() #############################################################################
# This function is used to order the rows in an indel table according to increasing (1 --> 9) indels in 
# the sorting_col column.
# CALLS HELPERS: NA
# CALLED IN: order_indels()
# ARGUMENTS:  numeric_table = the data frame containing an non-numeric indel column to be sorted
#             sorting_list = vector containing the numeric_table column values by which the table
#                            is to be sorted in increasing (1 --> 9) order
# OUTPUT: numerically sorted table (by sorting list)
numeric_order <- function(numeric_table, sorting_list){
  return(numeric_table[order(as.numeric(sorting_list), decreasing = FALSE),])
}


###### non_numeric_order() #############################################################################
# This function is used to order the rows in an indel table according to decreasing (Z --> A) indels in 
# the sorting_col column.
# CALLS HELPERS: NA
# CALLED IN: order_indels()
# ARGUMENTS:  non_numeric_table = the data frame containing an non-numeric indel column to be sorted
#             sorting_list = vector containing the non_numeric_table column values by which the table
#                            is to be sorted in decreasing (Z --> A) order
# OUTPUT: alphabetically sorted table (by sorting list)
non_numeric_order <- function(non_numeric_table, sorting_list){
  return(non_numeric_table[order(sorting_list, decreasing = FALSE),])
}


###### save_to_csv() #############################################################################
# This function is used to save the desired data frame to a .csv file (name chosen by user) in the chosen
# directory. The function is non-fruitful.
#
# CALLS HELPERS: NA
# CALLED IN: get_CRISPRessoBatch_allele_tbs()
# ARGUMENTS:  table_to_save = the data frame containing the data to save to the .csv file
#             save_directory = the directory to which to save the .csv file
#             save_name = the name to which to save the desired data table (table_to_save) as .csv
# OUTPUT: none
save_to_csv <- function(table_to_save, save_directory, save_name){
  csv_name <- paste(save_directory, paste(save_name,".csv", sep = ""), sep = "/")
  write.csv(table_to_save, file = csv_name, row.names = FALSE)
}


###### adjust_percent_frequency() #############################################################################
# This function recalculates percent frequencies of alleles so that they add to 100% for 
# each sample
#
# CALLS HELPERS: NA
# CALLED IN: get_CRISPRessoPooled_allele_tbs()
#
# ARGUMENTS:  final_allele_table = data frame allele frequency data
#             dataID = identifier for columns containing frequency data
# OUTPUT: final_allele_table = final_allele_table with adjusted frequencies such that all frequencies for
#                              a CRISPResso run/sample sum to 100%
adjust_percent_frequency <- function(final_allele_table, dataID){
  
  data_cols <- grep("reads", grep(dataID, names(final_allele_table), value = TRUE), 
                    invert = TRUE, 
                    value = TRUE)
  new_allele_total_perc <- colSums(select(final_allele_table, data_cols), na.rm = TRUE)
  
  for(n in data_cols){
    final_allele_table[[n]] <- final_allele_table[[n]] * 100/new_allele_total_perc[[n]] 
  }
  
  return(final_allele_table)
}
