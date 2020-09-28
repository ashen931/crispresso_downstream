#CRISPResso_allele_table_functs.R
# conda environment: crispresso_downstream_env
# last modified: 2020_08_19 Anne Shen

# for use with CRISPResso version 2.0.40
#
###  Dependencies:
# library(tidyselect)
# library(tidyverse)
# library(RColorBrewer)
# library(grid)
# options(scipen=999) #turn off scientific notation
#
### Functions:
# generate_updated_allele_tables()
#   get_run_table_by_crisp_mode()
#   extract_crispresso_run_info_and_allele_summaries()
#     regmatches_or_na()
#   generate_sample_allele_table_plots()
#     get_substitution_matrix_conversion_tb()
#     get_deletion_matrix_conversion_tb()
#     get_insertion_matrix_conversion_tb()
#     add_mutations_to_allele_matrix()
#     add_insertions_to_allele_matrix()
#     format_matrix_reference_rows
#     get_main_allele_tb()
#     get_freq_reads_tb()
#     get_legend_plot ()
#     get_allele_plot_composite_figure
#
### Sources:
#    CRISPRessoBatch_alleles_multiguide_multiallele_functs.R
#      get_batch_run_info()
#    CRISPRessoPooled_alleles_multiguide_functs.R
#      get_pooled_run_info()
#      reverse_complement() 
#    Summarize_off-target_editing_functs.R
#      save_composite_plot()


# tester 
#setwd("/Users/anneshen/Documents/local_working/local_Jing_BE/2020_1620_BE_NatureMed/20200810_BE_1620_rhAMPSeq_v2_0_40")
# setwd("/Users/anneshen/Documents/local_working/local_Jing_BE/2020_1620_BE_NatureMed/20191206_1620_input_rhAMPSeq_integrated3")
# 
# percent_freq_cutoff <- 0
# ref_nucleotide <- "C"
# target_nucleotides <- "T"
# crispressoPooled <- TRUE
# mode <- "BE_OT"
# 
# conversion_nuc_from <- ref_nucleotide
# conversion_nuc_to <- target_nucleotides
#
# allele_tb_csv <- "1620_BE_fig9_test.csv"

generate_updated_allele_tables <- function(updated_allele_tb_csv, mode, percent_freq_cutoff, 
                                           conversion_nuc_from, conversion_nuc_to, crispressoPooled){
  
  cat("CRISPResso run detailed allele table figure generation.\n",
      paste(Sys.time(), "\n", sep = ""),
      paste(getwd(), "\n", sep = ""),
      "percent_freq_cutoff: ", percent_freq_cutoff, "\n",
      "analysis_mode: ", mode, "\n",
      "allele_tb_csv_name: ", updated_allele_tb_csv, "\n",
      "ref_nucleotide: ", conversion_nuc_from, "\n",
      "target_nucleotide(s): ", conversion_nuc_to, "\n",
      "crispressoPooled :", crispressoPooled, "\n",
      "\n")
  
  #set date for saving figures
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  
  #read in allele sample table for allele table generation
  cat("Reading allele_tb_csv...", "\n")
  allele_tb_samples <- read.csv(updated_allele_tb_csv, stringsAsFactors = FALSE) %>%
    mutate(CRISPResso_dir_folder = mapply(function(x) ifelse(crispressoPooled, 
                                                             paste("CRISPRessoPooled_on_", x, sep = ""),
                                                             paste("CRISPRessoBatch_on_", x, sep = "")), 
                                          CRISPResso_dir_name))
  
  #get all summary file names
  if(grepl("BE", mode)){
    
    #if running OT analysis on base editing data, use BE summary tables
    conversion <- paste(conversion_nuc_from, "to", conversion_nuc_to, ".csv", sep = "")
    summary_file_suffix <- paste("BE_summary_", conversion, sep = "")
    
    list_summary_files <- list.files(pattern = summary_file_suffix)
    
  }else if(grepl("collapse", mode)){
    #get collapsed file suffix
    summary_file_suffix <- paste("collapsed_", percent_freq_cutoff, ".csv", sep = "")
    
    #if not running OT analysis on base editing data, use collapsed allele tables
    list_summary_files <- list.files(pattern = summary_file_suffix)
    
  }else{
    stop("No collapsed allele tables or base editing summary tables available.", call.=FALSE)
  }
  
  
  
  #extract CRISPResso run information for relevant samples
  list_crispresso_run_data_objects <- extract_crispresso_run_info_and_allele_summaries(allele_tb_samples, 
                                                                                       list_summary_files, 
                                                                                       crispressoPooled)
  all_info_tb <- list_crispresso_run_data_objects[[1]]
  all_samples_allele_tb_list <- list_crispresso_run_data_objects[[2]]
  
  
  #merge all_info_tb and allele_tb_samples
  #  1. calculate ref_cut_site, plot_first_index, plot_last_index
  #  2. generate the plot_ref_seq (get the reference amplicon sequence within the plot window)
  all_samples_info_tb <- all_info_tb %>%
    mutate(guide_aligned_ref = mapply(function(guide_seq, ref_seq) ifelse(regexpr(guide_seq, ref_seq)[1] <= 0, 
                                                                          reverse_complement(as.character(ref_seq)),
                                                                          as.character(ref_seq)),
                                      guide_seq, ref_seq),
           ref_cut_site = mapply(function(guide_seq, ref_seq, qw_center) regexpr(guide_seq, ref_seq)[1] + str_length(guide_seq) + qw_center,
                                 guide_seq, ref_seq, qw_center)) %>%
    mutate(plot_ref_seq = mapply(function(guide_aligned_ref, ref_cut_site, plot_window_size) substring(guide_aligned_ref,
                                                                                                       ref_cut_site-plot_window_size/2, 
                                                                                                       ref_cut_site+plot_window_size/2 - 1),
                                 guide_aligned_ref, ref_cut_site, plot_window_size),
           plot_first_index = ref_cut_site-plot_window_size/2,
           plot_last_index = ref_cut_site+plot_window_size/2 - 1,
           guide_first_index = mapply(function(guide_seq, ref_seq) regexpr(guide_seq, ref_seq)[1],
                                      guide_seq, ref_seq)) %>%
    right_join(allele_tb_samples, by = c("CRISPResso_dir_folder", "guide_seq" = "aligned_guide_seq")) %>%
    filter(!is.na(plot_window_size)) %>%
    select(-c("ref_seq", "guide_aligned_ref", "amplicon_sequence", "CRISPResso_dir_folder"))
  
  row.names(all_samples_info_tb) <- gsub("[[:punct:]]", "_", 
                                         paste(all_samples_info_tb$CRISPResso_dir_name, all_samples_info_tb$CRISPResso_run_name, sep = "__"))
  
  
  #generate and save (.pdf and .png) new allele tables
  generate_sample_allele_table_plots(all_samples_info_tb, all_samples_allele_tb_list, date)
  
}



  









get_run_table_by_crisp_mode <- function(crispressoPooled, run_log_string){
  
  if(crispressoPooled){ #if CRISPRessoPooled
    
    crispresso_run_info <- get_pooled_run_info(run_log_string)
    
  }else{ #if CRISPRessoPooled
    
    crispresso_run_info <- get_batch_run_info(run_log_string)
  }
  
  names(crispresso_run_info) <- c("guide_seq", "ref_seq", "ref_amplicon_name", "qw_size", "qw_center", "plot_window_size")
  
  return(crispresso_run_info)
}



extract_crispresso_run_info_and_allele_summaries <- function(allele_tb_samples, list_summary_files, crispressoPooled){
  
  cat("Extracting CRISPResso_on_ run data from CRISPResso_RUNNING_LOG.txt \n")
  
  #initialize empty data frame to store run summary information
  all_info_tb <- data.frame(guide_seq = c(), 
                            ref_seq = c(), 
                            ref_amplicon_name = c(),
                            qw_size = c(), 
                            qw_center = c(), 
                            plot_window_size = c(),
                            CRISPResso_dir_folder = c(),
                            stringsAsFactors = FALSE)
  
  #initialize list to store collapsed allele tables for desired CRISPResso runs
  all_samples_allele_tb_list <- list()
  
  #get CRISPResso run parameters
  for(crisp_outer_dir in unique(allele_tb_samples$CRISPResso_dir_name)){
    
    #get row indexes of allele_tb_samples that contain information on this crisp_outer_dir
    crisp_outer_dir_row_idx <- which(allele_tb_samples$CRISPResso_dir_name == crisp_outer_dir)
    
    #get CRISPResso batch/pool directory name
    crisp_outer_dir_folder <- allele_tb_samples$CRISPResso_dir_folder[crisp_outer_dir_row_idx[1]]
    cat("\nOpening ", crisp_outer_dir_folder, "\n")
    
    #get list of CRISPResso_on_ folders corresponding to desired CRISPResso runs
    list_crispresso_runs <- allele_tb_samples$CRISPResso_run_name[crisp_outer_dir_row_idx]
    
    #read in the collapsed_allele_tb for the crisp_outer_dir
    collapsed_tb_csv_name <- grep(paste(crisp_outer_dir, "_", sep = ""), list_summary_files, value = TRUE)
    collapsed_allele_tb <- read.csv(collapsed_tb_csv_name, stringsAsFactors = FALSE)
    #remove "X" from the beginning of columns beginning with a numerical name
    names(collapsed_allele_tb) <- gsub("^X(?=[0-9])", "", names(collapsed_allele_tb), perl = TRUE)
    
    
    #loop through crispresso runs in crisp_outer_dir_folder to extract run parameters
    for(crispresso_run in list_crispresso_runs){
      
      cat("Opening CRISPResso_on_", crispresso_run, "\n")
      
      #reconstruct run_log_file path 
      run_log_file <- paste(crisp_outer_dir_folder, 
                            paste("CRISPResso_on_", crispresso_run, sep = ""),
                            "CRISPResso_RUNNING_LOG.txt", sep = "/")
      
      cat("Reading ", run_log_file, "\n")
      
      #check that the run log exists
      if(file.exists(run_log_file)){
        
        run_log_string <- read.delim(run_log_file , sep = "\n", as.is = TRUE)
        crispresso_run_info <- get_run_table_by_crisp_mode(crispressoPooled, run_log_string[1]) %>%
          mutate(CRISPResso_dir_folder = crisp_outer_dir_folder)
        
        if(grepl("ERROR", run_log_string)){
          cat("CRISPRessoRun error, no results generated.\n")
          
        }else if(typeof(crispresso_run_info) == "character"){
          #If the desired allele is not found among the reference sequences, print "XXX allele not found."
          cat(crispresso_run_info, "\n")
          
        }else{
          
          #bind crispresso_run_info tables to generate composite table containing CRISPResso parameters
          # for all desired CRISPResso runs
          all_info_tb <- rbind(all_info_tb, crispresso_run_info)
          
          #select the "indel" and sample frequency/read columns and calculate indel lengths, etc.
          sample_allele_tb <- collapsed_allele_tb %>%
            select(c("indel", grep(paste(crispresso_run, "_reads", sep = ""), names(collapsed_allele_tb), value = TRUE)),
                   grep(paste(crispresso_run, "[__[ATCG]]?", sep = ""), names(collapsed_allele_tb), value = TRUE))
          
          names(sample_allele_tb) <- c("indel", "reads", "frequency")
          
          sample_allele_tb <- rbind(data.frame(indel = c("Reference", " ", "  "),
                                               frequency = c(0, 0, 0),
                                               reads = c(0, 0, 0)),
                                    sample_allele_tb)
          
          #calculate mutation indexes
          sample_allele_tb <- sample_allele_tb %>%
            mutate(substitution = mapply(function(indel) regmatches_or_na(indel, 
                                                                          gregexpr("[CTAG][0-9]{1,}[CTAG]", indel)), 
                                         indel),
                   in_del = mapply(function(indel) regmatches_or_na(indel, 
                                                                    regexpr("[+,-][0-9]{1,} ([[:print:]]{1,})", indel)), 
                                   indel)) %>%
            mutate(substitution_idx = mapply(function(substitution) as.numeric(regmatches_or_na(substitution, gregexpr("[-]?[0-9]{1,}", substitution))), 
                                             substitution),
                   substitution_to = mapply(function(substitution) regmatches_or_na(substitution, gregexpr("[ATCG]{1}$", substitution)), 
                                            substitution),
                   in_length = mapply(function(in_del) as.numeric(trimws(regmatches_or_na(in_del, 
                                                                                          regexpr("[+]{1}[0-9]{1,} ", in_del)))), 
                                      in_del),
                   in_seq = mapply(function(in_del) trimws(regmatches_or_na(in_del, regexpr("[ATCGN]{1,}", in_del))), 
                                   in_del),
                   in_start = mapply(function(in_del) as.numeric(trimws(gsub("[-\\)]$", "",(regmatches_or_na(in_del, 
                                                                                                             regexpr(" [-]{,2}[0-9]{1,}[-\\)]{1}", in_del)))))), 
                                     in_del),
                   del_length = mapply(function(in_del) as.numeric(gsub("^-", "", trimws(regmatches_or_na(in_del, 
                                                                                                          regexpr("[-]{1}[0-9]{1,} ", in_del))))), 
                                       in_del),
                   del_start = mapply(function(in_del) as.numeric(gsub("-$", "", gsub("[\\(\\)]", "", regmatches_or_na(in_del, 
                                                                                                                       regexpr("\\([-]{,2}[0-9]{1,}[-\\)]{1}", in_del))))), 
                                      in_del))
          sample_allele_tb$del_end <- sample_allele_tb$del_start + sample_allele_tb$del_length
          sample_allele_tb$in_end <- sample_allele_tb$in_start + sample_allele_tb$in_length
          
          #filter for indels represented in the sample
          sample_allele_tb <- filter(sample_allele_tb, !is.na(sample_allele_tb[,"frequency"]))
          
          #sort so that read counts/frequencies are in decreasing order
          sample_allele_tb <- sample_allele_tb[c(1:3, 
                                                 c(4:nrow(sample_allele_tb))[order(sample_allele_tb[4:nrow(sample_allele_tb),
                                                                                                    "frequency"], 
                                                                                   decreasing = TRUE)]),]
          
          #added the sample allele table to the all_samples_allele_tb_list
          all_samples_allele_tb_list[[gsub("[[:punct:]]", "_", paste(crisp_outer_dir, crispresso_run, sep = "__"))]] <- sample_allele_tb
        }
        
      }else{
        cat("*** CRISPRessoRun error, no results generated.\n")
      }
    }
  }
  
  return(list(all_info_tb, all_samples_allele_tb_list))
}



regmatches_or_na <- function(vector, match_data){
  
  match <- regmatches(vector, match_data)
  
  if(typeof(match) == "character" & length(match) < 1){ #if match_data generated by regexpr()
    return(NA)
    
  }else if(typeof(match) == "list" & length(match[[1]]) < 1){#if match_data generated by gregexpr()
    return(NA)
    
  }else{
    return(unlist(match))
  }
}


generate_sample_allele_table_plots <- function(all_samples_info_tb, all_samples_allele_tb_list, date){
  
  #loop through CRISPresso runs in all_samples_info_tb
  # 1. filter collapsed allele table to restrict table size and calculate matrix bp indices
  # 2. generate allele plot matrices for tracking mutations
  # 3. extract mutation index and conversion tables
  # 4. add mutations to allele plot matrices
  # 5. join collapsed allele tables and plot matrices & generate plotting tables
  # 6. calculate plot and font sizes
  # 7. generate updated CRISPResso Fig. 9 (allele tables)
  
  #initialize vector to track names of pdf files
  pdf_names <- c()
  
  for(run in row.names(all_samples_info_tb)){
    
    ####### 1. filter collapsed allele table to restrict table size and calculate matrix bp indices
    table_percent_cutoff <- all_samples_info_tb[run, "min_frequency_alleles_around_cut_to_plot"]
    table_row_cutoff <- all_samples_info_tb[run, "max_rows_alleles_around_cut_to_plot"] + 3
    
    run_summary_tb <- all_samples_allele_tb_list[[run]]
    sample_run_info <- all_samples_info_tb[run,]
    
    #filter allele table to restrict size (if too many alleles)
    if(nrow(run_summary_tb) > table_row_cutoff){
      
      #filter each allele table by table_percent_cutoff (skip first 3 reference rows)
      run_summary_tb <- run_summary_tb[c(1:3, c(4:nrow(run_summary_tb))[which(run_summary_tb[4:nrow(run_summary_tb), "frequency"] >= table_percent_cutoff)]),]
      
      if(nrow(run_summary_tb) > table_row_cutoff){
        #cut off each allele table to table_row_cutoff
        run_summary_tb <- run_summary_tb[table_row_cutoff,]
      }
    }
    
    
    #centered at the cut site, but index 1 is the first bp of the guide sequence
    bp_index_sequence <- (-sample_run_info$plot_window_size/2 - 
                            sample_run_info$qw_center +1):(sample_run_info$plot_window_size/2 - sample_run_info$qw_center)
    bp_index_start <- bp_index_sequence[1]
    bp_index_end <- sample_run_info$plot_window_size/2 - sample_run_info$qw_center
    
    #add this to the indel indexes to get matrix column indexes
    bp_idx_adj <- 1 - bp_index_start
    
    #adjust indel ranges to fall within plot window
    run_summary_tb$del_start <- mapply(function(x) ifelse(x < bp_index_start, bp_index_start, x), run_summary_tb$del_start)
    run_summary_tb$del_end <- mapply(function(x) ifelse(x > bp_index_end, bp_index_end, x), run_summary_tb$del_end) - 1
    run_summary_tb$in_start <- mapply(function(x) ifelse(x < bp_index_start, bp_index_start, x), run_summary_tb$in_start)
    run_summary_tb$in_end <- mapply(function(x) ifelse(x > bp_index_end, bp_index_end, x), run_summary_tb$in_end) - 1
    
    
    
    ####### 2. generate allele plot matrices for tracking mutations
    
    #generate initial allele matrix for tracking mutations
    allele_matrix <- as.matrix(data.frame("bp_1" = rep(sample_run_info$plot_ref_seq, times = nrow(run_summary_tb)),
                                          stringsAsFactors = FALSE) %>%
                                 separate(col = bp_1, into = paste("bp", bp_index_sequence, sep = "_"), 
                                          sep = 1:str_length(sample_run_info$plot_ref_seq)))
    
    
    
    ####### 3. extract mutation index and conversion tables
    
    #extract data frames of mutations to enter into the allele_matrix
    substitution_df <- get_substitution_matrix_conversion_tb(run_summary_tb, bp_idx_adj)
    deletion_df <- get_deletion_matrix_conversion_tb(run_summary_tb, bp_idx_adj)
    insertion_df <- get_insertion_matrix_conversion_tb(run_summary_tb, bp_idx_adj)
    
    #convert 0 frequency and reads to NA for plotting in run_summary_tb (do not do earlier as NAs are filtered out)
    run_summary_tb$frequency[1:3] <- NA
    run_summary_tb$reads[1:3] <- NA
    
    
    ####### 4. add mutations to allele plot matrices
    
    #add mutations to allele_matrix
    allele_matrix_sub <- add_mutations_to_allele_matrix( substitution_df, allele_matrix)
    allele_matrix_sub_del <- add_mutations_to_allele_matrix( deletion_df, allele_matrix_sub)
    allele_matrix_full <- add_insertions_to_allele_matrix(run_summary_tb, insertion_df, allele_matrix_sub_del, bp_idx_adj, bp_index_end)
    
    #format reference and sgRNA rows in allele_matrix_full
    allele_matrix_full_format <- format_matrix_reference_rows(sample_run_info, allele_matrix_full)
    
    #save run_summary_tb back into list with its modifications
    all_samples_allele_tb_list[[run]] <- run_summary_tb
    
    
    
    ####### 5. join collapsed allele tables and plot matrices & generate plotting tables
    
    joined_mat_tb <- cbind(run_summary_tb[,c("indel", "frequency", "reads")], allele_matrix_full_format) %>%
      transform(frequency = mapply(function(x) ifelse(!is.na(x), paste(as.character(round(x,2)), "%", sep = ""), NA), 
                                   frequency),
                reads = mapply(function(x) ifelse(!is.na(x), paste("(", as.character(x), " reads)", sep = ""), NA), 
                               reads))
    
    #generate plotting_tb for main allele table
    plotting_tb <- joined_mat_tb %>%
      gather(key = "bp", value = "nucleotide",vars_select(names(joined_mat_tb), contains("bp"))) %>%
      ## use pivot_longer() ONLY once conda r-tidyr is updated to >=1.0.0
      ## the following line performs the same function as the line using gather() above
      #pivot_longer(cols = vars_select(names(joined_mat_tb), contains("bp")), names_to = "bp", values_to = "nucleotide") %>%
      mutate(reference = mapply(function(x) ifelse(grepl("[PR]{1}", x), as.character(x), NA), nucleotide),
             nucleotide = gsub("sgRNA", NA, gsub("PAM", NA, nucleotide))) %>%
      mutate(insertion = mapply(function(x) ifelse(grepl("[ATCGN]{1}[\\*]{2}$", x), "I", 
                                                   ifelse(grepl("[0-9]{1,}", x), NA, sub("\\*", "", as.character(x)))), 
                                nucleotide),
             substitution = mapply(function(x) ifelse(grepl("[ATCGN]{1}\\*$", x), 
                                                      sub("\\*", "", as.character(x)), NA), nucleotide),
             nucleotide = gsub("[0-9]{1,}", NA, gsub("\\*", "", nucleotide)))
    
    plotting_tb$indel <- factor(plotting_tb$indel, levels = unique(plotting_tb$indel))
    plotting_tb$bp <- factor(plotting_tb$bp, levels = unique(plotting_tb$bp))
    
    
    #generate the data table containing frequency and read data for the allele table
    label_tb <- joined_mat_tb[,c("indel", "frequency", "reads")] %>%
      mutate(frequency = mapply(function(x) ifelse(is.na(x), x, paste(x, "  ", sep = "")), frequency)) %>%
      gather(key = "data", value = "value",  c("frequency", "reads"))
      ## use pivot_longer() ONLY once conda r-tidyr is updated to >=1.0.0
      ## the following line performs the same function as the line using gather() above
      # pivot_longer( cols = c("frequency", "reads"), names_to = "data", values_to = "value")
    label_tb$indel <- factor(label_tb$indel, levels = unique(label_tb$indel))
    
    
    ####### 6. calculate plot and font sizes
    
    #calculate composite plot relative widths and heights depending on the size of the allele table
    allele_tb_rel_width <- (max(str_length(run_summary_tb$indel), na.rm = TRUE) * 0.25 + sample_run_info$plot_window_size +
                              max(str_length(run_summary_tb$reads), na.rm = TRUE) * 0.25) %/% 10
    allele_to_read_freq_rel_width <- c(allele_tb_rel_width, 1)
    
    allele_tb_rel_height <- 1.5 + (nrow(run_summary_tb)-3) %/% 4
    allele_to_legend_rel_heights <- c(allele_tb_rel_height, 1.5)
    
    composite_plot_width_in <- allele_tb_rel_width * 3 + 2
    composite_plot_height_in <- allele_tb_rel_height + 2
    
    #calculate font sizes
    tile_font_size <- sample_run_info$plot_window_size * 0.05 + 1
    read_freq_font_size <- tile_font_size * 1.25
    indel_font_size <- nrow(run_summary_tb) * 0.05 + 10
    
    
    
    ####### 7. generate updated CRISPResso Fig. 9 (allele tables)
    
    #generate png plot name (also displayed at upper left of each image)
    composite_allele_plot_png_name <- paste(date, "__", run, "__allele_tb", sep = "")
    
    #generate individual plots
    cut_index <- sample_run_info$plot_window_size/2
    main_plot <- get_main_allele_tb(plotting_tb, cut_index, tile_font_size, indel_font_size)
    label_plot <- get_freq_reads_tb(label_tb, read_freq_font_size)
    legend_plot <- get_legend_plot()
    
    #generate composite allele plot object
    composite_grob_full <- get_allele_plot_composite_figure(main_plot, label_plot, legend_plot, 
                                           allele_to_read_freq_rel_width, allele_to_legend_rel_heights,
                                           composite_allele_plot_png_name)
    
    #save .pngs as individual files
    save_composite_plot(composite_allele_plot_png_name, composite_grob_full, 
                        plot_width_in = composite_plot_width_in, plot_height_in = composite_plot_height_in)
    
    
    
    #save .pdfs as individual files to preserve unique image sizes (will merge later)
    cairo_pdf(filename = paste(composite_allele_plot_png_name, ".pdf", sep = ""), onefile = TRUE, family = "Arial",
              width = composite_plot_width_in , height = composite_plot_height_in)
    
    print(composite_grob_full)
    
    #track pdf file names for combining and removing once all the pdfs are generated
    pdf_names <- c(pdf_names, paste(composite_allele_plot_png_name, ".pdf", sep = ""))
    
    #disconnect device for pdf printing
    dev.off()
    
  }
  
  #combine pdfs into single file
  pdf_combine(pdf_names, output = paste(date, "allele_tb.pdf", sep = "_"))
  
  #remove individual pdf files
  for(pdf_name in pdf_names){
    system(paste("rm ", pdf_name, sep = ""))
  }
  
  cat("\n")
  
}



get_substitution_matrix_conversion_tb <- function(run_summary_tb, bp_idx_adj){
  
  subs_tb_row_idx <- which(!is.na(run_summary_tb$substitution))
  
  if(length(subs_tb_row_idx) > 0){
    
    substitution_df <- data.frame(matrix_row_idx = unlist(mapply(function(x) rep(x, times = length(run_summary_tb$substitution_idx[[x]])), 
                                                                 subs_tb_row_idx)),
                                  matrix_col_idx = unlist(run_summary_tb$substitution_idx[subs_tb_row_idx]) + bp_idx_adj,
                                  matrix_new_symbol = paste(as.character(unlist(run_summary_tb$substitution_to[subs_tb_row_idx][[1]])),
                                                            "*", sep = ""),
                                  matrix_subs = unlist(run_summary_tb$substitution[subs_tb_row_idx]),
                                  stringsAsFactors = FALSE)
  }else{
    substitution_df <- NULL
  }

  return(substitution_df)
}


get_deletion_matrix_conversion_tb <- function(run_summary_tb, bp_idx_adj){
  
  deletion_tb_row_idx <- which(!is.na(run_summary_tb$del_length))
  
  if(length(deletion_tb_row_idx) > 0){
    
    deletion_df <- data.frame(matrix_row_idx = unlist(mapply(function(x) rep(x, times = (run_summary_tb$del_end[x] - 
                                                                                         run_summary_tb$del_start[x]) + 1), 
                                                             deletion_tb_row_idx)),
                              matrix_col_idx = unlist(mapply(function(x) run_summary_tb$del_start[x]:run_summary_tb$del_end[x], 
                                                             deletion_tb_row_idx)) + bp_idx_adj,
                              stringsAsFactors = FALSE)
    deletion_df <- filter(deletion_df, complete.cases(deletion_df))
    deletion_df$matrix_new_symbol <- "-"
    
  }else{
    deletion_df <- NULL
  }
  
  return(deletion_df)
}


get_insertion_matrix_conversion_tb <- function(run_summary_tb, bp_idx_adj){
  
  insertion_tb_row_idx <- which(!is.na(run_summary_tb$in_seq))
  
  if(length(insertion_tb_row_idx) > 0){
    
    insertion_df <- data.frame(matrix_row_idx = unlist(mapply(function(x) rep(x, times = str_length(run_summary_tb$in_seq[x])), 
                                                              insertion_tb_row_idx)),
                               matrix_col_idx = unlist(mapply(function(x) run_summary_tb$in_start[x]:run_summary_tb$in_end[x], 
                                                              insertion_tb_row_idx)) + bp_idx_adj,
                               matrix_new_symbol = paste(unlist(mapply(function(x) str_split(run_summary_tb$in_seq[x], ""), 
                                                                       insertion_tb_row_idx)), "**", sep = ""),
                               stringsAsFactors = FALSE)
  }else{
    insertion_df <- NULL
  }

  return(insertion_df)
}


add_mutations_to_allele_matrix <- function(mutation_df, plot_matrix){
  
  if(!is.null(mutation_df)){
    #loop through mutation_df and add mutations to the plot_matrix row by row
    for(n in 1:nrow(mutation_df)){
      plot_matrix[mutation_df$matrix_row_idx[n], mutation_df$matrix_col_idx[n]] <-  mutation_df$matrix_new_symbol[n]
    }
  }

  return(plot_matrix)
}


add_insertions_to_allele_matrix <- function(run_summary_tb, insertion_df, plot_matrix, bp_idx_adj, bp_index_end){
  
  for(idx in unique(insertion_df$matrix_row_idx)){
    
    #calculate the new plot_matrix column index of the shifted sequence (to the right)
    first_shifted_idx <- run_summary_tb$in_end[idx] + 1 + bp_idx_adj
    #calculate the last plot_matrix column index that will be kept in shifted sequence (within the plot window)
    last_kept_idx <- bp_index_end - run_summary_tb$in_length[idx] + bp_idx_adj
    
    #calculate the adjusted plot_matrix column indexes of the insertion start and plot window end
    adj_in_start <- run_summary_tb$in_start[idx] + bp_idx_adj
    adj_bp_index_end <- bp_index_end + bp_idx_adj
    
    insertion_shift <- plot_matrix[idx, adj_in_start:last_kept_idx]
    plot_matrix[idx, first_shifted_idx:adj_bp_index_end] <- insertion_shift
  }
  
  allele_matrix_sub_del_ins <- add_mutations_to_allele_matrix(insertion_df, plot_matrix)
  
  return(allele_matrix_sub_del_ins)
  
}


format_matrix_reference_rows <- function(sample_run_info, allele_matrix_full){
  matrix_guide_first_index <- sample_run_info$guide_first_index - sample_run_info$plot_first_index + 1
  matrix_guide_last_index <- matrix_guide_first_index + str_length(sample_run_info$guide_seq) -1
  matrix_pam_first_index <- matrix_guide_last_index  + 1
  matrix_pam_last_index <- matrix_pam_first_index + str_length(sample_run_info$pam) -1
  
  allele_matrix_full[3, ] <- NA
  allele_matrix_full[2, 1:(matrix_guide_first_index-1)] <- NA
  allele_matrix_full[2, matrix_guide_first_index:matrix_guide_last_index] <- "sgRNA"
  allele_matrix_full[2, matrix_pam_first_index:matrix_pam_last_index] <- "PAM"
  allele_matrix_full[2, (matrix_pam_last_index+1):ncol(allele_matrix_full)] <- NA
  
  return(allele_matrix_full)
}


get_main_allele_tb <- function(plotting_tb, cut_index, tile_font_size, indel_font_size){
  
  main_plot <- ggplot(data = plotting_tb,
                      aes(x = bp,
                          y = indel)) +
    geom_tile(aes(fill = reference),
              height = 0.5,
              width = 1, 
              alpha = 0.4) +
    geom_tile(aes(fill = nucleotide,
                  color = insertion),
              size = 0.55,
              height = 0.93,
              width = 0.93,
              alpha = 0.4) +
    geom_vline(xintercept = cut_index + 0.5,
               linetype = "dashed",
               color = "gray50") +
    xlab("") +
    ylab("") +
    geom_text(aes(label = nucleotide), 
              size = tile_font_size,
              family = "Arial") +
    geom_text(aes(label = substitution), 
              size = tile_font_size,
              fontface = "bold",
              family = "Arial Black") + #ADDED
    scale_x_discrete(limits = levels(plotting_tb$bp)) +
    scale_y_discrete(limits = rev(levels(plotting_tb$indel)),
                     breaks = rev(levels(plotting_tb$indel)),
                     labels = rev(levels(plotting_tb$indel))) +
    scale_fill_manual(name = "",
                      values = c("I" = "white",
                                 "A" = "#7fc97f", #light_green
                                 "T" = "#beaed4", #light_purple
                                 "C" = "#fdc086", #light_orange
                                 "G" = "#ffff99", #light_yellow
                                 "N" = "#c8c8c8", #light_gray
                                 "-" = "gray80",
                                 "PAM" = "steelblue1",
                                 "sgRNA" = "gray85"),
                      labels = c("I" = "insertion",
                                 "A" = "A",
                                 "T" = "T",
                                 "C" = "C",
                                 "G" = "G",
                                 "N" = "N",
                                 "-" = "deletion",
                                 "PAM" = "PAM",
                                 "sgRNA" = "sgRNA")) +
    scale_color_manual(name = "",
                       values = c("I" = "red", #"#c18072",
                                  "A" = "#7fc97f", #light_green
                                  "T" = "#beaed4", #light_purple
                                  "C" = "#fdc086", #light_orange
                                  "G" = "#ffff99", #light_yellow
                                  "N" = "#c8c8c8", #light_gray
                                  "-" = "gray80",
                                  "PAM" = "white",
                                  "sgRNA" = "white"),
                       labels = c("I" = "insertion",
                                  "A" = "A",
                                  "T" = "T",
                                  "C" = "C",
                                  "G" = "G",
                                  "N" = "N",
                                  "-" = "deletion",
                                  "PAM" = "PAM",
                                  "sgRNA" = "sgRNA")) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_text(size = indel_font_size, color = "black"),
          legend.position = "none")
  
  return(main_plot)
}


get_freq_reads_tb <- function(label_tb, tile_font_size){
  
  label_plot <- ggplot(data = label_tb,
                       aes(y = indel,
                           x = data)) +
    geom_tile(color = "white",
              fill = "white",
              height = 0.9,
              width = 1.5) +
    xlab("") +
    ylab("") +
    geom_text(aes(label = value), 
              size = tile_font_size,
              hjust = 1) +
    scale_y_discrete(limits = rev(levels(label_tb$indel)),
                     breaks = rev(levels(label_tb$indel)),
                     labels = rev(levels(label_tb$indel))) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          #text = element_text(family = "Helvitica"),
          legend.position = "none")
  
  return(label_plot)
}



get_legend_plot <- function(){
  
  #set legend for all plots
  legend_tb <- data.frame(row = rep(3:1, each = 6, times = 2),
                          column = c(rep(c(1:6), times = 3), rep(c(7:12), times = 3)),
                          key_lab = c("-", "deletion", rep(NA, 4), NA, "insertion", rep(NA, 4), 
                                      NA, "substitution", rep(NA, 4), NA, "sgRNA", rep(NA, 4), NA, "PAM", rep(NA, 10)),
                          bold_fill = c(rep(NA, 12), "bold", rep(NA, 23)),
                          color_fill = c("deletion", "w", rep("w", 16), "sgRNA", rep("w", 5), "PAM", "w", rep("w", 10)),
                          outline = c(rep(NA, 6), "insertion", rep(NA, 29)))
  
  legend_plot <- ggplot(data = legend_tb,
                        aes(x = column,
                            y = row)) +
    geom_tile(aes(fill = color_fill,
                  color = outline),
              size = 0.55,
              height = 0.93,
              width = 0.93, 
              alpha = 0.4) +
    geom_text(aes(label = key_lab), 
              size = 3,
              hjust = 0) +
    geom_text(aes(label = bold_fill), 
              size = 3,
              fontface = "bold") +
    xlab("") +
    ylab("") + 
    scale_fill_manual(values = c("w" = "white",
                                 "deletion" = "gray80",
                                 "PAM" = "steelblue1",
                                 "sgRNA" = "gray85")) +
    scale_color_manual(values = c("insertion" = "#c18072"), #"red"),
                       na.value = "white") +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          panel.spacing = unit(c(1,1,1,1), "mm"),
          #text = element_text(family = "Helvitica"),
          plot.background = element_rect(colour = "black", size = 1))
  
  return(legend_plot)
}



get_allele_plot_composite_figure <- function(main_plot, label_plot, legend_plot, 
                            allele_to_read_freq_rel_width, allele_to_legend_rel_heights,
                            composite_allele_plot_png_name){
  
  g_main <- ggplotGrob(main_plot)
  g_label <- ggplotGrob(label_plot)
  g_legend <- ggplotGrob(legend_plot)
  
  #generate title to go on the top left of the figure
  title_grob <- ggdraw() + draw_label(composite_allele_plot_png_name, hjust = 0, x = 0.01, y = 0.5)

  #generate white plot as filler
  filler <- ggplotGrob(ggplot() + theme(plot.background = element_rect(fill = "white"),
                                        panel.background = element_rect(fill = "white")))
  
  grid.newpage()
  
  composite_grob_allele_tb <- plot_grid(g_main, g_label, ncol = 2, align = "h", axis = "bt", 
                                        rel_widths = allele_to_read_freq_rel_width)
  
  composite_grob_legend <- plot_grid(filler, g_legend, filler, filler, filler, filler, ncol = 3, align = "h", axis = "bt", 
                                     rel_widths = rep(c(2, 1, 2), 2),
                                     rel_heights = rep(c(25, 1), each = 3))
  
  composite_grob_full <- plot_grid(title_grob, composite_grob_allele_tb, composite_grob_legend, ncol = 1, align = "h", axis = "bt", 
                                   rel_heights = c(0.5, allele_to_legend_rel_heights))
  
  return(composite_grob_full)
}




