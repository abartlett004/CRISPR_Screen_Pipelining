##--------------------------------------------------
# generates files of counts for each sample group
# files to be used as inputs for MAGeCK
# --------------------------------------------------

library(phylotools)
library(readr)

# function for reading in tables of counts and storing as list of dataframes
read_into_df_list <- function(filename){
  df <- read.table(filename, sep = '\t')
  names(df) <- c('Seq', 'Count')
  return(df)
  
}

# read in count tables specified in snakemake input
count_tables <- lapply(snakemake@input, read_into_df_list)
names(count_tables) <- snakemake@input

# store other snakemake parameters
sample_sheet_f <- snakemake@params[["sample_sheet"]]
outdir <- paste0(dirname(snakemake@output[[1]]), "/")
groups <- snakemake@params[["groups"]]
vs_plasmid_lib <- snakemake@params[["vs_plasmid_lib"]]


# read in reference file of library sequences, order alphabetically by seq
mutant_table <- read.fasta(snakemake@params[["mutant_table"]])
mutant_table <- mutant_table[order(mutant_table$seq.text),]

# create column of HGSE gene symbols from seq names, necessary formatting for MAGeCK
for (name in mutant_table$seq.name){
  symbol <- unlist(strsplit(name, '_'))[1]
  mutant_table$symbol[mutant_table$seq.name == name] <- symbol
}

# read in sample sheet
if(length(sample_sheet_f) >= 1) {
  if(endsWith(sample_sheet_f, ".csv")) {
    sample_sheet <- read_csv(sample_sheet_f)
  } else if(endsWith(sample_sheet_f, ".tsv")) {
    sample_sheet <- read_tsv(sample_sheet_f)
  } else {
    sample_sheet <- read_excel(sample_sheet_f)
  }
  
  # make flexible for capitalization
  colnames(sample_sheet) <- tolower(colnames(sample_sheet))
  
}
sample_sheet$group <- as.character(sample_sheet$group)
sample_sheet$control <- as.character(sample_sheet$control)

# if samples are to be normalized against abundance of seqs in plasmid library
if (length(vs_plasmid_lib) > 0){

  # read in provided file of seq abundances in plasmid lib, sort alphabetically by sequence
  plasmid_counts <- read.csv('/home/ec2-user/arrakis-internal-data/reference/libraries/200424_plasmid_norms_RCed.csv', fileEncoding="UTF-8-BOM")
  plasmid_counts <- plasmid_counts[order(plasmid_counts$Seq), ]

  # for each group of samples to be averaged and/or normalized together
  for (group in groups){

    # begin building MAGeCK input sheet using columns for seq name, HGSE symbol, and count in plasmid library
    input_sheet <- as.data.frame(cbind(mutant_table$seq.name, mutant_table$symbol))
    names(input_sheet) <- c('name', 'symbol')
    input_sheet$plasmid_counts <- plasmid_counts$Count

    # scrape samples in current group from sample sheet
    samples <- sample_sheet$name[grepl(group, sample_sheet$group) & grepl('No', sample_sheet$control)]

    # for each sample in current group
    for (sample in samples){

      # isolate sample's count table from list, order alphabetically by sequence
      count_table <- count_tables[grepl(sample, names(count_tables))]
      count_table <- as.data.frame(count_table)
      names(count_table) <- c('Seq', 'Count')

      count_table$Seq <- as.character(count_table$Seq)
      count_table$Count <- as.numeric(count_table$Count)
      count_table <- count_table[order(count_table$Seq),]

      # add column to input sheet containing seq counts for current sample
      input_sheet[sample] <- count_table$Count
    }
    
    # write input sheet containing plasmid counts and counts from all samples in group to file
    write_tsv(input_sheet, paste0(outdir, group, '.txt'))

  }
# if samples are to be normalized against unsorted samples 
} else {

  # for each group of samples to be averaged and/or normalized together
  for (group in groups){
    
    # determine which samples in group are controls and which are experimental samples
    controls <- sample_sheet$name[grepl(group, sample_sheet$group) & grepl('Yes', sample_sheet$control)]
    samples <- sample_sheet$name[grepl(group, sample_sheet$group) & grepl('No', sample_sheet$control)]
    
    # begin building MAGeCK input sheet using columns for seq name and HGSE symbol
    input_sheet <- as.data.frame(cbind(mutant_table$seq.name, mutant_table$symbol))
    names(input_sheet) <- c('name', 'symbol')
    
    # for each control sample in list of control samples 
    for (control in controls){

      # isolate control sample's count table from list, order alphabetically by sequence
      count_table <- count_tables[grepl(control, names(count_tables))]
      count_table <- as.data.frame(count_table)
      names(count_table) <- c('Seq', 'Count')
      
      count_table$Seq <- as.character(count_table$Seq)
      count_table$Count <- as.numeric(count_table$Count)
      count_table <- count_table[order(count_table$Seq),]

      # add column containing seq counts for current control sample to input sheet 
      input_sheet[control] <- count_table$Count
    }

    # for each experimental sample in list of experimental samples 
    for (sample in samples){

      # isolate experimental sample's count table from list, order alphabetically by sequence
      count_table <- count_tables[grepl(sample, names(count_tables))]
      count_table <- as.data.frame(count_table)
      names(count_table) <- c('Seq', 'Count')

      count_table$Seq <- as.character(count_table$Seq)
      count_table$Count <- as.numeric(count_table$Count)
      
      # add column containing seq counts for current experimental sample to input sheet
      count_table <- count_table[order(count_table$Seq),]
      input_sheet[sample] <- count_table$Count
    }
    
    # write input sheet containing all samples in group to file
    write_tsv(input_sheet, paste0(outdir, group, '.txt'))
  }
}
