## Calculate the protein abundance by summing up its peptide intensity
library(reshape2)
library(ggplot2)
library(stringr)
library(dplyr)
rm(list = ls())

## get the intensity of each TMT channel
get_Ion_cols = function(input_df, label) {
  if (label == "TMT10") {
    input_df = input_df[, c(
      'fraction_scan',
      'Ion_126.128',
      'Ion_127.125',
      'Ion_127.131',
      'Ion_128.128',
      'Ion_128.134',
      'Ion_129.131',
      'Ion_129.138',
      'Ion_130.135',
      'Ion_130.141',
      'Ion_131.138'
    )]
  }
  if (label == "TMT11") {
    input_df = input_df[, c(
      'fraction_scan',
      'Ion_126.128',
      'Ion_127.125',
      'Ion_127.131',
      'Ion_128.128',
      'Ion_128.134',
      'Ion_129.131',
      'Ion_129.138',
      'Ion_130.135',
      'Ion_130.141',
      'Ion_131.138',
      'Ion_131.144'
    )]
  }
  return(input_df)
}

combine_msgf_masic_df = function(msgf_result_df,
                                 masic_result_dir,
                                 masic_result_filename,
                                 label) {
  masic_result_path = paste(masic_result_dir, masic_result_filename, sep='/')
  masic_out_df = read.csv(masic_result_path, sep = '\t')
  SpecFile = substring(masic_result_filename, 1, nchar(masic_result_filename) -
                         17)
  masic_out_df$fraction_scan = paste(SpecFile, masic_out_df$ScanNumber)
  masic_out_df = get_Ion_cols(masic_out_df, label)
  
  combined_data = merge(msgf_result_df,
                        masic_out_df,
                        by.x = 'fraction_scan',
                        by.y = 'fraction_scan',)
  
  print(combined_data)
  
  return(combined_data)
}

## The combined database search result of MS-GF+
## Here we take the result of one fraction of TW experiment as example
msgf_result_path = './TW.microbiome.proteins_scan.F1_01.txt'
## Put all the Masic output files of each cohort in one folder
masic_result_dir = './masic_result_dir_TW'
masic_result_files = list.files(masic_result_dir, pattern = '_ReporterIons.txt')


msgf_result_df = read.csv(msgf_result_path, sep = '\t')

for (masic_result_filename in masic_result_files) {
  print(paste('Processing', masic_result_filename))
  
  combined_data_ratio_this_fraction = combine_msgf_masic_df(msgf_result_df, 
                                                            masic_result_dir, 
                                                            masic_result_filename, 
                                                            'TMT10')
  
  if (exists('combined_data_ratio_all')) {
    combined_data_ratio_all = rbind(combined_data_ratio_all,
                                    combined_data_ratio_this_fraction)
  } else {
    combined_data_ratio_all = combined_data_ratio_this_fraction
  }
}

notsum_frame = combined_data_ratio_all

#filter out the proteins with the reference channel intensity = 0 
notsum_frame = notsum_frame %>% filter(Ion_126.128 != 0)

ion_names = c(
  'Ion_126.128',
  'Ion_127.125',
  'Ion_127.131',
  'Ion_128.128',
  'Ion_128.134',
  'Ion_129.131',
  'Ion_129.138',
  'Ion_130.135',
  'Ion_130.141',
  'Ion_131.138'
)

# Summing up the protein intensity in each sample
sum_frame = aggregate(notsum_frame[, ion_names],
                      by = list(notsum_frame$proteinId, notsum_frame$experiment),
                      sum)

colnames(sum_frame) = c('proteinId', 'experiment', ion_names)

## Save the protein channel intensity file to 02.protein_channel_intensity
write.table(sum_frame,
            '../02.protein_channel_intensity/TW.microbiome.F1_01.proteins_intensity.test.txt',
            row.names = F,
            quote = F,
            sep = '\t')
