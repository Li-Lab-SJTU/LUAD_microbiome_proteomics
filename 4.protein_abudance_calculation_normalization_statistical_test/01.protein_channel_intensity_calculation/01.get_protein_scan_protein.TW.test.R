rm(list = ls())

library(dplyr)

## The path of MS-GF+ database search result, 
## Here we take one fraction of the TW result as an example
msgf_df_path = './20170207_LC_TMTB1_prot_F1_01.mzML.s3.smalldb.tsv.filtered.txt'

msgf_df = read.csv(msgf_df_path, sep = '\t')

msgf_df$fraction_scan = sapply(strsplit(msgf_df$X.SpecFile, '\\.'), function(x)
  x[[1]])
msgf_df$fraction_scan = paste(msgf_df$fraction_scan, msgf_df$ScanNum)

msgf_micro_df = msgf_df %>% filter(Class == 'Microbiome')

msgf_micro_df$proteinId = sapply(strsplit(msgf_micro_df$Protein, '\\|'), function(x) x[[2]])

msgf_out_df = subset(msgf_micro_df,select = c(fraction_scan, proteinId, experiment))
colnames(msgf_out_df) = c('fraction_scan', 'proteinId','experiment')

## Get the protein ID in each spectrum
out_path = './TW.microbiome.proteins_scan.F1_01.txt'
write.table(msgf_out_df, out_path, sep = '\t', quote = F, row.names = F)

