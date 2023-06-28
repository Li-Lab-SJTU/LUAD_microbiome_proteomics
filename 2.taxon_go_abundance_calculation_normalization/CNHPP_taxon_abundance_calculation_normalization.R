# calculate the taxon abundance and normalize the taxon abundance matrix

# Load libraries
library(openxlsx)
library(stringr)
library(getopt)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)

# load the ggplot theme file
source("../ggplot_theme/mytheme_normalization_plot.R")

# The text file that stores the information on the intensities of peptides
pep_intensity_path = './peptide_intensity_data/CNHPP_peptide_intensity.txt'
# The the unipept output file that stores the taxon information of peptides
unipept_out_path = './unipept_output_data/CNHPP_unipept_taxon.csv'

pep_intensity_df = read.csv(pep_intensity_path, sep = '\t')
unipept_out_df = read.csv(unipept_out_path)
# Unipept equates "I" and "L" in the analysis. 
# We use "@" to denote both amino acid "I" and "L"
unipept_out_df$peptide_ilsame = str_replace_all(unipept_out_df$peptide, '[IL]', '@')

# Add postfixes to the taxon name
unipept_out_df$superkingdom = paste(unipept_out_df$superkingdom, '_k', sep =
                                      '')
unipept_out_df$phylum = paste(unipept_out_df$phylum, '_p', sep = '')
unipept_out_df$class = paste(unipept_out_df$class, '_c', sep = '')
unipept_out_df$order = paste(unipept_out_df$order, '_o', sep = '')
unipept_out_df$family = paste(unipept_out_df$family, '_f', sep = '')
unipept_out_df$genus = paste(unipept_out_df$genus, '_g', sep = '')
unipept_out_df$species = paste(unipept_out_df$species, '_s', sep = '')


peptide_taxonomy_intensity_df = merge(pep_intensity_df,
                                      unipept_out_df,
                                      by.x = 'peptide',
                                      by.y = 'peptide_ilsame',
)

# Calculate the abundances of taxa at different levels
intensity_add = function(peptide_taxonomy_intensity_df,
                         taxa_name) {
  x = aggregate(
    peptide_taxonomy_intensity_df$intensity,
    by = list(
      peptide_taxonomy_intensity_df$sample,
      peptide_taxonomy_intensity_df[, taxa_name]
    ),
    sum
  )
  return(x)
}

sample_taxonomy_intensity_df = rbind(
  intensity_add(peptide_taxonomy_intensity_df, 'superkingdom'),
  intensity_add(peptide_taxonomy_intensity_df,  'phylum'),
  intensity_add(peptide_taxonomy_intensity_df, 'class'),
  intensity_add(peptide_taxonomy_intensity_df,  'order'),
  intensity_add(peptide_taxonomy_intensity_df,  'family'),
  intensity_add(peptide_taxonomy_intensity_df,  'genus'),
  intensity_add(peptide_taxonomy_intensity_df, 'species')
)
colnames(sample_taxonomy_intensity_df) = c('sample', 'taxonomy', 'intensity')

# Reshape the data to a taxon_sample_abundance matrix
out_df = acast(sample_taxonomy_intensity_df, taxonomy ~ sample, value.var = 'intensity', drop = F, fun.aggregate = mean)

# Delete the unassigned taxa in each level
out_df = out_df[!rownames(out_df) %in% c('_k',
                                         '_p',
                                         '_c',
                                         '_f',
                                         '_o',
                                         '_g',
                                         '_s'), ]
# Sample types: NAT or Tumor
sample_type = sapply(
  str_split(colnames(out_df), '_'),
  FUN = function(x)
    x[[2]]
)
sample_id = colnames(out_df)

out_df_header = rbind(sample_type, sample_id)

# Remove the high missing ratio (> 0.8) taxon in both tumor and NAT samples
data_matrix_frame = out_df
sample_type = out_df_header['sample_type',]

is_na_NAT_data_matrix_frame = is.na(data_matrix_frame[, sample_type == "NAT"])
is_na_Tumor_data_matrix_frame = is.na(data_matrix_frame[, sample_type ==
                                                          "Tumor"])

row_NAT_na_count_data_matrix_frame = apply(is_na_NAT_data_matrix_frame, 1, sum)
row_Tumor_na_count_data_matrix_frame = apply(is_na_Tumor_data_matrix_frame, 1, sum)

is_NAT_na_ratio_smallerthreshold_matrix_frame = row_NAT_na_count_data_matrix_frame /
                                                  ncol(is_na_NAT_data_matrix_frame) < 0.8
is_Tumor_na_ratio_smallerthreshold_matrix_frame = row_Tumor_na_count_data_matrix_frame / 
                                                  ncol(is_na_Tumor_data_matrix_frame) < 0.8
is_all_na_ratio_smallerthreshold_matrix_frame = is_NAT_na_ratio_smallerthreshold_matrix_frame |
                                                  is_Tumor_na_ratio_smallerthreshold_matrix_frame

data_matrix_NAsmallerthreshold_frame = data_matrix_frame[is_all_na_ratio_smallerthreshold_matrix_frame,]

# Distribution before normalization
raw_frame = data_matrix_NAsmallerthreshold_frame
raw_frame[is.na(raw_frame)] = 1
# Select 14 samples to check the distribution of taxa abundances  
raw_frame = log2(raw_frame[,c(1,2,11,12,21,22,31,32,41,42,51,52,61,62)])
plt_sample_names = colnames(raw_frame)
colnames(raw_frame) = paste('s', 1:14,sep='')
raw_frame_melt = melt(raw_frame)
# Boxplot
p1_1 = ggplot(data=raw_frame_melt %>% filter(value != 0), aes(x=Var2, y=value)) + 
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  ylab('log2(intensity)') +
  scale_y_continuous(limits=c(20,45))

# Density plot
p1_2 = ggplot(data=raw_frame_melt %>% filter(value != 0), aes(x=value, color=Var2)) + 
  geom_density(size=1) + 
  mytheme1 + 
  xlab('log2(intensity)') + 
  scale_color_manual(values=rev(mycolors)) +
  labs(col='Sample') +
  scale_x_continuous(limits=c(20,45)) +
  scale_y_continuous(labels=label_number(accuracy = 0.01))

ggsave('./normalization_plots/taxon/CNHPP_taxonomy_boxplot_before_norm.png',p1_1,dpi=600,width = 10, height = 5)
ggsave('./normalization_plots/taxon/CNHPP_taxonomy_density_before_norm.png',p1_2,dpi=600,width = 10, height = 5)

# sample loading normalization
# sl : sample loading
raw_frame = data_matrix_NAsmallerthreshold_frame
col_sums = colSums(raw_frame,na.rm = T)
target = mean(col_sums)
norm_facts = target / col_sums
sl_frame = sweep(raw_frame, 2, norm_facts, FUN='*')

plot_frame = sl_frame
plot_frame[is.na(plot_frame)] = 1
plot_frame = log2(plot_frame[,plt_sample_names])
colnames(plot_frame) = paste('s', 1:14,sep='')
plot_frame_melt = melt(plot_frame)
p2_1 = ggplot(data=plot_frame_melt %>% filter(value != 0), aes(x=Var2, y=value)) + 
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  ylab('log2(intensity)')+
  scale_y_continuous(limits=c(20,45))


p2_2 = ggplot(data=plot_frame_melt %>% filter(value != 0), aes(x=value, color=Var2)) + 
  geom_density(size=1) + 
  mytheme1 + 
  xlab('log2(intensity)') + 
  scale_color_manual(values=rev(mycolors)) +
  labs(col='Sample') +
  scale_x_continuous(limits=c(20,45)) +
  scale_y_continuous(labels=label_number(accuracy = 0.01))


ggsave('./normalization_plots/taxon/CNHPP_taxonomy_boxplot_sl.png',p2_1,dpi=600,width = 10,height = 5)
ggsave('./normalization_plots/taxon/CNHPP_taxonomy_density_sl.png',p2_2,dpi=600,width = 10,height = 5)

##write the data after sample loading normalization
write_frame = sl_frame

sample_type = sapply(
  str_split(colnames(write_frame), '_'),
  FUN = function(x)
    x[[2]]
)
sample_id = colnames(write_frame)
out_df_header = rbind(sample_type, sample_id)

write.table(
  write_frame,
  './output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.noheader',
  col.names = F,
  sep = '\t',
  quote = F
)

write.table(
  out_df_header,
  './output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.header',
  col.names = F,
  sep = '\t',
  quote = F
)
