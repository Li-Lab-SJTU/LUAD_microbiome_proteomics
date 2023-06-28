# Calculate the taxon abundance and normalize the taxon abundance matrix

# Load libraries
library(openxlsx)
library(stringr)
library(getopt)
library(dplyr)
library(tidyverse)
library(reshape2)
library(scales)

# load the ggplot theme file
source("../ggplot_theme/mytheme_normalization_plot.R")

# 10 plex TMT labels 
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

# The text file that stores the information on the intensities of peptides
pep_intensity_path = './peptide_intensity_data/TW_scan_peptide_intensity.txt'
# The the unipept output file that stores the taxon information of peptides
unipept_out_path = './unipept_output_data/TW_unipept_taxon.csv'
# The meta file that stores the sample name of each plex in TMT experiments
sample_info_path = './MS_experiment_sample_info/TW_MS_experiment_sample_info.xlsx'

pep_intensity_df = read.csv(pep_intensity_path, sep = '\t')
sample_info_df = read.xlsx(sample_info_path, sep = '\t')

# Add experiment name postfixes to the taxon name
batch_name = sapply(strsplit(sample_info_df$experiment,'-'),function(x) x[1])
for (i in 4:13) {
  sample_info_df[, i] = paste(sample_info_df[, i], 
                              batch_name, 
                              sep = '_')
}

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
                         ion_names,
                         taxa_name) {
  x = aggregate(
    peptide_taxonomy_intensity_df[, ion_names],
    by = list(
      peptide_taxonomy_intensity_df$experiment,
      peptide_taxonomy_intensity_df[, taxa_name]
    ),
    sum
  )
  return(x)
}

experiment_taxonomy_intensity_df = rbind(
  intensity_add(peptide_taxonomy_intensity_df, ion_names, 'superkingdom'),
  intensity_add(peptide_taxonomy_intensity_df, ion_names, 'phylum'),
  intensity_add(peptide_taxonomy_intensity_df, ion_names, 'class'),
  intensity_add(peptide_taxonomy_intensity_df, ion_names, 'order'),
  intensity_add(peptide_taxonomy_intensity_df, ion_names, 'family'),
  intensity_add(peptide_taxonomy_intensity_df, ion_names, 'genus'),
  intensity_add(peptide_taxonomy_intensity_df, ion_names, 'species')
)

colnames(experiment_taxonomy_intensity_df)[1:2] = c('experiment', 'taxonomy')

# Get the sample name according to the TMT channel in each MS experiments
x = melt(
  experiment_taxonomy_intensity_df,
  id.vars = c('experiment', 'taxonomy'),
  measure.vars = ion_names,
  variable.name = 'channel',
  value.name = 'intensity'
)
y = melt(
  sample_info_df,
  id.vars = c('experiment'),
  measure.vars = ion_names,
  variable.name = 'channel',
  value.name = 'sample'
)

z = merge(x, y, by = c('experiment', 'channel'))

# Removal of samples that are not suitable for further analysis
# Disqualified samples in the experiment
z = z %>% filter(!str_detect(sample, 'Disqualified'))
# QC2: pool of Late Stage Tumor
z = z %>% filter(!str_detect(sample, 'LateStageTumor'))

# Reshape the data to a taxon_sample_abundance matrix
out_df = acast(
  z,
  taxonomy ~ sample,
  value.var = 'intensity',
  drop = F,
  fun.aggregate = sum
)

out_df[out_df == 0] = NA

# Delete the unassigned taxa in each level
out_df = out_df[!rownames(out_df) %in% c('_k',
                                         '_p',
                                         '_c',
                                         '_f',
                                         '_o',
                                         '_g',
                                         '_s'),]
# Sample types: NAT or Tumor
sample_type = sapply(
  str_split(colnames(out_df), '_'),
  FUN = function(x)
    x[[2]]
)
sample_id = colnames(out_df)

out_df_header = rbind(sample_type, sample_id)
out_df_header[out_df_header == 'Normal Adjacent Tissue'] = 'NAT'

# Remove the high missing ratio (> 0.8) taxon in both tumor and NAT samples
data_matrix_frame = out_df
sample_type = out_df_header['sample_type', ]

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

data_matrix_NAsmallerthreshold_frame = data_matrix_frame[is_all_na_ratio_smallerthreshold_matrix_frame, ]

# Distribution before normalization
raw_frame = data_matrix_NAsmallerthreshold_frame
raw_frame[is.na(raw_frame)] = 1
# Select 14 samples to check the distribution of taxa abundances  
raw_frame = log2(raw_frame[,c(1,2,3,4,41,42,43,44,81,82,83,84,201,203)])
plt_sample_names = colnames(raw_frame)
colnames(raw_frame) = paste('s', 1:14,sep='')
raw_frame_melt = melt(raw_frame)
# Boxplot
p1_1 = ggplot(data=raw_frame_melt %>% filter(value != 0), aes(x=Var2, y=value)) + 
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  ylab('log2(intensity)') +
  scale_y_continuous(limits=c(10,30))

# Density plot
p1_2 = ggplot(data=raw_frame_melt %>% filter(value != 0), aes(x=value, color=Var2)) + 
  geom_density(size=1) + 
  mytheme1 + 
  xlab('log2(intensity)') + 
  scale_color_manual(values=rev(mycolors)) +
  labs(col='Sample',
       ) +
  scale_x_continuous(limits=c(10,30)) +
  scale_y_continuous(limits=c(0,0.3), labels=label_number(accuracy = 0.01))

ggsave('./normalization_plots/taxon/TW_taxonomy_boxplot_before_norm.png',p1_1,dpi=600,width = 10, height = 5)
ggsave('./normalization_plots/taxon/TW_taxonomy_density_before_norm.png',p1_2,dpi=600,width = 10, height = 5)

# internal reference scaling normalization (IRS)
raw_frame = data_matrix_NAsmallerthreshold_frame
ir_sample = data.frame(raw_frame[,str_detect(colnames(raw_frame),'QC')])
ir_sample_names = colnames(ir_sample)
sum(ir_sample==0,na.rm=T) ## there are some ref-channels that have 0
ir_sample[ir_sample==0] = NA
sum(ir_sample==0,na.rm=T)
ir_sample$average = apply(ir_sample, 1, function(x) exp(mean(log(x), na.rm=T)))
ir_fracs = ir_sample$average / ir_sample[,1:27]

data_exp = raw_frame[,str_detect(colnames(raw_frame), '_Proteome_Batch001$')]
ir_frac_exp = ir_fracs[, str_detect(colnames(ir_fracs), '_Proteome_Batch001$')]
irs_frame = data_exp * ir_frac_exp

ir_sample_postfix = sapply(str_split(ir_sample_names,'.Reported_'),function(x) x[2])[2:27]

for (exp_name in ir_sample_postfix) {
  data_exp = raw_frame[,str_detect(colnames(raw_frame), exp_name)]
  ir_frac_exp = ir_fracs[, str_detect(colnames(ir_fracs), exp_name)]
  irs_frame = cbind(irs_frame, data_exp * ir_frac_exp)
}

# Distribution after IRS
plot_frame = irs_frame
plot_frame[is.na(plot_frame)] = 1
plot_frame = log2(plot_frame[,plt_sample_names])
colnames(plot_frame) = paste('s', 1:14,sep='')
plot_frame_melt = melt(plot_frame)
p2_1 = ggplot(data=plot_frame_melt %>% filter(value != 0), aes(x=Var2, y=value)) + 
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  ylab('log2(intensity)') +
  scale_y_continuous(limits=c(10,30))

p2_2 = ggplot(data=plot_frame_melt %>% filter(value != 0), aes(x=value, color=Var2)) + 
  geom_density(size=1) + 
  mytheme1 + 
  xlab('log2(intensity)') + 
  scale_color_manual(values=rev(mycolors)) +
  labs(col='Sample') +
  scale_x_continuous(limits=c(10,30)) +
  scale_y_continuous(labels=label_number(accuracy = 0.01))

ggsave('./normalization_plots/taxon/TW_taxonomy_boxplot_IRS.png',p2_1,dpi=600,width = 10, height = 5)
ggsave('./normalization_plots/taxon/TW_taxonomy_density_IRS.png',p2_2,dpi=600,width = 10, height = 5)

# sample loading normalization (SL)
raw_frame = irs_frame
col_sums = colSums(raw_frame,na.rm = T)
target = mean(col_sums)
norm_facts = target / col_sums
irs_sl_frame = sweep(raw_frame, 2, norm_facts, FUN='*')

# Distribution after IRS + SL
plot_frame = irs_sl_frame
plot_frame[is.na(plot_frame)] = 1
plot_frame = log2(plot_frame[,plt_sample_names])
colnames(plot_frame) = paste('s', 1:14,sep='')
plot_frame_melt = melt(plot_frame)
p3_1 = ggplot(data=plot_frame_melt %>% filter(value != 0), aes(x=Var2, y=value)) + 
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  ylab('log2(intensity)') +
  scale_y_continuous(limits=c(10,30))

p3_2 = ggplot(data=plot_frame_melt %>% filter(value != 0), aes(x=value, color=Var2)) + 
  geom_density(size=1) + 
  mytheme1 + 
  xlab('log2(intensity)') + 
  scale_color_manual(values=rev(mycolors)) +
  labs(col='Sample') +
  scale_x_continuous(limits=c(10,30)) +
  scale_y_continuous(labels=label_number(accuracy = 0.01))

ggsave('./normalization_plots/taxon/TW_taxonomy_boxplot_IRS_SL.png',p3_1,dpi=600,width = 10, height = 5)
ggsave('./normalization_plots/taxon/TW_taxonomy_density_IRS_SL.png',p3_2,dpi=600,width = 10, height = 5)


# Write the data after IRS + SL normalization
write_frame = irs_sl_frame
write_frame = write_frame[,!str_detect(colnames(write_frame), 'QC')]

sample_type = sapply(
  str_split(colnames(write_frame), '_'),
  FUN = function(x)
    x[[2]]
)
sample_id = colnames(write_frame)

out_df_header = rbind(sample_type, sample_id)
out_df_header[out_df_header == 'Normal Adjacent Tissue'] = 'NAT'

write.table(
  write_frame,
  './output_taxon_abundance_data/TW.taxon_abundance.irs_sl.NAfiltered.txt.noheader',
  col.names = F,
  sep = '\t',
  quote = F
)

write.table(
  out_df_header,
  './output_taxon_abundance_data/TW.taxon_abundance.irs_sl.NAfiltered.txt.header',
  col.names = F,
  sep = '\t',
  quote = F
)