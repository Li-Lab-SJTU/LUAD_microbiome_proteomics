# Calculate the microbial protein abundance and normalize the protein abundance matrix

rm(list = ls())
library(openxlsx)
library(stringr)
library(getopt)
library(dplyr)
library(reshape2)
library(tidyverse)
library(limma)
library(edgeR)

# load the ggplot theme file
source("../ggplot_theme/mytheme_normalization_plot.R")

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

##########combine intensity and sample information##########
msgf_masic_combined_file_path = './01.protein_channel_intensity/TW.microbiome.all.proteins_channel_intensity.txt'
sample_info_path = '../2.taxon_go_abundance_calculation_normalization/MS_experiment_sample_info/TW_MS_experiment_sample_info.xlsx'

data_agged_df = read.csv(msgf_masic_combined_file_path, sep = '\t')
sample_info_df = read.xlsx(sample_info_path)
batch_name = sapply(strsplit(sample_info_df$experiment,'-'),function(x) x[1])

for (i in 4:13) {
  sample_info_df[, i] = paste(sample_info_df[, i], batch_name, sep =
                                '_')
}

x = melt(
  data_agged_df,
  id.vars = c('experiment', 'proteinId'),
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

z$sample = gsub('\n', '_', z$sample)

z = z %>% filter(!str_detect(sample, 'Taiwanese'))
z = z %>% filter(!str_detect(sample, 'Disqualified'))
z = z %>% filter(!str_detect(sample, 'LateStageTumor'))

out_df = acast(
  z,
  proteinId ~ sample,
  value.var = 'intensity',
  drop = F,
  fun.aggregate = sum
)
out_df[out_df == 0] = NA

sample_type = sapply(
  str_split(colnames(out_df), '_'),
  FUN = function(x)
    x[[2]]
)
sample_id = colnames(out_df)

out_df_header = rbind(sample_type, sample_id)
out_df_header[out_df_header == 'Normal Adjacent Tissue'] = 'NAT'

##########move the high missing ratio row##########
data_matrix_frame = out_df
sample_type = out_df_header['sample_type',]

is_na_NAT_data_matrix_frame = is.na(data_matrix_frame[, sample_type == "NAT"])
is_na_Tumor_data_matrix_frame = is.na(data_matrix_frame[, sample_type ==
                                                          "Tumor"])

row_NAT_na_count_data_matrix_frame = apply(is_na_NAT_data_matrix_frame, 1, sum)
row_Tumor_na_count_data_matrix_frame = apply(is_na_Tumor_data_matrix_frame, 1, sum)

is_NAT_na_ratio_smallerthreshold_matrix_frame = row_NAT_na_count_data_matrix_frame / ncol(is_na_NAT_data_matrix_frame) < 0.9
is_Tumor_na_ratio_smallerthreshold_matrix_frame = row_Tumor_na_count_data_matrix_frame / ncol(is_na_Tumor_data_matrix_frame) < 0.9
is_all_na_ratio_smallerthreshold_matrix_frame = is_NAT_na_ratio_smallerthreshold_matrix_frame |
  is_Tumor_na_ratio_smallerthreshold_matrix_frame

data_matrix_NAsmallerthreshold_frame = data_matrix_frame[is_all_na_ratio_smallerthreshold_matrix_frame,]

############distribution before normalization############
raw_frame = data_matrix_NAsmallerthreshold_frame
raw_frame[is.na(raw_frame)] = 1
raw_frame = log2(raw_frame[, c(1, 2, 3, 4, 41, 42, 43, 44, 81, 82, 83, 84, 201, 203)])
plt_sample_names = colnames(raw_frame)
colnames(raw_frame) = paste('s', 1:14, sep = '')
raw_frame_melt = melt(raw_frame)
p1_1 = ggplot(data = raw_frame_melt %>% filter(value != 0), aes(x = Var2, y =
                                                                  value)) +
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  scale_y_continuous(limits = c(10, 30))

p1_2 = ggplot(data = raw_frame_melt %>% filter(value != 0), aes(x = value, color =
                                                                  Var2)) +
  geom_density(size = 1) +
  mytheme1 +
  xlab('log2(intensity)') +
  scale_color_manual(values = rev(mycolors)) +
  labs(col = 'Sample',) +
  scale_x_continuous(limits = c(10, 30))

ggsave(
  './normalization_plots/TW/TW_protein_boxplot_before_norm.jpg',
  p1_1,
  dpi = 600,
  width = 18,
  height = 9
)
ggsave(
  './normalization_plots/TW/TW_protein_density_before_norm.jpg',
  p1_2,
  dpi = 600,
  width = 18,
  height = 9
)

############internal reference normalization############
raw_frame = data_matrix_NAsmallerthreshold_frame
ir_sample = data.frame(raw_frame[,str_detect(colnames(raw_frame),'QC')])
ir_sample_names = colnames(ir_sample)
sum(ir_sample==0,na.rm=T) ## there are some ref-channels that have 0
ir_sample[ir_sample==0] = 1
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

##plot the distribution after IRS
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
  scale_x_continuous(limits=c(10,30))

ggsave('./normalization_plots/TW/TW_protein_boxplot_irs.jpg',p2_1,dpi=600,width = 18, height = 9)
ggsave('./normalization_plots/TW/TW_protein_density_irs.jpg',p2_2,dpi=600,width = 18, height = 9)


write_frame = irs_frame[,!str_detect(colnames(irs_frame), 'QC')]

sample_type = sapply(
  str_split(colnames(irs_frame), '_'),
  FUN = function(x)
    x[[2]]
)
sample_id = colnames(irs_frame)

out_df_header = rbind(sample_type, sample_id)
out_df_header[out_df_header == 'Normal Adjacent Tissue'] = 'NAT'

############sample loading normalization############
raw_frame = irs_frame
col_sums = colSums(raw_frame, na.rm = T)
target = mean(col_sums)
norm_facts = target / col_sums
irs_sl_frame = sweep(raw_frame, 2, norm_facts, FUN = '*')

##plot the data after IRS_SL
plot_frame = irs_sl_frame
plot_frame[is.na(plot_frame)] = 1
plot_frame = log2(plot_frame[, plt_sample_names])
colnames(plot_frame) = paste('s', 1:14, sep = '')
plot_frame_melt = melt(plot_frame)
p3_1 = ggplot(data = plot_frame_melt %>% filter(value != 0), aes(x = Var2, y =
                                                                   value)) +
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  scale_y_continuous(limits = c(10, 30))

p3_2 = ggplot(data = plot_frame_melt %>% filter(value != 0), aes(x = value, color =
                                                                   Var2)) +
  geom_density(size = 1) +
  mytheme1 +
  xlab('log2(intensity)') +
  scale_color_manual(values = rev(mycolors)) +
  labs(col = 'Sample') +
  scale_x_continuous(limits = c(10, 30))

ggsave(
  './normalization_plots/TW/TW_protein_boxplot_irs_sl.jpg',
  p3_1,
  dpi = 600,
  width = 18,
  height = 9
)
ggsave(
  './normalization_plots/TW/TW_protein_density_irs_sl.jpg',
  p3_2,
  dpi = 600,
  width = 18,
  height = 9
)

##write the data of sample loading normalization
write_frame = irs_sl_frame[, !str_detect(colnames(irs_sl_frame), 'QC')]

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
  './03.normalized_protein_abundance/TW.all.protein_abundance.irs_sl.NAfiltered.txt.noheader',
  col.names = F,
  sep = '\t',
  quote = F
)

write.table(
  out_df_header,
  './03.normalized_protein_abundance/TW.all.protein_abundance.irs_sl.NAfiltered.txt.header',
  col.names = F,
  sep = '\t',
  quote = F
)
