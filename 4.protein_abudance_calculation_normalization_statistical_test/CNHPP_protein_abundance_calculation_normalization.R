# Calculate the microbial protein abundance and normalize the protein abundance matrix

rm(list = ls())
library(openxlsx)
library(stringr)
library(getopt)
library(dplyr)
library(reshape2)
library(ggplot2)
library(limma)

# load the ggplot theme file
source("../ggplot_theme/mytheme_normalization_plot.R")

################Calculate the abundances of proteins###################
protein_channel_intensity_path = './01.protein_channel_intensity/CNHPP.microbiome.all.proteins_channel_intensity.txt'

protein_channel_intensity_df = read.csv(protein_channel_intensity_path, sep = '\t')

colnames(protein_channel_intensity_df)[2] = 'sample'

out_df = acast(protein_channel_intensity_df, proteinId ~ sample, value.var = 'intensity', drop = F, fun.aggregate = mean)

sample_type = sapply(
  str_split(colnames(out_df), '_'),
  FUN = function(x)
    x[[2]]
)
sample_id = colnames(out_df)

out_df_header = rbind(sample_type, sample_id)

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
raw_frame = log2(raw_frame[,c(1,2,11,12,21,22,31,32,41,42,51,52,61,62)])
plt_sample_names = colnames(raw_frame)
colnames(raw_frame) = paste('s', 1:14,sep='')
raw_frame_melt = melt(raw_frame)
p1_1 = ggplot(data=raw_frame_melt %>% filter(value != 0), aes(x=Var2, y=value)) + 
  geom_boxplot() +
  mytheme1 +
  xlab('Sample') +
  ylab('log2(intensity)') +
  scale_y_continuous(limits=c(20,45))

p1_2 = ggplot(data=raw_frame_melt %>% filter(value != 0), aes(x=value, color=Var2)) + 
  geom_density(size=1) + 
  mytheme1 + 
  xlab('log2(intensity)') + 
  scale_color_manual(values=rev(mycolors)) +
  labs(col='Sample') +
  scale_x_continuous(limits=c(20,45))

ggsave('./normalization_plots/CNHPP/CNHPP_protein_boxplot_before_norm.jpg',p1_1,dpi=600,width = 18, height = 9)
ggsave('./normalization_plots/CNHPP/CNHPP_protein_density_before_norm.jpg',p1_2,dpi=600,width = 18, height = 9)

## sample loading normalization
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
  labs(col='Sample')+
  scale_x_continuous(limits=c(20,45))


ggsave('./normalization_plots/CNHPP/10samplesandref_taxonomy_boxplot_sl.jpg',p2_1,dpi=600,width = 18,height = 9)
ggsave('./normalization_plots/CNHPP/10samplesandref_taxonomy_density_sl.jpg',p2_2,dpi=600,width = 18,height = 9)

##write the data of sample loading normalization
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
  './03.normalized_protein_abundance/CNHPP.protein_abundance.sl.NAfiltered.txt.noheader',
  col.names = F,
  sep = '\t',
  quote = F
)

write.table(
  out_df_header,
  './03.normalized_protein_abundance/CNHPP.protein_abundance.sl.NAfiltered.txt.header',
  col.names = F,
  sep = '\t',
  quote = F
)
