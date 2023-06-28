library(dplyr)
library(stringr)
library(openxlsx)
library(reshape2)
library(ggplot2)
source("../../ggplot_theme/mytheme_basis.R")
                      
plot_df = read.xlsx('./taxa_abundance_ratio_staked_column_plot.xlsx')
plot_df$label = factor(plot_df$label, levels = rev(plot_df$label))

plot_df_2 = melt(plot_df, id.vars = 'label', measure.vars = c('TW_NAT','TW_Tumor',
                                                              'CPTAC_NAT','CPTAC_Tumor',
                                                              'CNHPP_NAT','CNHPP_Tumor'))
p1 = ggplot(plot_df_2,
            aes(variable, value)) +
  geom_col(aes(fill = label), position = 'fill', width = 0.6) +
  mytheme1 +
  labs(x = '', y = 'Relative Abundance(%)', fill = 'Order')  +
  scale_fill_manual(values=rev(plot_df$col)) +
  scale_y_continuous(expand = c(0, 0.01)) +
  coord_flip() +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

ggsave(
  "all_project_microbiome_distribution_order_phylum.png",
  p1,
  dpi = 300,
  width = 20,
  height = 10,
)