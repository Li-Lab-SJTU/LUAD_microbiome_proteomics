library(openxlsx)
library(patchwork)
library(dplyr)
library(tidyverse)
library(GO.db)
library(ggplot2)
library(ggtext)
library(reshape2)
source("../../ggplot_theme/mytheme_basis.R")

## get the significant up- or down- regulated GOs
twdf = read.xlsx("../../3.two_part_wilcoxon_test/test_result/test_result_TW_go.xlsx")
cptacdf = read.xlsx("../../3.two_part_wilcoxon_test/test_result/test_result_CPTAC_go.xlsx")
cnhppdf = read.xlsx("../../3.two_part_wilcoxon_test/test_result/test_result_CNHPP_go.xlsx")

getLabel = function(df) {
  ispval = df$p.adjust < 0.05
  islogfc = abs(df$log2fc) > 0.5
  isup = df$direction == 'Up'
  isdown = df$direction == 'Down'
  
  df[!ispval, 'label'] = 'NS'
  df[ispval & !islogfc, 'label'] = 'Sig'
  df[ispval & islogfc & isup, 'label'] = 'Up'
  df[ispval & islogfc & isdown, 'label'] = 'Down'
  return(df)
}

twdf = getLabel(twdf)
cptacdf = getLabel(cptacdf)
cnhppdf = getLabel(cnhppdf)

twdf = twdf %>% filter(label %in% c('Up', 'Down'))
cptacdf = cptacdf %>% filter(label %in% c('Up', 'Down'))
cnhppdf = cnhppdf %>% filter(label %in% c('Up', 'Down'))

## Annotation of go terms
godb = read.csv("go_database.txt", sep = '\t')

twdf = merge(twdf, godb, by.x = 'feature', by.y = 'GOID')
cptacdf = merge(cptacdf, godb, by.x = 'feature', by.y = 'GOID')
cnhppdf = merge(cnhppdf, godb, by.x = 'feature', by.y = 'GOID')
twdf$TERM = paste(twdf$TERM, ' [',twdf$ONTOLOGY,']',sep='')
cptacdf$TERM = paste(cptacdf$TERM, ' [',cptacdf$ONTOLOGY,']',sep='')
cnhppdf$TERM = paste(cnhppdf$TERM, ' [',cnhppdf$ONTOLOGY,']',sep='')
twdf$GoTerms = paste(twdf$feature, ' [',twdf$ONTOLOGY,']',sep='')
cptacdf$GoTerms = paste(cptacdf$feature, ' [',cptacdf$ONTOLOGY,']',sep='')
cnhppdf$GoTerms = paste(cnhppdf$feature, ' [',cnhppdf$ONTOLOGY,']',sep='')

# Only plot GO terms identified in more than one cohort
plotlist = read.table('go_morethan1cohort.txt')
plotlist = c(plotlist$V1)

get_plot = function(go_df, go_all_intensity_df, meta_df, cohort){
  morethan1_df = go_df %>% filter(feature %in% plotlist)
  go_all_intensity_df$feature = row.names(go_all_intensity_df)
  
  morethan1_df = merge(morethan1_df, go_all_intensity_df, by = 'feature')
  
  plot_df = melt(
    morethan1_df,
    id.vars = c('TERM', 'GoTerms', 'direction'),
    measure.vars = meta_df['sample_id', ],
    variable.name = 'sample',
    value.name = 'intensity'
  )
  
  if (cohort == 'TW'){
    plot_df$type = sapply(str_split(plot_df$sample, '_'), function(x)
      x[2])
    plot_df$type = gsub('Normal Adjacent Tissue', 'NAT', plot_df$type)
    plot_df$type = factor(plot_df$type, levels = c('Tumor','NAT'))
  } else if (cohort == 'CPTAC'){
    plot_df$type = sapply(str_split(plot_df$sample, '_'), function(x)
      x[2])
    plot_df$type = gsub('Primary Tumor', 'Tumor', plot_df$type)
    plot_df$type = gsub('Solid Tissue Normal', 'NAT', plot_df$type)
    plot_df$type = factor(plot_df$type, levels = c('Tumor','NAT'))
  } else if (cohort == 'CNHPP'){
    plot_df$type = sapply(str_split(plot_df$sample, '_'), function(x)
      x[2])
    plot_df$type = factor(plot_df$type, levels = c('Tumor','NAT'))
  }
  
  if (cohort == 'TW'){
    p1 = ggplot(plot_df, aes(
      x = GoTerms,
      y = log2(intensity),
      fill = type
    )) + 
      geom_boxplot(na.rm = T) + facet_grid(.~direction, scales = 'free_x', space = "free_x") + 
      mytheme1 + 
      theme(plot.title = element_text(hjust = 0.5),legend.position="bottom", axis.text.x = element_text(size = 12,color = "black", angle=45, vjust=1, hjust=1)) + 
      labs(y='log2(intensity)', x='', fill='', title = cohort) +
      scale_y_continuous(limits=c(15,23),breaks=round(seq(15,23,length.out=5),1),expand = c(0.03,0.03))+
      scale_fill_manual(values = c('#f23355', '#339dff'))
  } else if (cohort == 'CPTAC'){
    p1 = ggplot(plot_df, aes(
      x = GoTerms,
      y = log2(intensity),
      fill = type
    )) + 
      geom_boxplot(na.rm = T) + facet_grid(.~direction, scales = 'free_x', space = "free_x") + 
      mytheme1 + 
      theme(plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x = element_text(size = 12,color = "black", angle=45, vjust=1, hjust=1)) + 
      labs(x = '', y='log2(intensity)', fill='', title = cohort) +
      scale_y_continuous(limits=c(15,23),breaks=round(seq(15,23,length.out=5),1),expand = c(0.03,0.03)) +
      scale_fill_manual(values = c('#f23355', '#339dff'))
  } else if (cohort == 'CNHPP'){
    p1 = ggplot(plot_df, aes(
      x = GoTerms,
      y = log2(intensity),
      fill = type
    )) + 
      geom_boxplot(na.rm = T) + facet_grid(.~direction, scales = 'free_x', space = "free_x") + 
      mytheme1 + 
      theme(plot.title = element_text(hjust = 0.5), legend.position="none", axis.text.x = element_text(size = 12,color = "black", angle=45, vjust=1, hjust=1)) + 
      labs(x = '', y='log2(intensity)', fill='', title = cohort) +
      scale_y_continuous(limits=c(24,40),breaks=round(seq(24,40,length.out=5),1),expand = c(0.03,0.03)) +
      scale_fill_manual(values = c('#f23355', '#339dff'))
  }
  return(p1)
}

tw_go_all_intensity_df = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/TW.go_abundance.irs_sl.NAfiltered.txt.noheader',
  sep = '\t',
  row.names = 1
)
tw_meta_df = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/TW.go_abundance.irs_sl.NAfiltered.txt.header',
  sep = '\t',
  row.names = 1
)
colnames(tw_go_all_intensity_df) = tw_meta_df['sample_id', ]
p1 = get_plot(twdf, tw_go_all_intensity_df, tw_meta_df, cohort = 'TW')

cptac_go_all_intensity_df = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/CPTAC.go_abundance.irs_sl.NAfiltered.txt.noheader',
  sep = '\t',
  row.names = 1
)
cptac_meta_df = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/CPTAC.go_abundance.irs_sl.NAfiltered.txt.header',
  sep = '\t',
  row.names = 1
)
colnames(cptac_go_all_intensity_df) = cptac_meta_df['sample_id', ]
p2 = get_plot(cptacdf, cptac_go_all_intensity_df, cptac_meta_df, cohort = 'CPTAC')

cnhpp_go_all_intensity_df = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/CNHPP.go_abundance.sl.NAfiltered.txt.noheader',
  sep = '\t',
  row.names = 1
)
cnhpp_meta_df = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/CNHPP.go_abundance.sl.NAfiltered.txt.header',
  sep = '\t',
  row.names = 1
)
colnames(cnhpp_go_all_intensity_df) = cnhpp_meta_df['sample_id', ]
p3 = get_plot(cnhppdf, cnhpp_go_all_intensity_df, cnhpp_meta_df, cohort = 'CNHPP')

ggsave('./boxplot.go.morthan1cohort.png', dpi=600, p3+p2+p1, width =20, height = 10)