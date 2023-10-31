library(dplyr)
library(vegan)
library(reshape)
library(ggplot2)
library(scales)
source('../../ggplot_theme/mytheme_taxa_boxplot.R')

cnhppdata = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.noheader',
  sep = '\t',
  row.names = 1
)
cnhppinfo = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.header',
  sep = '\t',
  row.names = 1
)
cptacdata = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CPTAC.taxon_abundance.irs_sl.NAfiltered.txt.noheader',
  sep = '\t',
  row.names = 1
)
cptacinfo = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CPTAC.taxon_abundance.irs_sl.NAfiltered.txt.header',
  sep = '\t',
  row.names = 1
)
twdata = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/TW.taxon_abundance.irs_sl.NAfiltered.txt.noheader',
  sep = '\t',
  row.names = 1
)
twinfo = read.table(
  '../../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/TW.taxon_abundance.irs_sl.NAfiltered.txt.header',
  sep = '\t',
  row.names = 1
)

dataSelect = function(data, info, tax, project) {
  data = data.frame(t(data[rownames(data) == tax,]))
  colnames(data)[1] = 'Intensity'
  data$Project = project
  data$Type = c(t(info['sample_type',]))
  return(data)
}

#################Trichoderma_g#################
cnhppx = dataSelect(cnhppdata, cnhppinfo, 'Trichoderma_g', 'CNHPP')
cptacx = dataSelect(cptacdata, cptacinfo, 'Trichoderma_g', 'CPTAC')
twx = dataSelect(twdata, twinfo, 'Trichoderma_g', 'TW')

plotdata = rbind(cptacx, twx, cnhppx)
plotdata$Type = factor(plotdata$Type, level = c('Tumor', 'NAT'))

data_text =
  data.frame(
    label = c("p.adj=9e-4", "p.adj=7e-6", "p.adj=0.02"),
    Project = c('TW', 'CPTAC', 'CNHPP'),
    Type = c('NAT', 'Tumor', 'NAT'),
    y = c(6.30, 5.60, 7.50)
  )

plotdata$Intensity[plotdata$Intensity == 0] = NA
plotdata$Intensity = log10(plotdata$Intensity)
p1 = ggplot(plotdata, aes(y = Intensity, x = Type, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f23355', '#339dff')) +
  mytheme1 +
  theme(axis.title.x = element_blank(), legend.position = 'None') +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  ggtitle('Trichoderma_g') +
  labs(y = 'Log10(intensity)') +
  facet_grid(Project ~ ., scales = "free_y") +
  geom_text(data = data_text,
            mapping = aes(y = y, label = label),
            size=5)

ggsave('boxplot_Trichoderma_g.jpg', p1, width = 4, height = 6, dpi=1200)

#################Staphylococcus_g#################
cptacx = dataSelect(cptacdata, cptacinfo, 'Staphylococcus_g', 'CPTAC')
twx = dataSelect(twdata, twinfo, 'Staphylococcus_g', 'TW')

plotdata = rbind(cptacx, twx)
plotdata$Type = factor(plotdata$Type, level = c('Tumor', 'NAT'))

data_text =
  data.frame(
    label = c("p.adj=0.03", "p.adj=0.02"),
    Project = c('TW', 'CPTAC'),
    Type = c('NAT', 'NAT'),
    y = c(5.75, 5.75)
  )

plotdata$Intensity[plotdata$Intensity == 0] = NA
plotdata$Intensity = log(plotdata$Intensity, base = 10)
p1 = ggplot(plotdata, aes(y = Intensity, x = Type, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f23355', '#339dff')) +
  mytheme1 +
  theme(axis.title.x = element_blank(), legend.position = 'None') +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  ggtitle('Staphylococcus_g') +
  labs(y = 'Log10(intensity)') +
  facet_grid(Project ~ ., scales = "free_y") +
  geom_text(data = data_text,
            mapping = aes(y = y, label = label),
            size=5)

ggsave('boxplot_Staphylococcus_g.jpg', p1, width = 4, height = 6, dpi=1200)

#################Debaryomycetaceae_f#################
cnhppx = dataSelect(cnhppdata, cnhppinfo, 'Debaryomycetaceae_f', 'CNHPP')
cptacx = dataSelect(cptacdata, cptacinfo, 'Debaryomycetaceae_f', 'CPTAC')

plotdata = rbind(cptacx, cnhppx)
plotdata$Type = factor(plotdata$Type, level = c('Tumor', 'NAT'))

data_text =
  data.frame(
    label = c("p.adj=3e-5", "p.adj=0.03"),
    Project = c('CPTAC', 'CNHPP'),
    Type = c('Tumor', 'Tumor'),
    y = c(5.8, 11.0)
  )

plotdata$Intensity[plotdata$Intensity == 0] = NA
plotdata$Intensity = log(plotdata$Intensity, base = 10)
p1 = ggplot(plotdata, aes(y = Intensity, x = Type, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f23355', '#339dff')) +
  mytheme1 +
  theme(axis.title.x = element_blank(), legend.position = 'None') +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  ggtitle('Debaryomycetaceae_f') +
  labs(y = 'Log10(intensity)') +
  facet_grid(Project ~ ., scales = "free_y") +
  geom_text(data = data_text,
            mapping = aes(y = y, label = label),
            size=5)

ggsave('boxplot_Debaryomycetaceae_f.jpg', p1, width = 4, height = 6, dpi=1200)

#################Pasteurellaceae_f#################
cptacx = dataSelect(cptacdata, cptacinfo, 'Pasteurellaceae_f', 'CPTAC')
twx = dataSelect(twdata, twinfo, 'Pasteurellaceae_f', 'TW')

plotdata = rbind(cptacx, twx)
plotdata$Type = factor(plotdata$Type, level = c('Tumor', 'NAT'))

data_text =
  data.frame(
    label = c("p.adj=2e-3", "p.adj=0.04"),
    Project = c('TW', 'CPTAC'),
    Type = c('NAT', 'NAT'),
    y = c(5.75, 5.15)
  )

plotdata$Intensity[plotdata$Intensity == 0] = NA
plotdata$Intensity = log(plotdata$Intensity, base = 10)
p1 = ggplot(plotdata, aes(y = Intensity, x = Type, fill = Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f23355', '#339dff')) +
  mytheme1 +
  theme(axis.title.x = element_blank(), legend.position = 'None') +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  ggtitle('Pasteurellaceae_f') +
  labs(y = 'Log10(intensity)') +
  facet_grid(Project ~ ., scales = "free_y") +
  geom_text(data = data_text,
            mapping = aes(y = y, label = label),
            size=5)

ggsave('boxplot_Pasteurellaceae_f.jpg', p1, width = 4, height = 6, dpi=1200)