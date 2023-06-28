# survival analysis between the disease free survival (DFS) / overall survival (OS) of CNHPP LUAD patients 
# and microbial taxon abundance in tumors.

# Load libraries
library(dplyr)
library(vegan)
library(reshape)
library(openxlsx)
library(ggplot2)
library(survival)
library(survminer)

# ggplot theme for survival plot
source("../ggplot_theme/mytheme_survival_plot.R")

# Load CNHPP taxon abundance data matrix and the corresponding sample infor mation   
data = read.table(
  "../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.noheader",
  sep = '\t',
  row.names = 1
)

info = read.table(
  "../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.header",
  sep = '\t',
  row.names = 1
)

# Load the CNHPP clinical information
clinicalInfo = read.xlsx(
  "./CNHPP_clinical_mmc1.xlsx",
  sheet = 2
)

# Try to merge the survive sample info
info['ID',] = sapply(strsplit(t(info['sample_id',]), '_'), function(x) x[1])
info['ID',] = as.numeric(info['ID',])
clinicalInfo$ID = sapply(strsplit(clinicalInfo$Sample.ID, '_'), function(x) x[2])
clinicalInfo$ID = as.numeric(clinicalInfo$ID)

x = rbind(info, data)
x = data.frame(t(x))
x = merge(x,
          clinicalInfo[, c('ID', 'DFS', 'DFS_month', 'OS', 'OS_month')],
          by = 'ID',
          all.x = F,
          all.y = F)
x = data.frame(t(x))
test_df = data.frame(t(x))

# Only analyze microbial abundance in tumor samples
test_df = test_df %>% filter(sample_type == 'Tumor')
test_df$DFS_month = as.numeric(test_df$DFS_month)
test_df$DFS = as.numeric(test_df$DFS)
test_df$OS_month = as.numeric(test_df$OS_month)
test_df$OS = as.numeric(test_df$OS)

DFS_result = data.frame()
OS_result = data.frame()

# Univariate Cox proportional hazards models analysis for DFS
for (i in 4:68) {
  taxa = colnames(test_df)[i]
  test_new = cbind(test_df[, c('ID', 'DFS', 'DFS_month', 'OS', 'OS_month')], test_df[, i])
  colnames(test_new)[6] = 'X'
  test_new$X = log2(as.numeric(test_new$X))
  res.cox = coxph(Surv(DFS_month, DFS) ~ X, data = test_new)
  res.sum = summary(res.cox)
  DFS_result = rbind(
    DFS_result,
    data.frame(
      taxa = taxa,
      numSample = res.sum$n,
      numEvent = res.sum$nevent,
      coef = res.sum$coefficients[1],
      HR = res.sum$conf.int[, "exp(coef)"],
      HR.95L = res.sum$conf.int[, "lower .95"],
      HR.95H = res.sum$conf.int[, "upper .95"],
      pvalue = res.sum$coefficients[, "Pr(>|z|)"]
    )
  )
}

# Univariate Cox proportional hazards models analysis for OS
for (i in 4:68) {
  taxa = colnames(test_df)[i]
  test_new = cbind(test_df[, c('ID', 'DFS', 'DFS_month', 'OS', 'OS_month')], test_df[, i])
  colnames(test_new)[6] = 'X'
  test_new$X = log2(as.numeric(test_new$X))
  res.cox = coxph(Surv(OS_month, OS) ~ X, data = test_new)
  res.sum = summary(res.cox)
  OS_result = rbind(
    OS_result,
    data.frame(
      taxa = taxa,
      coef = res.sum$coefficients[1],
      numSample = res.sum$n,
      numEvent = res.sum$nevent,
      HR = res.sum$conf.int[, "exp(coef)"],
      HR.95L = res.sum$conf.int[, "lower .95"],
      HR.95H = res.sum$conf.int[, "upper .95"],
      pvalue = res.sum$coefficients[, "Pr(>|z|)"]
    )
  )
}

# Save the result of Cox models
write.xlsx(
  DFS_result,
  "./CNHPP_DFS_result.xlsx"
)
write.xlsx(
  OS_result,
  "./CNHPP_OS_result.xlsx"
)

OS_result_filtered = OS_result %>% filter(pvalue < 0.05)
DFS_result_filtered = DFS_result %>% filter(pvalue < 0.05)

bioForest = function(coxFile = null,
                     forestFile = null,
                     height = null,
                     width = null) {
  rt = coxFile
  gene <- rt$taxa
  hr <- sprintf("%.3f", rt$"HR")
  hrLow  <- sprintf("%.3f", rt$"HR.95L")
  hrHigh <- sprintf("%.3f", rt$"HR.95H")
  Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
  pVal <-
    ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  jpeg(
    file = forestFile,
    width = width,
    height = height,
    units = "px",
    res = 600
  )
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
  
  xlim = c(0, 3)
  par(mar = c(4, 2, 2, 1))
  plot(
    1,
    xlim = xlim,
    ylim = ylim,
    type = "n",
    axes = F,
    xlab = "",
    ylab = ""
  )
  text.cex = 1
  text(0, n:1, gene, adj = 0, cex = text.cex - 0.1)
  text(1.8, n:1, pVal, adj = 1, cex = text.cex)
  text(1.8,
       n + 1,
       'pvalue',
       cex = text.cex,
       font = 2,
       adj = 1)
  text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
  text(
    3.1,
    n + 1,
    'Hazard ratio(95% CI)',
    cex = text.cex,
    font = 2,
    adj = 1
  )
  
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(
    1,
    xlim = xlim,
    ylim = ylim,
    type = "n",
    axes = F,
    ylab = "",
    xaxs = "i",
    xlab = "Hazard ratio"
  )
  arrows(
    as.numeric(hrLow),
    n:1,
    as.numeric(hrHigh),
    n:1,
    angle = 90,
    code = 3,
    length = 0.05,
    col = "darkblue",
    lwd = 3.5
  )
  abline(
    v = 1,
    col = "black",
    lty = 2,
    lwd = 3
  )
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr),
         n:1,
         pch = 15,
         col = boxcolor,
         cex = 3)
  axis(1)
  dev.off()
}

bioForest(
  coxFile = OS_result_filtered,
  forestFile = "./Cox_bioForest_plots/OS_survival.jpeg",
  width = 6000,
  height = 1900
)
bioForest(
  coxFile = DFS_result_filtered,
  forestFile = "./Cox_bioForest_plots/DFS_survival.jpeg",
  width = 6000,
  height = 1900
)

#OS surv KM plot
OS_surv_data = test_df[, c(
  "OS",
  "OS_month",
  "Cytomegalovirus_g",
  "Human.betaherpesvirus.5_s",
  "Schizosaccharomyces.pombe_s"
)]
OS_surv_data$Cytomegalovirus_g = as.numeric(OS_surv_data$Cytomegalovirus_g)
OS_surv_data$Human.betaherpesvirus.5_s = as.numeric(OS_surv_data$Human.betaherpesvirus.5_s)
OS_surv_data$Schizosaccharomyces.pombe_s = as.numeric(OS_surv_data$Schizosaccharomyces.pombe_s)

res.cut <-
  surv_cutpoint(
    OS_surv_data,
    time = "OS_month",
    event = "OS",
    variables = c(
      "Cytomegalovirus_g",
      "Human.betaherpesvirus.5_s",
      "Schizosaccharomyces.pombe_s"
    )
  )
res.cat <- surv_categorize(res.cut)

for (name in c("Cytomegalovirus_g",
               "Human.betaherpesvirus.5_s",
               "Schizosaccharomyces.pombe_s")) {
  fit <- survfit(Surv(OS_month, OS) ~ res.cat[, name],  
                 data = res.cat)
  p1 = ggsurvplot(
    fit,
    data = res.cat,
    pval = TRUE,
    size = 0.5,
    surv.median.line = "hv",
    xlab = "OS(month)",
    legend.title = "",
    font.x = c(10, "bold", "black"),
    font.y = c(10, "bold", "black"),
    font.legend = c(10, "bold", "black"),
    font.tickslab = 10,
    pval.size = 3,
    pval.method.size = 3,
    censor.size = 4,
    legend.labs = c(paste0(name, "=High"), paste0(name, "=Low")),
    break.x.by = 20,
    xlim = c(0, 80),
    risk.table = T
  )
  
  ggsave(
    paste0('./KM_plots/OS_', name, '.png'),
    p1$plot,
    width = 7,
    height = 3,
    dpi = 600
  )
  
  
}

#DFS surv KM plot
DFS_surv_data = test_df[, c("DFS", "DFS_month", "Firmicutes_p", "Gammaproteobacteria_c")]
DFS_surv_data$Firmicutes_p = as.numeric(DFS_surv_data$Firmicutes_p)
DFS_surv_data$Gammaproteobacteria_c = as.numeric(DFS_surv_data$Gammaproteobacteria_c)

res.cut <-
  surv_cutpoint(
    DFS_surv_data,
    time = "DFS_month",
    event = "DFS",
    variables = c("Firmicutes_p", "Gammaproteobacteria_c")
  )
res.cat <- surv_categorize(res.cut)

for (name in c("Firmicutes_p", "Gammaproteobacteria_c")) {
  fit <- survfit(Surv(DFS_month, DFS) ~ res.cat[, name],  
                 data = res.cat)
  p1 = ggsurvplot(
    fit,
    data = res.cat,
    pval = TRUE,
    size = 0.5,
    surv.median.line = "hv",
    xlab = "DFS(month)",
    legend.title = "",
    font.x = c(10, "bold", "black"),
    font.y = c(10, "bold", "black"),
    font.legend = c(10, "bold", "black"),
    font.tickslab = 10,
    pval.size = 3,
    pval.method.size = 3,
    censor.size = 4,
    legend.labs = c(paste0(name, "=High"), paste0(name, "=Low")),
    xlim = c(0, 80),
    break.x.by = 20,
    risk.table = T
  )
  
  ggsave(
    paste0('./KM_plots/DFS_', name, '.png'),
    p1$plot,
    width = 7,
    height = 3,
    dpi = 600
  )
}
