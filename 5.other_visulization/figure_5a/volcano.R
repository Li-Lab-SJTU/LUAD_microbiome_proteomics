library(EnhancedVolcano)
library(openxlsx)
library(patchwork)
library(ggplot2)
source("../../ggplot_theme/mytheme_basis.R")

getColKeyvals = function(res) {
  keyvals <- ifelse(
    res$log2fc < -0.5 & res$p.adjust < 0.05,
    '#008000',
    ifelse(res$log2fc > 0.5 & res$p.adjust < 0.05, '#ff0000',
           'black')
  )
  keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == '#ff0000'] <- 'Up-regulated'
    names(keyvals)[keyvals == 'black'] <- 'Not-sig'
    names(keyvals)[keyvals == '#008000'] <- 'Down-regulated'
  return(keyvals)
}


res1 = read.xlsx(
  "../../3.two_part_wilcoxon_test/test_result/test_result_TW_go.xlsx"
)

res1Kv = getColKeyvals(res1)

p1 = EnhancedVolcano(
  res1,
  lab = res1$feature,
  x = 'log2fc',
  y = 'p.adjust',
  xlim = c(-2, 2),
  title = 'TW',
  subtitle = '',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  xlab = bquote(~ Log[2] ~ 'FC'),
  ylab = bquote(~ -Log[10] ~ 'P.adj'),
  colCustom = res1Kv,
  colAlpha = 1,
  cutoffLineType = 'twodash',
  cutoffLineWidth = 0.8,
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 2,
  caption = "",
  labSize  = 4.0
) +
  mytheme1 +
  theme(
    plot.title = element_text(hjust = 0.5, size = ),
    legend.position = "bottom",
    panel.border = element_blank(),
  ) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = round(seq(0, 15, length.out = 5), 0),
    expand = c(0.03, 0.03)
  ) +
  labs(color = '')

res2 = read.xlsx(
  "../../3.two_part_wilcoxon_test/test_result/test_result_CPTAC_go.xlsx"
)
res2Kv = getColKeyvals(res2)
p2 = EnhancedVolcano(
  res2,
  lab = res2$feature,
  x = 'log2fc',
  y = 'p.adjust',
  xlim = c(-2, 2),
  title = 'CPTAC',
  subtitle = '',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  xlab = bquote(~ Log[2] ~ 'FC'),
  ylab = bquote(~ -Log[10] ~ 'P.adj'),
  colCustom = res2Kv,
  colAlpha = 1,
  cutoffLineType = 'twodash',
  cutoffLineWidth = 0.8,
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 2,
  caption = "",
  labSize  = 4.0
) +
  mytheme1 +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = round(seq(0, 20, length.out = 5), 0),
    expand = c(0.03, 0.03)
  )

res3 = read.xlsx(
  "../../3.two_part_wilcoxon_test/test_result/test_result_CNHPP_go.xlsx"
)
res3Kv = getColKeyvals(res3)
p3 = EnhancedVolcano(
  res3,
  lab = res3$feature,
  x = 'log2fc',
  y = 'p.adjust',
  xlim = c(-2, 2),
  title = 'CNHPP',
  subtitle = '',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  xlab = bquote(~ Log[2] ~ 'FC'),
  ylab = bquote(~ -Log[10] ~ 'P.adj'),
  colCustom = res3Kv,
  colAlpha = 1,
  cutoffLineType = 'twodash',
  cutoffLineWidth = 0.8,
  legendPosition = 'bottom',
  legendLabSize = 10,
  legendIconSize = 2,
  caption = "",
  labSize  = 4.0
)  +
  mytheme1 +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = round(seq(0, 20, length.out = 5), 0),
    expand = c(0.03, 0.03)
  )


ggsave(
  './volcano.png',
  dpi = 600,
  p3 + p2 + p1,
  width = 20,
  height = 10
)