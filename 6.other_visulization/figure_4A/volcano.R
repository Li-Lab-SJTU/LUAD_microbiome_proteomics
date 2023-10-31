library(EnhancedVolcano)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(stringr)
source("../../ggplot_theme/mytheme_basis.R")

# plot theme
mytheme1 <- theme_bw() +
  theme(
    axis.text = element_text(
      size = 18,
      family = "Arial",
      color = "black"
    ),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.border = element_rect(color = "white"),
    axis.line = element_line(color = "black", size = 2.5),
    axis.ticks = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(2, "mm"),
    axis.title.x = element_text(
      size = 18,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 18,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      family = "Arial",
      face = "bold",
      color = "black"
    )
  )

# Custom Color Gradient
custom_palette = c("#253392",
                   "#3fb3e7",
                   "#1c980b",
                   "#e4bf3a",
                   "#ff0000",
                   "#c00000"
)


################################Volcano Plot for CNHPP##############################
res1 = read.xlsx(
  "../../3.two_part_wilcoxon_test/test_result/test_result_CNHPP_go.xlsx"
)

res1$p.adjust = res1$p.adjust + 1e-20
res1[res1$log2fc > 2, "log2fc"] = 2
res1[res1$log2fc < -2, "log2fc"] = -2


# 绘制火山图
p1 = ggplot(res1, aes(
  x = log2fc,
  y = -log10(p.adjust),
  fill = -log10(p.adjust),
  size = -log10(p.adjust)
)) +
  geom_point(shape = 21,
             stroke = 1,
             color = "black") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  geom_vline(
    xintercept = -0.5,
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  labs(x = "log2(FC)", y = "-log10(P.adj)") +
  theme_minimal() +
  scale_fill_gradientn(colors = custom_palette, values = c(0, 0.03, 0.06, 0.1, 0.5, 1), limits = c(0, 20)) +
  mytheme1 +
  theme(legend.position = "bottom") +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = round(seq(0, 20, length.out = 5), 0),
    expand = c(0.03, 0.03)
  ) +
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = round(seq(-2, 2, length.out = 5), 0),
    expand = c(0.03, 0.03)
  ) +
  guides(size = FALSE,color = FALSE) +
  labs(fill = "-log10(P.adj)") +
  scale_color_gradientn(colors = custom_palette, values = c(0, 0.03, 0.06, 0.1, 0.5, 1), limits = c(0, 20))

ggsave("CNHPP.volcano.1019.png",p1,dpi=600,height=12,width=9)

################################Volcano Plot for CPTAC##############################
res1 = read.xlsx(
  "../../3.two_part_wilcoxon_test/test_result/test_result_CPTAC_go.xlsx"
  )

res1$p.adjust = res1$p.adjust + 1e-20
res1[res1$log2fc > 2, "log2fc"] = 2
res1[res1$log2fc < -2, "log2fc"] = -2

# 绘制火山图
p1 = ggplot(res1, aes(
  x = log2fc,
  y = -log10(p.adjust),
  fill = -log10(p.adjust),
  size = -log10(p.adjust)
)) +
  geom_point(shape = 21,
             stroke = 1,
             color = "black") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  geom_vline(
    xintercept = -0.5,
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  labs(x = "log2(FC)", y = "-log10(P.adj)") +
  theme_minimal() +
  scale_fill_gradientn(colors = custom_palette, values = c(0, 0.03, 0.06, 0.1, 0.5, 1), limits = c(0, 20)) +
  mytheme1 +
  theme(legend.position = "bottom") +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = round(seq(0, 20, length.out = 5), 0),
    expand = c(0.03, 0.03)
  ) +
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = round(seq(-2, 2, length.out = 5), 0),
    expand = c(0.03, 0.03)
  ) +
  guides(size = FALSE,color = FALSE) +
  labs(fill = "-log10(P.adj)") +
  scale_color_gradientn(colors = custom_palette, values = c(0, 0.03, 0.06, 0.1, 0.5, 1), limits = c(0, 20))

ggsave("CPTAC.volcano.png",p1,dpi=600,height=12,width=9)

################################Volcano Plot for TW##############################
res1 = read.xlsx(
  "D:/lilab_data/sy.03.run_smalldb/00.all_results_table/10.go_sample_statistical_test/TW_go_irs_sl_statistical_test_result_median.inputR1R2combinelcacf.r1r2sum.nofiltration.xlsx"
)

res1$p.adjust = res1$p.adjust + 1e-20
res1[res1$log2fc > 2, "log2fc"] = 2
res1[res1$log2fc < -2, "log2fc"] = -2

# 绘制火山图
p1 = ggplot(res1, aes(
  x = log2fc,
  y = -log10(p.adjust),
  fill = -log10(p.adjust),
  size = -log10(p.adjust)
)) +
  geom_point(shape = 21,
             stroke = 1,
             color = "black") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  geom_vline(
    xintercept = -0.5,
    linetype = "dashed",
    color = "grey",
    size = 1.2
  ) +
  labs(x = "log2(FC)", y = "-log10(P.adj)") +
  scale_fill_gradientn(colors = custom_palette, values = c(0, 0.03, 0.06, 0.1, 0.5, 1), limits = c(0, max(-log10(res1$p.adjust)))) +
  mytheme1 +
  theme(legend.position = "bottom") +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = round(seq(0, 20, length.out = 5), 0),
    expand = c(0.03, 0.03)
  ) +
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = round(seq(-2, 2, length.out = 5), 0),
    expand = c(0.03, 0.03)
  ) +
  guides(size = FALSE, color=FALSE) +
  labs(fill = "-log10(P.adj)") + 
  scale_color_gradientn(colors = custom_palette, values = c(0, 0.03, 0.06, 0.1, 0.5, 1), limits = c(0, max(-log10(res1$p.adjust))))

ggsave("TW.volcano.png",p1,dpi=600,height=12,width=9)