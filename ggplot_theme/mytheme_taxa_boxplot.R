mytheme1 = theme_bw() +
  theme(
    axis.text = element_text(
      size = rel(1),
      family = "Arial",
      color = "black"
    ),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.border = element_rect(color = "black", size = 1),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.title.x = element_text(
      size = 12,
      family = "Arial",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 16,
      family = "Arial",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 16,
      family = "Arial",
      color = "black",
      angle=90, 
      hjust=0.5
    ),
    axis.text.y = element_text(
      size = 16,
      family = "Arial",
      color = "black"
    ),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    strip.text = element_text(size = 12),
  )

##14ɫɫֵ
mycolors = c('#4dbbd5', '#00a087', '#3c5488', 
             '#e64b35', '#f39b7f', '#8491b4', 
             '#91d1c2', '#FF0000', '#7e6148', 
             '#b09c85', '#668b8b', '#FFD39B',
             '#EE6AA7', '#25670F')
