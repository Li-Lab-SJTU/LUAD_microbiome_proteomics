mytheme1 <- theme_bw() +
  theme(
    axis.text = element_text(
      size = 20,
      family = "Arial",
      color = "black"
    ),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.border = element_rect(color = "white"),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1),
    axis.title.x = element_text(
      size = 22,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 22,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(size=18),
    legend.text = element_text(size=18),
    plot.title = element_text(hjust = 0.5,
                              size = 20,
                              family = "Arial",
                              face = "bold",
                              color = "black"
    )
    #legend.position = c(1,0.5)
  )
