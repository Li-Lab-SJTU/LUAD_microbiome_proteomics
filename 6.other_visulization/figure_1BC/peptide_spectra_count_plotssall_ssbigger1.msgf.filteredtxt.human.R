library(stringr)
library(getopt)
library(ggplot2)
library(dplyr)

# plot theme
mytheme1 <- theme_bw() +
  theme(
    axis.text = element_text(
      size = rel(1),
      family = "Arial",
      color = "black"
    ),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.border = element_rect(color = "white"),
    axis.line = element_line(color = "black", size = 2),
    axis.ticks = element_line(color = "black", size = 2),
    axis.ticks.length = unit(2, "mm"),
    axis.title.x = element_text(
      size = 16,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 16,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16),
    plot.title = element_text(hjust = 0.5,
                              size = 16,
                              family = "Arial",
                              face = "bold",
                              color = "black"
    )
  )

# plot function
generate_plot = function(plot_frame){
  p = ggplot() + geom_density(
    data = plot_frame,
    aes(x = negative_lg_SpecEValue, fill = Class),
    color = "black",
    alpha = 0.5
  )
  p = p + mytheme1 + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_continuous(limits = c(5,30)) + 
    xlab('-lg(SpecEValue)') +
    theme(legend.position = "bottom")
  return(p)
}

# Load MS-GF+ .tsv output file and generate psm score distribution plot
msgf_output_handle = function(msgf_output_path, plot_path) {
  print(paste('Handling', msgf_output_path))
  plot_frame = read.csv(msgf_output_path,
                      sep = '\t')
  p = generate_plot(plot_frame)
  ggsave(plot_path,p,dpi=600,height=4,width=5)
  return(p)
}

# Each subplot in figure 2 takes the combined file of all the result MS-GF+ .tsv output 
#    files of one project as the input file. It was very large () that we cannot share it.
# Instead, we will use 100000 lines in combined TW result file as an example.
msgf_output_path = './TW.PSMs.score_distribution.test.tsv'
plot_path = './TW.PSMs.score_distribution.test.legendbottom.png'
p=msgf_output_handle(msgf_output_path, plot_path)
