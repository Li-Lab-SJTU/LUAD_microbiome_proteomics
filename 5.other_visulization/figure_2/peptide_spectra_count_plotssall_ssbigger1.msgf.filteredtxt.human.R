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
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1),
    axis.title.x = element_text(
      size = 14,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 14,
      family = "Arial",
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    plot.title = element_text(hjust = 0.5,
                              size = 14,
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
  p = p + mytheme1 + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(limits = c(5,30)) + xlab('-lg(SpecEValue)')
  return(p)
}

# Load MS-GF+ .tsv output file and generate psm score distribution plot
msgf_output_handle = function(msgf_output_path, plot_path) {
  print(paste('Handling', msgf_output_path))
  plot_frame = read.csv(msgf_output_path,
                      sep = '\t')
  p = generate_plot(plot_frame)
  ggsave(plot_path,p,dpi=600) 
}

# Each subplot in figure 2 takes the combined file of all the result MS-GF+ .tsv output 
#    files of one project as the input file. It was very large () that we cannot share it.
# Instead, we will use 100000 lines in combined TW result file as an example.
msgf_output_path = './TW.PSMs.score_distribution.test.tsv'
plot_path = './TW.PSMs.score_distribution.test.png'
msgf_output_handle(msgf_output_path, plot_path)
