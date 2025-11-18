#' Get scatter plot from input dataframe
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom cowplot insert_xaxis_grob
#' @importFrom cowplot insert_xaxis_grob
#' @importFrom cowplot ggdraw
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 theme_classic
#' @importFrom grid unit
#' @importFrom ggthemes theme_few
#'

plot_scatter_box <- function(dt_sb,
                             var_x,
                             var_y,
                             col_g,
                             xlab,
                             ylab,
                             title_lab) {
  colors_fill = c(Reference = "#2f5c85", Query = "red")
  # colors_fill = c(Reference= "#2f5c85")
  pmain <- ggplot(dt_sb, ggplot2::aes_string(x = var_x, y = var_y, color = col_g)) +
    geom_point() +
    scale_color_manual(values = colors_fill) +
    theme_few() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    labs(title = title_lab, x = xlab, y = ylab)
  
  # pmain <- pmain + geom_text(data = subset(dt_sb, dt_sb[[col_g]] == "Query"),
  #                            ggplot2::aes_string(label = "batch", x = var_x, y = var_y),
  #                            nudge_y = 0.5, # 根据需要调整这个值来避免标签和点重叠
  #                            color = "red")
  
  xplot <- ggplot(dt_sb, ggplot2::aes_string(x = col_g, y = var_x, colour = col_g)) +
    geom_boxplot() +
    scale_color_manual(values = colors_fill) +
    coord_flip() +
    theme_classic()
  
  yplot <- ggplot(dt_sb, ggplot2::aes_string(x = col_g, y = var_y, colour = col_g)) +
    geom_boxplot() +
    scale_color_manual(values = colors_fill) +
    theme_classic()
  
  p1 <- cowplot::insert_xaxis_grob(pmain, xplot, grid::unit(.2, "null"), position = "top")
  p2 <- cowplot::insert_yaxis_grob(p1, yplot, grid::unit(.2, "null"), position = "right")
  pt_sb <- ggdraw(p2)
  return(pt_sb)
}
