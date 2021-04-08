constructGraphingUtils <- function() {
  
  public <- list()
  
  public$LARGE = "L"
  public$MEDIUM = "M"
  public$SMALL = "S"
  
  # ==========
  # ggplot
  # ==========
  
  public$ggplot <- function(..., size = public$LARGE, axes = TRUE) {
    textSize = if (size == public$LARGE) 20 else if (size == public$MEDIUM) 15 else 10
    plot <- ggplot(...) + theme_minimal()
    if (axes) {
      plot <- plot + 
        theme(text = element_text(size = textSize), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line.x = element_line(color="grey50"),
              axis.line.y = element_line(color="grey50"))
    } else {
      plot <- plot + 
        theme(text = element_text(size = textSize), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }
    plot
  }
  
  public$tiltX <- function(angle = 60, hjust = 1, vjust = 0.3) {
    theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust))
  }
  
  public$ggplot_continousColorGradient <- function() {
    scale_color_gradient2(low = "black", mid = "yellow", high = "white")
  }
  
  public$ggplot_tintBackground <- function() {
    theme(panel.background = element_rect(fill = "gray95", color = "white"))
  }
  
  # ==========
  # pheatmap
  # ==========
  
  public$heatmap <- function(datamatrix, decimalNums = 0, size = public$LARGE, ...) {
    sizeNum <- if (size == public$LARGE) 20 else if (size == public$MEDIUM) 13 else 8
    datamatrix %>% 
      pheatmap(display_numbers = decimalNums > 0,
               number_format = paste0("%.", as.character(decimalNums), "f"),
               fontsize = sizeNum,
               fontsize_number = sizeNum - 2,
               colorRampPalette(c("black", "yellow"))(100),
               ...)
  }
  
  return(public)
  
}