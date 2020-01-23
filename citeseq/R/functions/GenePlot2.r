library(MASS)
library(viridis)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



GenePlot2 = function (object, gene1, gene2, pch.use = 20,  
                      use.imputed = FALSE, use.scaled = FALSE,  use.raw = FALSE, 
                      plot.title = NULL, ...) {
  cell.ids <- object@cell.names
  data.use <- as.data.frame(x = FetchData(object = object, 
                                          vars.all = c(gene1, gene2), cells.use = cell.ids, use.imputed = use.imputed, 
                                          use.scaled = use.scaled, use.raw = use.raw))
  data.plot <- data.use[cell.ids, c(gene1, gene2)]
  names(x = data.plot) <- c("x", "y")
  data.plot = data.plot %>% dplyr::mutate(density = get_density(data.plot$x, data.plot$y, n = 100))
  
  
  p <- ggplot2::ggplot(data = data.plot, mapping = aes(x = x, y = y))
  
  p <- p + geom_point(aes(x,y, shape = pch.use, color = density), size = 0.3) + 
    scale_color_viridis() +  
    scale_shape_identity() 
  p <- p + labs(title = plot.title, x = gene1, y = gene2) + theme_light()
  return(p)     
}

