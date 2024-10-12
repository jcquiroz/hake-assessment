armo_talla <- function(x, breaks) {
  t = hist(x, breaks, plot = FALSE)
  tibble(
    bind = head(t[[1]],-1), mids = t[[4]], counts = t[[2]], 
    dens = t[[3]], freq = t[[2]]/sum(t[[2]])
  )
}

mitheme <- function(bz = 8){
  theme_bw(base_size = bz) + theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", colour = "#3A3333"), #family = "sans"  
    
    #axis.line = element_line(colour = "#3A3333"),
    axis.title.y = element_text(size = rel(1.1), face = "bold", colour = "#3A3333", vjust = 2.65),
    axis.title.x = element_text(size = rel(1.1), face = "bold", colour = "#3A3333"),
    
    plot.title = element_text(hjust = 0.0, size = rel(1.3)),
    plot.background = element_rect(fill = "#FEFFFE", colour = "#FEFFFE"),
    plot.margin = margin(0.5, 1, 1, 0.5, "cm"),
    
    strip.text = element_text(size = rel(0.8), face = "bold", colour = "#3A3333"),
    strip.placement.x = "inside", #strip.placement.y = "inside",
    strip.background = element_rect(colour = "white", fill = "white")
  )
}


mitheme2 <- function(bz = 8){
  theme_bw(base_size = bz) + theme(
    legend.position = "bottom",
    
    axis.title.y = element_text(vjust = 2.65),
    
    strip.text = element_text(size = rel(0.8), face = "bold", colour = "#3A3333"),
    strip.placement.x = "inside", #strip.placement.y = "inside",
    strip.background = element_rect(colour = "white", fill = "white")
  )
}

mycolor <- viridisLite::viridis(12, option = "turbo", begin = 0, end = 1)
mycolor2 <- viridisLite::viridis(7, option = "turbo", begin = 0, end = 1)


