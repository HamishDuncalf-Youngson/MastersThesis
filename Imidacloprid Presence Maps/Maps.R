rm(list=ls())


install.packages('maps')
library(ggplot2)
library(maps)
library(patchwork)


worldmap = map_data('world')
knitr::kable(head(worldmap, 20))

###### Making basemap of the UK ########
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group))

ggplot() + geom_polygon(data = worldmap, 
                        aes(x = long, 
                            y = lat, 
                            group = group)) + 
  coord_fixed(xlim = c(-10,3), 
              ylim = c(50.3, 59))

ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, 
                   y = lat, 
                   group = group)) + 
  coord_fixed(ratio = 1.3, 
              xlim = c(-10,3), 
              ylim = c(50, 59))

ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, 
                   group = group), 
               fill = 'gray90', 
               color = 'black') + 
  coord_fixed(ratio = 1.3, 
              xlim = c(-10,3), 
              ylim = c(50, 59))

ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, 
                   group = group), 
               fill = 'gray90', 
               color = 'black') + 
  coord_fixed(ratio = 1.3, 
              xlim = c(-10,3), 
              ylim = c(50, 59)) + 
  theme_void()

################ Putting in coordinate data ###########

S <- read.csv("Sites .csv")
library(sf)
points <- st_as_sf(S, coords = c("Long", "Lat"))


All <- ggplot() +
  geom_polygon(
    data = worldmap,
    aes(x = long, y = lat,
        group = group),
    fill = 'gray97',
    color = 'black'
  ) + geom_point(data = S, 
             mapping = aes(x = Long, y = Lat), shape = 21, fill = "blue", size= 1.2) +
  
  coord_fixed(ratio = 1.3,
              xlim = c(-10, 3),
              ylim = c(50, 59)) +
  theme_void() + labs(title = "Sites sampled")

All

#### Subset with rivers that have imidacloprid present

Imidarivers <- subset(S, Name %in% c('Wandle','Crane', 'Mildenhall Trout Farm' ) )

Imida <- ggplot() +
  geom_polygon(
    data = worldmap,
    aes(x = long, y = lat,
        group = group),
    fill = 'gray97',
    color = 'black'
  ) + geom_point(data = Imidarivers, 
                 mapping = aes(x = Long, y = Lat), shape = 21, fill = "red", size= 1.2) + coord_fixed(ratio = 1.3,
              xlim = c(-10, 3),
              ylim = c(50, 59)) +
  theme_void() + labs(title = "Imidacloprid detected")

Imida

Both <- All + Imida

Both

ggsave("/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Maps/Both.pdf", Both, dpi=1000)





    