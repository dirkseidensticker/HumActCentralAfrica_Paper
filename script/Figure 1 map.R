#### FIGURE 1 MAP ####

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

c14 <- data.table::fread("input/c14.csv", 
                         encoding = "UTF-8")

c14 <- sf::st_as_sf(c14, 
                    coords = c("LONG", "LAT"), 
                    remove = F, 
                    crs = 4326, 
                    na.fail = F)

sites	<- data.table::fread("input/sites.csv",
                           encoding = "UTF-8")

sites <- sf::st_as_sf(sites, 
                      coords = c("LONG", "LAT"), 
                      remove = F, 
                      crs = 4326, 
                      na.fail = F)

region.labels <- data.table::fread("input/region labels.csv", 
                                   encoding = "UTF-8")

region.labels <- sf::st_as_sf(region.labels, 
                              coords = c("LON", "LAT"), 
                              remove = F, 
                              crs = 4326, 
                              na.fail = F)

# ............................................................
# 2. LISTS OF SITES ----
# ............................................................

  # ............................................................
  # 2.1. Frequency of pottery group per site ----
  # ............................................................

sites.pottery.freq <- sites %>%
  dplyr::group_by(SITE, REGION) %>% 
  dplyr::summarise(Freq = length(POTTERY))

sites.pottery.freq <- st_cast(sites.pottery.freq, "POINT") # revertig MULTIPPOINTS back

  # ............................................................
  # 2.2. Only sites with dated unclassified potter ----
  # ............................................................

c14.pottery.indet <- dplyr::filter(c14, grepl("indet", POTTERY))

nrow(unique(c14.pottery.indet[,c("SITE", "LAT", "LONG")])) # number of sites with indet pottery

# ............................................................
# 3. FIGURE ----
# ............................................................

  # ............................................................
  # 3.1. Legend ----
  # ............................................................

legend <- sf::st_as_sf(data.frame(LAT = c(-6.6,-7.3,-8.0,-8.7,-9.4), 
                                  LONG = c(7.0,7.0,7.0,7.0,7.0), 
                                  Freq = c(20,1,5,10,20), 
                                  SIZE = c(1,2,2,2,2)),
                       coords = c("LONG", "LAT"), 
                       crs = 4326)

  # ............................................................
  # 3.2. Map ----
  # ............................................................

plt <- basemap() +
  geom_sf(data = sites.pottery.freq, 
          aes(alpha = Freq, 
              color = REGION), 
          size = 3) + 
  geom_sf(data = c14.pottery.indet, 
          aes(color = REGION), 
          size = 1) + 
  geom_sf_label(data = region.labels, 
                aes(label = REGION, 
                    fill = REGION)) + 
  geom_sf(data = legend[1,], 
          aes(alpha = Freq), 
          size = 1) + 
  geom_sf(data = legend[2:5,], 
          aes(alpha = Freq), 
          size = 3) + 
  coord_sf(xlim = c(7.5, 29), 
           ylim = c(-9.2, 6)) + 
  # styling & legend ----
  scale_fill_manual(values = region.labels$col) + 
  scale_color_manual(values = region.labels$col) + 
  annotate("text", x = 7.3, y = -6.6, label = paste("sites with unclassified pottery"),hjust=0,cex=2.5) +
  annotate("text", x = 7.3, y = -7.3, label = paste("sites with 1 pottery group"),hjust=0,cex=2.5) +
  annotate("text", x = 7.3, y = -8.0, label = paste("sites with 2-5 pottery groups"),hjust=0,cex=2.5) +
  annotate("text", x = 7.3, y = -8.7, label = paste("sites with 6-10 pottery groups"),hjust=0,cex=2.5) +
  annotate("text", x = 7.3, y = -9.4, label = paste("sites with >10 pottery groups"),hjust=0,cex=2.5) +
  scale_x_continuous(breaks = seq(8, 32, 2)) + 
  theme_few() + 
  theme(legend.position = "none", 
        axis.title = element_blank(),
        plot.background = element_rect(color = NA, 
                                       fill = NA))

windows() ; plt

ggsave("output/Figure 1 map.pdf", plt, 
       width = 8, height = 5.5)
ggsave("output/Figure 1 map.jpg", plt, 
       width = 8, height = 5.5)