#### SUPPLEMENTARY FIGURE 5 EFFECTIVE POPOLTATION SIZE ####

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

rcarbonA 	<- read.csv("output/rcarbonTest_timewindows_a-h.csv", 
                      encoding = "UTF-8")

coords <- data.table::fread("input/TabS3_coords.csv", 
                            encoding = "UTF-8")

coords <- sf::st_as_sf(coords, 
                       coords = c("Longitude",
                                  "Latitude"),
                       remove = F, 
                       crs = 4326, 
                       na.fail = F)

  # ............................................................
  # 1.1. READ INPUT FROM GENETIC RECONSTRUCTIO ----
  # ............................................................

f <- list.files("input/outputs_Ne_WC_Africa/")

for(i in 1:length(f)){
  d <- read.table(paste0("input/outputs_Ne_WC_Africa/", f[i]), 
                  sep="\t", header = T)
  
  d$CTR <- sub("\\_.*", "", f[i])
  d$POP <- str_extract(f[i], "_([^_]*)_")
  d$POP <- sub("_", "", d$POP)
  d$POP <- sub("_", "", d$POP)
  
  d$ID <- paste(d$POP, d$CTR, sep = " - ")
  
  if(i==1){
    data <- d
  }else{
    data <- rbind(data, d)
  }
}

data$calBCAD <- 2000 - data$GEN * 30 # calculate cal BC/AD age

  # ............................................................
  # 1.2. change prefixes in population/language names ----
  # ............................................................

dict = list(Akele = "Kele", 
            Bakota = "Kota", 
            Bapunu = "Punu",
            Bateke = "Teke",
            Eshira = "Shira",
            Eviya = "Viya",
            Makina = "Shiwe (Fang-Makina)",
            Obamba = "Bamba",
            Orungu = "Rungu")

for(i in 1:length(dict)){
  data$POP[data$POP == names(dict[i])] <- dict[[i]]
}

# remove generations <4 ----
data <- dplyr::filter(data, GEN >= 4)

# ............................................................
# 2. FILTERING INDIVIDUAL GROUPS ----
# ............................................................

selectionA <- dplyr::filter(data, POP %in% c("Kele",
                                             "Kota", 
                                             "Benga", 
                                             "Viya", 
                                             "Fang", 
                                             "Galoa", 
                                             "Shiwe (Fang-Makina)",
                                             "Shake", 
                                             "Tsogo"))

selectionB <- dplyr::filter(data, POP %in% c("Teke", 
                                             "Duma", 
                                             "Shira",
                                             "Bamba", 
                                             "Punu", 
                                             "Ndumu",
                                             "Nzebi"))

coordsA <- subset(coords, Community %in% unique(selectionA$POP))
coordsA$GRP <- "A"
coordsB <- subset(coords, Community %in% unique(selectionB$POP))
coordsB$GRP <- "B"

coordsGRP <- rbind(coordsA, coordsB)

# ............................................................
# 3. SUPPLEMENTARY FIGURE 5 ----
# ............................................................

  # ............................................................
  # 2.1 Effective Polulation Sizes (Ne) ----	
  # ............................................................

plt.ne <- ggplot() + 
  geom_rect(data = rcarbonA, 
            aes(xmin = FROM, 
                xmax = TO, 
                ymin = 0, 
                ymax = Inf, 
                fill = rcarbon), 
            alpha = .1) + 
  geom_line(data = selectionA, 
            aes(x = calBCAD, 
                y = NE, 
                group = ID), 
            color = "#006400") + 
  geom_line(data = selectionB, 
            aes(x = calBCAD, 
                y = NE, 
                group = ID), 
            color = "#cc5500") + 
  scale_x_continuous("cal BC/AD", 
    breaks = c(seq(-2000, 1500, 500), 1950),
    expand = c(0, 0), 
    sec.axis = sec_axis(~ (2000 - .) /30, 
                        name = "Generations ago",
                        breaks = seq(125, 0, -25))) + 
  coord_cartesian(xlim = c(-2100, 2000),
                  ylim = c(1e+3, 1e+11)) + 
  scale_y_log10("Effective population size, log10(Ne)", 
                expand = c(0, 0))+#, 
  annotation_logticks(sides = "l") + 
  guides(color = guide_legend(ncol = 2)) + 
  theme_classic() + 
  theme(legend.position = "none")

  # ............................................................
  # 2.2 panel B: Map ----	
  # ............................................................

map <- basemap() + 
  geom_label_repel(data = coordsGRP,
                   aes(x = Longitude, 
                       y = Latitude,
                       fill = GRP,
                       label = Community), 
                   size = 2, 
                   color = "white",
                   min.segment.length = 0,
                   segment.color = 'black') + 
  labs(fill = "") + 
  guides(fill = guide_legend(override.aes = aes(label = ""), 
                             ncol = 1)) + 
  scale_fill_manual(labels = c("North-Western Bantu", 
                               "West-Western Bantu"), 
                    values = c("#006400",
                               "#cc5500")) + 
  ggsn::north(coordsGRP, 
              anchor = c(x = 9.1, y = -2.8), 
              scale = .2) + 
  ggsn::scalebar(coordsGRP,
                 location  = "topright",
                 anchor = c(x = 10.2, y = -3.8),
                 dist = 100, dist_unit = "km",
                 transform = TRUE, model = "WGS84", 
                 height = .02, st.dist = .05, 
                 border.size = .1, st.size = 3) + 
  coord_sf(xlim = c(8.5, 14.5), 
           ylim = c(-4, 2.5), 
           label_graticule = "SE") + 
  theme_few() + 
  theme(axis.title = element_blank(), 
        legend.position = "bottom")

  # ............................................................
  # 2.3 Legend Element ----	
  # ............................................................

plt.legend1 <- ggplot() + 
  geom_rect(data = rcarbonA, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .3) + 
  scale_fill_manual("",values = c("#00bfc4", "#f8766d"), 
                    labels = c("more intense human activity", 
                               "less intense human activity")  ) +
  theme(legend.position = "top")
legend1 <- cowplot::get_legend(plt.legend1)

  # ............................................................
  # 3.3 Combine Elements ----	
  # ............................................................

plt1 <- cowplot::plot_grid(legend1, 
                           plt.ne,
                           ncol = 1, 
                           rel_heights = c(1, 10))

plt <- cowplot::plot_grid(plt1,
                          map, 
                          ncol = 2, 
                          rel_widths = c(1.5,1), 
                          labels = "auto")

ggsave("output/Figure 4 - Effective population size.pdf", 
       plt, 
       width = 8, height = 6)
ggsave("output/Figure 4 - Effective population size.jpg", 
       plt, 
       width = 8, height = 6)
