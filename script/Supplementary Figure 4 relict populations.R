#### SUPPLEMENTARY FIGURE 4 MAP OF RELICT POPULATIONS ####

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

c14 <- data.table::fread("input/c14.csv", 
                         encoding = "UTF-8")
c14.orig <- c14

c14 <- sf::st_as_sf(c14, 
                    coords = c("LONG", "LAT"), 
                    remove = F, 
                    crs = 4326, 
                    na.fail = F)

sites	<- data.table::fread("input/sites.csv",
                           encoding = "UTF-8")

rcarbon <- data.table::fread("input/rcarbonTest_timewindows_a-h.csv", 
                              encoding = "UTF-8")

phases <- data.table::fread("output/rcarbon_phases.csv", 
                            encoding = "UTF-8")

# ............................................................
# 2. SELECT MIA HIATUS PHASE ----
# ............................................................

mia.hiatus <- phases %>%
  dplyr::filter(PHASE == "low activity 2")

# ............................................................
# 3. CALIBRATE & SELECT DATES ----
# ............................................................

c14 <- dplyr::filter(c14, 
                     C14AGE > 0 &
                     C14STD > 0 & 
                     CLASS %in% c("Ia","Ib","Ic", "Id"))

cal <- rcarbon::calibrate(x = c14$C14AGE,
                          errors = c14$C14STD,
                          calCurves = 'intcal13', 
                          ncores = ncores, 
                          normalised = FALSE) #running calibration over 3 cores

cal.median <- rcarbon::medCal(cal) # return individual medians
cal.median <- 1950 - cal.median

c14$calBCAD <- cal.median

c14.hiatus <- dplyr::filter(c14, 
                              calBCAD > mia.hiatus$FROM & 
                              calBCAD < mia.hiatus$TO)

# ............................................................
# 4. SITES ----
# ............................................................

c14.sites <- sf::st_as_sf(unique(as.data.frame(c14[,c("SITE","LAT","LONG")])), 
                          coords = c("LONG", "LAT"), 
                          crs = 4326, 
                          remove = F)

hiatus.sites <- sf::st_as_sf(unique(as.data.frame(c14.hiatus[,c("SITE","LAT","LONG")])), 
                             coords = c("LONG", "LAT"), 
                             crs = 4326, 
                             remove = F)

# sites + c14.sites - hiatus.sites
# --------------------------------

sites <- unique(sites[,c("SITE", "LAT", "LONG")]) %>%
  sf::st_as_sf(coords = c("LONG", "LAT"), 
               remove = F, 
               crs = 4326, 
               na.fail = F)

sites.all <- rbind(sites, c14.sites)
sf::st_geometry(sites.all) <- NULL

sites.all <- sf::st_as_sf(unique(sites.all[,c("SITE","LAT","LONG")]), 
                          coords = c("LONG", "LAT"), 
                          crs = 4326, 
                          remove = F)
sf::st_geometry(sites.all) <- NULL

# sites from the hiatus 
h <- hiatus.sites
sf::st_geometry(h) <- NULL
h <- data.table(h)

not.hiatus.sites <- sites.all[!h, on=.(SITE, LAT, LONG)]

not.hiatus.sites <- sf::st_as_sf(unique(not.hiatus.sites[,c("SITE","LAT","LONG")]), 
                                 coords = c("LONG", "LAT"), 
                                 crs = 4326, 
                                 remove = F)

# how many sites?
sf::st_geometry(sites) <- NULL
sf::st_geometry(c14.sites) <- NULL

sites.doubl <- merge(x = sites, y = c14.sites, 
                     by = c("SITE","LAT","LONG"))

cat("",nrow(sites), "sites with known pottery groups \n", 
    nrow(c14.sites), "sites with 14C dates of class I \n",
    nrow(sites.doubl), "are in both lists! \n", 
    nrow(sites.all), "sites with described pottery or class I dates")

# small check of numbers: should be TRUE ;)
nrow(sites.all) - nrow(hiatus.sites) == nrow(not.hiatus.sites)
sites.all <- sf::st_as_sf(unique(sites.all[,c("SITE","LAT","LONG")]), 
                          coords = c("LONG", "LAT"), 
                          crs = 4326, 
                          remove = F)

# ............................................................
# 5. SUPPLEMENTARY FIGURE 4 ----
# ............................................................

# ............................................................
# 3.1. Minimap (insert) ----
# ............................................................


world <- ne_countries(scale = "medium", returnclass = "sf")

minimap <- ggplot(data = world) +
  geom_sf(color = NA, fill = "grey") + 
  geom_rect(xmin = 7, xmax = 30, 
            ymin = -10, ymax = 6.5, 
            fill = NA, color = "black") + 
  coord_sf(xlim = c(-15, 50), 
           ylim = c(-35, 35)) + 
  theme_void() + 
  theme(panel.border = element_rect(colour = "darkgrey", 
                                    fill = NA, size = .5))

  # ............................................................
  # 5.1. Map ----
  # ............................................................

plt.main <- basemap() +
  geom_sf(data = sites.all, 
          fill = "#8a8a8a",
          color = "black",
          shape = 21,
          size = 2) + 
  geom_sf(data = hiatus.sites, 
          fill = "#8a8a8a",
          color = "black",
          shape = 21,
          size = 2) + 
  geom_sf(data = hiatus.sites, 
          color = "red", 
          size = 1) + 
  ggsn::north(sites.all, 
              anchor = c(x = 29.5, y = 6.5)) + 
  ggsn::scalebar(sites.all,
                 location  = "topright",
                 anchor = c(x = 27, y = 6),
                 dist = 250, dist_unit = "km",
                 transform = TRUE, model = "WGS84", 
                 height = .01, st.dist = .025, 
                 border.size = .1, st.size = 3) + 
  geom_point(aes(x = 7, y = -7), shape = 21, size = 2, color = "black", fill = "#8a8a8a") + 
  annotate("text", x = 7.3, y = -7, label = paste0("all sites (Class I)"),hjust=0,cex=2.5) +
  annotate("text", x = 7.3, y = -7.5, label = paste0("2000 BC to AD 1900 (n = ",nrow(sites.all),")"), hjust = 0, cex = 2.5) +
  geom_point(aes(x = 7, y = -8.2), shape = 21, size = 2, color = "black", fill = "#8a8a8a") + 
  geom_point(aes(x = 7, y = -8.2), size = 1, color = "red") + # #f8766d
  annotate("text", x = 7.3, y = -8.2, label = paste0("sites occupied during"),hjust=0,cex=2.5) +
  annotate("text", x = 7.3, y = -8.7, label = paste0("AD 600-1000 (n = ", nrow(hiatus.sites), ")"),hjust=0,cex=2.5) +
  scale_x_continuous(breaks = seq(8, 32, 2)) + 
  coord_sf(xlim = c(7.5, 29), 
           ylim = c(-9.2, 6)) + 
  theme_few() + 
  theme(legend.position = "none", 
        axis.title = element_blank())

# ............................................................
# 3.3. Combine Map & Insert ----
# ............................................................

plt <- ggdraw() +
  draw_plot(plt.main) +
  draw_plot(minimap, 
            x = .05, y = .275, width = .15, height = .15)

windows() ; plt

ggsave("output/Supplementary Figure 4 relict_populations.pdf", plt, 
       width = 8, height = 5.5)
ggsave("output/Supplementary Figure 4 relict_populations.jpg", plt, 
       width = 8, height = 5.5)