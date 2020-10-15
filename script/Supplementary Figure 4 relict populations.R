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

sites <- sf::st_as_sf(sites, 
                      coords = c("LONG", "LAT"), 
                      remove = F, 
                      crs = 4326, 
                      na.fail = F)

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

# TODO ----
# - sites and c14-sites merge by unique name & lat & lon

c14.sites <- sf::st_as_sf(unique(c14[,c("SITE","LAT","LONG")]), 
#c14.sites <- sf::st_as_sf(unique(c14.orig[,c("SITE","LAT","LONG")]), 
                          coords = c("LONG", "LAT"), 
                          crs = 4326, 
                          remove = F)

hiatus.sites <- sf::st_as_sf(unique(c14.hiatus[,c("SITE","LAT","LONG")]), 
                             coords = c("LONG", "LAT"), 
                             crs = 4326, 
                             remove = F)

# sites + c14.sites - hiatus.sites

sites <- sites[,c("SITE", "LAT", "LONG")]

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
  # 5.1. Map ----
  # ............................................................

plt <- basemap() +
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
  #scale_colour_gradient2(low = "green", mid = "red", high = "blue",
  #                       midpoint = spd.hiatus.min, na.value = NA) + 
  # annotation
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

windows() ; plt

ggsave("output/Supplementary Figure 4 relict_populations.pdf", plt, 
       width = 8, height = 5.5)
ggsave("output/Supplementary Figure 4 relict_populations.jpg", plt, 
       width = 8, height = 5.5)