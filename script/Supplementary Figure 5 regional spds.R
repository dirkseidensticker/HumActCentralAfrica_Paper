#### SUPPLEMENTARY FIGURE 3 REGIONAL SPDs ####

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

c14 <- data.table::fread("input/c14.csv", 
                         encoding = "UTF-8")
region.labels <- data.table::fread("input/region labels.csv", 
                                   encoding = "UTF-8")
rcarbonA <- data.table::fread("output/rcarbonTest_timewindows_a-h.csv", 
                               encoding = "UTF-8")
rcarbonB <- data.table::fread("output/rcarbonTest_timewindows_i-k.csv", 
                              encoding = "UTF-8")

# ............................................................
# 2. SETUP ----
# ............................................................

timeRange = c(4000, 0)

n.bins <- 100 # time window for binning dates (in years) within the same site	
n.simulations <- 100 # nr of simulations

runm <- ceiling(mean(c14$C14STD)/10)*10 # time window for running mean (in years)	based mean standard error

# ............................................................
# 3. CALIBRATE ----
# ............................................................

c14 <- c14 %>%
  dplyr::filter(
    C14AGE > 0 & 
    C14STD > 0 & 
    CLASS %in% c("Ia","Ib","Ic", "Id") & # class filtering
    REGION %in% LETTERS[seq(from = 1, to = 11)] # regional filtering
  )

# calibrate & binning

cal <- rcarbon::calibrate(
  x = c14$C14AGE, 			
  errors = c14$C14STD,
  calCurves = 'intcal20',
  ncores = ncores, 
  normalised = FALSE)

bins <- binPrep(sites = c14$SITE,
                ages = c14$C14AGE,
                h = n.bins # xx years cut off value
) 

# run permTest ----
perm <- permTest(x = cal,
                 marks = c14$REGION,
                 timeRange = timeRange,
                 bins = bins,
                 nsim = n.simulations,
                 runm = runm)

summary(perm)

# extract raw results from permTest ----
index <- unique(c14$REGION)
reg.lst <- list()
bb.lst <- list()
for(i in 1:length(index)){
  #i = 1
  # extract SPD
  spd.regional <- perm$observed[[i]]
  spd.regional$REGION <- index[i]
  
  # extract & cbind CI
  spd.ci.regional <- cbind(spd.regional, 
                           as.data.frame(
                           perm$envelope[[i]])
  )
  
  spd.ci.regional$calBCAD <- 1950 - spd.ci.regional$calBP # convert to BC/AD
  
  reg.lst[[i]] <- spd.ci.regional
}
reg.res <- do.call(rbind, reg.lst)

# ............................................................
# 4. SUPPLEMENTARY FIGURE 3 - regional analysis in separate panels ----
# ............................................................

add <- dplyr::filter(c14, CLASS %in% c("Ia","Ib","Ic", "Id"))			
add <- as.data.frame(table(add$REGION))
region.labels <- merge(x = region.labels, by.x = "REGION",
                       y = add, by.y = "Var1")

colourCount = length(c(unique(reg.res$REGION)))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

region.max <- reg.res %>%
  dplyr::group_by(REGION) %>%
  dplyr::summarize(max = max(PrDens))

region.labels <- merge(x = region.labels,
                       y = region.max,
                       by = "REGION")

# rcarbon ----
rcarbon <- rcarbonA

rcarbon$REGION <- NA
for(i in 1:8){
  a <- rcarbonA
  a$REGION <- LETTERS[i]
  rcarbon <- rbind(rcarbon, a)
}
for(i in 9:11){
  a <- rcarbonB
  a$REGION <- LETTERS[i]
  rcarbon <- rbind(rcarbon, a)
}
rcarbon <- dplyr::filter(rcarbon, 
                         FROM > -2000 & TO < 1900 &
                          !is.na(REGION))

# dummy figure for legend
plt.legend <- ggplot() + 
  geom_rect(data = rcarbon, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .3) + 
  scale_fill_manual("",values = c("#00bfc4", "#f8766d"), 
                    labels = c("more intense human activity", 
                               "less intense human activity")  ) +
  theme(legend.position = "top")
legend <- cowplot::get_legend(plt.legend)

# plot ----
plt <- ggplot() + 
  geom_rect(data = rcarbon, 
            aes(xmin = FROM, xmax = TO, ymin = 0, ymax = Inf, fill = rcarbon), 
            alpha = .1) + 
  scale_fill_discrete("", 
                      labels = c("less intense human activity", 
                                 "more intense human activity")) + 
  geom_ribbon(data = reg.res, 
              aes(x = calBCAD, 
                  y = PrDens,
                  ymin = V1,
                  ymax = V2), 
              alpha = .2) + 
  geom_line(data = reg.res, 
            aes(x = calBCAD, 
                y = PrDens, 
                color = REGION)) + 
  scale_color_manual(values = region.labels$col,
                     guide = FALSE) + 
  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500), 1950),
                     limits = c(-2000, 2000), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Summed probability", expand = c(0, 0), 
                     breaks = scales::pretty_breaks(2), 
                     #limits = c(0, NA),
                     position = "left") + 
  geom_text(data = region.labels, 
            aes(x = -1950, y = .5*max, 
                label = paste0(region.labels$REGION, ") ", region.labels$LABEL, " (n = ", region.labels$Freq , ")"), color = REGION), hjust = 0, size = 3) + 
  facet_grid(REGION ~ ., 
             scales = "free", 
             space = "free") + 
  theme_classic() + 
  theme(legend.position = "none",
        #axis.line.y = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        #axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE))

  # ............................................................
  # 4.2. Individual Plot Elements ----
  # ............................................................

axis <- ggplot() +
  annotate(geom = "segment", 
           x = c(seq(-2000, 1900, 100), 1950), 
           xend = c(seq(-2000, 1900, 100), 1950),
           y = 0, yend = -1) + 
  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500), 1950),
                     limits = c(-2000, 2000), 
                     expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() + 
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "lines")
  )

  # ............................................................
  # 4.2. Combine Elements & Save Figure ----
  # ............................................................

plt.label <- cowplot::plot_grid(NULL, NULL, 
                                nrow = 2, 
                                rel_heights = c(16,3),
                                labels = "auto")

plt <- cowplot::plot_grid(plt, 
                   axis,
                   ncol = 1, 
                   rel_heights = c(20,1),
                   align = "v", axis = "lr")

plt <- cowplot::plot_grid(plt.label, plt, 
                          ncol = 2, 
                          rel_widths = c(1,30))


plt <- cowplot::plot_grid(legend, 
                          plt, 
                          ncol = 1, 
                          rel_heights = c(1,20))

windows() ; plt

ggsave("output/Supplementary Figure 5 regional spds - separate.pdf", 
       plt, 
       width = 8, height = 9)
ggsave("output/Supplementary Figure 5 regional spds - separate.jpg", 
       plt,
       width = 8, height = 9)