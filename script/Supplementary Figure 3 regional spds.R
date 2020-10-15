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

min.14C.age <- 0 # rcarbonsum fails on 0-dates
timeRange <- c(4000, 0)
  
# ............................................................
# 3. CALIBRATE REGIONAL SPDs ----
# ............................................................

c14.sel <- dplyr::filter(c14, 
                         C14AGE > min.14C.age & # age filtering
                           C14STD > 0 & 
                           CLASS %in% c("Ia","Ib","Ic", "Id") & # class filtering
                           REGION %in% LETTERS[seq(from = 1, to = 11)] # regional filtering
)


kernel <- ceiling(mean(c14.sel$C14STD)/10)*10 # time window for running mean (in years)	based mean standard error		

id.lst <- unique(c14.sel$REGION)
dates.list  <- list()
for (i in 1:length(id.lst)) {
  a <- c14.sel[c14.sel$REGION %in% id.lst[i],]
  
  res <- rcarbon.spd(a, 
                     timeRange, 
                     kernel)
  
  res[["dat"]]$REGION <- id.lst[i]
  res[["dat"]]$n <- nrow(a)
  dates.list[[i]] <- res[["dat"]]
}

spd.sel <- do.call(rbind, dates.list)

# ............................................................
# 4. SUPPLEMENTARY FIGURE 3 - regional analysis in separate panels ----
# ............................................................

add <- dplyr::filter(c14, CLASS %in% c("Ia","Ib","Ic", "Id"))			
add <- as.data.frame(table(add$REGION))
region.labels <- merge(x = region.labels, by.x = "REGION",
                       y = add, by.y = "Var1")

colourCount = length(c(unique(spd.sel$REGION)))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

region.max <- spd.sel %>%
  dplyr::group_by(REGION) %>%
  dplyr::summarize(max = max(grid.PrDens))

region.labels <- merge(x = region.labels,
                       y = region.max,
                       by = "REGION")

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

  # ............................................................
  # 4.1. Main Plot ----
  # ............................................................

plt <- ggplot() + 
  #geom_vline(xintercept = seq(-2000, 2000, 100), color = "#f5f5f5") + 
  geom_rect(data = rcarbon, 
            aes(xmin = FROM, xmax = TO, ymin = 0, ymax = Inf, fill = rcarbon), 
            alpha = .1) + 
  scale_fill_discrete("", 
                      labels = c("less intense human activity", 
                                 "more intense human activity")) + 
  geom_line(data = spd.sel, 
            aes(x = calBCAD, y = grid.PrDens, color = REGION)) + 
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
             space = "free", 
             scales="free_y") + 
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
                   #ggdraw(cowplot::get_x_axis(axis, position = "bottom")), 
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

ggsave("output/Supplementary Figure 3 regional spds - separate.pdf", 
       plt, 
       width = 8, height = 9)
ggsave("output/Supplementary Figure 3 regional spds - separate.jpg", 
       plt,
       width = 8, height = 9)