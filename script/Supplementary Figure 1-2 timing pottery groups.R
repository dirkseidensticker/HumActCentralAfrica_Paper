
#### SUPPLEMENTARY FIGURE 1 - TIMEFRAME PER REGION PER DATED POTTERY POTTERY  ####

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

c14 <- data.table::fread("input/c14.csv", 
                         encoding = "UTF-8")
po.gr <- data.table::fread("input/potterygroups.csv", 
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

c14.sel <- dplyr::filter(c14, 
                         C14AGE > min.14C.age & # age filtering
                         C14STD > 0 & 
                         CLASS %in% c("Ia","Ib","Ic", "Id") & # class filtering
                         REGION %in% LETTERS[seq(from = 1, to = 11)] # regional filtering
)

# ............................................................
# 3. CALIBRATE SPD PER POTTERY GROUP  ----
# ............................................................

id <- po.gr$POTTERY
res.lst <- list()
for(i in 1:length(id)){
  d <- dplyr::filter(c14.sel, 
                     grepl(id[i], 
                           c14.sel$POTTERY)) # filter for dates related to style
  d <- dplyr::filter(d, 
                    !grepl(paste0("\\(" , id[i], "\\)"), 
                                  d$POTTERY)) # remove cases in parantheses
  
  if(nrow(d)!=0){ # only if 14C-dates exist
    res <- rcarbon.spd(d, 
                       timeRange, 
                       kernel = ncores, 
                       oxcalnorm = T,
                       median = T)
    
    res[["dat"]]$median <- list(res[["median"]])
    res[["dat"]]$po.gr <- id[i]
    res[["dat"]]$n <- nrow(d)
    
    res.lst[[i]] <- res[["dat"]]
  }
}
po.prob <- do.call(rbind, res.lst)

# merge with po.gr
po.prob <- merge(x = po.prob, by.x = "po.gr", 
                 y = po.gr[,-c("DESCRIPTION")], by.y = "POTTERY")

# extract median
po.prob.med <- as.data.frame(tidyr::unnest(po.prob, median))
po.prob.med <- unique(po.prob.med[c("po.gr", "median", "REGION")])


# dummy figure for legend ----

plt.legend <- ggplot() + 
  geom_rect(data = rcarbonA, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .3) + 
  scale_fill_manual("",values = c("#00bfc4", "#f8766d"), 
                    labels = c("more intense human activity", 
                               "less intense human activity")  ) +
  theme(legend.position = "top")
legend <- cowplot::get_legend(plt.legend)

# ............................................................
# 4. SUPPLEMENTARY FIGURE 1  ----
# ............................................................

reg.labs <- po.prob %>% 
  dplyr::group_by(REGION) %>% 
  dplyr::slice(which.max(FROM))

reg.labs <- merge(x = reg.labs, 
                  y = region.labels, 
                  by = "REGION")

reg.labs$lab <- paste0(reg.labs$REGION, ") ", reg.labs$LABEL)

axis <- ggplot() +
  annotate(geom = "segment", 
           x = c(seq(-1200, 1900, 100), 1950), 
           xend = c(seq(-1200, 1900, 100), 1950),
           y = 0, yend = -1) + 
  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500), 1950),
                     limits = c(-1200, 2000), 
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
# 4.1 Regions a-h  ----
# ............................................................

rcarbonA$FROM 	<- replace(rcarbonA[,c("FROM")],rcarbonA[,c("FROM")]<(-1199),-1199)
rcarbonA$TO 		<- replace(rcarbonA[,c("TO")],rcarbonA[,c("TO")]>1999,1999)

plt.A <- ggplot() + 
  geom_line(data = dplyr::filter(po.prob, 
                          REGION %in% LETTERS[seq(from = 1, to = 8)]), 
            aes(x = calBCAD,
                y = reorder(po.gr, FROM),
                color = grid.PrDens), 
            size = 2) + 
  geom_rect(data = rcarbonA, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .1) +
  scale_fill_discrete("", 
                      labels = c("less intense human activity", 
                                 "more intense human activity")) + 
  geom_point(data = filter(po.prob.med, 
                           REGION %in% LETTERS[seq(from = 1, to = 8)]), 
             aes(x = median, y = po.gr), 
             color = "black", fill = "white", shape = 21, size = 1) +   
  scale_colour_gradient(low = "white", 
                        high = "black", 
                        guide = "none") + 
  geom_text(data = filter(reg.labs, 
                          REGION %in% LETTERS[seq(from = 1, to = 8)]), 
            aes(x = -1150, 
                y = po.gr, 
                label = lab), 
            vjust = .5, hjust = 0, 
            size = 2.5) + 
  scale_x_continuous("cal BC/AD", 
                     limits = c(-1200,2000), 
                     breaks = c(seq(-1500,1500,500), 1950), 
                     expand = c(0,0)) + 
  scale_y_discrete(position = "right") + 
  facet_grid(REGION ~ ., 
             scales = "free", 
             space = "free", 
             switch="both") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none", 
        strip.text.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE))

# ............................................................
# 4.2 Regions i-k  ----
# ............................................................

plt.B <- ggplot() + 
  geom_line(data = dplyr::filter(po.prob,
                                 REGION %in% LETTERS[seq(from = 9, to = 11)]), 
            aes(x = calBCAD,
                y = reorder(po.gr, FROM),
                color = grid.PrDens), 
            size = 2) + 
  geom_rect(data = rcarbonB, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .1) +
  scale_fill_discrete("", 
                      labels = c("less intense human activity", 
                                 "more intense activity")) + 
  geom_point(data = filter(po.prob.med, 
                           REGION %in% LETTERS[seq(from = 9, to = 11)]), 
             aes(x = median, y = po.gr), 
             color = "black", fill = "white", shape = 21, size = 1) +   
  
  scale_colour_gradient(low = "white", 
                        high = "black", 
                        guide = "none") + 
  geom_text(data = filter(reg.labs, 
                          REGION %in% LETTERS[seq(from = 9, to = 11)]), 
            aes(x = -1150, 
                y = po.gr, 
                label = lab), 
            vjust = .5, hjust = 0, 
            size = 2.5) + 
  scale_x_continuous("cal BC/AD", 
                     limits = c(-1200,2000), 
                     breaks = c(seq(-1500,1500,500), 1950), 
                     expand = c(0,0)) + 
  scale_y_discrete(position = "right") + 
  facet_grid(REGION ~ ., 
             scales = "free", 
             space = "free", 
             switch="both") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none", 
        strip.text.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE))


print(paste(length(unique(po.prob$po.gr)), "directly dated pottery groups"))

# ............................................................
# 4.3 Joint plot & Save Figure ----
# ............................................................

plt.row1 <- cowplot::plot_grid(NULL, NULL,
                               ncol = 1, 
                               rel_heights = c(14, 3.75),
                               labels = "auto")

plt.row2 <- cowplot::plot_grid(plt.A, plt.B, axis,
                               ncol = 1, 
                               rel_heights = c(14*1.25, 3*1.25, 1),
                               align = "v", axis = "lr")


plt.base <- cowplot::plot_grid(plt.row1, plt.row2, 
                          ncol = 2, 
                          rel_widths = c(1,30))

plt <- cowplot::plot_grid(legend, 
                          plt.base, 
                          ncol = 1, 
                          rel_heights = c(1,30))


#windows() ; plt

ggsave("output/Supplementary Figure 1 overview dated_stylegroups.pdf", 
       plt, 
       width = 8, height = 10)
ggsave("output/Supplementary Figure 1 overview dated_stylegroups.jpg", 
       plt,
       width = 8, height = 10)

# ............................................................
# 5. REPORT NUMBERS  ----
# ............................................................

length(unique(po.prob$po.gr)) # total nr of directly dated styles

sum(as.data.table(table(po.prob.med$po.gr))$N) # total nr of 14C dates associated to styles in Supp Fig 1

# total nr of 14C dates associated to styles in adrac; this number is lower than the one above, as some dates were used for two or more styles
nrow(dplyr::filter(c14, 
                   POTTERY != '-' &
                     POTTERY != 'indet' &
                     CLASS == "Ia"))

#### SUPPLEMENTARY FIGURE 2 - TIMEFRAME PER REGION PER UNDATED POTTERY POTTERY  ####

# ............................................................
# 6 SUBSET POTTERY GROUPS WITHOUT 14C DATES ----
# ............................................................

# groups with no 14C-dates associated to them
po.gr.14c <- unique(po.prob$po.gr)
po.gr.all <- unique(po.gr$POTTERY)
po.gr.diff <- po.gr.all[which(!po.gr.all %in% po.gr.14c)]

po.gr.rel <- base::subset(po.gr[,-c("DESCRIPTION")], POTTERY %in% po.gr.diff) %>%
  filter(REGION != '')

# ............................................................
# 7 PLOT ----
# ............................................................

# ............................................................
# 7.1 Labels ----
# ............................................................

reg.labs.rel <- po.gr.rel %>%
  dplyr::arrange(desc(POTTERY)) %>% # reverse order to match ggplot sorting
  dplyr::group_by(REGION) %>% 
  dplyr::slice(which.max(FROM))

reg.labs.rel <- merge(x = reg.labs.rel, 
                      y = region.labels, 
                      by = "REGION")

reg.labs.rel$lab <- paste0(reg.labs.rel$REGION, ") ", reg.labs.rel$LABEL)

# ............................................................
# 7.2 Figure ----
# ............................................................

plt.graph <- ggplot() + 
  geom_segment(data = dplyr::filter(po.gr.rel,
                                    REGION %in% LETTERS[seq(from = 1, to = 8)]), 
               aes(x = FROM, 
                   y = reorder(POTTERY, FROM), 
                   xend = TO, 
                   yend = POTTERY), 
               size = 3, alpha = 0.5) + #, linetype = "11") + 
  geom_rect(data = rcarbonA, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .1) +
  scale_fill_discrete("", 
                      labels = c("lower human activity", 
                                 "higher human activity")) + 
  geom_text(data = filter(reg.labs.rel, 
                          REGION %in% LETTERS[seq(from = 1, to = 8)]), 
            aes(x = -1150, 
                y = POTTERY, 
                label = lab), 
            vjust = .5, hjust = 0, 
            size = 2.5) + 
  scale_y_discrete(position = "right") + 
  scale_x_continuous("cal BC/AD", 
                     limits = c(-1200,2000), 
                     breaks = c(seq(-1500,1500,500), 1950), 
                     expand = c(0, 0)) +
  facet_grid(REGION ~ ., 
             scales = "free", 
             space = "free", 
             switch="both") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none", 
        strip.text.y = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE))


plt.axis <- cowplot::plot_grid(plt.graph,
                               axis,
                               ncol = 1,
                               rel_heights = c(20, 1.6), 
                               align = "v", axis = "lr")

plt <- cowplot::plot_grid(legend, 
                          plt.axis,
                          ncol = 1,
                          rel_heights = c(1, 20))

windows() ; plt

ggsave("output/Supplementary Figure 2 overview undated_stylegroups.pdf", 
       plt, 
       width = 8, height = 6)
ggsave("output/Supplementary Figure 2 overview undated_stylegroups.jpg", 
       plt,
       width = 8, height = 6)

