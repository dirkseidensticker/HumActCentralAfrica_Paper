#### FIGURE 2 OVERALL SPD ####

setwd("C:/Wannes/UGent Data/4 POSTDOC KMMA/R/radiocarbon - human activity/CentralAfricaHumanActivity v8_DS/")

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

c14 <- data.table::fread("input/c14.csv", 
                         encoding = "UTF-8")

# ............................................................
# 2. SETTINGS ----
# ............................................................

min.14C.age <- 50 # rcarbonTest fails on very small dates
n.bins <- 100 # time window for binning dates (in years) within the same site	
n.simulations <- 100 # nr of simulations
timeRange = c(4000, 0)# set the timerange of analysis in calBP, older date first
breaks = seq(4000, 0, -200) # 200 year blocks

xlim <- c(-2000,1900)

# ............................................................
# 4. DEFINE RCARBON MODELTEST FUNCTION & RETURN RAW VALUES ----
# ............................................................

rcarbonModelTestValues <- function(dates, 
                                   nullhypothesis){

  # 4.1 calibrate & binning ----
  cal = rcarbon::calibrate(x = dates$C14AGE, 			
                           errors = dates$C14STD,
                           calCurves = 'intcal13', 
                           ncores = ncores, 
                           normalised = FALSE)
  
  bins <- binPrep(sites = dates$SITE,
                  ages = dates$C14AGE,
                  h = n.bins # xx years cut off value
  ) 

  # 4.2 run modelTest ---
  
  print(paste(nullhypothesis, "nullhypothesis modelTest"))
  model <- modelTest(cal, 
                     errors = dates$C14STD, 
                     bins = bins, 
                     nsim = n.simulations, 
                     ncores = ncores, 
                     timeRange = timeRange, 
                     model = nullhypothesis, 
                     runm = kernel, 
                     raw = TRUE)
  
  # 4.6 extract raw results from modelTest ----
  
  # 4.6.1 spd & ci ----
  
  model.res <- data.frame(model$result)
  model.res$calBCAD <- 1950 - model.res$calBP
  model.res$model <- nullhypothesis
  
  # 4.6.2 booms&bust ----

  plot(get("model")) # NEVER TO BE REMOVED; PLOT NEEDS TO BE DRAWN!!!
  x = plot(get("model"), 
           bbty = 'n', 
           main = "model")
  
  lst <- list()
  cnt <- 1
  # 4.6.2.1 booms ----
  if(length(x$booms) != 0) {
    for(i in 1:length(x$booms)){
      FROM <- max(x$booms[[i]][[2]])
      TO <- min(x$booms[[i]][[2]])
      rcarbon <- "positive"
      bb <- data.frame(FROM = FROM,
                          TO = TO,
                          rcarbon = rcarbon)
      lst[[cnt]] <- bb
      cnt <- cnt + 1
    }
  }
  # 4.6.2.2 busts ----
  if(length(x$busts) != 0) {
    for(i in 1:length(x$busts)){
      FROM <- max(x$busts[[i]][[2]])
      TO <- min(x$busts[[i]][[2]])
      rcarbon <- "negative"
      bb <- data.frame(FROM = FROM,
                          TO = TO,
                          rcarbon = rcarbon)
      lst[[cnt]] <- bb
      cnt <- cnt + 1
    }
  }
  bb.res <- do.call(rbind, lst)
  if(length(bb.res) != 0){
    bb.res$MODEL <- nullhypothesis

    # 4.7 Convert to cal. BC/AD
    
    bb.res$FROM <- 1950 - bb.res$FROM
    bb.res$TO <- 1950 - bb.res$TO
    # cut excess data: 
    bb.res$FROM[bb.res$FROM <= xlim[1]] <- xlim[1]+1
    bb.res$TO[bb.res$TO >= xlim[2]] <- xlim[2]-1
  }

  # 4.4 Return list of results
  res <- list(
    ModelRes = model.res,
    BoomsBusts = bb.res
  )
  
  return(res)
}

# ............................................................
# 5. FIGURE 2 PANEL A - human activity record regions A-H (apply rcarbonTest function) ----
# ............................................................

  # ............................................................
  # 5.1 pottery record (based on ModelTest)----
  # ............................................................

c14.sel.A <- dplyr::filter(c14,
                           C14AGE > min.14C.age & # age filtering
                           C14STD > 0 & 
                           CLASS %in% c("Ia","Ib","Ic", "Id") & 
                           REGION %in% LETTERS[seq(from = 1, to = 8)] 	# regional filtering
                           )

kernel <- ceiling(mean(c14.sel.A$C14STD)/10)*10 # time window for running mean (in years)	based mean standard error		

# apply individual nullhypotheses
c14.sel.A.unf <- rcarbonModelTestValues(c14.sel.A, "uniform")
c14.sel.A.lin <- rcarbonModelTestValues(c14.sel.A, "linear")
c14.sel.A.exp <- rcarbonModelTestValues(c14.sel.A, "exponential")

# combine boom/bust phases
c14.sel.A.bb <- rbind(c14.sel.A.unf[["BoomsBusts"]], 
                      c14.sel.A.lin[["BoomsBusts"]], 
                      c14.sel.A.exp[["BoomsBusts"]])

  # ............................................................
  # 5.2 Refine high and low activity periods ----
  # ............................................................

c14.sel.A.bb_refined <- c14.sel.A.bb
c14.sel.A.bb_refined <- c14.sel.A.bb_refined [order(c14.sel.A.bb_refined$FROM),]
c14.sel.A.bb_refined$length <- abs(c14.sel.A.bb_refined$FROM - c14.sel.A.bb_refined$TO)
c14.sel.A.bb_refined <- c14.sel.A.bb_refined[c14.sel.A.bb_refined[,c("length")]>100,]
c14.sel.A.bb_refined$length <- NULL

  # ............................................................
  # 5.4 Define expansion and collapse periods ----
  # ............................................................

timeRange = c(4000, 0)

cal_pottery <- rcarbon::calibrate(x = c14.sel.A $C14AGE, 			
                         errors = c14.sel.A $C14STD,
                         calCurves = 'intcal13', 
                         ncores = ncores, 
                         normalised = FALSE)
spd_pottery <- spd(cal_pottery, 
          			timeRange = timeRange,
          			runm = kernel) #; windows() ; plot(spd_pottery)

change_method <- "breaks"
if(change_method == "breaks"){
	breaks = seq(4000, 0, -100) # 100 year blocks
	rate.of.change <- spd2rc(spd_pottery, breaks = breaks) #; windows() ; plot(rate.of.change)
	rc <- data.frame(calBP = rate.of.change$breaks[3:length(breaks)-1], rc = rate.of.change$roca)
	rc$rc.mean 	<- rc$rc
	rc$calBCAD 	<- 1950 - rc$calBP }
if(change_method == "backsight"){
	rate.of.change <- spd2rc(spd_pottery,backsight = 100)  #; windows() ; plot(rate.of.change)
	rc 		<- data.frame(calBP = rate.of.change$timeSequence, rc = rate.of.change$roca) 
	rc$rc.mean 	<- zoo::rollmean(rc$rc, 500, na.pad = T, align = "center")
	rc$calBCAD 	<- 1950 - rc$calBP }

rc_pos <- rc[rc[,c("rc.mean")]>0,]								# NEW method: two thresholds (cfr this line and the 4 lines below)
threshold_pos 	<- as.numeric(round(quantile(rc_pos$rc.mean,na.rm=T)[2],3))
rc_neg <- rc[rc[,c("rc.mean")]<0,]
threshold_neg 	<- as.numeric(round(quantile(rc_neg$rc.mean,na.rm=T)[2],3))
windows() ; plot(rc$calBCAD , rc$rc.mean) ; abline(h=0) ; polygon(x=c(-3000,3000,3000,-3000),y=c(threshold_pos,threshold_pos,threshold_neg,threshold_neg),col="grey50")

expansion1 <- rc[(rc[,c("calBCAD")] > (-1000)) & (rc[,c("calBCAD")] < 500),] 
expansion1 <- expansion1[expansion1[,c("rc.mean")] > threshold_pos ,] 				
expansion1 <- c(round(min(expansion1$calBCAD,na.rm=T)/100,0)*100 , round(max(expansion1$calBCAD,na.rm=T)/100,0)*100  ) ; expansion1

collapse <- rc[(rc[,c("calBCAD")] > 0) & (rc[,c("calBCAD")] < 1500),] 
collapse <- collapse[collapse[,c("rc.mean")] < threshold_neg,] 
collapse <- c(round(min(collapse$calBCAD,na.rm=T)/100,0)*100 , round(max(collapse$calBCAD,na.rm=T)/100,0)*100  ) ; collapse

expansion2 <- rc[ rc[,c("calBCAD")] > 500 & (rc[,c("calBCAD")] < 1900),] 
expansion2 <- expansion2[expansion2[,c("rc.mean")] > threshold_pos ,] 				
expansion2 <- c( round(min(expansion2$calBCAD,na.rm=T)/100,0)*100 , floor(max(expansion2$calBCAD,na.rm=T)/100)*100  ) ; expansion2

phases <- as.data.frame(c("low activity 1","expansion 1","high activity 1","collapse","low activity 2","expansion 2","high activity 2")) ; names(phases)<-"PHASE"
phases$FROM <- c(-1800,expansion1[1],expansion1[2],collapse[1],collapse[2],expansion2[1],expansion2[2])
phases$TO <- c(expansion1[1],expansion1[2],collapse[1],collapse[2],expansion2[1],expansion2[2],1900)

  # ............................................................
  # 5.5 FIG. 2.A ----
  # ............................................................

ymax <- as.data.frame(c14.sel.A.exp[["ModelRes"]]) ; ymax <- max(ymax$PrDens)			

plt.A <- ggplot() + 
  geom_ribbon(data = c14.sel.A.unf[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  geom_ribbon(data = c14.sel.A.lin[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  geom_ribbon(data = c14.sel.A.exp[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  #geom_line(data = spd_lithics, 
  #          aes(x = calBCAD, y = grid.PrDens), 
  #          linetype = 2) + 
  geom_line(data = c14.sel.A.exp[["ModelRes"]], # spd from one of the models
            aes(x = calBCAD, y = PrDens)) + 
  geom_rect(data = c14.sel.A.bb_refined,  
            aes(xmin = FROM, xmax = TO, ymin = 0, ymax = .65, fill = rcarbon), alpha = .1) + 
  geom_segment(aes(x = -3000, y = .1, xend = -2500, yend = .1, linetype = "dashed")) + 
  geom_segment(aes(x = -3000, y = .1, xend = -2500, yend = .1, linetype = "solid")) + 

  annotate("text", x=(-1950), y=ymax*1.53, label = paste("Congo Basin rainforest"),hjust=0,vjust=0.5,cex=4,col="black") +	

  annotate("text", x=mean(c(expansion1[1],collapse[1])) , y=ymax*1.53, label = paste("EARLY IRON AGE"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(expansion1[1],collapse[1])) , y=ymax*1.48, label = paste("AD",expansion1[1],"-",collapse[1]),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(expansion1[1],collapse[1]),y=c(ymax*1.45,ymax*1.45)),col="black",lwd=1)+

  annotate("text", x=mean(c(expansion2[1],1900)) , y=ymax*1.53, label = paste("LATE IRON AGE"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(expansion2[1],1900)) , y=ymax*1.48, label = paste("AD",expansion2[1],"- 1900"),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(expansion2[1],1900),y=c(ymax*1.45,ymax*1.45)),col="black",lwd=1)+

  annotate("text", x=mean(c(-1800,expansion1[1])) , y=ymax*1.03, label = paste("LOW ACTIVITY"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(-1800,expansion1[1])) , y=ymax*0.98, label = paste("-1800 BC - AD",expansion1[1]),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(-1800,expansion1[1]),y=c(ymax*0.95,ymax*0.95)),col="black",lwd=1)+

  annotate("text", x=mean(expansion1) , y=ymax*1.18, label = paste("EXPANSION"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(expansion1) , y=ymax*1.13, label = paste(expansion1[1],"BC - AD",expansion1[2]),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=expansion1,y=c(ymax*1.10,ymax*1.10)),col="black",lwd=1)+

  annotate("text", x=mean(c(expansion1[2],collapse[1])) , y=ymax*1.33, label = paste("HIGH ACTIVITY"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(expansion1[2],collapse[1])) , y=ymax*1.28, label = paste("AD",expansion1[2],"-",collapse[1]),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(expansion1[2],collapse[1]),y=c(ymax*1.25,ymax*1.25)),col="black",lwd=1)+

  annotate("text", x=mean(collapse) , y=ymax*1.18, label = paste("COLLAPSE"),hjust=0.5,vjust=0,cex=3,col="black") +	
  annotate("text", x=mean(collapse) , y=ymax*1.13, label = paste("AD",collapse[1],"-",collapse[2]),hjust=0.5,vjust=0,cex=3,col="black") +	
  geom_line(aes(x=collapse,y=c(ymax*1.10,ymax*1.10)),col="black",lwd=1)+

  annotate("text", x=mean(c(collapse[2],expansion2[1])) , y=ymax*1.03, label = paste("LOW ACTIVITY"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(collapse[2],expansion2[1])) , y=ymax*0.98, label = paste("AD",collapse[2],"-",expansion2[1]),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(collapse[2],expansion2[1]),y=c(ymax*0.95,ymax*0.95)),col="black",lwd=1)+

  annotate("text", x=mean(expansion2) , y=ymax*1.18, label = paste("EXPANSION"),hjust=0.5,vjust=0,cex=3,col="black") +	
  annotate("text", x=mean(expansion2) , y=ymax*1.13, label = paste("AD",expansion2[1],"-",expansion2[2]),hjust=0.5,vjust=0,cex=3,col="black") +		
  geom_line(aes(x=expansion2,y=c(ymax*1.1,ymax*1.1)),col="black",lwd=1)+

  annotate("text", x=mean(c(expansion2[2],1700)) , y=ymax*1.33, label = paste("HIGH ACTIVITY"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(expansion2[2],1700)) , y=ymax*1.28, label = paste("AD",expansion2[2],"-",1900),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(expansion2[2],1900),y=c(ymax*1.25,ymax*1.25)),col="black",lwd=1)+

  #scale_linetype_manual("Congo Basin rainforest",
  #                      values = c(dashed = "dashed", solid = "solid"), 
  #                      labels = c(
  #					  #paste("Pottery-related inventory (n =", nrow(c14.sel.A),"dates)"), 
  #                                #paste("Lithics-only inventory (n =", nrow(c14.sel.A.li),"dates)")), 
  #					   paste("Pottery-related"), 
  #                                 paste("Lithics-only")), 
  #                      breaks = c("solid", "dashed")) + 
  #scale_linetype_manual("Congo Basin rainforest",values=c("",""),labels=c("",""), breaks = c("solid", "dashed"))+
  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500)), 
                     limits = c(xlim[[1]], xlim[[2]]),
                     expand = c(0,0)) + 
  scale_y_continuous("Summed probability", 
                     limits = c(-.01,ymax*1.57),
			 expand = c(0, 0)) + 
  scale_fill_manual(values = c("#00bfc4", "#f8766d"), guide = FALSE) + 

  annotate(geom = "segment", x = seq(-2000, 2000, 100), xend = seq(-2000, 2000, 100),
           y = 0, yend = -.01) + 
  geom_hline(yintercept = 0) + 

  theme_classic() + 
  theme(#legend.position = c(.01, 1.02), 
		legend.position = "none",
        #legend.justification = c("left", "top"), 
        #legend.background = element_blank(), 
        #axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        #axis.ticks.x = element_blank(), 
        axis.title.x = element_blank())

windows(6.7,4) ; plt.A

# ............................................................
# 6. FIGURE 2 PANEL B - human activity record regions I-K (apply rcarbonTest function) ----
# ............................................................

  # ............................................................
  # 6.1 lithics record (not via rcarbonTest function) ----
  # ............................................................

c14.sel.B.li <- dplyr::filter(c14, 
                              CLASS == "IIa" & 
                                REGION %in% LETTERS[seq(from = 9, to = 11)])

cal = rcarbon::calibrate(x = c14.sel.B.li$C14AGE,
                         errors = c14.sel.B.li$C14STD,
                         calCurves = 'intcal13', 
                         ncores = ncores, 
                         normalised = FALSE)

spd_lithics_B <- spd(cal, 
           timeRange = timeRange, 
           runm = kernel)
spd_lithics_B <- as.data.frame(spd_lithics_B[2])
spd_lithics_B$calBCAD <- 1950 - spd_lithics_B$grid.calBP
spd_lithics_B <- filter(spd_lithics_B, grid.PrDens > 0)

  # ............................................................
  # 6.2 pottery record ----
  # ............................................................

c14.sel.B <- dplyr::filter(c14,
                           C14AGE > min.14C.age & # age filtering
                           C14STD > 0 & 
                           CLASS %in% c("Ia","Ib","Ic", "Id") & # class filtering
                           REGION %in% LETTERS[seq(from = 9, to = 11)] # regional filtering
)

# apply individual nullhypotheses
c14.sel.B.unf <- rcarbonModelTestValues(c14.sel.B, "uniform")
c14.sel.B.lin <- rcarbonModelTestValues(c14.sel.B, "linear")
c14.sel.B.exp <- rcarbonModelTestValues(c14.sel.B, "exponential")

# combine boom/bust phases
c14.sel.B.bb <- rbind(c14.sel.B.unf[["BoomsBusts"]], 
                      c14.sel.B.lin[["BoomsBusts"]], 
                      c14.sel.B.exp[["BoomsBusts"]])

  # ............................................................
  # 6.3 Refine high and low activity periods ----
  # ............................................................

c14.sel.B.bb_refined <- c14.sel.B.bb
c14.sel.B.bb_refined <- c14.sel.B.bb_refined [order(c14.sel.B.bb_refined$FROM),]
c14.sel.B.bb_refined$length <- abs(c14.sel.B.bb_refined$FROM - c14.sel.B.bb_refined$TO)
c14.sel.B.bb_refined <- c14.sel.B.bb_refined[c14.sel.B.bb_refined[,c("length")]>100,]
c14.sel.B.bb_refined$length <- NULL

  # ............................................................
  # 6.4 FIG. 2.B ----
  # ............................................................

ymax2 <- as.data.frame(c14.sel.B.exp[["ModelRes"]]) ; ymax2 <- max(ymax2$PrDens)			

plt.B <- ggplot() + 
  geom_ribbon(data = c14.sel.B.unf[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  geom_ribbon(data = c14.sel.B.lin[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  geom_ribbon(data = c14.sel.B.exp[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  #geom_line(data = spd_lithics_B,
  #          aes(x = calBCAD, y = grid.PrDens), 
  #          linetype = 2) + 
  geom_line(data = c14.sel.B.unf[["ModelRes"]], # spd from one of the models
            aes(x = calBCAD, y = PrDens)) + 
  geom_rect(data = c14.sel.B.bb_refined, 
            aes(xmin = FROM, xmax = TO, ymin = 0, ymax = Inf, fill = rcarbon), alpha = .1) + 

  annotate("text", x=(-1950), y=ymax2, label = paste("Congo Basin woodland and Bioko"),hjust=0,vjust=0.5,cex=4,col="black") +	

  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500)), 
                     limits = c(xlim[[1]], xlim[[2]]),
                     expand = c(0,0)) + 
  scale_y_continuous("", expand = c(0, 0)) + 
  coord_cartesian(xlim = c(xlim[[1]], xlim[[2]]), clip = "off") + 
  scale_fill_manual(values = c("#00bfc4", "#f8766d"), guide = FALSE) + 
  geom_segment(aes(x = -3000, y = .01, xend = -2500, yend = .01, linetype = "dashed")) + 
  geom_segment(aes(x = -3000, y = .01, xend = -2500, yend = .01, linetype = "solid")) + 
  #scale_linetype_manual("Periphery of the rainforest",
  #                      values = c(dashed = "dashed", solid = "solid"), 
  #                      labels = c(paste("Pottery inventory n =", nrow(c14.sel.B)), 
  #                                 paste("Lithics-only n =", nrow(c14.sel.B.li))), 
  #                      breaks = c("solid", "dashed")) + 
  annotate(geom = "segment", x = seq(-2000, 2000, 100), xend = seq(-2000, 2000, 100),
           y = 0, yend = -.0075) + 
  geom_hline(yintercept = 0) + 
  theme_classic() + 
  theme(#legend.position = c(.001, 1), 
		legend.position = "none",
        #legend.justification = c("left", "top"), 
        #legend.background = element_blank(), 
        axis.line.x = element_blank())
 
windows(6.7,2) ; plt.B

# ............................................................
# 7. PLOT ----
# ............................................................

# dummy figure for legend ----
plt.legend <- ggplot() + 
  geom_rect(data = c14.sel.A.bb, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .3) + 
  scale_fill_manual("",values = c("#00bfc4", "#f8766d"), 
                    labels = c("more intense human activity", 
                               "less intense human activity")  ) +
  theme(legend.position = "top")
legend <- cowplot::get_legend(plt.legend)

plt.graph <- cowplot::plot_grid(plt.A, 
                                plt.B, 
                                ncol = 1, 
                                labels = "auto",  align = "v", axis = "lr", 
                                rel_heights = c(2, 0.5))
plt <- cowplot::plot_grid(legend, 
                          plt.graph, 
                          ncol = 1,
                          rel_heights = c(1, 10))

windows() ; plt

# ............................................................
# 8. WRITE ----
# ............................................................

ggsave("output/Figure 2 overall spd - main analysis.pdf", 
       plt, 
       width = 8, height = 5.5)
ggsave("output/Figure 2 overall spd - main analysis.jpg", 
       plt,
       width = 8, height = 5.5)

write.csv(c14.sel.A.bb_refined, 
          "output/rcarbonTest_timewindows_a-h.csv", 
          fileEncoding = "UTF-8",
          row.names = F)

write.csv(c14.sel.B.bb_refined, 
          "output/rcarbonTest_timewindows_i-k.csv", 
          fileEncoding = "UTF-8",
          row.names = F)

write.csv(phases, 
          "output/rcarbon_phases.csv", 
          fileEncoding = "UTF-8",
          row.names = F)