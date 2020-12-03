#### FIGURE 2 OVERALL SPD ####

# setwd("C:/Wannes/UGent Data/4 POSTDOC KMMA/R/radiocarbon - human activity/CentralAfricaHumanActivity v8_DS/")

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

min.14C.age <- 0 # rcarbonTest fails on 0 bp dates
n.bins <- 100 # time window for binning dates (in years) within the same site	
n.simulations <- 100 # nr of simulations
timeRange = c(4000, 0) # set the timerange of analysis in calBP, older date first
breaks = seq(4000, 0, -200) # 200 year blocks

xlim <- c(-2000,1900)

nsim <- 500 # number of simulations for KDE
kernelbw <- 100 # kernel bandwidth for KDE

cal.norm <- FALSE
calCurve <- 'intcal20'

SensitivityAnalysis <- FALSE

# ............................................................
# 3. DEFINE RCARBON MODELTEST FUNCTION & RETURN RAW VALUES ----
# ............................................................

rcarbonModelTestValues <- function(dates, 
                                   nullhypothesis){

  # 3.1 calibrate & binning ----
  cal = rcarbon::calibrate(x = dates$C14AGE, 			
                           errors = dates$C14STD,
                           calCurves = calCurve, 
                           ncores = ncores, 
                           normalised = cal.norm)
  
  bins <- binPrep(sites = dates$SITE,
                  ages = dates$C14AGE,
                  h = n.bins # xx years cut off value
  ) 

  print(paste("based on", nrow(dates), "dates and", length(unique(bins)), "bins"))
  
  # 3.2 run modelTest ----
  print(paste(nullhypothesis, "nullhypothesis modelTest"))
  
  if(nullhypothesis == "logistic"){ # 3.2.1 logistic model ----
    
    # logistic growth model
    # see https://cran.r-project.org/web/packages/rcarbon/vignettes/rcarbon.html

      spd = rcarbon::spd(
      cal,
      timeRange = timeRange,
      bins = bins,
      runm = kernel
    )

    # define xmid for nls function/inflection point of the curve as max of SPD
    logFit.mid <- plyr::round_any(
      spd$grid[which.max(spd$grid$PrDens),"calBP"],
      500, 
      f = ceiling
      )
    
    logFit <- stats::nls(
      PrDens ~ SSlogis(calBP, 
                       Asym, 
                       xmid, 
                       scale
      ),
      data = spd$grid,
      control = stats::nls.control(
        maxiter = 200
      ),
      start = list(Asym = .3, # TODO: default .2 causes error with intcal20 !!!
                   xmid = logFit.mid, # start value
                   scale = -100
      )
    )

    logFitDens = data.frame(
      calBP = spd$grid$calBP,
      PrDens = stats::SSlogis(
        input = spd$grid$calBP,
        Asym = coefficients(logFit)[1],
        xmid = coefficients(logFit)[2],
        scal = coefficients(logFit)[3]
        )
      )
    
    model <- rcarbon::modelTest(
      cal, 
      errors = dates$C14STD, 
      bins = bins,
      nsim = n.simulations,
      timeRange = timeRange, 
      model = "custom",
      predgrid = logFitDens, 
      runm = kernel,
      raw = TRUE
      )

  } else { # 3.2.2 build in models ----
    
    model <- modelTest(
      cal, 
      errors = dates$C14STD, 
      bins = bins, 
      nsim = n.simulations, 
      ncores = ncores, 
      timeRange = timeRange, 
      model = nullhypothesis, 
      runm = kernel, 
      raw = TRUE
      )
    
  }
  
  # 3.3 extract raw results from modelTest ----
  
  # 3.3.1 spd & ci ----
  model.res <- data.frame(model$result)
  model.res$calBCAD <- 1950 - model.res$calBP
  model.res$model <- nullhypothesis
  
  # 3.3.2 booms&bust ----
  plot(get("model")) # NEVER TO BE REMOVED; PLOT NEEDS TO BE DRAWN!!!
  x = plot(get("model"), 
           bbty = 'n', 
           main = "model")
  
  lst <- list()
  cnt <- 1
  # 3.3.2.1 booms ----
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
  # 3.3.2.2 busts ----
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

    # 3.4 Convert to cal. BC/AD
    
    bb.res$FROM <- 1950 - bb.res$FROM
    bb.res$TO <- 1950 - bb.res$TO
    # cut excess data: 
    bb.res$FROM[bb.res$FROM <= xlim[1]] <- xlim[1]+1
    bb.res$TO[bb.res$TO >= xlim[2]] <- xlim[2]-1
  }

  # 3.5 Remove boom/bust phases earlier than modeled phases ----
  model.start <- head(model.res[model.res$lo > 0,], 1)$calBCAD
  bb.res <- bb.res[bb.res$FROM > model.start, ]
  
  # 3.6 Return list of results
  res <- list(
    ModelRes = model.res,
    BoomsBusts = bb.res
  )
  
  return(res)
}

# ............................................................
# 4 DEFINE RCARBON CKDE FUNCTION & RETURN RAW VALUES ----
# ............................................................

rcarbonKDE <- function(dates){

  # 4.1 optimal kernel (cf. McLaughlin 2018) ----
  kernel <- density(dates$C14STD, 
                    kernel = "gaussian") 
  
  kernel.df <- data.frame(x = kernel$x, 
                          y = kernel$y)
  
  kernel.val <- round(kernel.df[which.max(kernel.df$y),"x"])
  
  print(paste("kernel =", kernel.val))
  
  # 4.2 calibrate & binning ----
  cal <- rcarbon::calibrate(
    x = dates$C14AGE, 			
    errors = dates$C14STD,
    calCurves = calCurve, 
    ncores = ncores, 
    normalised = cal.norm)
  
  bins <- rcarbon::binPrep(
    sites = dates$SITE,
    ages = dates$C14AGE,
    h = n.bins)

  # 4.3 sample dates 6 KDE ----
  rand <- rcarbon::sampleDates(
    cal, 
    bins = bins, 
    nsim = nsim)
  
  print(paste(length(rand$weight), "randomly sampled calibrated 14C dates"))
  
  ckde = rcarbon::ckde(
    rand,
    timeRange = timeRange,
    bw = kernel.val,
    normalised = cal.norm)

  # 4.4 extract data ----
  kde <- as.data.frame(ckde$res.matrix)
  kde$calBP <- as.numeric(
    seq(timeRange[1], timeRange[2], -1)
  )
  kde <- tibble::column_to_rownames(kde, var = "calBP")
  
  # 4.5 correct KDE probability to match SPD (expected value) ----
  corr.fct <- length(unique(bins))
  
  print(corr.fct)
  
  kde <- tibble::rownames_to_column(kde, var = "calBP")
  kde$calBCAD <- 1950 - as.numeric(kde$calBP)
  kde.melt <- reshape2::melt(kde, id.vars = c("calBCAD", "calBP"))
  kde.melt$value <- kde.melt$value * corr.fct # apply correction factor
  
  return(kde.melt)
}

# ............................................................
# 5. FIGURE 2 PANEL A - human activity record regions A-H (apply rcarbonTest function) ----
# ............................................................

  # ............................................................
  # 5.1 pottery record (based on ModelTest)----
  # ............................................................

c14.sel.A <- dplyr::filter(
  c14,
  C14AGE > min.14C.age & # age filtering
  C14STD > 0 & 
  REGION %in% LETTERS[seq(from = 1, to = 8)] 	# regional filtering
)

# filter CLASS within main analysis
if(SensitivityAnalysis != TRUE){
  c14.sel.A <- dplyr::filter(
    c14.sel.A,
    CLASS %in% c("Ia","Ib","Ic", "Id"))
}
nrow(c14.sel.A)

kernel <- ceiling(mean(c14.sel.A$C14STD)/10)*10 # time window for running mean (in years)	based mean standard error		

# apply individual nullhypotheses
c14.sel.A.unf <- rcarbonModelTestValues(c14.sel.A, "uniform")
c14.sel.A.lin <- rcarbonModelTestValues(c14.sel.A, "linear")
c14.sel.A.exp <- rcarbonModelTestValues(c14.sel.A, "exponential")
c14.sel.A.log <- rcarbonModelTestValues(c14.sel.A, "logistic")

# combine boom/bust phases
c14.sel.A.bb <- rbind(c14.sel.A.unf[["BoomsBusts"]], 
                      c14.sel.A.lin[["BoomsBusts"]], 
                      c14.sel.A.exp[["BoomsBusts"]],
                      c14.sel.A.log[["BoomsBusts"]])

  # ............................................................
  # 5.2 Refine high and low activity periods ----
  # ............................................................

c14.sel.A.bb_refined <- c14.sel.A.bb
c14.sel.A.bb_refined <- c14.sel.A.bb_refined [order(c14.sel.A.bb_refined$FROM),]

  # ............................................................
  # 5.3 Define expansion and collapse periods ----
  # ............................................................

cal_pottery <- rcarbon::calibrate(x = c14.sel.A$C14AGE, 			
                         errors = c14.sel.A$C14STD,
                         calCurves = calCurve, 
                         ncores = ncores, 
                         normalised = cal.norm)

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

rc_pos <- rc[rc[,c("rc.mean")]>0,]
threshold_pos <- as.numeric(round(quantile(rc_pos$rc.mean,na.rm=T)[3],5))
rc_neg <- rc[rc[,c("rc.mean")]<0,]
threshold_neg <- as.numeric(round(quantile(rc_neg$rc.mean,na.rm=T)[3],5))
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
  # 5.4 Run KDE ----
  # ............................................................

c14.sel.A.kde <- rcarbonKDE(c14.sel.A)

  # ............................................................
  # 5.5 FIG. 2.A ----
  # ............................................................

ymax <- as.data.frame(c14.sel.A.exp[["ModelRes"]]) ; ymax <- max(ymax$PrDens)			

if(SensitivityAnalysis == TRUE){
  bb.rect.ymax <- .8
} else {
  bb.rect.ymax <- .65
}

plt.A <- ggplot() + 
  geom_ribbon(data = c14.sel.A.log[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .3) + 
  geom_line(data = c14.sel.A.exp[["ModelRes"]], # spd from one of the models
            aes(x = calBCAD, y = PrDens)) + 
  geom_rect(data = c14.sel.A.bb_refined,  
            aes(xmin = FROM, xmax = TO, ymin = 0, ymax = bb.rect.ymax, fill = rcarbon), alpha = .1) + 
  geom_line(data = c14.sel.A.kde, aes(x = calBCAD, y = value, group = variable), alpha = .01) + 
  geom_segment(aes(x = -3000, y = .1, xend = -2500, yend = .1, linetype = "dashed")) + 
  geom_segment(aes(x = -3000, y = .1, xend = -2500, yend = .1, linetype = "solid")) + 

  annotate("text", x=(-1950), y=ymax*1.53, label = paste("Congo Basin rainforest"),hjust=0,vjust=0.5,cex=4,col="black") +	

  annotate("text", x=mean(c(expansion1[1],collapse[1])) , y=ymax*1.53, label = paste("EARLY IRON AGE"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(expansion1[1],collapse[1])) , y=ymax*1.48, label = paste(abs(expansion1[1]),"BC - AD",collapse[1]),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(expansion1[1],collapse[1]),y=c(ymax*1.45,ymax*1.45)),col="black",lwd=1)+

  annotate("text", x=mean(c(expansion2[1],1900)) , y=ymax*1.53, label = paste("LATE IRON AGE"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(expansion2[1],1900)) , y=ymax*1.48, label = paste("AD",expansion2[1],"- 1900"),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
  geom_line(aes(x=c(expansion2[1],1900),y=c(ymax*1.45,ymax*1.45)),col="black",lwd=1)+

  annotate("text", x=mean(c(-1800,expansion1[1])) , y=ymax*1.03, label = paste("LOW ACTIVITY"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(c(-1800,expansion1[1])) , y=ymax*0.98, label = paste("1800 -",abs(expansion1[1]), "BC"),hjust=0.5,vjust=0,cex=3,col="black") + 
  geom_line(aes(x=c(-1800,expansion1[1]),y=c(ymax*0.95,ymax*0.95)),col="black",lwd=1)+

  annotate("text", x=mean(expansion1) , y=ymax*1.18, label = paste("EXPANSION"),hjust=0.5,vjust=0,cex=3,col="black") +			
  annotate("text", x=mean(expansion1) , y=ymax*1.13, label = paste(abs(expansion1[1]),"BC - AD",expansion1[2]),hjust=0.5,vjust=0,cex=3,col="black") +			# added by WH
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

  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500)), 
                     limits = c(xlim[[1]], xlim[[2]]),
                     expand = c(0,0)) + 
  scale_y_continuous("Summed probability", 
                     limits = c(-.01,ymax*1.57),
			 expand = c(0, 0)) + 
  scale_fill_manual(values = c("#f8766d", "#00bfc4"), guide = FALSE) + 

  annotate(geom = "segment", x = seq(-2000, 2000, 100), xend = seq(-2000, 2000, 100),
           y = 0, yend = -.01) + 
  geom_hline(yintercept = 0) + 

  theme_classic() + 
  theme(
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.title.x = element_blank())

windows(6.7,4) ; plt.A

# ............................................................
# 6. FIGURE 2 PANEL B - human activity record regions I-K (apply rcarbonTest function) ----
# ............................................................

  # ............................................................
  # 6.1 pottery record ----
  # ............................................................

c14.sel.B <- dplyr::filter(
  c14,
  C14AGE > min.14C.age & # age filtering
  C14STD > 0 & 
  REGION %in% LETTERS[seq(from = 9, to = 11)] # regional filtering
)

# filter CLASS within main analysis
if(SensitivityAnalysis != TRUE){
  c14.sel.B <- dplyr::filter(
    c14.sel.B,
    CLASS %in% c("Ia","Ib","Ic", "Id"))
}
nrow(c14.sel.B)

# apply individual nullhypotheses
c14.sel.B.unf <- rcarbonModelTestValues(c14.sel.B, "uniform")
c14.sel.B.lin <- rcarbonModelTestValues(c14.sel.B, "linear")
c14.sel.B.exp <- rcarbonModelTestValues(c14.sel.B, "exponential")
c14.sel.B.log <- rcarbonModelTestValues(c14.sel.B, "logistic")

# combine boom/bust phases
c14.sel.B.bb <- rbind(c14.sel.B.unf[["BoomsBusts"]], 
                      c14.sel.B.lin[["BoomsBusts"]], 
                      c14.sel.B.exp[["BoomsBusts"]],
                      c14.sel.B.log[["BoomsBusts"]])

  # ............................................................
  # 6.2 Refine high and low activity periods ----
  # ............................................................

c14.sel.B.bb_refined <- c14.sel.B.bb
c14.sel.B.bb_refined <- c14.sel.B.bb_refined [order(c14.sel.B.bb_refined$FROM),]

  # ............................................................
  # 6.3 FIG. 2.B ----
  # ............................................................

ymax2 <- as.data.frame(c14.sel.B.exp[["ModelRes"]]) ; ymax2 <- max(ymax2$PrDens)			

plt.B <- ggplot() + 
  geom_ribbon(data = c14.sel.B.log[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .3) + 
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
  scale_fill_manual(values = c("#f8766d", "#00bfc4"), guide = FALSE) + 
  geom_segment(aes(x = -3000, y = .01, xend = -2500, yend = .01, linetype = "dashed")) + 
  geom_segment(aes(x = -3000, y = .01, xend = -2500, yend = .01, linetype = "solid")) + 
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

# 7.1 dummy figure for legend ----
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

if(SensitivityAnalysis == TRUE){
  mainfig.path <- "output/Supplementary Figure 6 overall spd - sensitivity analysis"
  supplfig.path <- "output/Supplementary Figure 6 individual models - sensitivity analysis"
} else {
  mainfig.path <- "output/Figure 2 overall spd - main analysis"
  supplfig.path <- "output/Supplementary Figure 3 individual models - main analysis"
}

ggsave(paste0(mainfig.path, ".pdf"), 
       plt, 
       width = 8, height = 5.5)
ggsave(paste0(mainfig.path, ".jpg"), 
       plt,
       width = 8, height = 5.5)

if(SensitivityAnalysis != TRUE){ # only export in main analysis
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
}

# ............................................................
# 9. SUPPLEMENTAL FIG. 5 ----
# ............................................................

# 9.1 comine results of individual models ----
c14.sel.A.unf[["ModelRes"]]$REGION <- "A-H"
c14.sel.A.unf[["ModelRes"]]$GRP <- paste(c14.sel.A.unf[["ModelRes"]]$REGION, c14.sel.A.unf[["ModelRes"]]$model)

c14.sel.A.lin[["ModelRes"]]$REGION <- "A-H"
c14.sel.A.lin[["ModelRes"]]$GRP <- paste(c14.sel.A.lin[["ModelRes"]]$REGION, c14.sel.A.lin[["ModelRes"]]$model)

c14.sel.A.exp[["ModelRes"]]$REGION <- "A-H"
c14.sel.A.exp[["ModelRes"]]$GRP <- paste(c14.sel.A.exp[["ModelRes"]]$REGION, c14.sel.A.exp[["ModelRes"]]$model)

c14.sel.A.log[["ModelRes"]]$REGION <- "A-H"
c14.sel.A.log[["ModelRes"]]$GRP <- paste(c14.sel.A.log[["ModelRes"]]$REGION, c14.sel.A.log[["ModelRes"]]$model)

c14.sel.B.unf[["ModelRes"]]$REGION <- "I-K"
c14.sel.B.unf[["ModelRes"]]$GRP <- paste(c14.sel.B.unf[["ModelRes"]]$REGION, c14.sel.B.unf[["ModelRes"]]$model)

c14.sel.B.lin[["ModelRes"]]$REGION <- "I-K"
c14.sel.B.lin[["ModelRes"]]$GRP <- paste(c14.sel.B.lin[["ModelRes"]]$REGION, c14.sel.B.lin[["ModelRes"]]$model)


c14.sel.B.exp[["ModelRes"]]$REGION <- "I-K"
c14.sel.B.exp[["ModelRes"]]$GRP <- paste(c14.sel.B.exp[["ModelRes"]]$REGION, c14.sel.B.exp[["ModelRes"]]$model)

c14.sel.B.log[["ModelRes"]]$REGION <- "I-K"
c14.sel.B.log[["ModelRes"]]$GRP <- paste(c14.sel.B.log[["ModelRes"]]$REGION, c14.sel.B.log[["ModelRes"]]$model)

models.all <- rbind(c14.sel.A.unf[["ModelRes"]],
                    c14.sel.A.lin[["ModelRes"]], 
                    c14.sel.A.exp[["ModelRes"]],
                    c14.sel.A.log[["ModelRes"]],
                    c14.sel.B.unf[["ModelRes"]],
                    c14.sel.B.lin[["ModelRes"]],
                    c14.sel.B.exp[["ModelRes"]],
                    c14.sel.B.log[["ModelRes"]])

# 9.2 combinbe boom/bust phases ----
c14.sel.A.bb_refined$REGION <- "A-H"
c14.sel.B.bb_refined$REGION <- "I-K"

c14.sel.bb_refined <- rbind(c14.sel.A.bb_refined, 
                            c14.sel.B.bb_refined)

names(c14.sel.bb_refined)[names(c14.sel.bb_refined) == 'MODEL'] <- 'model'

# 9.3. giver order of models ----
models.all$model <- ordered(models.all$model, levels = c("uniform", "linear", "exponential", "logistic"))
c14.sel.bb_refined$model <- ordered(c14.sel.bb_refined$model, levels = c("uniform", "linear", "exponential", "logistic"))

model.labs <- models.all %>%
  dplyr::group_by(REGION, model) %>%
  dplyr::summarize(max = max(PrDens))

# 9.4 generate plot ----
plt.suppl <- ggplot() + 
  geom_rect(data = c14.sel.bb_refined,  
            aes(xmin = FROM, xmax = TO, 
                ymin = 0, ymax = Inf, 
                fill = rcarbon), 
            alpha = .2) + 
  scale_fill_discrete("", 
                      labels = c("less intense human activity", 
                                 "more intense human activity")) + 
  geom_ribbon(data = models.all, 
              aes(x = calBCAD, y = PrDens, 
                  ymin = lo, ymax = hi, 
                  group = GRP), 
              alpha = .3) + 
  geom_line(data = models.all, 
            aes(x = calBCAD, y = PrDens, 
                group = GRP)) + 
  geom_text(data = model.labs, 
            aes(x = -1950, y = .75*max, 
                label = paste(model, "growth")), 
            hjust = 0, size = 3, fontface = "bold") + 
  facet_grid(REGION + model ~ ., 
             scales = "free", 
             space = "free") + 
  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500), 1950),
                     limits = c(xlim[[1]], xlim[[2]]), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Summed probability", 
                     breaks = seq(0, max(models.all$PrDens), .1)) + 
  theme_classic() + 
  theme(legend.position = "top", 
        strip.text.y = element_text(size = 0, color = "white"),
        #strip.text.y.left = element_blank(),# element_text(angle = 0, face = "bold"), 
        strip.background = element_blank(), 
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE))

# ............................................................
# 9.5. Individual Plot Elements ----
# ............................................................

axis <- ggplot() +
  annotate(geom = "segment", 
           x = c(seq(-2000, 1900, 100), 1950), 
           xend = c(seq(-2000, 1900, 100), 1950),
           y = 0, yend = -1) + 
  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1500, 500), 1950),
                     limits = c(xlim[[1]], xlim[[2]]), 
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
# 9.6. Combine Elements & Save Figure ----
# ............................................................

plt.label <- cowplot::plot_grid(NULL, NULL, NULL, 
                                nrow = 3, 
                                rel_heights = c(1, 15, 3.25),
                                labels = c("", "a", "b"))

plt.suppl.1 <- cowplot::plot_grid(plt.suppl, axis,
                          ncol = 1, 
                          rel_heights = c(22.5,1),
                          align = "v", axis = "lr")

plt.suppl.2 <- cowplot::plot_grid(plt.label, plt.suppl.1, 
                          ncol = 2, 
                          rel_widths = c(1,30))

ggsave(paste0(supplfig.path, ".pdf"),
       plt.suppl.2,
       width = 8, height = 10)
ggsave(paste0(supplfig.path, ".jpg"),
       plt.suppl.2,
       width = 8, height = 10)
