
#### FIGURE 3 POTTERY STYLES ####

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

  # ............................................................
  # 1.1 geo data ----
  # ............................................................

rivers10 <- ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")
lakes10 <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass="sf")
coast10 <- ne_download(scale = 10, type = "coastline", category = "physical", returnclass="sf")
land10 <- ne_download(scale = 10, type = "land", category = "physical", returnclass="sf")
boundary_lines_land10 <- ne_download(scale = 10, type = "boundary_lines_land", category = "cultural", returnclass="sf")

  # ............................................................
  # 1.2 domain specific data ----
  # ............................................................

sites		<- data.table::fread("input/sites.csv",encoding = "UTF-8")
pottery 	<- data.table::fread("input/potterygroups.csv", encoding = "UTF-8") 

rcarbon 	<- read.csv("output/rcarbonTest_timewindows_a-h.csv", encoding = "UTF-8")			
phases 	<- read.csv("output/rcarbon_phases.csv", encoding = "UTF-8")			

names(sites) <-  c("REGION","SITE","LAT","LONG","STYLE","SOURCE" )
names(pottery) <-  c("ID","STYLE","FROM","TO","DESCRIPTION","REGION")

# ............................................................
# 2. CHOICES ----
# ............................................................

merge_cameroon 	<- "FALSE"
concave 		<- "TRUE"
ttest_on_bins 	<- "TRUE"
binyears 		<- 100 
x_limits 		<- c(-1000,1800)	

regions <- c("A","B","C","D","E","F","G","H")

# ............................................................
# 2. PREPARE DATASETS ----
# ............................................................

nrow(sites)
sites <- filter(sites,REGION %in% regions)				
nrow(sites)

nrow(pottery)
pottery <- filter(pottery,REGION %in% regions)
nrow(pottery)

if(merge_cameroon == "TRUE"){
	tmp <- as.data.frame(pottery) ; tmp <- tmp[tmp[,c("STYLE")] %in% c("Obobogo","Malongo/Bissiang","Bwamb?/Mpoengu","Dibamba E","Okala/Epona/Yindo"),]
	pottery[which(pottery[,c("STYLE")] == "Okala/Epona/Yindo"),c("FROM")] <- mean(tmp$FROM)
	pottery[which(pottery[,c("STYLE")] == "Okala/Epona/Yindo"),c("TO")] <- mean(tmp$TO)
	pottery <- pottery[-(which(pottery[,c("STYLE")] == "Obobogo")),]
	pottery <- pottery[-(which(pottery[,c("STYLE")] == "Malongo/Bissiang")),]
	pottery <- pottery[-(which(pottery[,c("STYLE")] == "Bwamb?/Mpoengu")),]
	pottery <- pottery[-(which(pottery[,c("STYLE")] == "Dibamba E")),]
	sites$STYLE <- replace(sites[,c("STYLE")],sites[,c("STYLE")]=="Obobogo","Okala/Epona/Yindo")
	sites$STYLE <- replace(sites[,c("STYLE")],sites[,c("STYLE")]=="Malongo/Bissiang","Okala/Epona/Yindo")
	sites$STYLE <- replace(sites[,c("STYLE")],sites[,c("STYLE")]=="Bwamb?/Mpoengu","Okala/Epona/Yindo")
	sites$STYLE <- replace(sites[,c("STYLE")],sites[,c("STYLE")]=="Dibamba E","Okala/Epona/Yindo")   }

sites <- st_as_sf(sites, 
                  coords = c("LONG", "LAT"), 
                  remove = F, 
                  crs = 4326, 
                  na.fail = F)

add <- pottery ; add $DESCRIPTION<-NULL ; add $REGION<-NULL
sites <- merge(x = sites, 
		           y = add,
               by = "STYLE", 
               all.x = T)

sites <- dplyr::filter(sites, 
                       REGION %in% LETTERS[seq(from = 1, to = 11)]) # delete those outside study region

rcarbon$FROM 	<- replace(rcarbon[,c("FROM")],rcarbon[,c("FROM")]<(-999),-999)
rcarbon$TO 	<- replace(rcarbon[,c("TO")],rcarbon[,c("TO")]>x_limits[2],x_limits[2])

# ............................................................
# 3. BINS ----
# ............................................................

breaks <- seq(-1000, 2000, binyears)
class <- seq(1,length(breaks), 1)		#class <- LETTERS[seq( from = 1, to = length(breaks) )]
breaks <- data.frame(breaks, class)
for(i in 1:nrow(breaks)){breaks[i, "labels"] <- paste0(breaks[i,"class"], ": ", breaks[i,"breaks"], "/", breaks[i+1,"breaks"])}

# ............................................................
# 4. ANALYSIS ----
# ............................................................

  # ............................................................
  # 4.1. Frequency of sites per pottery group ----
  # ............................................................

pottery.sites.freq <- as.data.frame(stats::aggregate(SITE ~ STYLE,data = sites, FUN = length)) 	

  # ............................................................
  # 4.2. Area per pottery group (Concave/Convex hull method) ----
  # ............................................................

id <- dplyr::filter(pottery.sites.freq, SITE > 2)

if(concave == "TRUE"){ 								# concave hull ----  # see https://github.com/joelgombin/concaveman
  pottery.sites.area <- sf::st_multipolygon()
  pottery.sites.area <- st_sf(polygons = st_sfc(st_polygon()))
  sf::st_crs(pottery.sites.area) <- (4326)
  pottery.sites.area$Type <- NA
  for(i in 1:nrow(id)){
    sites.f <- dplyr::filter(sites, STYLE == id[i,1])
    conc.hull <- concaveman(sites.f)
    conc.hull$Type = id[i, "STYLE"]
    pottery.sites.area <- rbind(pottery.sites.area, conc.hull)
  }

}else{ 										# convex hull ----
  lst <- lapply(1:nrow(id), function(x){
    sites.f <- dplyr::filter(sites, STYLE == id[x,1])
    s <- sites.f
    sf::st_geometry(sites.f) <- NULL
    ch <- grDevices::chull(sites.f[,c("LONG", "LAT")]) # return pts on hull in clockwise order
    coords <- sites.f[c(ch, ch[1]), ] # closed polygon
    sf::st_polygon(list(as.matrix(coords[,c("LONG", "LAT")])))
  })
  pottery.sites.area <- sf::st_sf(Type = id[, "STYLE"], st_sfc(lst))
  sf::st_crs(pottery.sites.area) <- (4326)
}

pottery.sites.area$AREA <- sf::st_area(pottery.sites.area)
pottery.sites.area$AREA <- as.numeric(pottery.sites.area$AREA)/1E9 # convert m2 into k(ilo) km2
names(pottery.sites.area)[names(pottery.sites.area) == 'Type'] <- 'STYLE'

  # ............................................................
  # 4.3. Frequency of pottery groups per bin ----
  # ............................................................

pottery.cent <- data.frame(matrix(ncol = ncol(pottery)+1, nrow = 0))
x <- c(names(pottery), "CLASS")
colnames(pottery.cent) <- x

for (i in 1:length(pottery$STYLE)){
  for (j in 1:(nrow(breaks)-1)) {
    if(pottery[i,"TO"] > breaks[j,"breaks"] & 
       pottery[i,"FROM"] < breaks[j+1,"breaks"]){
      l <- pottery[i,]
      l$CLASS <- breaks[j,"labels"]
      pottery.cent <- rbind(pottery.cent, as.data.frame(l))
    }

  }
}
pottery.cent$AGE <- (as.numeric(sub("/.*", "", sub(".*? ", "", pottery.cent$CLASS))) + as.numeric(sub(".*/", "", sub(".*? ", "", pottery.cent$CLASS)))) / 2
pottery.cent$AGE.jitter 	<- jitter(pottery.cent$AGE, 2)

  # ............................................................
  # 4.4. merge into meta tables  ----
  # ............................................................

pottery.cent.meta <- pottery.cent %>%
  dplyr::select(-DESCRIPTION) %>%
  dplyr::left_join(pottery.sites.freq, by = "STYLE") %>%
  dplyr::left_join(pottery.sites.area, by = "STYLE")

sites.cent <- merge(x = sites, 							# merge sites per style with class (200-year century list)
                    y = dplyr::select(pottery.cent, -DESCRIPTION), 
                    by = "STYLE", 
                    allow.cartesian = TRUE)

  # ............................................................
  # 4.5. extremes  ----
  # ............................................................

#head(pottery.cent.meta)

#hist(pottery.cent.meta$FROM,breaks= (max(pottery.cent.meta$FROM)-min(pottery.cent.meta$FROM))/100   )
#hist(pottery.cent.meta$TO,breaks= (max(pottery.cent.meta$TO)-min(pottery.cent.meta$TO))/100   )

pottery.cent.freq		<- pottery.cent.freq[pottery.cent.freq[,c("Var1")] < 1800  ,] 	# filter out last class (historical time)
nrow(pottery.cent.meta)
pottery.cent.meta		<- pottery.cent.meta[pottery.cent.meta[,c("AGE")] < 1800  ,] 	# filter out last class (historical time)
nrow(pottery.cent.meta)

  # ............................................................
  # 4.6. t-test  ----
  # ............................................................

sites_per_stylegroup 		<- pottery.cent.meta 									; nrow(sites_per_stylegroup)

tmp1 <- sites_per_stylegroup[ sites_per_stylegroup[,c("AGE")]> phases[2,c("FROM")] & sites_per_stylegroup[,c("AGE")]< phases[2,c("TO")] ,] 	; nrow(tmp1)
tmp2 <- sites_per_stylegroup[ sites_per_stylegroup[,c("AGE")]> phases[3,c("FROM")] & sites_per_stylegroup[,c("AGE")]< phases[3,c("TO")] ,]	; nrow(tmp2)
tmp3 <- sites_per_stylegroup[ sites_per_stylegroup[,c("AGE")]> phases[5,c("FROM")] & sites_per_stylegroup[,c("AGE")]< phases[5,c("TO")] ,]	; nrow(tmp3)
tmp4 <- sites_per_stylegroup[ sites_per_stylegroup[,c("AGE")]> phases[6,c("FROM")] & sites_per_stylegroup[,c("AGE")]< phases[6,c("TO")] ,] 	; nrow(tmp4)	
tmp5 <- sites_per_stylegroup[ sites_per_stylegroup[,c("AGE")]> phases[7,c("FROM")] & sites_per_stylegroup[,c("AGE")]< phases[7,c("TO")] ,] 	; nrow(tmp5)	

pvalue1 <- t.test(tmp1$SITE ,tmp2$SITE,alternative="greater")$p.value ; pvalue1 
pvalue2 <- t.test(tmp4$SITE ,tmp5$SITE,alternative="greater")$p.value ; pvalue2

tmp1 <- unique(tmp1[,c("STYLE","SITE","AREA","FROM","TO")])
tmp2 <- unique(tmp2[,c("STYLE","SITE","AREA","FROM","TO")])
tmp3 <- unique(tmp3[,c("STYLE","SITE","AREA","FROM","TO")])
tmp4 <- unique(tmp4[,c("STYLE","SITE","AREA","FROM","TO")])
tmp5 <- unique(tmp5[,c("STYLE","SITE","AREA","FROM","TO")])

phases[2,c("n_STYLES")] <- nrow(tmp1) ; phases[2,c("mean_SITES")] <- round(mean(tmp1[,c("SITE")],na.rm=T),0) ; phases[2,c("mean_AREA")] <- round(mean(tmp1[,c("AREA")],na.rm=T),3) 
phases[3,c("n_STYLES")] <- nrow(tmp2) ; phases[3,c("mean_SITES")] <- round(mean(tmp2[,c("SITE")],na.rm=T),0) ; phases[3,c("mean_AREA")] <- round(mean(tmp2[,c("AREA")],na.rm=T),3)  
phases[5,c("n_STYLES")] <- nrow(tmp3) ; phases[5,c("mean_SITES")] <- round(mean(tmp3[,c("SITE")],na.rm=T),0) ; phases[5,c("mean_AREA")] <- round(mean(tmp3[,c("AREA")],na.rm=T),3)  
phases[6,c("n_STYLES")] <- nrow(tmp4) ; phases[6,c("mean_SITES")] <- round(mean(tmp4[,c("SITE")],na.rm=T),0) ; phases[6,c("mean_AREA")] <- round(mean(tmp4[,c("AREA")],na.rm=T),3)  
phases[7,c("n_STYLES")] <- nrow(tmp5) ; phases[7,c("mean_SITES")] <- round(mean(tmp5[,c("SITE")],na.rm=T),0) ; phases[7,c("mean_AREA")] <- round(mean(tmp5[,c("AREA")],na.rm=T),3)  
phases


tmp <- tmp1
for (i in 1:nrow(tmp1)){ if(tmp1[i,c("STYLE")] %in% tmp2[,c("STYLE")]){
		if( mean(tmp1[i,c("FROM")] + tmp1[i,c("TO")]) < phases[1,c("TO")] ){ 	tmp2 <- tmp2[-which(tmp2[,c("STYLE")]==tmp1[i,c("STYLE")]) , ]  
		}else{											tmp <- tmp[-which(tmp[,c("STYLE")]==tmp1[i,c("STYLE")]) , ]} }} ; tmp1<-tmp 
tmp <- tmp4
for (i in 1:nrow(tmp4)){ if(tmp4[i,c("STYLE")] %in% tmp5[,c("STYLE")]){
		if( mean(tmp4[i,c("FROM")] + tmp4[i,c("TO")]) < phases[4,c("TO")] ){ 	tmp5 <- tmp5[-which(tmp5[,c("STYLE")]==tmp4[i,c("STYLE")]) , ]  
		}else{											tmp <- tmp[-which(tmp[,c("STYLE")]==tmp4[i,c("STYLE")]) , ]} }} ; tmp4<-tmp 
tmp <- tmp3
for (i in 1:nrow(tmp3)){ if(tmp3[i,c("STYLE")] %in% c(tmp1[,c("STYLE")],tmp2[,c("STYLE")],tmp4[,c("STYLE")],tmp5[,c("STYLE")])){tmp <- tmp[-which(tmp[,c("STYLE")]==tmp3[i,c("STYLE")]) , ]}};tmp3 <- tmp


if(ttest_on_bins == "FALSE"){
	pvalue1 <- t.test(tmp1$SITE ,tmp2$SITE,alternative="greater")$p.value ; pvalue1 
	pvalue2 <- t.test(tmp4$SITE ,tmp5$SITE,alternative="greater")$p.value ; pvalue2 }

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

area_per_stylegroup 		<- pottery.cent.meta 									; nrow(area_per_stylegroup)
area_per_stylegroup 		<- area_per_stylegroup[!is.na(area_per_stylegroup[,c("AREA")]),]  			; nrow(area_per_stylegroup)
area_per_stylegroup 		<- area_per_stylegroup[area_per_stylegroup[,c("AREA")]>0,]  			; nrow(area_per_stylegroup)
	
tmp1 <- area_per_stylegroup[ area_per_stylegroup[,c("AGE")]> phases[2,c("FROM")] & area_per_stylegroup[,c("AGE")]< phases[2,c("TO")] ,] 	; nrow(tmp1) 											; nrow(tmp1)
tmp2 <- area_per_stylegroup[ area_per_stylegroup[,c("AGE")]> phases[3,c("FROM")] & area_per_stylegroup[,c("AGE")]< phases[3,c("TO")] ,]	; nrow(tmp2)
tmp3 <- area_per_stylegroup[ area_per_stylegroup[,c("AGE")]> phases[5,c("FROM")] & area_per_stylegroup[,c("AGE")]< phases[5,c("TO")] ,]	; nrow(tmp3)
tmp4 <- area_per_stylegroup[ area_per_stylegroup[,c("AGE")]> phases[6,c("FROM")] & area_per_stylegroup[,c("AGE")]< phases[6,c("TO")] ,]	; nrow(tmp4)	
tmp5 <- area_per_stylegroup[ area_per_stylegroup[,c("AGE")]> phases[7,c("FROM")] & area_per_stylegroup[,c("AGE")]< phases[7,c("TO")] ,]	; nrow(tmp5)	

pvalue3 <- t.test(tmp1$AREA ,tmp2$AREA,alternative="greater")$p.value ; pvalue3
pvalue4 <- t.test(tmp4$AREA ,tmp5$AREA,alternative="greater")$p.value ; pvalue4 

tmp1 <- unique(tmp1[,c("STYLE","SITE","AREA","FROM","TO")])
tmp2 <- unique(tmp2[,c("STYLE","SITE","AREA","FROM","TO")])
tmp3 <- unique(tmp3[,c("STYLE","SITE","AREA","FROM","TO")])
tmp4 <- unique(tmp4[,c("STYLE","SITE","AREA","FROM","TO")])
tmp5 <- unique(tmp5[,c("STYLE","SITE","AREA","FROM","TO")])

tmp <- tmp1
for (i in 1:nrow(tmp1)){ if(tmp1[i,c("STYLE")] %in% tmp2[,c("STYLE")]){
		if( mean(tmp1[i,c("FROM")] + tmp1[i,c("TO")]) < phases[1,c("TO")] ){ 	tmp2 <- tmp2[-which(tmp2[,c("STYLE")]==tmp1[i,c("STYLE")]) , ]  
		}else{											tmp <- tmp[-which(tmp[,c("STYLE")]==tmp1[i,c("STYLE")]) , ]} }} ; tmp1<-tmp 
tmp <- tmp4
for (i in 1:nrow(tmp4)){ if(tmp4[i,c("STYLE")] %in% tmp5[,c("STYLE")]){
		if( mean(tmp4[i,c("FROM")] + tmp4[i,c("TO")]) < phases[4,c("TO")] ){ 	tmp5 <- tmp5[-which(tmp5[,c("STYLE")]==tmp4[i,c("STYLE")]) , ]  
		}else{											tmp <- tmp[-which(tmp[,c("STYLE")]==tmp4[i,c("STYLE")]) , ]} }} ; tmp4<-tmp 

if(ttest_on_bins == "FALSE"){
	pvalue3 <- t.test(tmp1$AREA ,tmp2$AREA,alternative="greater")$p.value ; pvalue3
	pvalue4 <- t.test(tmp4$AREA ,tmp5$AREA,alternative="greater")$p.value ; pvalue4 }

# ............................................................
# 5. Figure 4 ----
# ............................................................

col_expansion 	<- "gray50"
col_regionalization <- "gray30"
col_collapse 		<- "black"

# set filter for pottery groups highlighted ----
gab1 <- data.frame(STYLE = c("Okala/Epona/Yindo"), 
                   col = c("#39b600"))
gab2 <- data.frame(STYLE = c("Okanda", "Otoumbi", "Oveng"), 
                   col = c("#00c08d", "#95a900", "#cf9400"))
icb1 <- data.frame(STYLE = c("Imbonga"), 
                   col = c("#f8766d"))
icb2 <- data.frame(STYLE = c("Monkoto", "Lokondola", "Yete"), 
                   col = c("#ff71c9", "#ed8141", "#cf78ff"))
icb3 <- data.frame(STYLE = c("Bondongo"), 
                   col = c("#00abfd"))
icb4 <- data.frame(STYLE = c("Besongo", "Bolondo", "Malelembe"), 
                   col = c("#00bbdb", "#00c1aa", "#9590ff"))
filt <- rbind(gab1, gab2, icb1, icb2, icb3, icb4)
filt <- dplyr::arrange(filt, as.character(STYLE))

  # ............................................................
  # 5.1. Subfigure A: Frequency of pottery groups per 200 years ----
  # ............................................................

pottery.cent.freq <- as.data.frame(table(pottery.cent$AGE))
pottery.cent.freq$Var1 <- as.numeric(as.character(pottery.cent.freq$Var1))

ymax1 <- max(pottery.cent.freq$Freq)

plt.freq <- ggplot() +
  geom_rect(data = filter(rcarbon, FROM > -2000 & TO < 1900), 										 # added by WH
		aes(xmin = FROM, xmax = TO, ymin = 0, ymax = ymax1*1.40, fill = rcarbon), alpha = .1) + 				 # added by WH
  scale_fill_discrete("",labels = c("period of less intense human activity   ", "period of more intense human activity      ")) +
  geom_bar(data = filter(pottery.cent.freq, Var1 < 1800), 
           aes(x = Var1, weight = Freq), 
           fill = "white", color = "#333333", width = binyears*3/4) + 

  geom_bar(data = filter(pottery.cent.freq, Var1 > phases[2,c("FROM")], Var1 < phases[2,c("TO")]), aes(x = Var1 , weight = Freq),  # added by WH
			fill = col_expansion, color = "#333333", width = binyears*3/4) + 
  geom_bar(data = filter(pottery.cent.freq, Var1 > phases[3,c("FROM")], Var1 < phases[3,c("TO")]), aes(x = Var1 , weight = Freq),  # added by WH
			fill = col_regionalization, color = "#333333", width = binyears*3/4) + 
  geom_bar(data = filter(pottery.cent.freq, Var1 > phases[6,c("FROM")], Var1 < phases[6,c("TO")]), aes(x = Var1 , weight = Freq),  # added by WH
			fill = col_expansion, color = "#333333", width = binyears*3/4) + 
  geom_bar(data = filter(pottery.cent.freq, Var1 > phases[7,c("FROM")], Var1 < 1800), aes(x = Var1 , weight = Freq),
			fill = col_regionalization, color = "#333333", width = binyears*3/4) + 

  coord_cartesian(xlim = c(x_limits[1], 1850), 
                  #ylim = c(0, ymax1*1.3)) + 
                  ylim = c(0, ymax1*1.43)) + 

  annotate("text", x = (phases[2,c("FROM")]/2+phases[3,c("TO")]/2), y = ymax1*1.38, label = paste("EARLY IRON AGE"),hjust=0.5,vjust=1,cex=3,col="black") +	
  geom_line(aes(x=c(phases[2,c("FROM")]+20,phases[3,c("TO")]-20),y=c(ymax1*1.28,ymax1*1.28)),col="black",lwd=1)+

  annotate("text", x = (phases[6,c("FROM")]/2+phases[7,c("TO")]/2), y = ymax1*1.38, label = paste("LATE IRON AGE"),hjust=0.5,vjust=1,cex=3,col="black") +	
  geom_line(aes(x=c(phases[6,c("FROM")]+20,phases[7,c("TO")]-20),y=c(ymax1*1.28,ymax1*1.28)),col="black",lwd=1)+

  annotate("text", x = (phases[2,c("FROM")]/2+phases[2,c("TO")]/2), y = ymax1*1.20, label = paste("EXPANSION"),hjust=0.5,vjust=1,cex=3,col=col_expansion) +	
  geom_line(aes(x=c(phases[2,c("FROM")]+20,phases[2,c("TO")]-20),y=c(ymax1*1.1,ymax1*1.1)),col=col_expansion,lwd=1)+

  annotate("text", x = (phases[3,c("FROM")]/2+phases[3,c("TO")]/2), y = ymax1*1.20, label = paste("HIGH ACTIVITY"),hjust=0.5,vjust=1,cex=3,col=col_regionalization) +
  geom_line(aes(x=c(phases[3,c("FROM")]+20,phases[3,c("TO")]-20),y=c(ymax1*1.1,ymax1*1.1)),col=col_regionalization,lwd=1)+

  annotate("text", x = (phases[4,c("FROM")]/2+phases[4,c("TO")]/2), y = ymax1*1.0, label = paste("COLLAPSE"),hjust=0.5,vjust=1,cex=3,col=col_collapse) +
  geom_line(aes(x=c(phases[4,c("FROM")]+20,phases[4,c("TO")]-20),y=c(ymax1*0.9,ymax1*0.9)),col=col_collapse,lwd=1)+

  annotate("text", x = (phases[5,c("FROM")]/2+phases[5,c("TO")]/2), y = ymax1*0.85, label = paste("LOW ACTIVITY"),hjust=0.5,vjust=1,cex=3,col=col_collapse) +
  geom_line(aes(x=c(phases[5,c("FROM")]+20,phases[5,c("TO")]-20),y=c(ymax1*0.75,ymax1*0.75)),col=col_collapse,lwd=1)+

  annotate("text", x = (phases[6,c("FROM")]/2+phases[6,c("TO")]/2), y = ymax1*1.20, label = paste("EXPANSION"),hjust=0.5,vjust=1,cex=3,col=col_expansion) +
  geom_line(aes(x=c(phases[6,c("FROM")]+20,phases[6,c("TO")]-20),y=c(ymax1*1.1,ymax1*1.1)),col=col_expansion,lwd=1)+

  annotate("text", x = (phases[7,c("FROM")]/2+phases[7,c("TO")]/2)-50, y = ymax1*1.20, label = paste("HIGH ACTIVITY"),hjust=0.5,vjust=1,cex=3,col=col_regionalization) +
  geom_line(aes(x=c(phases[7,c("FROM")]+20,phases[7,c("TO")]-20),y=c(ymax1*1.1,ymax1*1.1)),col=col_regionalization,lwd=1)+

  scale_x_continuous(breaks = c(seq(-2000, 1800, 200)), 
                     expand = c(0,0)) + 
  scale_y_continuous("Number of  \n pottery groups", 
                     expand = c(0, 0)) + #, 
  theme_classic() + 
  theme(text = element_text(size=12),			# added by WH
	    	axis.title.x = element_blank(),
      	legend.margin = margin(t = 0, r = 0, b = 0, l = 2, unit = "cm"),
      	legend.box.spacing = unit(c(0, 0, 0.1, 0), "cm"),
        legend.position = "none", 
        legend.title = element_blank())+
  guides(fill = guide_legend(reverse = TRUE))		# added by WH

  # ............................................................
  # 5.2. Subfigure B: Frequency of sites per pottery group ----
  # ............................................................

ymax2 <- max(pottery.cent.meta$SITE,na.rm=T)										# added by WH

plt.sites <- ggplot() + 
  geom_rect(data = filter(rcarbon, FROM > -2000 & TO < 1900), 										 # added by WH
		aes(xmin = FROM, xmax = TO, ymin = 0, ymax = ymax2*1.6, fill = rcarbon), alpha = .1) + 				 # added by WH
  scale_fill_discrete("",labels = c("period of less intense human activity   ", "period of more intense human activity   ")) +  # added by WH

  geom_boxplot(data = pottery.cent.meta, 
               aes(x = AGE, y = SITE, group = AGE), 
               width = binyears*3/4, outlier.shape = 19) + 

  geom_boxplot(data = filter(pottery.cent.meta, SITE> 0, AGE > -500, AGE < phases[2,c("TO")]), aes(x = AGE, y = SITE, group = AGE),
               width = binyears*3/4, fill=col_expansion,
               outlier.shape = NA) +
  geom_boxplot(data = filter(pottery.cent.meta, SITE> 0, AGE > phases[3,c("FROM")], AGE < phases[3,c("TO")]), aes(x = AGE, y = SITE, group = AGE), 
               width = binyears*3/4, fill=col_regionalization,
               outlier.shape = NA) +
  geom_boxplot(data = filter(pottery.cent.meta, SITE> 0, AGE > phases[6,c("FROM")], AGE < phases[6,c("TO")]), aes(x = AGE, y = SITE, group = AGE), 
               width = binyears*3/4, fill=col_expansion,
               outlier.shape = NA) +
  geom_boxplot(data = filter(pottery.cent.meta, SITE> 0, AGE > phases[7,c("FROM")], AGE < phases[7,c("TO")]), aes(x = AGE, y = SITE, group = AGE), 
               width = binyears*3/4, fill=col_regionalization,
               outlier.shape = NA) +

  annotate("text", x = (-500/2+phases[2,c("TO")]/2), y = ymax2*1.25, label = paste("HOMOGENEITY"),hjust=0.5,vjust=1,cex=3,col=col_expansion) +	
  annotate("text", x = (phases[3,c("FROM")]/2+phases[3,c("TO")]/2), y = ymax2*1.25, label = paste("REGIONALIZATION"),hjust=0.5,vjust=1,cex=3,col=col_regionalization) +	
  geom_line(aes(x=c(-500,phases[2,c("TO")]-20),y=c(ymax2*1.1,ymax2*1.1)),col=col_expansion,lwd=1)+
  geom_line(aes(x=c(phases[3,c("FROM")]+20,phases[3,c("TO")]-20),y=c(ymax2*1.1,ymax2*1.1)),col=col_regionalization,lwd=1)+
  annotate("text", x = (phases[6,c("FROM")]/2+phases[6,c("TO")]/2)-25, y = ymax2*1.25, label = paste("HOMOGENEITY"),hjust=0.5,vjust=1,cex=3,col=col_expansion) +	
  annotate("text", x = (phases[7,c("FROM")]/2+phases[7,c("TO")]/2)-30, y = ymax2*1.25, label = paste("REGIONALIZATION"),hjust=0.5,vjust=1,cex=3,col=col_regionalization) +
  geom_line(aes(x=c(phases[6,c("FROM")]+20,phases[6,c("TO")]-20),y=c(ymax2*1.1,ymax2*1.1)),col=col_expansion,lwd=1)+
  geom_line(aes(x=c(phases[7,c("FROM")]+20,phases[7,c("TO")]-20),y=c(ymax2*1.1,ymax2*1.1)),col=col_regionalization,lwd=1)+

  coord_cartesian(xlim = c(x_limits[1], 1850), 
                  ylim = c(0, ymax2*1.3)) + 			
  scale_color_manual(values = as.character(filt$col)) + 
  scale_x_continuous(breaks = c(seq(-2000, 1800, 200)), 
                     expand = c(0,0)) + 
  scale_y_sqrt("Number of sites\n per pottery group", 
               breaks = c(0,5,10,25,50,75,100),
               expand = c(0, 0)) + 
  theme_classic() + 
  theme(legend.position = "none",
	text = element_text(size=12),				
        axis.title.x = element_blank())

  # ............................................................
  # 5.3. Subfigure C: Mean Distribution Area ----
  # ............................................................

ymax3 <- max(area_per_stylegroup$AREA,na.rm=T)

plt.area <- ggplot() + 
  geom_rect(data = filter(rcarbon, FROM > -2000 & TO < 1900), 										 # added by WH
		aes(xmin = FROM, xmax = TO, ymin = 0, ymax = ymax3*1.6, fill = rcarbon), alpha = .1) + 				 # added by WH
  scale_fill_discrete("",labels = c("period of less intense human activity   ", "period of more intense human activity   ")) +  # added by WH

  geom_boxplot(data = filter(pottery.cent.meta, AREA > 0),											# changed by WH: formerly no filter
               aes(x = AGE, y = AREA, group = AGE), 
               width = binyears*3/4, outlier.shape = 19) + 

  geom_boxplot(data = filter(pottery.cent.meta, AREA > 0, AGE > -500, AGE < phases[2,c("TO")]), aes(x = AGE, y = AREA, group = AGE),   # added by WH
               width = binyears*3/4, fill=col_expansion,
               outlier.shape = NA) +
  geom_boxplot(data = filter(pottery.cent.meta, AREA > 0, AGE > phases[3,c("FROM")], AGE < phases[3,c("TO")]), aes(x = AGE, y = AREA, group = AGE), 
               width = binyears*3/4, fill=col_regionalization,
               outlier.shape = NA) +
  geom_boxplot(data = filter(pottery.cent.meta, AREA > 0, AGE > phases[6,c("FROM")], AGE < phases[6,c("TO")]), aes(x = AGE, y = AREA, group = AGE), 
               width = binyears*3/4, fill=col_expansion,
               outlier.shape = NA) +
  geom_boxplot(data = filter(pottery.cent.meta, AREA > 0, AGE > phases[7,c("FROM")], AGE < phases[7,c("TO")]), aes(x = AGE, y = AREA, group = AGE), 
               width = binyears*3/4, fill=col_regionalization,
               outlier.shape = NA) +

  annotate("text", x = (-500/2+phases[2,c("TO")]/2), y = ymax3*1.25, label = paste("HOMOGENEITY"),hjust=0.5,vjust=1,cex=3,col=col_expansion) +	
  annotate("text", x = (phases[3,c("FROM")]/2+phases[3,c("TO")]/2), y = ymax3*1.25, label = paste("REGIONALIZATION"),hjust=0.5,vjust=1,cex=3,col=col_regionalization) +
  geom_line(aes(x=c(-500,phases[2,c("TO")]-20),y=c(ymax3*1.1,ymax3*1.1)),col=col_expansion,lwd=1)+
  geom_line(aes(x=c(phases[3,c("FROM")]+20,phases[3,c("TO")]-20),y=c(ymax3*1.1,ymax3*1.1)),col=col_regionalization,lwd=1)+
  annotate("text", x = (phases[6,c("FROM")]/2+phases[6,c("TO")]/2)-25, y = ymax3*1.25, label = paste("HOMOGENEITY"),hjust=0.5,vjust=1,cex=3,col=col_expansion) +	
  annotate("text", x = (phases[7,c("FROM")]/2+phases[7,c("TO")]/2)-30, y = ymax3*1.25, label = paste("REGIONALIZATION"),hjust=0.5,vjust=1,cex=3,col=col_regionalization) +
  geom_line(aes(x=c(phases[6,c("FROM")]+20,phases[6,c("TO")]-20),y=c(ymax3*1.1,ymax3*1.1)),col=col_expansion,lwd=1)+
  geom_line(aes(x=c(phases[7,c("FROM")]+20,phases[7,c("TO")]-20),y=c(ymax3*1.1,ymax3*1.1)),col=col_regionalization,lwd=1)+

  coord_cartesian(xlim = c(x_limits[1], 1850), 
                  ylim = c(0, ymax3*1.3)) + 
  scale_color_manual(values = as.character(filt$col)) + 
  scale_x_continuous("cal BC/AD", 
                     breaks = c(seq(-2000, 1800, 200)), 
                     expand = c(0,0)) + 
  scale_y_sqrt("Distribution area of \n pottery groups (1000 km^2)",
               #breaks = c(0, (1 * 2^(0:7))),
               breaks = c(0,5,10,25,50,75,100),
               expand = c(0,0)) + 
  theme_classic() + 
  theme(legend.position = "none",
	 text = element_text(size=12))

  # ............................................................
  # 5.4. Subfigure D: Maps ----
  # ............................................................

# Base Map ----
plt.basemap <- ggplot() + 
  geom_sf(data = coast10, size = .5, color = '#808080') + 
  geom_sf(data = rivers10, size = .5, color = '#808080') + 
  geom_sf(data = lakes10, fill = '#808080', color = NA) + 
  geom_sf(data = boundary_lines_land10, size = .1, color = 'black') + 
  theme_few() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.caption = element_text(hjust = 0.5), 
        legend.position = "none", 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
		plot.margin = unit(c(0, 0, 0, 0), "cm"),			# added by WH
        plot.background = element_rect(color = NA, fill = NA))

# Map function ----
plt.map.fct <- function(x, y, z){
  class <- dplyr::filter(pottery.cent.meta, STYLE %in% x[[1]])$CLASS
  sites.class <- dplyr::filter(sites.cent, CLASS %in% class)
  plt.map <- plt.basemap + 
    geom_sf(data = sites.class, color = "grey", size = 1) + 
    geom_sf(data = dplyr::filter(pottery.sites.area, STYLE %in% x[[1]]), 
            fill = "black", color = NA, alpha = .1) + 
    geom_sf(data = dplyr::filter(sites, STYLE %in% x[[1]]), 
            aes(fill = STYLE), shape = 21, color = "white") + 
    scale_fill_manual(values = as.character(x$col)) + 
    labs(title = z)  + 
    theme(plot.title = element_text(size = 12, 
                                    colour= "white")	,					# changed by WH: size was 10
		plot.margin = unit(c(0.2, 0.2, 0.2,0.2), "cm"))

  if(y == "gab"){
    plt.map <- plt.map + coords.gab
  }else if(y == "icb"){
    plt.map <- plt.map + coords.icb
  }
  
  plt.txt <- ggplot() + 
    geom_text(data = x, 
              aes(x = 1, y = as.numeric(row.names(x)), 
                  label = STYLE, color = STYLE),
              size = 4.5, fontface = "bold") + 
    scale_color_manual(values = as.character(x$col)) + 
    scale_y_reverse(limit = c(4,0.8)) + 			# changed by WH: limit was c(4,0)
    theme_void() + 
    theme(legend.position = "none",
		plot.margin = unit(c(0, 0, 0, 0), "cm"),			# added by WH
		text = element_text(size=12) )		# line added by WH
  
  plt.map.comp <- cowplot::plot_grid(plt.map, 
                                     plt.txt, 
                                     ncol = 1, 
                                     rel_heights = c(2.5,1)) 
  
  return(plt.map.comp)
}

# Coords ----
coords.icb <- coord_sf(xlim = c(16, 23.5), 
                   ylim = c(-2, 2.1))
coords.gab <- coord_sf(xlim = c(8, 15.5), 
                       ylim = c(-2.2, 1.9))


# Perform functions ----

plt.map.gab1 <- plt.map.fct(gab1, "gab", "~200 BC")+	# changed by WH: was 400-200 BC
		theme(plot.background = element_rect(color = "white", size=4, fill = col_expansion))
plt.map.gab2 <- plt.map.fct(gab2, "gab", "AD ~200")+
		theme(plot.background = element_rect(color = "white", size=4, fill = col_regionalization))

plt.map.icb1 <- plt.map.fct(icb1, "icb", "~200 BC")+	# changed by WH: was 400-200 BC
		theme(plot.background = element_rect(color = "white", size=4, fill = col_expansion))
plt.map.icb2 <- plt.map.fct(icb2, "icb", "AD ~200")+
		theme(plot.background = element_rect(color = "white", size=4, fill = col_regionalization))
plt.map.icb3 <- plt.map.fct(icb3, "icb", "AD ~1100")+
		theme(plot.background = element_rect(color = "white", size=4, fill = col_expansion))
plt.map.icb4 <- plt.map.fct(icb4, "icb", "AD ~1700")+
		theme(plot.background = element_rect(color = "white", size=4, fill = col_regionalization))

plt.map.gab <- cowplot::plot_grid(plt.map.gab1, plt.map.gab2,
                                  nrow = 1, align = "h", axis = "t")
plt.map.icb <- cowplot::plot_grid(plt.map.icb1, plt.map.icb2, plt.map.icb3, plt.map.icb4,
                                  nrow = 1, align = "h", axis = "t",rel_widths = c(1,1,1,1))

plt.map <- plot_grid(plt.map.icb,
				plot_grid(plt.map.gab, NULL, rel_widths = c(1,1)),
                     ncol = 1)

# ............................................................
# 6. Composed Figures ----
# ............................................................

# dummy figure for legend ----
plt.legend <- ggplot() + 
  geom_rect(data = rcarbon, 
            aes(xmin = FROM, xmax = TO, ymin = -Inf, ymax = Inf, fill = rcarbon), alpha = .3) + 
  scale_fill_manual("",values = c("#00bfc4", "#f8766d"), 
                    labels = c("more intense human activity", 
                               "less intense human activity")  ) +
  theme(legend.position = "top")
legend <- cowplot::get_legend(plt.legend)

plt <- cowplot::plot_grid(plt.freq,
                          plt.sites, 
                          plt.area,
                          plt.map.icb,
                          ncol = 1, 
                          align = "v", axis = "lr", 
                          #labels = "auto", rel_heights = c(0.8, 1, 1.2, 0.8))
                          labels = "auto", rel_heights = c(1, 1, 1, 0.8))

plt <- cowplot::plot_grid(legend, 
                          plt, 
                          ncol = 1,
                          rel_heights = c(1, 30))

windows(14,10) ; plt
ggsave("output/Figure 3 pottery_groups_regionalization_icb.pdf", plt, width = 8, height = 10)
ggsave("output/Figure 3 pottery_groups_regionalization_icb.jpg", plt, width = 8, height = 10)