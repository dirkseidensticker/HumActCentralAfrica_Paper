#### TABLES ####

#### OVERVIEW OF CLASSES
#### OVERVIEW OF REGIONS
#### OVERVIEW OF RADIOCARBON DATES
#### OVERVIEW OF POTTERY STYLES
#### OVERVIEW OF SITES

source("script/pkg.R")
source("script/fct.R")

# ............................................................
# 1. READ INPUT FILES ----
# ............................................................

c14 <- data.table::fread("input/c14.csv", 
                         encoding = "UTF-8")
c14.orig <- c14

sites	<- data.table::fread("input/sites.csv",
                           encoding = "UTF-8")

pottery <- data.table::fread("input/potterygroups.csv", 
                             encoding = "UTF-8")

region.labels <- data.table::fread("input/region labels.csv", 
                                   encoding = "UTF-8")

references <- read.table("input/references.csv", 			
                         encoding = "UTF-8", 
                         header = T)

# ............................................................
# 2. CORRECT DOUBLE SITE NAMES in sites  ----							
#		problem: 	9 site names were used twice, in different regions
#				therefore, some sites are kicked out from site counts in Supplementary Table 2 and 4
#		solution:	add the coordinates behind the problematic site names (separated with "_") to make these site names unique 
#				the merger between sites and coords is undone when Supp Table 3 is finished, so it's not written to the output
# ............................................................

check <- unique(sites[,c("SITE","LAT","LONG")]) 
check <- aggregate(LAT~SITE,data=check,FUN=length) ; names(check) <- c("SITE","N_coords")
check <- check[ check[,c("N_coords")] >  1,] ; nrow(check) ; check 

sites$SITE_ori <- sites$SITE
sites$change <- rep("no")
for(i in 1:nrow(check)){ sites$change <- replace(sites[,c("change")],sites[,c("SITE")]==check[i,c("SITE")],"change" ) }
for(i in 1:nrow(sites)){ 
	if(sites[i,c("change")] == "change"){		
		sites[i,c("SITE")] <- paste(sites[i,c("SITE")],sites[i,c("LAT")],sites[i,c("LONG")],sep=" ")  }}
corrected_sitenames <- as.data.frame(sites) ;  corrected_sitenames <- corrected_sitenames[corrected_sitenames[,c("change")]=="change",]
corrected_sitenames <- unique(corrected_sitenames[,c("SITE_ori","SITE")])
sites$SITE_ori <- NULL
sites$change <- NULL

tmp <- as.data.frame(c14[,c("LABNR","SITE","LAT","LONG")])
tmp$change <- rep("no")
for(i in 1:nrow(check)){ tmp$change <- replace(tmp[,c("change")],tmp[,c("SITE")]==check[i,c("SITE")],"change" ) }
for(i in 1:nrow(tmp)){ 
	if(tmp[i,c("change")] == "change"){	
		tmp[i,c("SITE")] <- paste(tmp[i,c("SITE")],tmp[i,c("LAT")],tmp[i,c("LONG")],sep=" ")  }}
tmp$change <- NULL
c14$SITE<- tmp$SITE

check <- unique(sites[,c("SITE","LAT","LONG")]) 
check <- aggregate(LAT~SITE,data=check,FUN=length) ; names(check) <- c("SITE","N_coords")
check <- check[ check[,c("N_coords")] >  1,] 
nrow(check)  # THIS CHECK SHOULD BE ZERO

# ............................................................
# 3. OVERVIEW OF CLASSES ----
# ............................................................

table_classes <- stats::aggregate(LABNR ~ CLASS, 
                           data = c14,
                           FUN = length)

names(table_classes) <- c("Class","N_all_regions")

table_classes$Archaeological_association <- c(
  "relevant dates, strong archaeological context",
  "relevant dates, medium archaeological context",
  "relevant dates, poor archaeological context",
  "relevant dates, strong proxy for human activity (charred fruit remains but no artefacts)",
  "irrelevant dates, related to lithic artefacts only",
  "irrelevant dates, doubtful archaeological context (post-depositional mixing)",
  "irrelevant dates, lack of archaeological context",
  "unreliable dates due to presumed lab errors",
  "unreliable dates on shells (biased due to old-carbon age offsets)",
  "unreliable dates on bulk organic matter (from sediment cores)"
  )

# Regions A-H:
tmp1 <- dplyr::filter(c14, REGION %in% LETTERS[seq(from = 1, to = 8)])
tmp1 <- stats::aggregate(LABNR ~ CLASS, 
                         data = tmp1,
                         FUN = length) ; names(tmp1) <- c("Class","N_regions_A-H")


tmp2 <- dplyr::filter(c14, REGION %in% LETTERS[seq(from = 9, to = 11)])
tmp2 <- stats::aggregate(LABNR ~ CLASS, 
                         data = tmp2,
                         FUN = length) ; names(tmp2) <- c("Class","N_regions_I-K")

table_classes <- table_classes %>%
  dplyr::left_join(tmp1, by = "Class") %>%
  dplyr::left_join(tmp2, by = "Class") %>%
  dplyr::select(Class, 
                Archaeological_association, 
                N_all_regions, 
                `N_regions_A-H`, 
                `N_regions_I-K`)

table_classes[is.na(table_classes)] <- 0

# ............................................................
# 4. OVERVIEW OF RADIOCARBON DATES ----
# ............................................................

table_dates <- c14.orig

add <- table_classes[,c("Class","Archaeological_association")] ; names(add)[1]<-"CLASS"	
table_dates <- merge(table_dates,add,by="CLASS")
table_dates$CLASS <- paste(table_dates$CLASS,table_dates$Archaeological_association,sep=" ")

table_dates <- table_dates %>%
  dplyr::arrange(REGION,
                 CLASS,
                 SITE) %>%
  dplyr::select("LABNR",
                "C14AGE",
                "C14STD",
                "MATERIAL",
                "COUNTRY",
                "REGION",
                "SITE",
                #"FEATURE",		
                #"FEATURE_DESC",
                "LAT",
                "LONG",
                "CLASS",
                "LITHICS",
                "POTTERY",
                "IRON",
                "SOURCES")

# ............................................................
# 5. OVERVIEW OF POTTERY STYLES ----
# ............................................................

#c14$CLASS <- gsub("[^::A-Z::]","", c14$CLASS) # reduce to main classes

c14$CLASS <- sub("IIa", "II", c14$CLASS )	
c14$CLASS <- sub("IIb", "II", c14$CLASS )
c14$CLASS <- sub("IIc", "II", c14$CLASS )
c14$CLASS <- sub("IIIa", "III", c14$CLASS )	
c14$CLASS <- sub("IIIb", "III", c14$CLASS )
c14$CLASS <- sub("IIIc", "III", c14$CLASS )
sort(unique(c14$CLASS))

  # ............................................................
  # 5.1. DATES PER POTTERY STYLE FROM ADRAC ----
  # ............................................................

tmp1 <- c14 %>% 
  dplyr::filter(CLASS == "Ia") %>%
  dplyr::select("LABNR", "POTTERY")

tmp1.n_dates <- unique(tmp1$POTTERY)
tmp1.n_dates <- unlist(strsplit(as.character(tmp1.n_dates), "; ")) # split fields with multiple entries
tmp1.n_dates <- gsub("\\(", "", tmp1.n_dates) # remove parantheses
tmp1.n_dates <- gsub("\\)", "", tmp1.n_dates)
tmp1.n_dates <- data.frame(POTTERY = unique(tmp1.n_dates))
tmp1.n_dates$N_DATES <- NA

for(i in 1:nrow(tmp1.n_dates)){
  d <- dplyr::filter(tmp1, 
                     grepl(tmp1.n_dates[i,"POTTERY"], 
                           tmp1$POTTERY)) # filter for dates related to style
  d <- dplyr::filter(d, 
                     !grepl(paste0("\\(" , tmp1.n_dates[i,"POTTERY"], "\\)"), 
                            d$POTTERY)) # remove cases in parantheses  
  tmp1.n_dates[i,"N_DATES"] <- nrow(d)
}

  # ............................................................
  # 5.2. REGIONS PER POTTERY STYLE FROM SPACIALIST ----
  # ............................................................

tmp2.n_regions <- unique(sites[,c("POTTERY","REGION")])
tmp2.n_regions <- stats::aggregate(REGION ~ POTTERY, 
                                   data = tmp2.n_regions,
                                   FUN = length) ; names(tmp2.n_regions) <- c("POTTERY","N_REGIONS")

  # ............................................................
  # 5.3. SITES PER POTTERY STYLE FROM SPACIALIST ----
  # ............................................................

tmp3.n_sites <- unique(sites[,c("POTTERY","SITE")])
tmp3.n_sites <- stats::aggregate(SITE ~ POTTERY,
                                 data = tmp3.n_sites,
                                 FUN = length) ; names(tmp3.n_sites) <- c("POTTERY","N_SITES")
tmp3.n_sites$MAJOR_REGION <- NA

for(i in 1:nrow(tmp3.n_sites)){
  a <- dplyr::filter(sites, POTTERY == tmp3.n_sites[i,c("POTTERY")])
  a <- stats::aggregate(SITE ~ REGION, 
                         data = a,
                         FUN = length) 
  a <- a[  a[,c("SITE")]== max(a[,c("SITE")]),]
  tmp3.n_sites[i,c("MAJOR_REGION")] <- as.character(a[1,1])
}

  # ............................................................
  # 5.4. WRAP TABLE - POTTERY STYLES ----
  # ............................................................

table_pottery <- pottery %>%
  dplyr::left_join(tmp1.n_dates, by = "POTTERY") %>%
  dplyr::left_join(tmp2.n_regions, by = "POTTERY") %>%
  dplyr::left_join(tmp3.n_sites, by = "POTTERY") %>%
  replace(is.na(.), 0) %>% 				
  dplyr::select(#"ID", 	
                "POTTERY",
                "FROM",
                "TO",
                "MAJOR_REGION",
                "N_REGIONS",
                "N_DATES",
                "N_SITES",
                "DESCRIPTION")

  # ............................................................
  # 5.5. CHECK ----
  # ............................................................

sum(table_pottery[,c("N_DATES")], na.rm = T)	# total nr of 14C dates associated to styles in Supp Fig 1
nrow(c14[c14$CLASS == "Ia", ])			# total nr of 14C dates associated to styles in adrac ; this number is lower than the one above, as some dates were used for two or more styles

  # ............................................................
  # 5.6. FOOTNOTES ----
  # ............................................................

names(table_pottery)[which(names(table_pottery)=="N_DATES")] <- "N_DATES *"

double_dates <- nrow(c14[  grepl(";",c14$POTTERY) &  c14$CLASS == "Ia" ,])		
add <- table_pottery[1:2,] ; add[,2] <- as.character(add[,2]) ; for(j in 1:ncol(add)){ add[1,j] <- "" ; add[2,j] <- ""}
add[1,1] <- paste("*",double_dates,	" radiocarbon dates are attributed to two or three different pottery styles occurring in an intimately mixed assemblage.")
add[2,1] <- c(				"Therefore, the sum of radiocarbon dates of all pottery styles exceeds the number of class II radiocarbon dates in Supplementary Table 1.")
table_pottery <- rbind(table_pottery,add)

# ............................................................
# 6. OVERVIEW OF SITES ----
# ............................................................

  # ............................................................
  # 6.1. DATES PER CLASS PER SITE ----
  # ............................................................

table_sites <- reshape2::dcast(c14, SITE  ~ CLASS, 			# changed by WH: deleted "+ REGION + LAT + LONG"
                             value.var = "LABNR",
                             fun.aggregate = length)
colnames(table_sites)[2:7] <- paste("N_DATES_CLASS_", colnames(table_sites)[2:7], sep = "")		

tmp4.n_dates <- c14 %>% 
  dplyr::group_by(SITE) %>% 						# changed by WH: deleted ", LAT, LONG"
  dplyr::summarise(N_DATES = length(SITE))

  # ............................................................
  # 6.2. POTTERY STYLES PER SITE ----
  # ............................................................

tmp4.n_pottery <- stats::aggregate(POTTERY ~ SITE, 
                                   data = sites, 
                                   FUN = length)
colnames(tmp4.n_pottery) <- c("SITE", "N_POTTERYSTYLES")

# CHECK: c14 class Ia sites not in sites
check2 <- table_sites[(table_sites[,c("N_DATES_CLASS_Ia")] != 0 ),]
check2 <- check2[ !(check2[,c("SITE")] %in% tmp4.n_pottery[,c("SITE")]),] 
nrow(check2)	# THIS CHECK SHOULD BE ZERO	 

  # ............................................................
  # 6.3. METADATA ----	
  # ............................................................

meta1 <- unique(as.data.frame(sites[,c("SITE","REGION","LAT","LONG")]))
meta2 <- unique(as.data.frame(c14[,c("SITE","REGION","LAT","LONG")]) )
meta2 <- meta2[ !(meta2[,c("SITE")] %in% meta1[,c("SITE")]),] 	# sites have better coords, so prefer sites
meta <- unique(rbind(meta1,meta2)) #; nrow(meta)

#head(table_sites)
#head(tmp4.n_dates)
#head(tmp4.n_pottery)
#head(meta)

  # ............................................................
  # 6.3. WRAP TABLE - SITES ----
  # ............................................................

table_sites <- table_sites %>% 
  dplyr::full_join(tmp4.n_dates, by = "SITE") %>%		# corrected WH: full_join instead of left_join; deleted , "LAT", "LONG"
  dplyr::full_join(tmp4.n_pottery, by = "SITE") %>%		# corrected WH: full_join 
  dplyr::full_join(meta, by = "SITE") %>%			# added WH
  replace(is.na(.), 0) %>% 
  dplyr::select("SITE",
                "REGION",
                "LAT",
                "LONG",
                "N_POTTERYSTYLES",
                "N_DATES",
                "N_DATES_CLASS_Ia",
                "N_DATES_CLASS_Ib",
                "N_DATES_CLASS_Ic",
                "N_DATES_CLASS_Id",
                "N_DATES_CLASS_II",
                "N_DATES_CLASS_III") %>%			
  dplyr::ungroup() %>%
  dplyr::arrange(REGION,
                 SITE)

for(i in 1:nrow(corrected_sitenames)){
	table_sites$SITE <- replace(table_sites[,c("SITE")],table_sites[,c("SITE")]== corrected_sitenames[i,c("SITE")],corrected_sitenames[i,c("SITE_ori")]) }

  # ............................................................
  # 6.4. CHECK ----
  # ............................................................

check3 <- table_sites
check3$SUMS <- rowSums(check3[which(colnames(table_sites)=="N_DATES_CLASS_Ia"):which( colnames(table_sites)=="N_DATES_CLASS_III")])   
check3$DIFF <- check3$N_DATES - check3$SUMS
nrow(filter(check3, DIFF != 0))	# THIS CHECK SHOULD BE ZERO

# ............................................................
# 7. OVERVIEW OF REGIONS ----
# ............................................................

  # ............................................................
  # 7.1. DATES PER REGION ----
  # ............................................................

table_regions <- reshape2::dcast(c14, REGION ~ CLASS, 
                             value.var = "LABNR", 
                             fun.aggregate = length)

table_regions$N_DATES_TOTAL <- rowSums(table_regions[,c(which(colnames(table_regions)=="Ia"):
                                                which(colnames(table_regions)=="III"))])					

colnames(table_regions)[which(colnames(table_regions)=="Ia"):
                    which(colnames(table_regions)=="III")] <- 
  paste("N_DATES_CLASS_", 
        colnames(table_regions)[which(colnames(table_regions)=="Ia"):
                            which(colnames(table_regions)=="III")], 
        sep = "")

  # ............................................................
  # 7.2. SITES PER REGION ----
  # ............................................................

tmp6 <- sites %>%
  dplyr::distinct(SITE, LAT, LONG, REGION) %>%
  dplyr::group_by(REGION) %>%
  dplyr::summarise(N_SITES_potterygroups = length(SITE))

tmp7 <- c14 %>%
  dplyr::distinct(SITE, LAT, LONG, REGION) %>%
  dplyr::filter(!SITE %in% sites$SITE) %>%
  dplyr::group_by(REGION) %>%
  dplyr::summarise(N_SITES_other = length(SITE))

# CHECK: sites with coords from multiple regions
check1 <- unique(sites[,c("SITE","REGION")])
check1 <- aggregate(REGION ~ SITE,
                    data = check1,
                    FUN = length) ; names(check1) <- c("SITE","N_REGIONS")
check1 <- check1[ check1[,c("N_REGIONS")] >  1,]
nrow(check1)	# THIS CHECK SHOULD BE ZERO (as a result of the prep script)

  # ............................................................
  # 7.3. STYLES PER REGION ----
  # ............................................................

tmp8 <- c14 %>% 
  dplyr::filter(CLASS == 'Ia') %>%
  dplyr::select(REGION, POTTERY) %>% 
  mutate(POTTERY = strsplit(POTTERY, "; ")) %>% # separate multiple entries
  unnest(POTTERY)

tmp8 <- tmp8[!grepl("\\(", tmp8$POTTERY),] # remove cases in parantheses
tmp8 <- unique(tmp8)
tmp8 <- filter(tmp8, POTTERY != '-' & POTTERY != 'indet')

tmp8d <- aggregate(POTTERY ~ REGION,
                   data = tmp8,
                   FUN = length) ; names(tmp8d) <- c("REGION","N_POTTERYSTYLES_dated")

tmp8$pottery_region <- paste(tmp8$POTTERY, tmp8$REGION, sep="-")

tmp9 <- as.data.frame(unique(sites[,c("POTTERY","REGION")]))
tmp9$pottery_region <- paste(tmp9$POTTERY, tmp9$REGION, sep="-")
tmp9 <- tmp9[!(tmp9$pottery_region %in% tmp8$pottery_region),] # remove dated styles

tmp9f <- dplyr::filter(tmp9, POTTERY %in% tmp8$POTTERY) # styles that are dated in general/but in different regions not
tmp9f <- aggregate(POTTERY ~ REGION, # count of per-se dated styles in another region
                   data = tmp9f, 
                   FUN = length) ; names(tmp9f) <- c("REGION","N_POTTERYSTYLES_dated")

tmp8d <- aggregate(N_POTTERYSTYLES_dated ~ REGION, # combine dated sites in main region + other regions
                   data = rbind(tmp8d, tmp9f), 
                   FUN = sum)

tmp9b <- tmp9[!(tmp9$POTTERY %in% tmp8$POTTERY),] # remove dated groups
tmp9d <- aggregate(POTTERY ~ REGION,
                   data = tmp9b, 
                   FUN = length) ; names(tmp9d) <- c("REGION","N_POTTERYSTYLES_undated")

  # ............................................................
  # 7.4. WRAP TABLE - REGIONS ----
  # ............................................................

table_regions <- table_regions %>%
  dplyr::left_join(tmp6, by = "REGION") %>%
  dplyr::left_join(tmp7, by = "REGION") %>%
  dplyr::left_join(tmp8d, by = "REGION") %>%
  dplyr::left_join(tmp9d, by = "REGION") %>%
  dplyr::left_join(region.labels[,c("REGION", "LABEL")], by = "REGION") %>% 
  replace(is.na(.), 0) %>%
  dplyr::mutate(N_POTTERYSTYLES_TOTAL = N_POTTERYSTYLES_dated + N_POTTERYSTYLES_undated) %>%
  dplyr::mutate(N_SITES_TOTAL = N_SITES_potterygroups + N_SITES_other) %>%
  tidyr::unite(REGION, REGION, LABEL, sep = ") ") %>%
  dplyr::select("REGION",
                "N_DATES_TOTAL",
                "N_DATES_CLASS_Ia",
                "N_DATES_CLASS_Ib",
                "N_DATES_CLASS_Ic",
                "N_DATES_CLASS_Id",
                "N_DATES_CLASS_II",
                "N_DATES_CLASS_III",	
                "N_POTTERYSTYLES_TOTAL",
                "N_POTTERYSTYLES_dated",
                "N_POTTERYSTYLES_undated",
                "N_SITES_TOTAL",
                "N_SITES_potterygroups",
                "N_SITES_other")

  # ............................................................
  # 7.5. Table footnotes ----
  # ............................................................

names(table_regions)[which(names(table_regions)=="N_POTTERYSTYLES_TOTAL")] <- "N_POTTERYSTYLES_TOTAL *"
names(table_regions)[which(names(table_regions)=="N_SITES_TOTAL")] <- "N_SITES_TOTAL **"

double_styles <- nrow(table_pottery[table_pottery[,c("N_REGIONS")] > 1,])
add <- table_regions[1:2,] ; add[,1] <- as.character(add[,1]) ; for(j in 1:ncol(add)){ add[1,j] <- "" ; add[2,j] <- ""}
add[1,1] <- paste("*",double_styles," Pottery styles occur in two regions. Therefore, the sum of pottery styles of all regions exceeds the number of pottery styles in Supplementary Table 1.")
add[2,1] <- paste("**",			" Sites with both pottery groups as well as unclassified pottery assemblages are counted only as sites with pottery groups.")
table_regions <- rbind(table_regions,add)

# ............................................................
# 8. WRITE OUTPUT ----
# ............................................................

write.csv(table_classes,
          "output/Table S1 Classes.csv", 
          fileEncoding = "UTF-8", 
          row.names = F)

write.csv(table_regions,
          "output/Table S2 Regions.csv",
          fileEncoding = "UTF-8",
          row.names = F)

write.csv(table_dates,
          "output/Data S1 Radiocarbon dates.csv", 
          fileEncoding = "UTF-8",
          row.names = F)

write.csv(table_pottery,
          "output/Data S2 Pottery styles.csv", 
          fileEncoding = "UTF-8",
          row.names = F)

#write.csv(references,
#          "output/Data S3 References.csv",
#          fileEncoding = "UTF-8", 
#          row.names = F)

write.csv(table_sites,
          "output/Data S4 Sites.csv",
          fileEncoding = "UTF-8", 
          row.names = F)

#tab.lst <- list("Data S1 Radiocarbon dates" = table_dates,
#                "Data S2 Pottery styles" = table_pottery,
#                #"Data S3 References" = references,
#               "Data S4 Sites" = table_sites)
#write.xlsx(tab.lst, file = "output/Data S1-S4.xlsx")

# ............................................................
# 9. INTEGRITY CHECKS ----
# ............................................................

# CHECK NR OF DATES
sum(table_classes$N_all_regions)
nrow(table_dates)
sum(table_sites$N_DATES)
sum(as.numeric(table_regions$N_DATES_TOTAL),na.rm=T)

# CHECK NR OF REGIONS ASSOCIATED WITH STYLES
sum(as.numeric(table_pottery[,c("N_REGIONS")]),na.rm=T)
sum(as.numeric(table_regions[,c("N_POTTERYSTYLES_TOTAL *")]),na.rm=T)

sum(as.numeric(table_pottery[,c("N_SITES")]),na.rm=T)
sum(as.numeric(table_sites[,c("N_POTTERYSTYLES")]),na.rm=T)

sum(as.numeric(table_regions[,c("N_SITES_TOTAL **")]),na.rm=T)
nrow(table_sites)

sum(as.numeric(table_regions[,c("N_SITES_potterygroups")]),na.rm=T)
nrow(table_sites[ table_sites[,c("N_POTTERYSTYLES")] >0,])

# ad number for indet sites == 223