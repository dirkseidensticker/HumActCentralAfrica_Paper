graphics.off()
rm(list=ls(all=TRUE))

library(boot)
library(concaveman)
library(cowplot)
library(data.table)
library(dunn.test)
library(geojsonsf)
library(ggplot2)
library(ggrepel)
library(ggsn)
library(ggthemes)
library(cowplot)
library(dplyr)
library(lwgeom)
library(openxlsx)
library(raster)
library(rcarbon)
library(rgdal)
library(rnaturalearth)
library(rnaturalearthdata)
library(parallel)
library(RColorBrewer)
library(reshape2)
library(scales)
library(sf)
library(stringr)
library(tidyr)
library(viridis)
library(zoo)

ncores = (detectCores() - 1)

rm(list=ls(all=TRUE))
graphics.off()

quiet <- function(x) { sink(tempfile()) ; on.exit(sink()) ; invisible(force(x)) } 

