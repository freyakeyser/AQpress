#' @title Propagule pressure calculation

#' @description calcAQpress will calculate propagule pressure for rivers using information on the locations and use of aquaculture sites.
#'
#' Propagule pressure is calculated for each river, and is equal to the sum (for all aquaculture sites) of 1/(distance between the river and each stocked aquaculture site).
#'
#' All coordinates must be in decimal degrees, with E/W hemisphere specified using +/-, and all of your coordinates must be in the water.
#' When in doubt, assign river coordinates to the river's mouth.
#' Using an upstream location may cause the function to fail, as the raster resolution may place such locations on land.
#'
#' Two CSV files are created and saved to your specified working directory. prop.press_final.csv is the year-specific propagule pressure for each river. prop.press_avg.csv is the average propagule pressure for each river.
#'
#' @param AQsites Your aquaculture site dataframe. Must have columns: Site.ID, Lat, Long, prov
#' @param rivercoords Your river coordinate dataframe. Must have columns: River, Lat, Long
#' @param inventory Your inventory dataframe. Must have columns: Site.ID, Lat, Long, Year, prov, totalfish (Note: totalfish does not have to be perfectly accurate as it will be replaced by 1's to indicate that a site was stocked in a given year. Do not include 0's for unstocked (fallow) sites.)
#' @param dir The directory where you would like to write csv files
#' @import rgeos
#' @import plyr
#' @import raster
#' @importFrom reshape2 melt
#' @importFrom geosphere distm distHaversine
#' @importFrom gdistance transition geoCorrection costDistance
#' @importFrom utils write.csv
#' @importFrom sp SpatialPoints CRS
#' @export
#' @rdname calcAQpress


calcAQpress <- function(AQsites, rivercoords, inventory, dir, inputRast, rastName, saveRast, distType){

  options(scipen = 999)

  AQsites <- subset(AQsites, select=c(Site.ID, Lat, Long, prov))
  rivercoords <- subset(rivercoords, select=c(River, Lat, Long))
  inventory <- subset(inventory, select=c(Site.ID, Lat, Long, Year, prov, totalfish))

  AQsites$Site.ID <- as.character(AQsites$Site.ID)
  AQsites$Lat <- as.numeric(as.character(AQsites$Lat))
  AQsites$Long <- as.numeric(as.character(AQsites$Long))

  rivercoords$River <- as.character(rivercoords$River)
  rivercoords$Lat <- as.numeric(as.character(rivercoords$Lat))
  rivercoords$Long <- as.numeric(as.character(rivercoords$Long))

  inventory$Site.ID <- as.character(inventory$Site.ID)
  inventory$Lat <- as.numeric(as.character(inventory$Lat))
  inventory$Long <- as.numeric(as.character(inventory$Long))
  inventory$prov <- as.character(inventory$prov)
  inventory$totalfish <- as.numeric(as.character(inventory$totalfish))

  alllong <- c(AQsites$Long, rivercoords$Long)
  alllat <- c(AQsites$Lat, rivercoords$Lat)
  extentAQ <- c(min(alllong)-0.25, max(alllong)+0.25, min(alllat)-0.25, max(alllat)+0.25)

  if(inputRast==TRUE){
    rasterAQ <- rastName
  }
  else{
  # get shapefiles
  NAm <- crop(getData("GADM", country=c("CAN", "USA"), level=1), extent(extentAQ))

  print("(1/7) Found shapefile")

  horizdist <- distm(c(extentAQ[1], extentAQ[3]), c(extentAQ[2], extentAQ[3]), fun=distHaversine)/1000
  vertdist <- distm(c(extentAQ[1], extentAQ[3]), c(extentAQ[1], extentAQ[4]), fun=distHaversine)/1000

  ncolrast <- ifelse(horizdist < 1000 & vertdist < 1000, round(horizdist*10,0), round(horizdist*4,0))
  nrowrast <- ifelse(horizdist < 1000 & vertdist < 1000, round(vertdist*10,0), round(vertdist*4,0))

  r <- raster(ncols=ncolrast, nrows=nrowrast)
  extent(r) <- extentAQ
  rasterAQ <- rasterize(NAm, r, fun="first")

  # change land values to make them impassable
  rasterAQ@data@values <- ifelse(is.na(rasterAQ@data@values)=="FALSE", 1000000, 1)
  }
  if(saveRast=="TRUE"){
  rasterAQ <<- rasterAQ}
  else {rasterAQ <- rasterAQ}

  print("(2/7) Created raster from shapefile")

  # make sure all sites are "in water"
  vals <- data.frame(extract(x = rasterAQ, y = SpatialPoints(coords = cbind(alllong, alllat),
                                                          proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))))
  badvals <- data.frame(cbind(Long=alllong[vals>1], Lat=alllat[vals>1]))

  badvals <- join(badvals, AQsites, type="left")
  badvals <- join(badvals, rivercoords, type="left")

  if(length(badvals$Long) > 0){
    print(badvals)
  }

  if(length(badvals$Long) > 0)
  stop("Some points are on land. Use plot(rasterAQ) and click() to adjust placement of points
          in AQsites and rivercoords data.frames manually prior to running calcAQpress(). Check printed dataframe above for failed coordinates.")

  print("(3/7) All coordinates are in water")

  # make transitionlayers
  trAQ <- transition(1/rasterAQ, transitionFunction=mean, directions=8)
  trAQc <- geoCorrection(trAQ, type="c")

  dAQ <- costDistance(trAQc, fromCoords = cbind(as.numeric(AQsites$Long), as.numeric(AQsites$Lat)),
                      toCoords = cbind(as.numeric(rivercoords$Long), as.numeric(rivercoords$Lat)))

  print("(4/7) Calculated least-cost distances between AQ sites and specified rivers")

  dAQ <- as.data.frame(dAQ)
  colnames(dAQ) <- rivercoords$River
  dAQ <- cbind(dAQ, AQsites$Site.ID, as.character(AQsites$prov))
  colnames(dAQ)[dim(dAQ)[2] - 1] <- "Site.ID"
  colnames(dAQ)[dim(dAQ)[2]] <- "prov"

  if(distType=="distance"){

    distances <- melt(data=dAQ, id.vars = c("Site.ID","prov"), variable.name = "River")

    avgdistance <- ddply(.data=distances, .(River),
                         summarize,
                         avgdist = mean(value, na.rm=TRUE))

    save0 <- paste0(dir, "/avg_distance_final_", Sys.Date(), ".csv")

    write.csv(x = avgdistance, file=save0)
  }

  dist2AQ <- join(dAQ, inventory, type="left")
  dist2AQ <- subset(dist2AQ, is.na(totalfish)==FALSE)

  print("(5/7) Joined inventory data")

  dist2AQ[,1:length(rivercoords)] <- apply(dist2AQ[,1:length(rivercoords)], 2, function(x) as.character(x))
  dist2AQ[,1:length(rivercoords)] <- apply(dist2AQ[,1:length(rivercoords)], 2, function(x) as.numeric(x))

  dist2AQ[dist2AQ == 0] <- 1000

  # change totalfish to 1 (1 year)
  dist2AQ$totalfish <- 1

  # For each AQ site, caclulate total fish/distance from river
  fishdist <- c(dist2AQ$totalfish/dist2AQ[,1:length(rivercoords)])
  dist2AQtest <- cbind(Year=dist2AQ$Year, data.frame(fishdist))

  dist2AQmelt <- melt(data=dist2AQtest, id.vars = "Year", variable.name = "River")

  # For each year and river, calculate total propagule pressure
  prop.press <- ddply(.data=dist2AQmelt, .(Year, River),
                      summarize,
                      prop.press=sum(value, na.rm=T))

  print("(6/7) Calculated and saved propagule pressure by year")

  prop.press_final <- join(prop.press, rivercoords, type="left")

  save1 <- paste0(dir, "/prop.press_final_", Sys.Date(), ".csv")

  write.csv(x = prop.press_final, file = save1)

  prop.press_avg <- ddply(.data=prop.press_final, .(River, Lat, Long),
                           summarize,
                           averagepress = mean(prop.press),
                           nyears = length(unique(Year)))

  print("(7/7) Calculated and saved average propagule pressure for each river")

  save2 <- paste0(dir, "/prop.press_avg_", Sys.Date(), ".csv")

  write.csv(x = prop.press_avg, file = save2)

  }



