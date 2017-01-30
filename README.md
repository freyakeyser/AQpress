# AQpress
## Propagule pressure calculation

### Installation:
devtools::install_github("freyakeyser/AQpress")

### Input: 
AQsites: a dataframe of aquaculture sites, with Site.ID, Lat, and Long (decimal degrees)
    
rivercoords: a dataframe of river coordinates, with River, Lat, Long (decimal degrees), and prov (province or state). Make sure that the coordinates in this data frame are at the river mouth. 
    
inventory: a dataframe of aquaculture inventory, with Site.ID, Year, Lat, Long (decimal degrees), and approximate number of fish on site each year.
    
### Output:
~/prop.press_final_YYYY-MM-DD.csv: Propagule pressures by river and year

~/prop.press_avg_YYYY-MM-DD.csv: Average propagule pressure by river for all years
