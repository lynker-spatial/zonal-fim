# Load required libraries using pacman
pacman::p_load(
  AOI,        # For handling areas of interest
  sf,         # For spatial vector data manipulation
  terra,      # For handling raster data
  raster,     # Needed for fasterize
  zonal,      # For zonal statistics
  dplyr,      # For data manipulation
  fasterize,  # For fast rasterization of vector data
  arrow,      # For working with Apache Parquet files
  glue        # For string interpolation
)

# Define file paths:
# > In this directory, you need the CoastalBathy Grid, the bary_triangles.gpkg
base <- '/Users/mikejohnson/Downloads/coastal_data'


gpkg <- glue('{base}/ElementPolygons.gpkg')
weights_file <- glue('{base}/ElementPolygons.parquet')
d <- rast(glue('{base}/TBDEM_AtlanticGulf_Mosaic_NWM_3_Revised_v4_COG_4326.tif'))
mask_file <- glue("{base}/mask.tif")

# Read and reproject vector data
gpkg_layer <- 'ElementPolygons'
pg5070 <- st_transform(read_sf(gpkg, gpkg_layer), 5070)

# Create a grid based on the bounding box of the input geometry
tess <- bbox_get(pg5070) |> 
  st_make_grid() |> 
  st_as_sf()

# Compute the convex hull of the input geometry
hull <- st_as_sf(st_union(pg5070))

# Filter grid cells intersecting the hull
tess_small <- st_filter(tess, hull)

# Initialize a list to store weighted grid results
ll <- list()

# Iterate through each grid cell and compute weights
for (i in 1:nrow(tess_small)) {
  tmp <- st_filter(pg5070, tess_small[i, ])
  
  message("There are ", 
          nrow(tmp), 
          " features in tesselation ", 
          i, 
          "/", 
          nrow(tess_small))
  
  ll[[i]] <- weight_grid(d, geom = tmp, ID = 'pg_id')
}

# Masking out the ocean
if(!file.exists(mask_file)) {
  mask <- st_zm(read_sf(gpkg, "mask"))
  f <- fasterize(mask, raster(d))
  writeRaster(f, mask_file, overwrite = TRUE)
}

# Combine elevation data with mask values
tab <- data.frame(
  elevation = terra::values(d),
  mask = values(rast(mask_file)),
  cell = 1:ncell(d)
) |> 
  setNames(c("elevation", 'mask', "cell"))

# Combine all grid results and join with elevation/mask data
fin <-  left_join(bind_rows(ll), tab, by = 'cell')

# Preview the final dataset
head(fin)

# Save the final dataset to a Parquet file
unlink(weights_file)
write_parquet(fin, weights_file)