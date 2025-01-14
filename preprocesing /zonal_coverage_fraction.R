# --------------------------------------------------- Make chunks of DEM --------------------------------------------------- #

# Load the required library
library(terra)

# Load the raster file
dem_path <- "my_path/COG_4326.tif"
raster_file <- rast(dem_path)
# Define the maximum number of pixels per chunk
max_pixels <- 2.5e+07

# Calculate the number of chunks needed along rows and columns
nrow_chunks <- ceiling(nrow(raster_file) * ncol(raster_file) / max_pixels)
chunk_height <- ceiling(nrow(raster_file) / sqrt(nrow_chunks))
chunk_width <- ceiling(ncol(raster_file) / sqrt(nrow_chunks))

# chunk, crop, and save
chunk_counter <- 1
for (i in seq(1, nrow(raster_file), by = chunk_height)) {
  for (j in seq(1, ncol(raster_file), by = chunk_width)) {
    # Calculate the extent for the current chunk
    xmin <- xFromCol(raster_file, j)
    xmax <- xFromCol(raster_file, min(j + chunk_width - 1, ncol(raster_file)))
    ymin <- yFromRow(raster_file, min(i + chunk_height - 1, nrow(raster_file)))
    ymax <- yFromRow(raster_file, i)
    chunk_extent <- ext(xmin, xmax, ymin, ymax)
    
    # Crop 
    chunk <- crop(raster_file, chunk_extent)
    output_file <- paste0("my_path/chunk_", chunk_counter, ".tif")
    
    # Write the chunk to a GeoTIFF file
    writeRaster(chunk, output_file, overwrite = TRUE)
    
    # Increment the chunk counter
    chunk_counter <- chunk_counter + 1
    print(paste0("Completed ", output_file))
  }
}
cat("Splitting completed and files saved.")


# --------------------------------------------------- Zonal coverage fraction --------------------------------------------------- #
# Load  libraries
library(terra)
library(sf)
library(zonal)
library(parallel)
library(arrow)
library(raster)


# DEM chunks
chunk_folder <- "my_path/chunks/"  
triangle_file <- "my_path/bary_triangles.gpkg"

# Output folder
output_folder <- "my_path/cf/"  
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# List all chunk files
chunk_files <- list.files(chunk_folder, pattern = "\\.tif$", full.names = TRUE)

process_chunk <- function(chunk_file, triangle_file, output_folder) {
  gdf <- st_read(triangle_file, layer="triangles")
  chunk_raster <- rast(chunk_file)
  
  # Skip processing if the chunk contains only NaN values
  if (all(is.na(values(chunk_raster)))) {
    cat("Skipping empty chunk:", chunk_file, "\n")
    return(NULL)
  }
  
  # Run zonal::weight_grid
  result <- zonal::weight_grid(chunk_raster, gdf, "pg_id")
  # Define output file name
  chunk_name <- gsub(".*/|\\.tif$", "", chunk_file)
  output_file <- file.path(output_folder, paste0(chunk_name, "_result_local_indexed.parquet"))
  write_parquet(result, output_file)
  # To change: cleanup memory
  rm(gdf, chunk_raster, result)
  cat("Processed and saved chunk:", chunk_file, "\n")
}

# Set up parallel processing
num_cores <- 4
cl <- makeCluster(num_cores)
clusterExport(cl, c("chunk_files", "triangle_file", "output_folder", "process_chunk"))
clusterEvalQ(cl, library(terra))
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(zonal))
clusterEvalQ(cl, library(arrow))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(parallel))

# Run the function on all chunk files in parallel
parLapply(cl, chunk_files, function(chunk_file) {
  process_chunk(chunk_file, dem_path, triangle_file, output_folder)
})

# Stop the cluster
stopCluster(cl)
cat("Processing completed, and results have been saved in", output_folder, "\n")
