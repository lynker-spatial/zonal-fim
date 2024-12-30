# coastal-fim
A repository to compute SCHISM derived depth at 30m grids


## Progress 
1- implemented barycentric computations
2- implemented batching DEM and zonal (coverage fraction) parallel execution for large domains
3- next: write algorithm to re-index cell ids from batches to a global index matching original DEM 
4- next: Write tests for package


## Report
https://docs.google.com/document/d/1DoPeE0IRVHkjqabqTUaX5aWCnPZn9Mdv/edit?usp=sharing&ouid=110666552849114372265&rtpof=true&sd=true

## Test
1- Executed the entire process on Tampa: Pass
2- Executed everything other than coverage fraction interpolation on entire AtlGolf domain: Pass
3- Compare with Linear interpolation: Pass  