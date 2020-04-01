## ------------------------------------------------------------------------
library(sf)
library(raster)
library(DIFMR)
library(GWmodel)
library(reticulate)
np <- import("numpy")

set.seed(12345)

f_coords <- cbind(c(0, 400), c(0, 300))
pol <- st_as_sfc(st_bbox(st_multipoint(f_coords)))
field <- st_sf(pol, crs = 32616)
rst <- raster(field, res = 10)

trial_grd = rasterToPoints(rst, spatial = TRUE)
gridded(trial_grd) = TRUE
pts <- sp::coordinates(trial_grd)
dMat <- gw.dist(pts, pts)

bws = seq(5000,2000,-100)
bws = c(bws, seq(2000,1000,-50))
bws = c(bws, seq(1000,500,-25))
bws = c(bws, seq(500,100,-5))
bws = c(bws, seq(100,1,-1))


wwf = list()

i = 620
dist.vi <- dMat[, i]

i = 1
for(bw in bws){
  rst$W.i <- gw.weight(dist.vi, bw, kernel = "gaussian", adaptive = TRUE)
  wwf[[paste0('W', bw)]] = round(255 * rst$W.i)
}
wwr = stack(wwf)
wwm = as.array(wwr)

wfile = './data/weigths_zoom.npy'
np$save(wfile, wwm)


