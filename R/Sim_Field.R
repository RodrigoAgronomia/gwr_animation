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
trial <- aggregate(rst, c(4, 1))
trial$id <- 1:ncell(trial)
trial$pcol = 1:ncol(trial)
trial$prow = rep(1:nrow(trial), each = ncol(trial))

pts = rasterToPoints(rst, spatial = TRUE)
pts = st_as_sf(pts)
cnt = st_centroid(st_geometry(field))

pols <- rasterToPolygons(trial)
pols <- st_as_sf(pols)
pols <- get_block_ids(pols)

Treat_rates <- c(0, 50, 100, 150, 200)
pols_treat <- design_treat_graeco(pols, length(Treat_rates))
pols_treat$Treat <- Treat_rates[pols_treat$F1]

## ------------------------------------------------------------------------
trial$Treat <- pols_treat$Treat[pols_treat$id]
plot(trial$Treat)

rst$Treat <- resample(trial$Treat, rst, method = "ngb")
plot(rst$Treat)

## ------------------------------------------------------------------------
trial_grd <- rasterToPoints(rst, spatial = TRUE)
gridded(trial_grd) <- TRUE

m <- gstat::vgm(psill = 1, model = "Gau",
                range = 50,
                nugget = 0)

set.seed(12345)
g.dummy <- gstat::gstat( 
  formula = z ~ 1,
  dummy = TRUE, beta = 0,
  model = m, nmax = 10
)
rst_sim <- predict(g.dummy, trial_grd, nsim = 1)
rst$Yield_Ref <- scale(stack(rst_sim))
rst$Yield_Ref = 1e4 + 1e3 * rst$Yield_Ref
plot(rst$Yield_Ref)

m <- gstat::vgm(psill = 1, model = "Gau",
                range = 100,
                nugget = 0)
set.seed(12354)
g.dummy <- gstat::gstat(
  formula = z ~ 1,
  dummy = TRUE, beta = 0,
  model = m, nmax = 10
)
rst_sim <- predict(g.dummy, trial_grd, nsim = 1)
rst$NR_gwr = -5 * rst_sim$sim1

rst$Treat_Yield = rst$NR_gwr * rst$Treat
rst$Yield_Obs = rst$Yield_Ref + rst$Treat_Yield

treatm = as.array(rst)

treat_file = './data/treat.npy'
np$save(treat_file, treatm)



trial_grd = rasterToPoints(rst, spatial = TRUE)
gridded(trial_grd) = TRUE
pts <- sp::coordinates(trial_grd)
dMat <- gw.dist(pts, pts)
gwr.fml <- as.formula(Yield_Obs ~ poly(Treat, 1))

bw = 40
gwr.model <- gwr.basic(gwr.fml, trial_grd,
                       bw = bw, dMat = dMat,
                       kernel = "gaussian"
)
print(gwr.model)


wwf = list()

i = 1
for(i in 1:nrow(pts)){
  print(i)
  dist.vi <- dMat[, i]
  rst$W.i <- gw.weight(dist.vi, bw, kernel = "gaussian", adaptive = TRUE)
  wwf[[paste0('W', i)]] = round(255 * rst$W.i)
}
wwr = stack(wwf)
wwm = as.array(wwr)

wfile = './data/weigths.npy'
np$save(wfile, wwm)


# 
# 
# 
# gwr_r <- gwr.model$SDF
# gridded(gwr_r) = TRUE
# gwr_rst = stack(gwr_r)
# trial_grd$NR_est = gwr_rst$poly.Treat..1.[]
# 
# lm0 = lm(NR_gwr ~ NR_est, trial_grd)
# summary(lm0)
# trial_grd$NR_est = predict(lm0)
# 
# plot(trial_grd$NR_est, trial_grd$NR_gwr)
# 
# plot(rst$NR_gwr)
# plot(gwr_rst$Treat)



