# Multivariate and Spatiotemporal Geostatistics {#sec-stgeostatistics}

library(gstat)
library(sp)
library(sf)

data(meuse)
coordinates(meuse) = ~x+y

data(meuse.grid)
gridded(meuse.grid) = ~x+y

data("meuse.riv")
meuse.sr = SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)),"meuse.riv")))

spplot(meuse[,"zinc"], main="Zinc concentrations", 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))
spplot(meuse[,"lead"], main="Lead concentrations", 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))



# cokriging of the four heavy metal variables
meuse.g <- gstat(id="zn", formula=log(zinc)~1, data=meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cu", log(copper)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cd", log(cadmium)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "pb", log(lead)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.g, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.g)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid)

library(lattice)
pl1 <- spplot(z["zn.pred"], main="log-zinc predictions")
pl2 <- spplot(z["cu.pred"], main="log-copper predictions")
pl3 <- spplot(z["cd.pred"], main="log-cadmium predictions")
pl4 <- spplot(z["pb.pred"], main="log-lead predictions")
print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))

load("/Users/iprincetech/Downloads/2024-12-05 - Multivariate Geostatistics-20241205/lecture_08_data.RData")

sfc <- st_geometry(a2.sf)[match(colnames(aqsel),
                                a2.sf$station_european_code)] |>
  st_transform(crs)

library(xts)
library(stars)
st_as_stars(NO2 = as.matrix(aqsel)) |>
  st_set_dimensions(names = c("time", "station")) |>
  st_set_dimensions("time", index(aqsel)) |>
  st_set_dimensions("station", sfc) -> no2.st
no2.st

# From this, we can compute the spatiotemporal variogram using

load(file = "/Users/iprincetech/Downloads/2024-12-05 - Multivariate Geostatistics-20241205/vst.RData")
#library(gstat)
#v.st <- variogramST(NO2~1, no2.st[,1:(24*31)], tlags = 0:48, 
#					cores = getOption("mc.cores", 2))
#save(list = "v.st", file = "vst.RData")

v1 <- plot(v.st)
v2 <- plot(v.st, map = FALSE, legend = list())
print(v1, split = c(1,1,2,1), more = TRUE)
print(v2, split = c(2,1,2,1), more = FALSE)

plot(v.st, wireframe=TRUE)

# To this sample variogram, we can fit a variogram model. One relatively
# flexible model we try here is the product-sum model [@RJ-2016-014], fitted by

# product-sum
prodSumModel <- vgmST("productSum",
                      space = vgm(150, "Exp", 200000, 0),
                      time = vgm(20, "Sph", 6, 0),
                      k = 2)


#v.st$dist = v.st$dist / 1000
StAni <- estiStAni(v.st, c(0,200000))
(fitProdSumModel <- fit.StVariogram(v.st, prodSumModel,
                                    fit.method = 7, stAni = StAni, method = "L-BFGS-B",
                                    control = list(parscale = c(1,100000,1,1,0.1,1,10)),
                                    lower = rep(0.0001, 7)))



plot(v.st, fitProdSumModel, wireframe = FALSE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 150))
plot(v.st, model = fitProdSumModel, wireframe = TRUE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 195))


set.seed(1331)
pt <- st_sample(de, 2)
t <- st_get_dimension_values(no2.st, 1)
st_as_stars(list(pts = matrix(1, length(t), length(pt)))) |>
  st_set_dimensions(names = c("time", "station")) |>
  st_set_dimensions("time", t) |>
  st_set_dimensions("station", pt) -> new_pt

load("/Users/iprincetech/Downloads/2024-12-05 - Multivariate Geostatistics-20241205/new_ts.RData")
#no2.st <- st_transform(no2.st, crs)
#new_ts <- krigeST(NO2~1, data = no2.st["NO2"], newdata = new_pt,
#				  nmax = 50, stAni = StAni, modelList = fitProdSumModel,
#				  progress = FALSE)
#save(list = "new_ts", file = "data/new_ts.RData")

plot(as.xts(new_ts[2]))

st_bbox(de) |>
  st_as_stars(dx = 10000) |>
  st_crop(de) -> grd
d <- dim(grd)
t4 <- t[(1:4 - 0.5) * (3*24*30)]
st_as_stars(pts = array(1, c(d[1], d[2], time = length(t4)))) |>
  st_set_dimensions("time", t4) |>
  st_set_dimensions("x", st_get_dimension_values(grd, "x")) |>
  st_set_dimensions("y", st_get_dimension_values(grd, "y")) |>
  st_set_crs(crs) -> grd.st

load("/Users/iprincetech/Downloads/2024-12-05 - Multivariate Geostatistics-20241205/new_int.RData")
# new_int <- krigeST(NO2~1, data = no2.st["NO2"], newdata = grd.st,
# 				   nmax = 200, stAni = StAni, modelList = fitProdSumModel,
# 				   progress = FALSE)
# names(new_int)[2] = "NO2"

# save(list = "new_int", file = "data/new_int.RData")

library(viridis)
library(viridisLite)
library(ggplot2)
g <- ggplot() + coord_equal() +
  scale_fill_viridis() +
  theme_void() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0))
g + geom_stars(data = new_int, aes(fill = NO2, x = x, y = y)) + 
  facet_wrap(~as.Date(time), nrow = 1) +
  geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
  geom_sf(data = no2.sf, col = 'grey', cex = .5) + 
  coord_sf(lims_method = "geometry_bbox")

### Irregular space time data

library(sp)
library(spacetime)
library(gstat)
library(lattice)

# create n space-time points over [0,1] x [0,1] x [Now, Now+some days]
t0 = Sys.time() # now
n = 1000
set.seed(13131) # fix outcomes
x = runif(n)
y = runif(n)
t = t0 + 1e6 * runif(n)
z = rnorm(n)
stidf = STIDF(SpatialPoints(cbind(x,y)), sort(t), data.frame(z=z))

stplot(stidf, number=21, main="random spatio-temporal noise")

library(sftime)
sft = st_as_sftime(stidf)

# create a regular 20 x 20 x 10 grid of prediction locations:
grd = as(SpatialGrid(GridTopology(c(0.025,0.025), c(.05, .05), c(20,20))), "SpatialPixels")
tgrd = seq(min(t)+10000, max(t)-10000, length.out = 10)

stf = STF(grd, tgrd)
#stf = STFDF(grd, tgrd, data.frame(x=rep(0,400*10)))

library(stars)
st = st_as_stars(stf)

# define a variogram model
sumMetricModel <- vgmST("sumMetric",
                        space=vgm(1/6, "Sph", 0.25, 1/60),
                        time =vgm(2/6, "Exp",  1e5, 1/60),
                        joint=vgm(0.4, "Exp", 0.3, 0.1),
                        stAni=1/1e6)
attr(sumMetricModel, "temporal unit") <- "secs"

dg <- data.frame(spacelag=rep(c(0.001,1:10)/10,6), 
                 timelag=rep(0:5*50e3, each=11))
#wireframe(model~spacelag+timelag,
#          variogramSurface(sumMetricModel, dist_grid = dg),
#          scales=list(arrows=F),
#          drape=T, col.regions=bpy.colors(),
#          zlim=c(0,1.2),
#          main="imposed sum-metric model")

locKrig_sft <- krigeST(z~1, sft, st, sumMetricModel, nmax=20, computeVar = T)
locKrig <- krigeST(z~1, stidf, stf, sumMetricModel, nmax=20, computeVar = T)
stplot(locKrig[,,"var1.pred"], col.regions=bpy.colors(), scales=list(draw=T))
plot(locKrig_sft[1], col = sf.colors(), breaks = "equal")
stplot(locKrig[,,"var1.var"], col.regions=bpy.colors(), scales=list(draw=T))
plot(locKrig_sft[2], col = sf.colors(), breaks = "equal")

st$foo = 0
st_as_sf(st, long = TRUE) |> st_as_sftime() -> st.sftime
locKrig_sft <- krigeST(z~1, sft, st.sftime, sumMetricModel, nmax=20, computeVar = T)
plot(locKrig_sft["var1.pred"])




# Solutions Exercise 1, lecture 8 - MV and spatiotemporal geostats
library(gstat)
library(sp)
data(meuse)
coordinates(meuse) = ~x+y

data(meuse.grid)
gridded(meuse.grid) = ~x+y

data("meuse.riv")
meuse.sr = SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)),"meuse.riv")))

spplot(meuse[,"zinc"], main="Zinc concentrations", 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))
spplot(meuse[,"lead"], main="Lead concentrations", 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))

# cokriging of the four heavy metal variables
meuse.g <- gstat(id="zn", formula=log(zinc)~1, data=meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cu", log(copper)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cd", log(cadmium)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "pb", log(lead)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.g, cutoff=1000)
str(x,2)

plot(x)

meuse.fit = fit.lmc(x, meuse.g)

str(meuse.fit)
plot(x, model = meuse.fit)

z <- predict(meuse.fit, newdata = meuse.grid)

library(lattice)
pl1 <- spplot(z["zn.pred"], main="log-zinc predictions")
pl2 <- spplot(z["cu.pred"], main="log-copper predictions")
pl3 <- spplot(z["cd.pred"], main="log-cadmium predictions")
pl4 <- spplot(z["pb.pred"], main="log-lead predictions")
print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))


## # zinc separate
v_zinc <- variogram(log(zinc)~1, meuse)
fitv_zinc <- fit.variogram(v_zinc, vgm(1, "Sph", 900, 1))
plot(v_zinc, model=fitv_zinc)

pred_zinc <- krige(log(zinc)~1, meuse, meuse.grid, fitv_zinc)

z$zn.sep <- pred_zinc$var1.pred
z$zn.diff <- z$zn.sep - z$zn.pred

spplot(z[,c("zn.pred", "zn.sep")], 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))

spplot(z[,"zn.diff"], 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))

z$zn.pred.exp <- exp(z$zn.pred)
z$zn.sep.exp <- exp(z$zn.sep)
z$zn.diff.exp <- z$zn.sep.exp - z$zn.pred.exp

spplot(z[,c("zn.pred.exp", "zn.sep.exp")], 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))
spplot(z[,"zn.diff.exp"], 
       sp.layout=list(list("sp.polygons", meuse.sr, col="darkblue")))


# Exercise 2

library(stars)
library(gstat)

no2.daily  = aggregate(no2.st, "1 day", mean, na.rm = TRUE)

v.daily <- variogramST(NO2 ~ 1, no2.daily)

plot(v.daily, main = "Daily Mean Spatiotemporal Variogram")
plot(v.st, main = "Hourly Spatiotemporal Variogram")


# Exercise 3: Fit different spatiotemporal variogram models and choose the “best fit”

# Estimate spatiotemporal anisotropy
StAni <- estiStAni(v.st, c(0, 200000))

# Product-Sum Model 
prodSumModel <- vgmST("productSum",
                      space = vgm(150, "Exp", 200000, 0),
                      time = vgm(20, "Sph", 6, 0),
                      k = 2)

fitProdSumModel <- fit.StVariogram(v.st, prodSumModel,
                                   fit.method = 7, stAni = StAni, method = "L-BFGS-B",
                                   control = list(parscale = c(1, 100000, 1, 1, 0.1, 1, 10)),
                                   lower = rep(0.0001, 7))

# Plot Product-Sum Model
plot(v.st, fitProdSumModel, wireframe = FALSE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 150))
plot(v.st, model = fitProdSumModel, wireframe = TRUE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 195))


# Separable Model 
separableModel <- vgmST("separable",
                        space = vgm(150, "Exp", 200000, 0),
                        time = vgm(20, "Sph", 6, 0),
                        sill = 200)

fitSeparableModel <- fit.StVariogram(v.st, separableModel,
                                     fit.method = 7, stAni = StAni, method = "L-BFGS-B",
                                     control = list(parscale = c(1, 100000, 1, 1, 100)),
                                     lower = rep(0.0001, 5))

# Plot Separable Model
plot(v.st, fitSeparableModel, wireframe = FALSE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 150))
plot(v.st, model = fitSeparableModel, wireframe = TRUE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 195))


# ---- Metric Model ----
jointVgm <- vgm(150, "Exp", 200000, 0)  # Joint variogram

metricModel <- vgmST("metric",
                     joint = jointVgm,
                     stAni = StAni)

fitMetricModel <- fit.StVariogram(v.st, metricModel,
                                  fit.method = 7, stAni = StAni, method = "L-BFGS-B",
                                  control = list(parscale = c(1, 100000, 1, 100)),
                                  lower = rep(0.0001, 4))

# Plot Metric Model
plot(v.st, fitMetricModel, wireframe = FALSE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 150))
plot(v.st, model = fitMetricModel, wireframe = TRUE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 195))


#  Sum-Metric Model
spatialVgm <- vgm(150, "Exp", 200000, 0)
temporalVgm <- vgm(20, "Sph", 6, 0)

sumMetricModel <- vgmST("sumMetric",
                        space = spatialVgm,
                        time = temporalVgm,
                        joint = jointVgm,
                        stAni = StAni)

fitSumMetricModel <- fit.StVariogram(v.st, sumMetricModel,
                                     fit.method = 7, stAni = StAni, method = "L-BFGS-B",
                                     control = list(parscale = c(1, 100000, 1, 1, 0.1, 1, 1, 100000, 1, 1)),
                                     lower = rep(0.0001, 10))

# Plot Sum-Metric Model
plot(v.st, fitSumMetricModel, wireframe = FALSE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 150))
plot(v.st, model = fitSumMetricModel, wireframe = TRUE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 195))


set.seed(1331)
pt <- st_sample(de, 2)
t <- st_get_dimension_values(no2.st, 1) 
st_as_stars(list(pts = matrix(1, length(t), length(pt)))) |>
  st_set_dimensions(names = c("time", "station")) |>
  st_set_dimensions("time", t) |>
  st_set_dimensions("station", pt) -> new_pt


## Kriging for Product-Sum Model
new_ts_prodSum <- krigeST(NO2~1, data = no2.st["NO2"], newdata = new_pt,
                          nmax = 50, stAni = StAni, modelList = fitProdSumModel,
                          progress = TRUE)

## Separable Model(Kriging)
new_ts_separable <- krigeST(NO2~1, data = no2.st["NO2"], newdata = new_pt,
                            nmax = 50, stAni = StAni, modelList = fitSeparableModel,
                            progress = TRUE)

## Metric Model(Kriging)
new_ts_metric <- krigeST(NO2~1, data = no2.st["NO2"], newdata = new_pt,
                         nmax = 50, stAni = StAni, modelList = fitMetricModel,
                         progress = TRUE)

# Kriging for Sum-Metric Model
new_ts_sumMetric <- krigeST(NO2~1, data = no2.st["NO2"], newdata = new_pt,
                            nmax = 50, stAni = StAni, modelList = fitSumMetricModel,
                            progress = TRUE)



# Remove rows with missing values
v.st_clean <- na.omit(v.st)

# Or impute missing values (e.g., using the mean)
v.st$dist[is.na(v.st$dist)] <- mean(v.st$dist, na.rm = TRUE)
v.st$gamma[is.na(v.st$gamma)] <- mean(v.st$gamma, na.rm = TRUE)

# Convert `timelag` to numeric (hours)
v.st$timelag_numeric <- as.numeric(v.st$timelag, units = "hours")

# If you have spatial data in `sf` or `sp` format
if (inherits(v.st, "sf")) {
  spatial_coords <- st_coordinates(v.st)
  print(head(spatial_coords))
} else {
  print("Spatial data not in 'sf' format.")
}


#Exercise 4:Carry out a spatiotemporal interpolation for daily mean values for the days corresponding to those in the lecture, and compare the results.


library(viridis)
library(viridisLite)
library(ggplot2)

#setting up the Product-Sum model for spatiotemporal variogram modeling
prodSumModel <- vgmST("productSum",
                      space = vgm(50, "Exp", 200, 0),
                      time = vgm(20, "Sph", 40, 0),
                      k = 2)
# spatiotemporal variogram for daily mean concentrations have already been computed in exercise 2 (v.daily)

v.daily <- variogramST(NO2 ~ 1, no2.daily)

#estimating spatiotemporal anisotropy based on v.daily (already computed in exercise 2)
StAni = estiStAni(v.daily, c(0,20000))

#fitting the prodSumModel based on daily mean concentration (v.daily)
(fitProdSumModel <- fit.StVariogram(v.daily, prodSumModel, fit.method = 7,
                                    stAni = StAni, method = "L-BFGS-B",
                                    control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                    lower = rep(0.0001, 7)))

plot(v.daily, fitProdSumModel, wireframe = FALSE, all = TRUE, scales = list(arrows=FALSE), zlim = c(0,50))

plot(v.daily, model = fitProdSumModel, wireframe = TRUE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 195))

#spatiotemporal kriging using the fitted variogram model (fitProdSumModel_vdaily)

SpTempKr <- krigeST(NO2~1, data = no2.daily["NO2"], newdata = grd.st,
                    nmax = 100, stAni = StAni, modelList = fitProdSumModel,
                    progress = TRUE)
names(SpTempKr)[2] = "NO2"

#Visualization:
#spatiotemporal interpolation results (SpTempKr)

g <- ggplot() + coord_equal() +
  scale_fill_viridis() +
  theme_void() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0))
g + geom_stars(data = SpTempKr, aes(fill = NO2, x = x, y = y)) +
  facet_wrap(~as.Date(time)) +
  geom_sf(data = st_cast(de, "MULTILINESTRING")) +
  geom_sf(data = no2.sf, col = 'grey', cex = .5) +
  coord_sf(lims_method = "geometry_bbox")




