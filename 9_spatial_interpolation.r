#### 9.1 ####
# 図9.1
# (a) ガウス関数
plot(density(0, kernel="gaussian", bw=1), main="", xlab="", ylab="",
lwd=8, cex.axis=1.8)
# (b) イパネクニコフ関数
plot(density(0, kernel="epanechnikov", bw=1), main="", xlab="", ylab="",
lwd=8, cex.axis=1.8)
# (c) 四次関数
plot(density(0, kernel="biweight", bw=1), main="", xlab="", ylab="",
lwd=8, cex.axis=1.8)

# 図9.2 ssの数字をそれぞれ入れ替えて図を作成する((a) バンド幅=1)
ss <-1

gau.all <- function(x){
1/sqrt(2*pi*ss)*exp(-x^2/(2*ss^2)) +
1/sqrt(2*pi*ss)*exp(-(x+4)^2/(2*ss^2)) +
1/sqrt(2*pi*ss)*exp(-(x-2)^2/(2*ss^2)) +
1/sqrt(2*pi*ss)*exp(-(x-5)^2/(2*ss^2)) +
1/sqrt(2*pi*ss)*exp(-(x+1.5)^2/(2*ss^2))}
curve(gau.all, xlim=c(-10, 10), ylim=c(0, 0.8), main="", xlab="", ylab="", lwd=8, cex.axis=1.8)

gau1 <-  function(x){1/sqrt(2*pi*ss)*exp(-x^2/(2*ss^2))}
gau2 <-  function(x){1/sqrt(2*pi*ss)*exp(-(x+4)^2/(2*ss^2))}
gau3 <-  function(x){1/sqrt(2*pi*ss)*exp(-(x-2)^2/(2*ss^2))}
gau4 <-  function(x){1/sqrt(2*pi*ss)*exp(-(x-5)^2/(2*ss^2))}
gau5 <-  function(x){1/sqrt(2*pi*ss)*exp(-(x+1.5)^2/(2*ss^2))}

curve(gau1, lwd=4, add=TRUE)
curve(gau2, lwd=4, add=TRUE)
curve(gau3, lwd=4, add=TRUE)
curve(gau4, lwd=4, add=TRUE)
curve(gau5, lwd=4, add=TRUE)

# ２次元平面空間での密度関数
x <- runif(50)*10
y <- runif(50)*10
xy <- cbind(x, y)
poly0 <- cbind(c(0,10,10,0), c(0,0,10,10))
# 図9.3(a)
image(kernel2d(xy, poly0, h0=1, nx=100, ny=100), col=gray((0:20)/20),
cex.axis=1.8)
# 図9.3(b)
image(kernel2d(xy, poly0, h0=0.5, nx=100, ny=100), col=gray((0:20)/20),
cex.axis=1.5)
# 図9.3(c)
image(kernel2d(xy, poly0, h0=0.7, nx=100, ny=100), col=gray((0:20)/20),
cex.axis=1.5)
# 最小二乗誤差法によるバンド幅の決定
library(splancs)
poly0 <- cbind(c(0,1,1,0), c(0,0,1,1))
Mse2d <- mse2d(as.points(X), poly0, nsmse=50, range=1)
# 図9.5
plot(Mse2d$h[5:50],Mse2d$mse[5:50], type="l", ylab="MSE", xlab="バンド幅", lwd=6, cex.axis=1.3, cex.lab=1.2)
points(Mse2d$h[which.min(Mse2d$mse)], Mse2d$mse[which.min(Mse2d$mse)], pch=1, cex=3, lwd=2)
# 交差検証対数尤度関数によるバンド幅の決定
library(spatialkernel)
x <- runif(600)
y <- runif(600)
mks <- sample(c("a","b", "c"), 600, replace=TRUE)
pts <- cbind(x, y)
h <- seq(0.01, 1, by=0.01)
cv <- cvloglk(pts, mks, h=h)$cv
# 図9.6
plot(h, cv, type="l", ylab="交差検証対数尤度", xlab="バンド幅", lwd=6, cex.axis=1.3, cex.lab=1.2)
points(h[which.max(cv)], cv[which.max(cv)], pch=1, cex=3, lwd=2)

#### 9.2～　データの作成と表示 ####
# データの読み込み
spm.shp <- readShapePoints("tma_spm.shp")
spm <- cbind(spm.shp$ID, spm.shp$X, spm.shp$Y, spm.shp$SPM07)
colnames(spm) <- c("ID", "X", "Y", "SPM07")
spm <- as.data.frame(spm)
ward.shp <- readShapePoly("Ward.shp")
mesh.grid <- read.table("mesh.csv", header=TRUE, sep=",")
coordinates(mesh.grid) <- c("X", "Y")
mesh.grid <- as(mesh.grid, "SpatialPixelsDataFrame")

#### 9.2 空間内挿手法：逆距離加重法(IDW) ####
spm.idw1 <- idw(SPM07*1000~1, locations=~X+Y, data=spm, mesh.grid, idp=2)
spplot(spm.idw1["var1.pred"])

#### 9.3 ヴァリオグラム ####
## 探索的ヴァリオグラム分析 [#md53f8a0]
# coordinatesの設定
coordinates(spm) <- c("X","Y")
# 等方性モデル
# 定数項モデル
spm.var0 <- variogram(SPM07*1000~1, data=spm)
plot(spm.var0)
# 緯度経度によるトレンド
spm.var1 <- variogram(SPM07*1000~X+Y, data=spm)
# 図9.9
var.cld <- variogram(SPM07*1000~X+Y, data=spm, cloud=TRUE)
plot(var.cld$dist, var.cld$gamma, pch=19, cex=1.2, xlab="distance", ylab="gamma", cex.axis=1.3, cex.lab=1.2)
# 図9.10
plot(spm.var1$dist, spm.var1$gamma, pch=1, lwd=2, cex=1.5, ylim=c(0, 80), xlab="distance", ylab="gamma", cex.axis=1.3, cex.lab=1.2)
## ヴァリオグラム・モデル
# 図9.12(1)　指数モデル
plot(variogramLine(vgm(psill=25, model="Exp", range=28000, nugget=45), 70000, 100), type="l", cex=1.5, lwd=4, ylab="semivariance", xlab="distance", cex.axis=1.2, cex.lab=1.2, ylim=c(40,80))
points(spm.var1$dist, spm.var1$gamma, cex=1.5, lwd=2)
spm.model1 <- vgm(psill=25, model="Exp", range=28000, nugget=45)
plot(spm.var1, spm.model1, cex=1.5, lwd=4)
# 図9.12(2)　球形モデル
plot(variogramLine(vgm(psill=25, model="Sph", range=60000, nugget=45), 70000, 100), type="l", cex=1.5, lwd=4, ylab="semivariance", xlab="distance", cex.axis=1.2, cex.lab=1.2, ylim=c(40,80))
points(spm.var1$dist, spm.var1$gamma, cex=1.5, lwd=2)
# 図9.12(3)　線形モデル
plot(variogramLine(vgm(psill=25, model="Lin", range=56000, nugget=45), 70000, 100), type="l", cex=1.5, lwd=4, ylab="semivariance", xlab="distance", cex.axis=1.2, cex.lab=1.2, ylim=c(40,80))
points(spm.var1$dist, spm.var1$gamma, cex=1.5, lwd=2)
# 図9.12(4)　ガウスモデル
plot(variogramLine(vgm(psill=20, model="Gau", range=35000, nugget=50), 70000, 100), type="l", cex=1.5, lwd=4, ylab="semivariance", xlab="distance", cex.axis=1.2, cex.lab=1.2, ylim=c(40,80))
points(spm.var1$dist, spm.var1$gamma, cex=1.5, lwd=2)
# 図9.12(5)　ナゲット効果モデル
plot(variogramLine(vgm(psill=0, model="Nug", nugget=70), 70000, 100), type="l", cex=1.5, lwd=4, ylab="semivariance", xlab="distance", cex.axis=1.2, cex.lab=1.2, ylim=c(40,80))
points(spm.var1$dist, spm.var1$gamma, cex=1.5, lwd=2)
# 図9.12(6)　Maternモデル
plot(variogramLine(vgm(psill=25, model="Mat", range=30000, nugget=45), 70000, 100), type="l", cex=1.5, lwd=4, ylab="semivariance", xlab="distance", cex.axis=1.2, cex.lab=1.2, ylim=c(40,80))
points(spm.var1$dist, spm.var1$gamma, cex=1.5, lwd=2)

## ヴァリオグラム・モデルのfit：推定方法による違い
# 球形モデル
spm.model2 <- vgm(psill=25, model="Sph", range=60000, nugget=45)
spm.fit <- fit.variogram(spm.var1, spm.model2)
plot(spm.var1, spm.model2, cex=1.5, lwd=4)
# 球形モデル：WLS
fit.variogram(spm.var1, spm.model2, fit.method=7)
# 球形モデル：OLS
fit.variogram(spm.var1, spm.model2, fit.method=6)
# 制限付き最尤法
fit.variogram.reml(SPM07*1000~X+Y,data=spm, model=vgm(25, "Sph", 60000, 45))
# 図9.13
# WLS
plot(variogramLine(vgm(psill=25, model="Sph", range=60000, nugget=45), 70000, 100), type="l", cex=1.5, lwd=4, ylab="semivariance", xlab="distance", cex.axis=1.2, cex.lab=1.2, ylim=c(40,80))
# OLS
lines(variogramLine(vgm(psill=25.15494, model="Sph", range=93784.44, nugget=48.34058), 100000, 100), cex=1.5, lwd=4,lty=2)
# REML
lines(variogramLine(vgm(psill=28.94858, model="Sph", range=60000, nugget=51.51024), 100000, 100), cex=1.5, lwd=4,lty=3)
# バリオグラム
points(spm.var1$dist, spm.var1$gamma, cex=1.5, lwd=2)
# 凡例
legend("bottomright", legend=c("WLS", "OLS", "REML"), lty=c(1,2,3), lwd=c(4,4,4), cex=1.5)

## 9.3.2 異方性モデル
spm.var2 <- variogram(SPM07*1000~X+Y, data=spm, alpha=0:3*90)
plot(spm.var2)
spm.anis1 <- vgm(psill=25, model="Gau", range=35000, nugget=50, anis=c(0, 0.8))
plot(spm.var2, spm.anis1)

#### 9.4 クリギング ####
# 単純クリギング
spm.gs <- gstat(id="ID", formula=SPM07*1000~1, data=spm, model=spm.model1, beta=mean(spm$SPM07*1000))
spm.ps <- predict(spm.gs, mesh.grid)
spplot(spm.ps[1])
spplot(spm.ps[2])
# 通常クリギング
spm.go <- gstat(id="ID", formula=SPM07*1000~1, data=spm, model=spm.model1)
spm.po <- predict(spm.go, mesh.grid)
spplot(spm.po[1])
spplot(spm.po[2])
# 普遍クリギング
spm.gu <- gstat(id="ID", formula=SPM07*1000~X+Y, data=spm, model=spm.model1)
spm.pu <- predict(spm.gu, mesh.grid)
# 図9.17
spplot(spm.pu[1])
spplot(spm.pu[2])
# ブロック・クリギング
spm.bl <- krige(SPM07*1000~1, data=spm, mesh.grid, model=spm.model1, block=c(50,50))
spplot(spm.bl[1])
# シミュレーション：Sequential Gaussian simulation
spm.sim <- krige(SPM07*1000~1, data=spm, mesh.grid, spm.fit, nsim=6, nmax=100)
spplot(spm.sim)
# シミュレーション：Sequential Indicator simulation
spm.simI <- krige(SPM07*1000~1, data=spm, mesh.grid, spm.fit, indicators = TRUE ,nsim=6, nmax=100)
spplot(spm.simI)
# ベイジアン・クリギング

## クロス・ヴァリデーション
# 残差
mod.set <- sample(1:402,200)
spm.mod <- spm[mod.set,]
spm.val <- spm[-mod.set,]
spm.mod.var <- variogram(SPM07*1000~X+Y, data=spm.mod)
plot(spm.mod.var)
spm.mod.fit <- fit.variogram(spm.mod.var, spm.model1)
spm.val.pr <- krige(SPM07*1000~X+Y,spm.mod, spm.val, spm.mod.fit)
spm.res.kr <- spm.val$SPM07*1000 - spm.val.pr$var1.pred
spm.res.mean <- spm.val$SPM07*1000 - mean(spm.val$SPM07*1000)
R2 <- 1-sum(spm.res.kr^2)/sum(spm.res.mean^2)
R2
# Z値
spm.cv5 <- krige.cv(SPM07*1000~X+Y, data=spm, spm.fit, nfold=5)

## 参考
# クリギング：線形回帰モデル
# トレンドなし
spm.krg1 <- krige(SPM07*1000~1, locations=~X+Y, data=spm, mesh.grid, degree=2)
spplot(spm.krg1["var1.pred"], main = "ordinary kriging predictions")
# 緯度経度によるトレンド
spm.krg2 <- krige(SPM07*1000~X+Y, locations=~X+Y, data=spm, mesh.grid)
spplot(spm.krg2["var1.pred"], main = "ordinary kriging predictions")
