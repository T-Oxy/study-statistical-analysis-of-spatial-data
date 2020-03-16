# パッケージの読み込み
library(spdep) # 空間データの操作
library(maptools)
library(classInt)
library(spgwr) # 地理的加重回帰モデル
library(nlme) # マルチレベルモデル
# library(lme4)  # マルチレベルモデル
library(mgcv) # 一般化加法モデル
library(MASS)
# データの読み込み
lph <- read.table("lph.csv", sep=",", header=T)
summary(lph)
kanto <- readShapePoly("kanto_area.shp", IDvar="JCODE")

#### 10.1 回帰モデルと空間的自己相関 ####
# 図10.1
# 地価データの表示
pal1 <- gray.colors(n=4,start=1,end=0.3)
q_kanto <- classIntervals(round(kanto$LPH,1), n=4, style="quantile")
q_kanto_Col <- findColours(q_kanto,pal1)
plot(kanto,col=q_kanto_Col)
legend("bottomleft", fill=attr(q_kanto_Col,"palette"),
legend=names(attr(q_kanto_Col,"table")), cex=1.1, bty="n")
# 夜間人口密度
q_kanto_p <- classIntervals(round(kanto$POPD,1), n=4, style="quantile")
q_kanto_Col_p <- findColours(q_kanto_p,pal1)
plot(kanto,col=q_kanto_Col_p)
legend("topleft",fill=attr(q_kanto_Col_p,"palette"),
legend=names(attr(q_kanto_Col_p,"table")), cex=1.1, bty="n")
# 第三次産業従業人口密度
q_kanto_e <- classIntervals(round(kanto$EMP3D,1), n=4, style="quantile")
q_kanto_Col_e <- findColours(q_kanto_e,pal1)
plot(kanto,col=q_kanto_Col_e)
legend("topleft",fill=attr(q_kanto_Col_e,"palette"),
legend=names(attr(q_kanto_Col_e,"table")), cex=1.1, bty="n")
# 線形回帰モデル
# 通常最小二乗法
# GISデータの属性からモデル推定
lph.lm <- lm(LPH~POPD+EMP3D,data=kanto)
summary(lph.lm)
# 誤差の算出
lph.lm.resid <- resid(lph.lm)
kanto$lm.resid <- lph.lm.resid
# 図10.2
# 誤差の表示
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 主題図の作成
# 地価モデル（lph.lm）の誤差表示
q_lm.resid <- classIntervals(round(kanto$lm.resid,1), n=4, style="quantile")
q_lm.resid_Col <- findColours(q_lm.resid,pal1)
plot(kanto,col=q_lm.resid_Col)
legend("bottomleft",fill=attr(q_lm.resid_Col,"palette"),
legend=names(attr(q_lm.resid_Col,"table")), cex=1.1, bty="n")
# lphデータからモデル推定
# lph.lm <- lm(LPH~POPD+EMP3D,data=kanto)
# summary(lph.lm)
#ベイズ線形回帰モデル
# MCMCpackのMCMCregressを使う方法
# パッケージの呼び出し
library(MCMCpack)
lph.mcmc <- MCMCregress(LPH~POPD+EMP3D,data=kanto, b0=0, B0=1e-6, c0=1e-2, d0=1e-2, mcmc=100000, burnin=10000)
summary(lph.mcmc)
# 図10.3
plot(lph.mcmc)
# Jagsを使う方法
# パッケージの読み込み
library(R2jags)
# JAGS用データ作成
# 説明変数・被説明変数・誤差項
# 変数の指定
# 被説明変数
y <- as.vector(lph$LPH)
# 地域数
n <- length(y)
# 説明変数
x <- as.matrix(cbind(
lph$POPD,
lph$EMP3D
))
# 線形回帰モデルの誤差項（最小二乗法）
model.lm <- lm(y~x[,1]+x[,2])
summary(model.lm)
#
####
# JAGSコード(lm.txt)
model{
for(i in 1:n){
y[i] ~ dnorm(mu[i], tau)
mu[i] <- b0+b1*x[i,1]+b2*x[i,2]
}
b0 ~ dnorm(0,1.0E-6)
b1 ~ dnorm(0,1.0E-6)
b2 ~ dnorm(0,1.0E-6)
tau ~ dgamma(0.01, 0.01)
sigma <- 1/sqrt(tau)
}
#
####
# Rコード
# JAGS変数設定
# データ
data <- list("n", "y", "x")
# MCMC初期値（事前情報）
in1 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], tau=1)
in2 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], tau=1)
in3 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], tau=1)
inits <- list(in1,in2,in3)
# パラメータ
parameters <- c("b0", "b1", "b2", "tau", "sigma")
# モデルファイル
model.file <- system.file(package="R2jags", "model", "lm.txt")
# MCMC
lm.jags <- jags(data=data, inits=inits, parameters, n.iter=10000,
n.burnin=1000, n.chains=3, model.file=model.file)
print(lm.jags, digits=3)
# jpeg(filename="fig1.jpg")
plot(lm.jags)
# dev.off()
traceplot(lm.jags)
# jpeg(filename="fig2.jpg")
# traceplot(lm.jags, varname="b0")
# dev.off()
lm.fit <- update(lm.jags)
print(lm.fit, digits=3)
# 空間的自己相関
# 隣接行列の作成
coords <- matrix(0, nrow=length(kanto$LPH), ncol=2)
coords[,1] <- kanto$Easting
coords[,2] <- kanto$Northing
lph.tri.nb <- tri2nb(coords)
# Moran's I
moran.test(kanto$LPH, nb2listw (lph.tri.nb,style="W"))
moran.test(kanto$POPD, nb2listw (lph.tri.nb,style="W"))
moran.test(kanto$EMP3D, nb2listw (lph.tri.nb,style="W"))
# 誤差項の空間的自己相関
lm.morantest(lph.lm, nb2listw(lph.tri.nb, style="W"))
# moran.test(lph.lm$residuals, nb2listw (lph.tri.nb,style="W"))

#### 10.2 可変集計単位問題 ####

# 仮想データの作成
library(spatstat)
library(sp)
px <- rnorm(1000, mean=0, sd=2.5)
py <- rnorm(1000, mean=0, sd=2.5)
pnt <- ppp(px, py, c(-10,10), c(-10,10))
plot(pnt)
quad400 <- quadrats(pnt, 20, 20)
plot(quad400, main="")
points(pnt)
quad100 <- quadrats(pnt, 10, 10)
plot(quad100, main="")
points(pnt)
quad25 <- quadrats(pnt, 5, 5)
plot(quad25, main="")
points(pnt)
#
pvar1 <- rnorm(1000, mean=0, sd=1)
pvar2 <- rnorm(1000, mean=0, sd=1)
beta0 <- rnorm(1000, mean=0.8, sd=1.2)
beta1 <- rnorm(1000, mean=4, sd=1.8)
beta2 <- rnorm(1000, mean=1.2, sd=1.4)
pz <- beta0 + beta1 * pvar1 + beta2 * pvar2
lm.out1 <- lm(pz~pvar1+pvar2)
summary(lm.out1)

pnt_mat <- matrix(0, ncol=4, nrow=1000)
colnames(pnt_mat) <- c("grd_id", "pz", "pvar1", "pvar2")
pnt_df <- as.data.frame(pnt_mat)
pnt_df$grd_id <- (floor(px)+11)*100+(floor(py)+11)
pnt_df$pz <- pz
pnt_df$pvar1 <- pvar1
pnt_df$pvar2 <- pvar2
head(pnt_df)
grd400 <- cbind(
as.matrix(tapply(pnt_df$pz, factor(pnt_df$grd_id),sum)),
as.matrix(tapply(pnt_df$pvar1, factor(pnt_df$grd_id),sum)),
as.matrix(tapply(pnt_df$pvar2, factor(pnt_df$grd_id),sum)))
colnames(grd400) <- cbind("pz", "pvar1", "pvar2")
head(grd400)
grd400_df <- as.data.frame(grd400)
lm.out2 <- lm(pz~pvar1+pvar2, data=grd400_df)
summary(lm.out2)

pnt_df$grd_id2 <- ceiling((floor(px)+11)/2) *100 +ceiling((floor(py)+11)/2)
grd100 <- cbind(
as.matrix(tapply(pnt_df$pz, factor(pnt_df$grd_id2),sum)),
as.matrix(tapply(pnt_df$pvar1, factor(pnt_df$grd_id2),sum)),
as.matrix(tapply(pnt_df$pvar2, factor(pnt_df$grd_id2),sum)))
colnames(grd100) <- cbind("pz", "pvar1", "pvar2")
grd100 <- grd100/4
head(grd100)
grd100_df <- as.data.frame(grd100)
lm.out3 <- lm(pz~pvar1+pvar2, data=grd100_df)
summary(lm.out3)

pnt_df$grd_id3 <- ceiling((floor(px)+11)/4) *100 +ceiling((floor(py)+11)/4)
grd25 <- cbind(
as.matrix(tapply(pnt_df$pz, factor(pnt_df$grd_id3),sum)),
as.matrix(tapply(pnt_df$pvar1, factor(pnt_df$grd_id3),sum)),
as.matrix(tapply(pnt_df$pvar2, factor(pnt_df$grd_id3),sum)))
colnames(grd25) <- cbind("pz", "pvar1", "pvar2")
grd25 <- grd25/16
head(grd25)
grd25_df <- as.data.frame(grd25)
lm.out4 <- lm(pz~pvar1+pvar2, data=grd25_df)
summary(lm.out4)

lm.out <- as.matrix(cbind(
lm.out1$coefficient, lm.out2$coefficient, lm.out3$coefficient,
lm.out4$coefficient))
colnames(lm.out) <- cbind("lm.out1", "lm.out2", "lm.out3", "lm.out4")
matplot(t(lm.out), type="b", pch=1:4, axes=F, xlim=c(0.7,4),
ylim=c(min(lm.out),max(lm.out)),
xlab="モデル", ylab="パラメータ", col=1)
axis(1, 0:4)
# axis(2, min(lm.out):max(lm.out))
# legend("topleft", c("定数項", "説明変数1", "説明変数2"), col=1, lty=1:4)

#### 10.3 一般化モデル ####
# 一般化線型モデル [#x44c320d]
# glm()関数を使った一般化線形モデルの推定
lph.glm <- glm(LPH~POPD+EMP3D+offset(log(S)), data=lph)
summary(lph.glm)
# 一般化加法モデル [#cfbee64c]
# gam()関数を使った一般化加法モデルの推定→library(mgcv)
lph.gam1 <- gam(LPH~POPD+EMP3D+s(Easting, Northing), data=lph)
summary(lph.gam1)
lph.gam2 <- gam(LPH~POPD+EMP3D+offset(log(S))+s(Easting, Northing), data=lph)
summary(lph.gam2)
AIC(lph.lm, lph.glm, lph.gam1, lph.gam2)
# Moran固有ベクトルによる空間ラグの表現 [#wa010849]
# 空間フィルタリングによる空間ラグ
lph.SF <- SpatialFiltering(LPH~POPD+EMP3D, data=lph, nb=lph.tri.nb, style="W")
lph.lm.SF <- lm(LPH~POPD+EMP3D+fitted(lph.SF), data=lph)
summary(lph.lm.SF)
# 空間ラグを考慮した一般化線形回帰モデル
# 計算に時間を要する場合がある
lph.ME <- ME(LPH~POPD+EMP3D, data=lph, offset=log(S), listw=nb.w, alpha=0.5)
lph.glm.ME <- glm(LPH~POPD+EMP3D+offset(log(S))+fitted(lph.ME), data=lph)
summary(lph.glm.ME)

#### 10.4 自己回帰モデル ####
# 隣接行列の作成
coords <- matrix(0, nrow=length(kanto$LPH), ncol=2)
coords[,1] <- kanto$Easting
coords[,2] <- kanto$Northing
lph.tri.nb <- tri2nb(coords)
# 同時自己回帰モデル
# lph.sar <- spautolm(LPH~POPD+EMP3D, data=kanto,
# nb2listw(lph.tri.nb, style="W"), family="SAR", method="full")
lph.sar <- spautolm(LPH~POPD+EMP3D, data=kanto,
nb2listw(lph.tri.nb, style="W"), family="SAR", method="eigen",
verbose=TRUE)
summary(lph.sar)
# 誤差項の表示
lph.sar.resid <- resid(lph.sar)
kanto$sar.resid <- lph.sar.resid
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 図10.4
# 地価モデル（lph.lm）の誤差表示
q_lm.resid <- classIntervals(round(kanto$sar.resid,1), n=4, style="quantile")
q_lm.resid_Col <- findColours(q_lm.resid,pal1)
plot(kanto,col=q_lm.resid_Col)
legend("bottomleft",fill=attr(q_lm.resid_Col,"palette"),
legend=names(attr(q_lm.resid_Col,"table")), cex=1.1, bty="n")
# 条件付き自己回帰モデル [#z6104f58]
# lph.car <- spautolm(LPH~POPD+EMP3D, data=kanto,
# nb2listw(lph.tri.nb, style="W"), family="CAR", method="full")
lph.car <- spautolm(LPH~POPD+EMP3D, data=kanto,
nb2listw(lph.tri.nb, style="W"), family="CAR", method="eigen",
verbose=TRUE)
summary(lph.car)
# 誤差項の表示
lph.car.resid <- resid(lph.car)
kanto$car.resid <- lph.car.resid
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 図10.5
# 地価モデル（lph.lm）の誤差表示
q_lm.resid <- classIntervals(round(kanto$car.resid,1), n=4, style="quantile")
q_lm.resid_Col <- findColours(q_lm.resid,pal1)
plot(kanto,col=q_lm.resid_Col)
legend("bottomleft",fill=attr(q_lm.resid_Col,"palette"),
legend=names(attr(q_lm.resid_Col,"table")), cex=1.1, bty="n")

#### 10.5 空間的自己相関モデル ####
# 空間的従属性
LPH.lag <- lag.listw(nb2listw(lph.tri.nb, style="W"), kanto$LPH)
# 図10.6
plot(LPH.lag, kanto$LPH, ylab="y", xlab="Wy", cex=1.5, lwd=2,
cex.axis=1.3, cex.lab=1.2)
# 図10.7
plot(ecdf(kanto$LPH), main="", cex=1.5, cex.axis=1.3, cex.lab=1.2)
# 図10.8
plot(ecdf(LPH.lag), main="", cex=1.5, cex.axis=1.3, cex.lab=1.2)

## 10.5.1 空間的自己回帰モデル（空間同時自己回帰ラグモデル）
# 最尤法
lph.lag <- lagsarlm(LPH~POPD+EMP3D, data=kanto,
nb2listw(lph.tri.nb, style="W"))
# lph.lag <- lagsarlm(LPH~POPD+EMP3D, data=kanto,
# nb2listw(lph.tri.nb, style="W"), method="eigen", quiet=FALSE)
summary(lph.lag)
# 誤差項の表示
lph.lag.resid <- resid(lph.lag)
kanto$lag.resid <- lph.lag.resid
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 図10.9
# 地価モデル（lph.lm）の誤差表示
q_lm.resid <- classIntervals(round(kanto$lag.resid,1), n=4, style="quantile")
q_lm.resid_Col <- findColours(q_lm.resid,pal1)
plot(kanto,col=q_lm.resid_Col)
legend("bottomleft",fill=attr(q_lm.resid_Col,"palette"),
legend=names(attr(q_lm.resid_Col,"table")), cex=1.1, bty="n")
# 二段階最小二乗法
lph.stsls <- stsls(LPH~POPD+EMP3D, data=kanto, nb2listw(lph.tri.nb, style="W"))
summary(lph.stsls)
# 最小二乗法
# LPH.lag <- lag.listw(nb2listw(lph.tri.nb, style="W"), lph$LPH)
# plot(LPH.lag, lph$LPH)
# plot(ecdf(lph$LPH))
# plot(ecdf(LPH.lag))
lph.lag2 <- lm(lph$LPH~LPH.lag+lph$POPD+lph$EMP3D)
summary(lph.lag2)
# ベイズ法
# JAGS用データ作成
# 説明変数・被説明変数・誤差項
# 変数の指定
# 被説明変数
y <- as.vector(lph$LPH)
# 地域数
n <- length(y)
# 地域セグメント数
m <- max(lph$seg)
# 説明変数
x <- as.matrix(cbind(
lph$POPD,
lph$EMP3D,
lph$seg,
lph$ID
))
# 線形回帰モデルの誤差項（最小二乗法）
model.lm <- lm(y~x[,1]+x[,2])
summary(model.lm)
#
# 空間隣接行列の作成
coords <- as.matrix(cbind(lph$Easting,lph$Northing))
# 空間重み付け行列の作成
nb <- tri2nb(coords)
# nb.mat <- nb2mat(nb, style="B")
nb.mat <- nb2mat(nb, style="W")
# 空間ウェイト付き変数
yL <- nb.mat %*% y
xL <- nb.mat %*% x
#
####
# JAGSコード(slag1.txt)
# spatial lag model
model{
for (i in 1 : n) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- rho*yL[i,1]+b0+b1*x[i,1]+b2*x[i,2]
}
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
b2 ~ dnorm(0.0, 1.0E-6)
tau ~ dgamma(0.001, 0.001)
rho ~ dnorm(0.0, 1.0E-6)
sigma <- 1/sqrt(tau)
}
#
####
# Rコード
# JAGS変数設定
# データ
data <- list("n", "y", "x", "yL")
# MCMC初期値（事前情報）
in1 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], rho=0, tau=1)
in2 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], rho=0.5, tau=1)
in3 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], rho=1, tau=1)
inits <- list(in1,in2,in3)
# パラメータ
parameters <- c("b0", "b1", "b2", "rho", "tau", "sigma")
# モデルファイル
model.file <- system.file(package="R2jags", "model", "slag1.txt")
# MCMC
slag.jags <- jags(data=data, inits=inits, parameters, n.iter=10000,
n.burnin=1000, n.chains=3, model.file=model.file)
print(slag.jags, digits=3)
plot(slag.jags)
traceplot(slag.jags)
slag.fit <- update(slag.jags)
print(slag.fit, digits=3)

## 10.5.2 誤差項の空間的自己回帰モデル（空間同時自己回帰誤差モデル）
# 最尤法
lph.err <- errorsarlm(LPH~POPD+EMP3D, data=kanto,
nb2listw(lph.tri.nb, style="W"))
# lph.err <- errorsarlm(LPH~POPD+EMP3D, data=kanto,
# nb2listw(lph.tri.nb, style="W"), method="eigen", quiet=FALSE)
summary(lph.err)
# 誤差項の表示
lph.err.resid <- resid(lph.err)
kanto$err.resid <- lph.err.resid
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 図10.10
# 地価モデル（lph.lm）の誤差表示
q_lm.resid <- classIntervals(round(kanto$err.resid,1), n=4, style="quantile")
q_lm.resid_Col <- findColours(q_lm.resid,pal1)
plot(kanto,col=q_lm.resid_Col)
legend("bottomleft",fill=attr(q_lm.resid_Col,"palette"),
legend=names(attr(q_lm.resid_Col,"table")), cex=1.1, bty="n")
# Generalised Moments estimator
# vv <- lph.GMerr$vv
# lambda <- lph.GMerr$lambda
# s2 <- lph.GMerr$s2
lph.GMerr <- GMerrorsar(LPH~POPD+EMP3D, data=kanto,
nb2listw(lph.tri.nb, style="W"))
summary(lph.GMerr)
lseq <- seq(-0.1, 0.5, 0.01)
s2seq <- seq(65, 65.6, 0.01)
lsseq <- as.matrix(expand.grid(lseq, s2seq))
res <- numeric(nrow(lsseq))
for (i in seq(along=res)) res[i] <- spdep:::.kpgm(lsseq[i,,drop=TRUE],
v=lph.GMerr$vv, verbose=TRUE)
SGDF <- SpatialPixelsDataFrame(lsseq, data=data.frame(fn=res))
fullgrid(SGDF) <- TRUE
# 図10.11
image(SGDF, "fn", col=grey.colors(10, 1, 0.3, 2.2), axes=TRUE)
title(xlab=expression(lambda), ylab=expression(sigma^2))
contour(SGDF, "fn", add=TRUE)
points(c(lph.GMerr$lambda, lph.err$lambda), c(lph.GMerr$s2, lph.err$s2),
pch=c(4, 3), lwd=2)
text(c(lph.GMerr$lambda, lph.err$lambda), c(lph.GMerr$s2, lph.err$s2),
labels=c("GM", "ML"), pos=c(2, 4), offset=0.5)
# Spatial Error Model
####
# JAGSコード(sem1.txt)
# spatial error model
model{
for (i in 1 : n) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- lambda*yL[i,1]+b0*(1-lambda)+b1*x[i,1]+b2*x[i,2]
-b1*lambda*xL[i,1]-b2*lambda*xL[i,2]
}
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
b2 ~ dnorm(0.0, 1.0E-6)
tau ~ dgamma(0.001, 0.001)
lambda ~ dunif(0,1)
sigma <- 1/sqrt(tau)
}
#
####
# Rコード
# JAGS変数設定
# データ
data <- list("n", "y", "x", "yL", "xL")
# MCMC初期値（事前情報）
in1 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], lambda=0, tau=1)
in2 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], lambda=0.5, tau=1)
in3 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], lambda=1, tau=1)
inits <- list(in1,in2,in3)
# パラメータ
parameters <- c("b0", "b1", "b2", "lambda", "tau", "sigma")
# モデルファイル
model.file <- system.file(package="R2jags", "model", "sem1.txt")
# MCMC
sem.jags <- jags(data=data, inits=inits, parameters, n.iter=10000,
n.burnin=1000, n.chains=3, model.file=model.file)
print(sem.jags, digits=3)
plot(sem.jags)
traceplot(sem.jags)
sem.fit <- update(sem.jags)
print(sem.fit, digits=3)

## 10.5.3 Spatial Durbinモデル
# 最尤法
lph.durbin <- lagsarlm(LPH~POPD+EMP3D, data=kanto,
nb2listw(lph.tri.nb, style="W"), type="mixed")
#  lph.durbin <- lagsarlm(LPH~POPD+EMP3D, data=kanto,
# nb2listw(lph.tri.nb, style="W"), type="mixed", method="eigen", quiet=FALSE)
summary(lph.durbin)
# 誤差項の表示
lph.durbin.resid <- resid(lph.durbin)
kanto$durbin.resid <- lph.durbin.resid
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 図10.12
# 地価モデル（lph.lm）の誤差表示
q_lm.resid <- classIntervals(round(kanto$durbin.resid,1), n=4, style="quantile")
q_lm.resid_Col <- findColours(q_lm.resid,pal1)
plot(kanto,col=q_lm.resid_Col)
legend("bottomleft",fill=attr(q_lm.resid_Col,"palette"),
legend=names(attr(q_lm.resid_Col,"table")), cex=1.1, bty="n")
# ベイズ法
####
# JAGSコード(sdm1.txt)
# spatial durbin model
model{
for (i in 1 : n) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- rho*yL[i,1]+b0+b1*x[i,1]+b2*x[i,2]+g1*xL[i,1]+g2*xL[i,2]
}
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
b2 ~ dnorm(0.0, 1.0E-6)
g1 ~ dnorm(0.0, 1.0E-6)
g2 ~ dnorm(0.0, 1.0E-6)
tau ~ dgamma(0.001, 0.001)
rho ~ dnorm(0.0, 1.0E-6)
sigma <- 1/sqrt(tau)
}
#
####
# Rコード
# JAGS変数設定
# データ
data <- list("n", "y", "x", "yL", "xL")
# MCMC初期値（事前情報）
in1 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], g1=model.lm$coefficients[2],
g2=model.lm$coefficients[3], rho=0, tau=1)
in2 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], g1=model.lm$coefficients[2],
g2=model.lm$coefficients[3], rho=0.5, tau=1)
in3 <- list(b0=model.lm$coefficients[1],b1=model.lm$coefficients[2],
b2=model.lm$coefficients[3], g1=model.lm$coefficients[2],
g2=model.lm$coefficients[3], rho=1, tau=1)
inits <- list(in1,in2,in3)
# パラメータ
parameters <- c("b0", "b1", "b2", "g1", "g2", "rho", "tau", "sigma")
# モデルファイル
model.file <- system.file(package="R2jags", "model", "sdm1.txt")
# MCMC
sdm.jags <- jags(data=data, inits=inits, parameters, n.iter=10000,
n.burnin=1000, n.chains=3, model.file=model.file)
print(sdm.jags, digits=3)
plot(sdm.jags)
traceplot(sdm.jags)
sdm.fit <- update(sdm.jags)
print(sdm.fit, digits=3)

## 空間的自己相関モデルを用いた予測
predict(lph.err)
predict(lph.lag)
predict(lph.durbin)

## 10.5.5 空間従属性の検定
lm.LMtests(lph.lm, nb2listw(lph.tri.nb), test=c("LMerr", "LMlag", "RLMerr",
"RLMlag", "SARMA"))
# LMerr: 誤差項の空間依存性に関する検定
# LMlag: 空間ラグ依存変数がない場合の検定
# RLMerr: 空間ラグ依存変数がない場合の誤差項の空間依存性に関する検定
# RLMlag: 空間ラグ依存変数がある場合の誤差項の空間依存性に関する検定
# SARMA: ポルトマント検定
# lm.LMtests(lph.lm, nb2listw(lph.tri.nb))
# lm.LMtests(residuals(lph.lm), nb2listw(lph.tri.nb))

#### 10.6 マルチレベルモデル ####
# lmer()関数を使ったマルチレベルモデルの推定
# パッケージ
library(lme4)
# 固定効果：傾き、ランダム効果：切片
lph.lme1 <- lmer(LPH~POPD+EMP3D+(1|AREA), data=lph)
summary(lph.lme1)
ranef(lph.lme1)
# random.effects(lph.lme1)
# 固定効果：切片、ランダム効果：傾き
lph.lme2 <- lmer(LPH~1+(0+POPD+EMP3D|AREA), data=lph)
summary(lph.lme2)
# random.effects(lph.lme2)
ranef(lph.lme2)
# mcmcsamp()関数を使ったマルチレベルモデルのベイズ推定
lph.mcmc1 <- mcmcsamp(lph.lme1, n=1000)
print(lph.mcmc1)
densityplot(lph.mcmc1)
#
lph.mcmc2 <- mcmcsamp(lph.lme2, n=1000)
print(lph.mcmc2)
densityplot(lph.mcmc2)
# ベイズ推定
# bayesmパッケージのrhierLinearModel()関数を使う推定方法
# 市区町村単位でパラメータ推定しているため計算に時間を要します
library(bayesm)
# mcmcの設定
R=10000
keep=1
# 事前分布を設定
reg=levels(factor(lph$JCODE))
nreg=length(reg)
nvar=3 	#説明変数２つ＋定数項
# 変数を設定
regdata=NULL
for (j in 1:nreg) {
       y=lph$LPH[lph$JCODE==reg[j]]
       iota=c(rep(1,length(y)))
       X=cbind(iota, lph$POPD[lph$JCODE==reg[j]],
               lph$EMP3D[lph$JCODE==reg[j]])
       regdata[[j]]=list(y=y,X=X)}
Z=matrix(c(rep(1,nreg)),ncol=1)
Data1=list(regdata=regdata,Z=Z)
Mcmc1=list(R=R,keep=1)
set.seed(66)
# モデル推定
out=rhierLinearModel(Data=Data1,Mcmc=Mcmc1)
# 推定結果の表示
# summary(out$Deltadraw, burnin=1000)
# summary(out$Vbetadraw, burnin=1000)
summary(t(out$betadraw[1,,]), burnin=1000)
plot(out$betadraw, burnin=1000)
#
out_mat <- matrix(0,nrow=nreg,ncol=nvar+3)
colnames(out_mat) <- c("JCODE","POPD","EMP3D", "const","b.POPD","b.EMP3D")
out_mat <- as.data.frame(out_mat)
out_mat$JCODE <- lph$JCODE
out_mat$POPD  <- lph$POPD
out_mat$EMP3D <- lph$EMP3D
for(i in 1:nreg){
	for(j in 1:nvar){
		out_mat[i,j+3] <- mean(out$betadraw[i,j,1001:10000])
}}
# 推定結果の可視化
# パラメータ推定結果を地図属性に結合
ID.match <- match(kanto$JCODE, out_mat$JCODE)
out_mat1 <- out_mat[ID.match,]
row.names(out_mat1) <- kanto$JCODE
kanto_out <- spCbind(kanto, kanto1)
#
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 図10.13(a)
q_kanto <- classIntervals(round(kanto_out$b.POPD,2), n=4, style="quantile")
q_kanto_Col <- findColours(q_kanto,pal1)
plot(kanto_out,col=q_kanto_Col)
legend("bottomleft",fill=attr(q_kanto_Col,"palette"),
legend=names(attr(q_kanto_Col,"table")), cex=1.1, bty="n")
# 図10.13(b)
q_kanto <- classIntervals(round(kanto_out$b.EMP3D,2), n=4, style="quantile")
q_kanto_Col <- findColours(q_kanto,pal1)
plot(kanto_out,col=q_kanto_Col)
legend("bottomleft",fill=attr(q_kanto_Col,"palette"),
legend=names(attr(q_kanto_Col,"table")), cex=1.1, bty="n")

#### 10.7 地理的加重回帰モデル ####
# パッケージ
library(spgwr)
# バンド幅の計算
lph.bw <- gwr.sel(LPH~POPD+EMP3D, data=lph, coords=coords)
# lph.bw.csvはgwr.selの計算過程を出力したもの
lph.bw.dat <- read.table("lph.bw.csv", sep=",", header=TRUE)
# 図10.14
plot(lph.bw.dat, type="o", cex=2, lwd=3, ylab="CV", xlab="バンド幅",
cex.axis=1.3, cex.lab=1.2)
# 地理的加重回帰モデルの推定
lph.gwr <- gwr(LPH~POPD+EMP3D, data=lph, coords=coords, bandwidth=lph.bw,  hatmatrix=TRUE)
summary(lph.gwr$SDF)
# 主題図の作成
# 回帰係数POPDの表示
kanto$gwr.popd <- lph.gwr$SDF$POPD
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
q_gwr <- classIntervals(round(kanto$gwr.popd,2), n=4, style="quantile")
q_gwr_Col <- findColours(q_gwr,pal1)
plot(kanto,col=q_gwr_Col)
legend("bottomleft",fill=attr(q_gwr_Col,"palette"),
legend=names(attr(q_gwr_Col,"table")), cex=1.1, bty="n")
# 回帰係数EMP3Dの表示
kanto$gwr.emp3d <- lph.gwr$SDF$EMP3D
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
q_gwr <- classIntervals(round(kanto$gwr.emp3d,2), n=4, style="quantile")
q_gwr_Col <- findColours(q_gwr,pal1)
plot(kanto,col=q_gwr_Col)
legend("bottomleft",fill=attr(q_gwr_Col,"palette"),
legend=names(attr(q_gwr_Col,"table")), cex=1.1, bty="n")
# 標準誤差の表示
lph.gwr.resid <- resid(lph.gwr)
kanto$gwr.pred.se <- lph.gwr$SDF$pred.se
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
q_gwr.se <- classIntervals(round(kanto$gwr.pred.se,2), n=4, style="quantile")
q_gwr.se_Col <- findColours(q_gwr.se,pal1)
plot(kanto,col=q_gwr.se_Col)
legend("bottomleft",fill=attr(q_gwr.se_Col,"palette"),
legend=names(attr(q_gwr.se_Col,"table")), cex=1, bty="n")
# 回帰係数localR2の表示
kanto$gwr.localR2 <- lph.gwr$SDF$localR2
# 色パレットの作成
pal1 <- gray.colors(n=4,start=1,end=0.3)
# 地価モデル（lph.lm）の誤差表示
q_gwr <- classIntervals(round(kanto$gwr.localR2,2), n=4, style="quantile")
q_gwr_Col <- findColours(q_gwr,pal1)
plot(kanto,col=q_gwr_Col)
legend("bottomleft",fill=attr(q_gwr_Col,"palette"),
legend=names(attr(q_gwr_Col,"table")), cex=1.1, bty="n")
