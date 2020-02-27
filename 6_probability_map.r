library(spdep)
library(maptools)

#### 6.1 粗率 ####
# 地図データの読み込み
jpn_pref <- readShapePoly("jpn_pref.shp",IDvar="PREF_CODE")
pref_pnt <- readShapePoints("pref_gov.shp")
hd06 <- read.table("hd06.csv", sep=",", header=TRUE)
summary(hd06)
# 都道府県境界ポリゴンデータへの属性データhd06のマッチング
ID.match <- match(jpn_pref$PREF_CODE, hd06$PREF_CODE)
jpn_hd06 <- hd06[ID.match,]
jpn_pref_hd06 <- spCbind(jpn_pref, jpn_hd06)
summary(jpn_pref_hd06)
#確率地図の作成
# probmap()関数
jpn_hd06_pm <- probmap(jpn_pref_hd06$HD06, jpn_pref_hd06$POPJ06/100)
summary(jpn_hd06_pm)
# 主題図の作成
library(classInt)
brks1 <- c(0, 2000, 3000, 4000, 20000)
brks2 <- c(0, 125, 150, 175, 200)
brks3 <- c(0, 0.2, 0.6, 0.9, 1.0)
cols <- grey(4:1/4)
#図6.1
# 観測値
plot(jpn_pref, col=cols[findInterval(jpn_pref_hd06$HD06, brks1, all.inside=TRUE)])
legend("topleft", fill=cols, legend=leglabs(brks1), cex=1.5, bty="n")
# title(main="Observed number of heart disease death (2006)")
#図6.2
# 粗率（jpn_hd06_pm$raw）
plot(jpn_pref, col=cols[findInterval(jpn_hd06_pm$raw, brks2, all.inside=TRUE)])
legend("topleft", fill=cols, legend=leglabs(brks2), cex=1.5, bty="n")
title(main="Raw rate of heart disease per 100,000 (2006)")


#### 6.2 相対危険度 ####
#図6.3
# 期待値（jpn_hd06_pm$expCount）
plot(jpn_pref, col=cols[findInterval(jpn_hd06_pm$expCount, brks1, all.inside=TRUE)])
legend("topleft", fill=cols, legend=leglabs(brks1), cex=1.5, bty="n")
title(main="Expected number of heart disease death (2006)")
#図6.4
# 相対危険度（jpn_hd06_pm$relRisk）
plot(jpn_pref, col=cols[findInterval(jpn_hd06_pm$relRisk, brks2, all.inside=TRUE)])
legend("topleft", fill=cols, legend=leglabs(brks2), cex=1.5, bty="n")
title(main="Relevant risk of heart disease per 100,000 (2006)")

#### 6.3 ポアソン確率（jpn_hd06_pm$pmap） ####
# ポアソン分布
plot(dpois(0:50, lambda=1), type="l", xlim=c(0,20), lwd=2, ylab="確率",
xlab="x",cex.axis=1.2,cex.lab=1.2,font.lab=4)
text(3.5, 0.36, "Po(1)", cex=1.3)
lines(dpois(0:50, lambda=3), type="l", lwd=2)
text(5, 0.24, "Po(3)", cex=1.3)
lines(dpois(0:50, lambda=5), type="l", lwd=2)
text(7, 0.19, "Po(5)", cex=1.3)
lines(dpois(0:50, lambda=10), type="l", lwd=2)
text(11, 0.14, "Po(10)", cex=1.3)
# ポアソン確率
#図6.6
plot(jpn_pref, col=cols[findInterval(jpn_hd06_pm$pmap, brks3, all.inside=TRUE)])
legend("topleft", fill=cols, legend=leglabs(brks3), cex=1.5, bty="n")
title(main="Poisson probability of heart disease (2006)")
#図6.7
hist(jpn_hd06_pm$pmap, col="grey", ylim=c(0,40), main="", xlab="確率",
ylab="度数", cex.axis=1.3, cex.lab=1.2)

#### 6.4 相対リスクのベイズ推定 ####
# 経験ベイズ推定（Marshallのグローバルな経験ベイズ推定値）
jpn_hd06_ebg <- EBest(jpn_pref_hd06$HD06, jpn_pref_hd06$POPJ06/100)
#図6.8
plot(jpn_pref, col=cols[findInterval(jpn_hd06_ebg$estmm, brks2, all.inside=TRUE)])
legend("topleft", fill=cols, legend=leglabs(brks2), cex=1.5, bty="n")
title(main="Empirical Bayes estimates of heart disease death rate per 100,000 (2006)", cex.main=0.9)
# localな経験ベイズ推定（Marshallのローカルな経験ベイズ推定値）
pref_gov <- read.table("pref_gov.txt",",",header=TRUE, row.names=2)
coords <- matrix(0,nrow(pref_gov),2)
coords[,1] <- pref_gov$X
coords[,2] <- pref_gov$Y
pref.knn <- knearneigh(coords, k=4)
pref.knn.nb <- knn2nb(pref.knn, row.names=rownames(pref_gov))
jpn_hd06_ebl <- EBlocal(jpn_pref_hd06$HD06, jpn_pref_hd06$POPJ06/100, pref.knn.nb)
#図6.9
plot(jpn_pref, col=cols[findInterval(jpn_hd06_ebl$est, brks2, all.inside=TRUE)])
plot(pref.knn.nb, coords, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks2), cex=1.5, bty="n")
title(main="Local empirical Bayes estimates of heart disease death rate per 100,000 (2006)", cex.main=0.8)
# 経験ベイズ推定値によるMoran's Iの繰り返し検定
EBImoran.mc(jpn_pref_hd06$HD06, jpn_pref_hd06$POPJ06/100,
nb2listw(pref.knn.nb, style="B", zero.policy=TRUE),
nsim=999, zero.policy=TRUE)
# 通常のMoran's Iの繰り返し検定
moran.mc(jpn_pref_hd06$HD06/(jpn_pref_hd06$POPJ06/100),
nb2listw(pref.knn.nb, style="B", zero.policy=TRUE),
nsim=999, zero.policy=TRUE)
# Γ(ν, α)関数
x <- 1:1000/1000
plot(x, dgamma(x, 20,50), type="l", xlim=c(0,1), ylim=c(0,10),lwd=2,
ylab="確率",xlab="x",main=expression(Gamma(nu,alpha)),
cex.axis=1.2,cex.lab=1.2,font.lab=4)
text(0.44, 5, expression(Gamma(20,50)), cex=1.2)
lines(x, dgamma(x, 1,20), type="l", lwd=2)
text(0.14, 8, expression(Gamma(1,20)), cex=1.2)
lines(x, dgamma(x, 50,70), type="l", lwd=2)
text(0.7, 4.4, expression(Gamma(50,70)), cex=1.2)
# ポアソンーガンマモデル（経験ベイズ）
# epiRを使う方法
install.packages("epiR")
library(epiR)
jpn_hd06_sm <- epi.empbayes(jpn_hd06$HD06, jpn_hd06_pm$expCount)
jpn_hd06_sm
# SpatialEpiを使う方法
install.packages("SpatialEpi")
library(SpatialEpi)
jpn_hd06_eB <- eBayes(jpn_hd06$HD06, jpn_hd06_pm$expCount)
jpn_hd06_eB
# 対数正規モデル
# R2.15.0以前（DClusterのみ）
jpn_hd06_ln <- lognormalEB(jpn_hd06$HD06, jpn_hd06_pm$expCount)
jpn_hd06_ln
# 6.4.2
# ポアソン・ガンマモデル（階層ベイズと経験ベイズ：事例＝交通事故）R2.15.0以前
library(DCluster)
library(pscl)
library(R2jags)
data84 <- read.table("data84.csv", sep=",", header=TRUE)
summary(data84)
r <- sum(data84$TrAcc)/sum(data84$Pop)
data84$TAExpected <- data84$Pop * r
data84$TARR <- data84$TrAcc / data84$TAExpected
# 観測値のヒストグラム
# 図6.10(a)
hist(data84$TrAcc, 0:30000, xlim=c(0,1000), ylim=c(c,40),
main="", ylab="度数", xlab="", cex.axis=1.3, cex.lab=1.2)
# 期待値のヒストグラム
hist(data84$TAExpected, 0:30000, xlim=c(0,1000), ylim=c(c,20))
# 相対リスクのヒストグラム
#図6.10(b)
hist(data84$TARR, col="grey", ylim=c(0,1000), main="", ylab="度数", xlab="",
cex.axis=1.3, cex.lab=1.2)
# Marshallのグローバルな経験ベイズ推定量
data84$EB <- EBest(data84$TrAcc, data84$TAExpected)

# 経験ベイズ推定
eb <- eBayes(data84$TrAcc,data84$TAExpected)
data84$EBPG <- eb$SMR

#図6.12
hist(data84$EBPG, col="grey", main="", ylab="度数", xlab="", ylim=c(0, 1000),
cex.axis=1.3, cex.lab=1.2)

# 階層ベイズ推定
n <- nrow(data84)
acc.obs <- data84$TrAcc
acc.exp <- data84$TAExpected
data <- list("n","acc.obs","acc.exp")
inits <- list(nu=1,alpha=1)
parameters<- c("theta","nu","alpha")
model.file <- system.file(package="R2jags", "model", "poisson_gamma.txt")
model84_pg <- jags(data=data, inits=inits, parameters, n.iter=1000,
n.burnin=100, n.chains=1, model.file=model.file)
print(model84_pg, digits=3)
# 図6.13
plot(model84_pg)

# 対数正規モデル  R2.15.0以前
data84_ln <- lognormalEB(data84$TrAcc, data84$TAExpected)
data84_ln
#図6.14
hist(data84_ln$smthrr, col="grey", main="", xlab="", ylab="度数",
cex.axis=1.3, cex.lab=1.2)
