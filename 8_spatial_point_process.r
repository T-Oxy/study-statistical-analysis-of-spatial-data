#### データの作成と表示 ####
# パッケージの読み込み
library(spatstat)
# データの読み込み
X600 <- read.table("X600.txt", sep=",", header=TRUE)
m <- sample(c("c","m","y","b"), 600, replace=TRUE)
m <- factor(m, levels=c("c","m","y","b"))
X <- ppp(X600$x, X600$y, c(0,1), c(0,1), marks=m)
# データ作成
x <- runif(600)
y <- runif(600)
m <- sample(c("c","m","y","b"), 600, replace=TRUE)
m <- factor(m, levels=c("c","m","y","b"))
X <- ppp(x, y, c(0,1), c(0,1), marks=m)
# データ表示
# 図8.1
plot(split(X)$b, main="")
points(split(X)$m, pch=2)
points(split(X)$y, pch=3)
points(split(X)$c, pch=5)

#### 8.1 コドラート法によるランダム性の検定 (Complete Spatial Randomness: CSR) ####
quadratcount(X, nx=5, ny=5)
plot(X)
plot(quadratcount(X, nx=5, ny=5), cex=1.2, add=TRUE)
# ティーセン分割
plot(dirichlet(X))
# ドロネー三角網図
plot(delaunay(X))
# コドラート法によるランダム性のχ2検定
all points
quadrat.test(X, nx=5, ny=5)
# 図8.4
plot(X, cex=0.2, main="")
plot(quadrat.test(X, nx=5, ny=5),  col=1, cex=1.5, add=TRUE)
# 'cyan' points
quadrat.test(split(X)$c, nx=5, ny=5)
plot(split(X)$c)
plot(quadrat.test(split(X)$c, nx=5, ny=5), add=TRUE)

#### 8.3 コルモゴロフ・スミルノフ検定 ####
xcoord <- function(x,y) {x}
kstest(split(X)$c, xcoord)
# 図8.5
plot(kstest(split(X)$c, xcoord), lwd=2, cex.axis=1.3, cex.lab=1.2)

#### 8.4 ポアソンモデルの適用:　λ_θ(x,y)=θ_0+θ_1*x+θ_2*y ####
modelc <- ppm(split(X)$c, ~x+y)
# 図8.6
plot(modelc, how="image", se=FALSE)
plot(predict(modelc, type="trend", ngrid=512))

#### 8.5 距離に基づく関数を用いる方法 ####
## 8.5.1 境界効果
x <- runif(100)
y <- runif(100)
X1 <- ppp(x, y, c(0,1), c(0,1))
Y1 <- ppp(x,y,c(0.1,0.9),c(0.1,0.9))
# 距離マップ 図8.7
plot(distmap(Y1), main="")
points(Y1, col="white", pch=19)
# 境界効果 図8.8(a)
plot(X1, pch=19, main="")
rect(0.1,0.1,0.9,0.9, lty=2, lwd=3)
# 図8.8(b)
plot(distmap(X1), main="")
points(X1, col="white", pch=19)
rect(0.1,0.1,0.9,0.9, lty=2, lwd=3)
# 図8.8(c)
plot(distmap(Y1), xlim=c(0,1), ylim=c(0,1), main="")
points(Y1, pch=19)

## 8.5.2 F関数：一定半径以内  図8.9
plot(Fest(split(X)$c), theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)

## 8.5.3 G関数：最近隣距離  図8.10
plot(Gest(split(X)$c), theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)

## 8.5.4 K関数：ペアワイズ距離  図8.11
plot(Kest(split(X)$c), theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)

## 8.5.5 L関数：一般に用いられるK関数  図8.12
plot(Lest(X), theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)

## 8.5.6 J関数：{1-G}/{1-F}  図8.13
plot(Jest(X), km~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)

## 8.5.7 ペア相関関数：K'/πr^2  図8.14
plot(pcf(X), iso~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)

## 不均一なK関数 (inhomogeneous K-function)
predc <- predict(modelc, locations=X)
Ki <- Kinhom(X, predc)
plot(Ki)
Ki2 <- Kinhom(X)
plot(Ki2)

#### 8.6 マーク付き点過程 ####
# marked data
plot(split(X))
plot(density(split(X)))
plot(smooth.ppp(X))
# 要約統計量
M <- marktable(X, R=0.1)
# multi-type
# Gij関数
plot(Gcross(X, "c", "b"))
plot(alltypes(X, "G"))

# Kij関数
# 図8.15(a)
plot(Kcross(X, "c", "b"), iso~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)
# 図8.15(b)
plot(alltypes(X, "K"), iso~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)

# Lij関数
plot(Lcross(X, "c", "b"))
# Jij関数
plot(Jcross(X, "c", "b"))
plot(alltypes(X, "J"))
# i-to-any
# Gdot関数
plot(alltypes(X, "Gdot"))
# マーク相関
eqfun <- function(m1, m2) { m1 == m2}
M <- markcorr(X, eqfun, correction="translate", method="density", kernel="epanechnikov")
plot(M)

## ランダム性の検定
# 帰無仮説：マーク付き点過程がポアソン均一過程に従う
E <- envelope(X, Kcross, nsim=99, i="b", j="c")
plot(E)
# 要素の独立性
E <- envelope(X, Kcross, nsim=99, i="b", j="c", simulate=expression(rshift(X, radius=0.1)))
plot(E)

#### 8.7 包絡分析と適合度検定 ####
# ポイントワイズ包絡分析
E <- envelope(X, Kest, nsim=99)
plot(E)
# マーク別包絡分析
# シアン　F関数
Fenv.c <- envelope(split(X)$c, fun=Fest, nsim=100)
xx <- seq(from=0, to=max(Fenv.c$r),length=length(Fenv.c$lo))
# 図8.16(a)
plot(Fenv.c, theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)
polygon(c(xx,max(Fenv.c$r)-xx), c(Fenv.c$hi, sort(Fenv.c$lo, decreasing=TRUE)), col="grey", border="grey")
plot(Fenv.c, theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4, add=TRUE)
# シアン　G関数
Genv.c <- envelope(split(X)$c, fun=Gest, nsim=100)
xx <- seq(from=0, to=max(Genv.c$r),length=length(Genv.c$lo))
# 図8.16(b)
plot(Genv.c, theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4)
polygon(c(xx,max(Genv.c$r)-xx), c(Genv.c$hi, sort(Genv.c$lo, decreasing=TRUE)), col="grey", border="grey")
plot(Genv.c, theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3, font.lab=4, add=TRUE)
# シアン　K関数
Kenv.c <- envelope(split(X)$c, fun=Kest, nsim=100)
xx <- seq(from=0, to=max(Kenv.c$r),length=length(Kenv.c$lo))
# 図8.16(c)
plot(Kenv.c, theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3,font.lab=4)
polygon(c(xx,max(Kenv.c$r)-xx), c(Kenv.c$hi, sort(Kenv.c$lo, decreasing=TRUE)), col="grey", border="grey")
plot(Kenv.c, theo~r, main="", cex.axis=1.2, cex.lab=1.2, lwd=3,font.lab=4, add=TRUE)
# マゼンダ
Kenv.m <- envelope(split(X)$m, fun=Kest, nsim=100)
xx <- seq(from=0, to=max(Kenv.m$r),length=length(Kenv.m$lo))
plot(Kenv.m, theo~r, lwd=2, main="")
polygon(c(xx,max(Kenv.m$r)-xx), c(Kenv.m$hi, sort(Kenv.m$lo, decreasing=TRUE)), col="grey", border="grey")
plot(Kenv.m, theo~r, lwd=2, add=T, main="")
# イエロー
Kenv.y <- envelope(split(X)$y, fun=Kest, nsim=100)
xx <- seq(from=0, to=max(Kenv.y$r),length=length(Kenv.y$lo))
plot(Kenv.y, theo~r, lwd=2, main="")
polygon(c(xx,max(Kenv.y$r)-xx), c(Kenv.y$hi, sort(Kenv.y$lo, decreasing=TRUE)), col="grey", border="grey")
plot(Kenv.y, theo~r, lwd=2, add=T, main="")
# イエロー
Kenv.b <- envelope(split(X)$b, fun=Kest, nsim=100)
xx <- seq(from=0, to=max(Kenv.b$r),length=length(Kenv.b$lo))
plot(Kenv.b, theo~r, lwd=2, main="")
polygon(c(xx,max(Kenv.b$r)-xx), c(Kenv.b$hi, sort(Kenv.b$lo, decreasing=TRUE)), col="grey", border="grey")
plot(Kenv.b, theo~r, lwd=2, add=T, main="")
# モンテカルロ検定
E <- envelope(X, Kest, nsim=99, global=T)
plot(E)
