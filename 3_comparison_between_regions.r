library(sp)

#### 3.1 密度 ####
library(spatstat)
# 500人の居住地の座標を正規分布に則り決定
px <- rnorm(500, mean=0.5, sd=0.15)
py <- rnorm(500, mean=0.5, sd=0.15)
pz <- as.data.frame(rep(1, 500))
colnames(pz) <- c("pz")

# quadratcount(Divides window into quadrats and counts the numbers of points in each quadrat)用
pnt <- ppp(px, py, c(0,1), c(0,1))
plot(pnt, type="p")
plot(density(pnt), 0.1)
contour(density(pnt), add=T)
plot(pnt, type="p", add=T)

# overlay(Create a new Raster* object, based on two or more Raster* objects)用
pnt_xy <- cbind(px, py)
pnt_sp <- SpatialPoints(data.frame(px, py))
pnt_spdf <- SpatialPointsDataFrame(pnt_sp, pz)

# 正方領域で分割!
quadratcount(pnt, nx=4, ny=4)  #x,y方向共に4分割してできた計16地区について人口数計算
quadratcount(pnt, nx=4, ny=4)/(0.25^2)  #人口密度計算
plot(pnt, type="p")  #各地区での人口可視化
plot(quadratcount(pnt, nx=4, ny=4), add=T, col="red")  #各地区での人口密度可視化

#複雑な形状の地区で分割!
poly1 <- cbind(  #地区の頂点のx座標とy座標を指定(poly1~poly4)
c(0, 0.75, 0.75, 0.5, 0.5, 0.25, 0.25, 0, 0),
c(0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 0))
poly2 <- cbind(
c(0.25, 0.75, 0.75, 0.25, 0.25),
c(0.5, 0.5, 0.75, 0.75, 0.5))
poly3 <- cbind(
c(0.75, 1, 1, 0.75, 0.75, 0.5, 0.5, 0.75, 0.75),
c(0, 0, 0.75, 0.75, 0.5, 0.5, 0.25, 0.25, 0))
poly4 <- cbind(
c(0, 1, 1, 0, 0),
c(0.75, 0.75, 1, 1, 0.75))
poly1_pl <- Polygons(list(Polygon(poly1)), "poly1")
poly2_pl <- Polygons(list(Polygon(poly2)), "poly2")
poly3_pl <- Polygons(list(Polygon(poly3)), "poly3")
poly4_pl <- Polygons(list(Polygon(poly4)), "poly4")
poly_sp <- SpatialPolygons(list(poly1_pl, poly2_pl, poly3_pl, poly4_pl))
poly_spdf <- SpatialPolygonsDataFrame(poly_sp, data.frame(c(1:4),
row.names=c("poly1", "poly2", "poly3", "poly4")))  #地区データまとめ作成
plot(poly_spdf)
plot(pnt, type="p", add=T)
overlay(pnt_spdf, poly_spdf, fn=colSums)
area.poly_spdf <- c(0.25^2*6, 0.25^2*2, 0.25^2*4, 0.25^2*4)  #各地区面積
overlay(pnt_spdf, poly_spdf, fn=colSums)/area.poly_spdf


#### 3.2 属性値の基本統計量と標準化 ####
library(spdep)
library(maptools)
library(classInt)

#地価分布データを可視化しよう!
lph <- read.table("lph2010.csv", sep=",", header=TRUE)  #横浜市内の住宅地地価分布データの読み込み
x <- lph$Easting
y <- lph$Northing
z <- as.data.frame(cbind(lph$JCODE, lph$lph2010))  #座標データ設定
colnames(z) <- c("JCODE", "lph2010")
lph_sp <- SpatialPoints(data.frame(x, y))
lph_spdf <- SpatialPointsDataFrame(lph_sp, z)  #SpatialPointsクラス設定
yoko_spdf <- lph_spdf[lph_spdf$JCODE>14000 & lph_spdf$JCODE<14200,]  #範囲指定
pal0 <- c("grey80","black")   #色指定
q_lph <- classIntervals(round(yoko_spdf$lph2010/10000,1), n=4, style="fisher")
q_lph_Col <- findColours(q_lph,pal0)
plot(yoko_spdf,col=q_lph_Col, pch=20, cex=1.8, axes=TRUE, cex.axis=1.5)  #地価分布データ可視化
legend("topright",fill=attr(q_lph_Col,"palette"),
legend=names(attr(q_lph_Col,"table")), cex=1.3, bty="n")

# 平均・分散・標準編纂
mean(lph$lph2010)  # 標本平均
var(lph$lph2010) * (length(lph$lph2010)-1)/length(lph$lph2010)  # 標本分散
sqrt(var(lph$lph2010) * (length(lph$lph2010)-1)/length(lph$lph2010))  # 標本標準偏差
var(lph$lph2010)  # 不偏分散
sd(lph$lph2010)  # 不偏標準偏差

#標準化:平均0,標準偏差1の分布へ
lph$lph2010.scale <- scale(lph$lph2010)  # 標準化
summary(lph$lph2010.scale)
hist(lph$lph2010.scale, col="grey", xlim=c(-2,4), ylim=c(0,4000), breaks=48,
main="", xlab="標準化後の地価", ylab="頻度", cex.axis=1.3, cex.lab=1.2)
hist(scale(yoko_spdf$lph2010), col="grey", breaks=36,
xlim=c(-2,6), ylim=c(0,300), main="", xlab="標準化後の地価", ylab="頻度",
cex.lab=1.2, cex.axis=1.3)  #標準化後の地価分布表示

#歪度:分布の対称性 /尖度:分布の裾の重さ
library(e1071)
skewness(lph$lph2010.scale)  #歪度
kurtosis(lph$lph2010.scale)  #尖度
mean((lph$lph2010.scale-mean(lph$lph2010.scale))^3)/(sd(lph$lph2010.scale)^3)　 #歪度の計算
mean((lph$lph2010.scale-mean(lph$lph2010.scale))^4)/(sd(lph$lph2010.scale)^4)  #尖度の計算


#### 3.3 地域属性の差の比較 ####
## SPM(浮遊状物質)データを分析しよう! ##

library(RColorBrewer)
spm.shp <- readShapePoints("tma_spm.shp")

#地域内でのSPM分布可視化
# カラーパレットを指定する方法
pal0 <- c("grey","grey2")
q_spm <- classIntervals(spm.shp$SPM07, n=5, style="quantile")
q_spm_Col <- findColours(q_spm,pal0)
plot(spm.shp, col=q_spm_Col, pch=20, cex=1.5)
legend("topright",fill=attr(q_spm_Col,"palette"),
legend=names(attr(q_spm_Col,"table")),cex=1.2, bty="n")
# spplotを使う方法
spplot(spm.shp, c("SPM07"), cex=1.5, col.regions=brewer.pal(5, "Greys"))

tma_spm <- read.table("tma_spm.csv", sep=",", header=TRUE)
hist(tma_spm$SPM07)
plot(density(tma_spm$SPM07))

# 地区別集計
tma_spm$KCODE <- floor(tma_spm$ID/1000000)
tapply(tma_spm$SPM07, factor(tma_spm$KCODE), length)
tapply(tma_spm$SPM07, factor(tma_spm$KCODE), mean)
tapply(tma_spm$SPM07, factor(tma_spm$KCODE), var)
tapply(tma_spm$SPM07, factor(tma_spm$KCODE), sd)

# 都県別に公示地価地点データを抽出
spm_11 <- tma_spm[tma_spm$KCODE==11,]
spm_12 <- tma_spm[tma_spm$KCODE==12,]
spm_13 <- tma_spm[tma_spm$KCODE==13,]
spm_14 <- tma_spm[tma_spm$KCODE==14,]

# コルモゴロフ・スミルノフ検定:異なる2標本のデータの分布の一致性を調べる検定
ks.test(spm_11$SPM07, "pnorm", mean(spm_11$SPM07), sd(spm_11$SPM07))
ks.test(spm_12$SPM07, "pnorm", mean(spm_12$SPM07), sd(spm_12$SPM07))
ks.test(spm_13$SPM07, "pnorm", mean(spm_13$SPM07), sd(spm_13$SPM07))
ks.test(spm_14$SPM07, "pnorm", mean(spm_14$SPM07), sd(spm_14$SPM07))

# SPA観測値分布ヒストグラム作図
hist(spm_13$SPM07, col="grey", xlim=c(0.04, 0.10), ylim=c(0,25),
xlab="SPM(mg/m3)", ylab="頻度", main="", cex.lab=1.2, cex.axis=1.3)
# SPA観測値分布累積分布関数作図
plot(ecdf(spm_13$SPM07), do.point=FALSE, verticals=TRUE, main="",
lwd=2, cex.axis=1.3, cex.lab=1.2)
z <- seq(0.04, 0.09, by=0.001)
lines(z, pnorm(z, mean=mean(spm_13$SPM07), sd=sd(spm_13$SPM07)), lty=2, lwd=2)

# 分散の差の検定
var.test(spm_13$SPM07, spm_14$SPM07)
# 平均値の差の検定（等分散）
t.test(spm_13$SPM07, spm_14$SPM07, var.equal=T)
# 平均値の差の検定（不等分散）
t.test(spm_13$SPM07, spm_14$SPM07, var.equal=F)

# Wilcoxonの順位和検定:標本数が少なかったり正規分布に従わない場合に使うノンパラメトリック検定
wilcox.test(spm_13$SPM07, spm_14$SPM07)

## 地価データをベイズ的方法で比較 ##
## 市区町村のようなより小さい空間集計単位での比較には、
##  マルコフ連鎖モンテカルロ法で標本数が多い場合の事後分布を推定し,
##    平均値に差があるとした場合の確率を比較するノンパラメトリック手法を用いる


lph_14201 <- lph[lph$JCODE==14201,]  # 横須賀市
lph_14204 <- lph[lph$JCODE==14204,]  # 鎌倉市
lph_14208 <- lph[lph$JCODE==14208,]  # 逗子市
lph_14210 <- lph[lph$JCODE==14210,]  # 三浦市
lph_14301 <- lph[lph$JCODE==14301,]  # 葉山町
length(lph_14201$lph2010)
length(lph_14204$lph2010)
length(lph_14208$lph2010)
length(lph_14210$lph2010)
length(lph_14301$lph2010)

# 確率密度関数
plot(density(lph_14201$lph2010), xlim=c(0, 5e+5), ylim=c(0,3.0e-5),
main="", xlab="地価(円/㎡)", ylab="密度", lty=1, lwd=2,
cex.axis=1.3, cex.lab=1.2)
lines(density(lph_14204$lph2010), lty=2, lwd=2)
lines(density(lph_14208$lph2010), lty=3, lwd=2)
lines(density(lph_14210$lph2010), lty=4, lwd=2)
lines(density(lph_14301$lph2010), lty=5, lwd=2)
rect(320000, 2e-5, 500000, 3e-5,col="white")
segments(350000, 2.9e-5, 400000, 2.9e-5, lty=1, lwd=2)
text(450000,2.9e-5, "横須賀市", cex=1.3)
segments(350000, 2.7e-5, 400000, 2.7e-5, lty=2, lwd=2)
text(450000,2.7e-5, "鎌倉市", cex=1.3)
segments(350000, 2.5e-5, 400000, 2.5e-5, lty=3, lwd=2)
text(450000,2.5e-5, "逗子市", cex=1.3)
segments(350000, 2.3e-5, 400000, 2.3e-5, lty=4, lwd=2)
text(450000,2.3e-5, "三浦市", cex=1.3)
segments(350000, 2.1e-5, 400000, 2.1e-5, lty=5, lwd=2)
text(450000,2.1e-5, "葉山町", cex=1.3)
# 箱ひげ図
boxplot(lph_14201$lph2010, lph_14204$lph2010, lph_14208$lph2010,
lph_14210$lph2010, lph_14301$lph2010,
names=c("横須賀市","鎌倉市","逗子市","三浦市","葉山町"),
cex.axis=1.3, lwd=2)

# 2群間の平均
# 逗子市と葉山町
y1 <- lph_14208$lph2010
y2 <- lph_14301$lph2010
n1<-length(y1) ; n2<-length(y2)
# prior parameters
mu0<-mean(c(y1,y2)) ; g02<-var(y1)
del0<-0 ; t02<-var(y1)
s20<-15; nu0<-1
# starting values
mu<- ( mean(y1) + mean(y2) )/2
del<- ( mean(y1) - mean(y2) )/2
# Gibbs sampler
MU<-DEL<-S2<-NULL
Y12<-NULL
set.seed(1)
for(s in 1:10000)
{
 ##update s2
 s2<-1/rgamma(1,(nu0+n1+n2)/2,
       (nu0*s20+sum((y1-mu-del)^2)+sum((y2-mu+del)^2) )/2)
 ##update mu
 var.mu<-  1/(1/g02+ (n1+n2)/s2 )
 mean.mu<- var.mu*( mu0/g02 + sum(y1-del)/s2 + sum(y2+del)/s2 )
 mu<-rnorm(1,mean.mu,sqrt(var.mu))
 ##update del
 var.del<-  1/(1/t02+ (n1+n2)/s2 )
 mean.del<- var.del*( del0/t02 + sum(y1-mu)/s2 - sum(y2-mu)/s2 )
 del<-rnorm(1,mean.del,sqrt(var.del))
 ##save parameter values
 MU<-c(MU,mu) ; DEL<-c(DEL,del) ; S2<-c(S2,s2)
 Y12<-rbind(Y12,c(rnorm(2,mu+c(1,-1)*del,sqrt(s2) ) ) )
}

# 住宅地地価の事後分布グラフ作図
plot(density(Y12[,1]), col="black", lwd=3, cex.axis=1.3, cex.lab=1.2,
main="", xlab="地価（円/㎡）", ylab="密度")
lines(density(Y12[,2]), col="grey", lwd=3)
text(3e+5,1e-5,"逗子市", cex=1.3)
text(3e+5,9e-6,"葉山町", cex=1.3)
lines(x=c(2.5e+5,2.75e+5),y=c(1e-5,1e-5), lwd=3, col="black")
lines(x=c(2.5e+5,2.75e+5),y=c(9e-6,9e-6), lwd=3, col="grey")


#### 3.4 地域間格差 ####
## 3.4.1 ジニ係数:1寄りなら格差大、0寄りなら平等であることを示す尺度
library(reldist)
gini(lph$lph2010)
inc <- read.table("inc2006.csv", sep=",", header=TRUE)
library(reldist)
gini(inc$INC2006)
library(ineq)
inc.gini <- Lc(inc$INC2006)
inc.gini
# ローレンツ曲線作図:45度線と比較して格差調査
plot(inc.gini, cex.axis=1.3, cex.lab=1.2, cex.main=1.3, lwd=3)

## 3.4.2 地域変動係数:大きいほど格差大であることを示す尺度
Cv <- sd(lph$lph2010)/mean(lph$lph2010)
Cv

## 3.4.3 地域特化係数:1より大きければ、地域iが産業部門kに特化していることを示す尺度
CLpop05 <- read.table("CLpop05.csv", sep=",", header=TRUE)
CLpop05$pop05CL3 <- (CLpop05$POP3 / CLpop05$POP05) / (sum(CLpop05$POP3) /  sum(CLpop05$POP05))

jpn_pref <- readShapePoly("jpn_pref.shp",IDvar="PREF_CODE")
plot(jpn_pref, col="grey")
ID.match <- match(jpn_pref$PREF_CODE, CLpop05$PREFCODE)
ID.match
jpn_CL <-CLpop05[ID.match,]
summary(jpn_CL)
jpn_pref_CL<- spCbind(jpn_pref,jpn_CL)
names(jpn_pref_CL)
summary(jpn_pref_CL)

# カラーパレットの作成
library(maptools)
library(classInt)
pal0 <- c("white","grey","grey2")
# 等量分類
q_pref <- classIntervals(round(jpn_pref_CL$pop05CL3,2), n=5, style="quantile")
plot(q_pref, pal=pal0)
q_pref_Col <- findColours(q_pref,pal0)
plot(jpn_pref_CL,col=q_pref_Col)
#title("地域特化係数（65歳以上人口）")
legend("topleft", cex=1.5, fill=attr(q_pref_Col,"palette"),
legend=names(attr(q_pref_Col,"table")),bty="n")
