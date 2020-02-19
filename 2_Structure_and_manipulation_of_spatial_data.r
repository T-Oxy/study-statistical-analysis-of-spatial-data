##spdepパッケージ:空間統計の定番パッケージ
library(spdep)

#### 空間データ生成 ####
library(gpclib)
# 2.1.1 点データ(ポイントデータ)作成
px <- runif(30,0,10)
py <- runif(30,0,10)
pz <- as.data.frame(px+py)
colnames(pz) <- c("pz")
pnt_xy <- cbind(px, py)
pnt_sp <- SpatialPoints(data.frame(px, py))
pnt_spdf <- SpatialPointsDataFrame(pnt_sp, pz)
plot(pnt_spdf, cex=2, pch=19, lwd=3)

# 2.1.2 線データ(ラインデータ)作成
line1 <- cbind(runif(4,0,10), runif(4,0,10))
line2 <- cbind(runif(4,0,10), runif(4,0,10))
line3 <- cbind(runif(4,0,10), runif(4,0,10))
line1_ln <- Lines(list(Line(line1)), "line1")
line2_ln <- Lines(list(Line(line2)), "line2")
line3_ln <- Lines(list(Line(line3)), "line3")
line_sp <- SpatialLines(list(line1_ln, line2_ln, line3_ln))
line_spdf <- SpatialLinesDataFrame(line_sp, data.frame(c(1:3),
row.names=c("line1", "line2", "line3")))
plot(line_spdf, lty=1, lwd=3)

# 2.1.3 面データ(ポリゴンデータ)作成
yoko <- readShapePoly("yoko.shp")
yoko_coord <- coordinates(yoko)
plot(yoko, col="grey", border="black", lwd=3)

# 2.1.4 グリッドデータ(ラスターデータ,ピクセルデータ)作成:格子状に空間分割
grid_topo <- GridTopology(cellcentre.offset=c(1, 1),
cellsize=c(2, 2), cells.dim=c(5, 5))
grid_sgdf <- SpatialGridDataFrame(grid_topo,
data=as.data.frame(runif(25, min=1, max=50)))
# image(grid_sgdf, col=grey((0:50)/51))
# グリッドデータの面オブジェクトへの変換
sp_grd <- as.SpatialPolygons.GridTopology(grid_topo)

plot(sp_grd, lwd=3)
text(getSpPPolygonsLabptSlots(sp_grd),
getSpPPolygonsIDSlots(sp_grd))
sp_spdf<- SpatialPolygonsDataFrame(sp_grd,
data=data.frame(c(1:25),
row.names=sapply(slot(sp_grd, "polygons"),
function(i) slot(i, "ID"))))
spplot(sp_spdf)


#### 単一レイヤでの操作 ####
# 2.2.1 面積計算
poly3_gpc <- as(poly3, "gpc.poly")
area.poly(poly3_gpc)

# 2.2.2 セントロイド(ポリゴンデータの重心)
library(geosphere)
plot(yoko, col="white", border="grey", lwd=3)
points(centroid(yoko), lwd=3, pch=19)

# 2.2.3 距離計算
x <- c(5, 2, 6, 8,10)
y <- c(6, 4,10, 0, 7)
xy <- t(rbind(x,y)) #ポイントデータの座標
plot(xy,xlim=c(0,10),ylim=c(0.10),cex=2,pch=19)
abline(h=0:10,v=0:10,lty=2)
dist(xy,method="euclidean")  #直線距離
dist(xy,method="manhattan")  #マンハッタン距離:地点間の座標の差の絶対値

# 2.2.4 属性テーブルの結合...2テーブル間のIDわマッチングさせて新テーブル作成
library(maptools)
jpn_pref <- readShapePoly("jpn_pref.shp",IDvar="PREF_CODE")
jpn_COD <- read.table("COD.csv",sep=",",header=TRUE)
ID.match <- match(jpn_pref$PREF_CODE, jpn_COD$PREF_CODE)
ID.match
jpn_COD1 <- jpn_COD[ID.match,]
summary(jpn_COD1)
jpn_pref_COD <- spCbind(jpn_pref,jpn_COD1)
names(jpn_pref_COD)
summary(jpn_pref_COD)

# 2.2.5 ディゾルブ:同じ属性を持つオブジェクトやグリッドを単一のオブジェクトに統合
library(maptools)
jpn <- readShapeSpatial("japan_ver62.shp", IDvar="JCODE",
proj4string=CRS("+proj=tmerc +ellps=GRS80 +units=m +datum=JGD2000"))
library(gpclib)
gpclibPermit()
jpn1 <- jpn
pref1 <- unionSpatialPolygons(jpn1, round(jpn1$JCODE / 1000))

# 2.2.6 サブセットの抽出
yoko <- jpn[jpn$JCODE > 14100,]
yoko <- yoko[yoko$JCODE < 14120,]
# writeSpatialShape(yoko, "yoko")
# 図2.14
plot(jpn, xlim=c(139.4, 139.8), ylim=c(35.02, 35.8), col="grey95", border="grey50",lwd=2)
plot(yoko, col="grey", border="black", add=TRUE, lwd=3)

# 2.2.7 バッファリング例:任意の点データから一定距離の範囲内に圏域作成,オブジェクト近傍領域定義法
library(geosphere)
library(sp)
library(maptools)
library(spatstat)  #以上、パッケージの読み込み
x <- c(5, 2, 6, 8,10)
y <- c(6, 4,10, 0, 7)
temp <- as.data.frame(cbind(x,y))
colnames(temp) <- c("x", "y")
polys<-list()
for(i in 1:nrow(temp)) {
 discbuff<-disc(radius=1, centre=c(temp$x[i], temp$y[i]))
    discpoly<-Polygon(rbind(cbind(discbuff$bdry[[1]]$x,
y=discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1],
y=discbuff$bdry[[1]]$y[1])))
    polys<-c(polys, discpoly)
}
spolys<-list()
for(i in 1:length(polys)) {
 spolybuff<-Polygons(list(polys[[i]]), ID=row.names(temp)[i])
 spolys<-c(spolys, spolybuff)
}
polybuffs<-SpatialPolygons(spolys)

plot(polybuffs, col="grey", lwd=2)
points(temp, col="black", pch=20, cex=3)

# 2.2.8 ボロノイ分割:隣接する任意の2点を結ぶ線分の垂直二等分線で構成されるポリゴンで領域分割,オブジェクト最近隣領域定義法
x <- c(5, 2, 6, 8,10)
y <- c(6, 4,10, 0, 7)
plot(deldir(x,y), axes=FALSE, cex=3)

#### 複数レイヤでの操作 ####
# 2.3.1 ポリゴンデータ内のポイントデータの属性の集計
overlay(pnt_spdf, sp_spdf)
overlay(pnt_spdf, sp_spdf, fn=colSums)  #合計
overlay(pnt_spdf, sp_spdf, fn=mean)  #平均
overlay(pnt_spdf, sp_spdf, fn=sd)  #標準偏

# 2.3.2 面オブジェクトに対するオーバーレイ操作
library(gpclib)
##オーバーレイする面オブジェクト
p1 <- cbind(c(0, 5, 5, 0, 0), c(0, 0, 5, 5, 0))
p2 <- cbind(c(5, 7.5, 5, 2.5, 5), c(2.5, 5, 7.5, 5, 2.5))
##gpc.polyに変換
p1_gpc <- as(p1, "gpc.poly")
p2_gpc <- as(p2, "gpc.poly")
##作図
plot(p1_gpc, poly.args = list(lwd=8), axes=FALSE, xlab="", ylab="")
plot(p2_gpc, poly.args = list(lwd=8), axes=FALSE, xlab="", ylab="")
###
plot(append.poly(p1_gpc, p2_gpc), poly.args = list(lwd=8), axes=FALSE,
xlab="", ylab="")
### erase -> setdiff() :削除
poly_diff <- setdiff(p1_gpc, p2_gpc)
plot(append.poly(p1_gpc, p2_gpc), poly.args = list(lwd=8), axes=FALSE,
xlab="", ylab="")
plot(poly_diff, poly.args=list(col="grey", lwd=8), add=TRUE)
###union -> union() :論理和
poly_union <- union(p1_gpc, p2_gpc)
plot(append.poly(p1_gpc, p2_gpc), poly.args = list(lwd=8), axes=FALSE,
xlab="", ylab="")
plot(poly_union, poly.args=list(col="grey", lwd=8), add=TRUE)
###intersect -> intersect() :論理積
poly_intersect <- intersect(p1_gpc, p2_gpc)
plot(append.poly(p1_gpc, p2_gpc), poly.args = list(lwd=8), axes=FALSE,
xlab="", ylab="")
plot(poly_intersect, poly.args=list(col="grey", lwd=8), add=TRUE)
###identity -> append.poly() :順位づけの積
poly_append <- append.poly(p1_gpc, p2_gpc)
plot(poly_append, poly.args=list(col="grey", lwd=8), axes=FALSE,
xlab="", ylab="")
#plot(append.poly(p1_gpc, p2_gpc), poly.args=list(col="grey"))
###update ：更新
poly_update <- append.poly(setdiff(p1_gpc, p2_gpc), intersect(p1_gpc, p2_gpc))
plot(append.poly(p1_gpc, p2_gpc), poly.args=list(lwd=8), axes=FALSE,
xlab="", ylab="")
plot(poly_update, poly.args=list(col="grey", lwd=8), add=TRUE)
