# 全体を通じて使用するパッケージ
library(spdep)

#### 5.1　空間隣接行列 ####
# 例１:横浜市
# 地図データの読み込み
yoko<- readShapePoly("yoko.shp", IDvar="JCODE")
# 座標テーブルの作成
yoko_coords <- coordinates(yoko)
# ドロネー三角網図
yoko.tri.nb <- tri2nb(yoko_coords)
plot(yoko, border="white", col="grey")
plot(yoko.tri.nb, yoko_coords, add=TRUE)
# 再近隣k地点
yoko.knn4 <- knearneigh(yoko_coords, k=4)
yoko.knn4.nb <- knn2nb(yoko.knn4, row.names=rownames(yoko$JCODE))
plot(yoko, border="white", col="grey")
plot(yoko.knn4.nb, yoko_coords, add=TRUE)
# k=4とk=3の違い
yoko.knn3 <- knearneigh(yoko_coords, k=3)
yoko.knn3.nb <- knn2nb(yoko.knn3, row.names=rownames(yoko$JCODE))
diffs <- diffnb(yoko.knn3.nb, yoko.knn4.nb)
plot(yoko, border="white", col="grey")
plot(yoko.knn3.nb, yoko_coords, lwd=3, add=TRUE)
plot(diffs, yoko_coords, col="red", lty=2, lwd=2, add=TRUE)
# 距離により隣接関係を定義
yoko.r.nb <- dnearneigh(yoko_coords, 0, 0.06)
plot(yoko, border="white", col="grey")
plot(yoko.r.nb, yoko_coords, add=TRUE)
# 面オブジェクトの隣接関係
yoko.poly.nb <- poly2nb(yoko)
plot(yoko, border="white", col="grey")
plot(yoko.poly.nb, yoko_coords, add=TRUE)

# 例２:関東地方全域
# 地図データの読み込み
kanto <- readShapePoly("kanto_area.shp", IDvar="JCODE")
# 座標テーブルの作成
coords <- t(sapply(slot(kanto, "polygons"),function(x) slot(x, "labpt")))
# ドロネー三角網図
kanto.tri.nb <- tri2nb(coords, row.names=rownames(kanto$JCODE))
plot(kanto)
plot(kanto.tri.nb, coords, add=T)
# 再近隣k地点
kanto.knn4 <- knearneigh(coords, k=4)
kanto.knn4.nb <- knn2nb(kanto.knn4, row.names=rownames(kanto$JCODE))
plot(kanto)
plot(kanto.knn4.nb, coords, add=TRUE)
# k=4とk=3の違い
kanto.knn3 <- knearneigh(coords, k=3)
kanto.knn3.nb <- knn2nb(kanto.knn3, row.names=rownames(kanto$JCODE))
diffs <- diffnb(kanto.knn3.nb, kanto.knn4.nb)
plot(kanto)
plot(kanto.knn4.nb, coords, add=TRUE)
plot(diffs, coords, add=TRUE, col="red", lty=2)
# 距離により隣接関係を定義
kanto.r.nb <- dnearneigh(coords, 0, 0.05, row.names=rownames(kanto$JCODE))
plot(kanto)
plot(kanto.r.nb, coords, add=TRUE)
# 面オブジェクトの隣接関係
kanto <- readShapePoly("kanto_area.shp", IDvar="JCODE")
kanto_coords <- t(sapply(slot(kanto, "polygons"),function(x) slot(x, "labpt")))
kanto.poly.nb <- poly2nb(kanto)
plot(kanto)
plot(kanto.poly.nb, kanto_coords, add=TRUE)

#### 5.2　空間重み付け行列 ####

#### 5.3　空間的自己相関 ####
# 関東圏地価データの分析
lph <- read.table("lph.csv", sep=",", header=TRUE)
summary(lph)
kanto <- readShapePoly("kanto_area.shp", IDvar="JCODE")
# ドロネー三角網図
coords <- matrix(0, nrow=nrow(lph), ncol=2)
coords[,1] <- lph$Easting
coords[,2] <- lph$Northing
lph.tri.nb <- tri2nb(coords)
# Moran's I
moran.test(lph$LPH, nb2listw (lph.tri.nb,style="W"))
moran.test(lph$LPH, nb2listw (lph.tri.nb,style="B"))
moran.test(lph$LPH, nb2listw (lph.tri.nb,style="C"))
moran.test(lph$LPH, nb2listw (lph.tri.nb,style="S"))

moran.test(lph$LPH, nb2listw (lph.tri.nb,style="W"))
moran.test(lph$POPD, nb2listw (lph.tri.nb,style="W"))
moran.test(lph$EMP3D, nb2listw (lph.tri.nb,style="W"))
# Geary's C
geary.test(lph$LPH, nb2listw (lph.tri.nb,style="W"))
# Join Count統計量
lph.hi.low <- cut(lph$LPH,
breaks=c(0,mean(lph$LPH),max(lph$LPH)),
labels=c("low","high"))
names(lph.hi.low) <- lph$JCODE
joincount.multi(lph.hi.low, nb2listw(lph.tri.nb, style="B"))
joincount.test(lph.hi.low, nb2listw(lph.tri.nb, style="B"))
by(card(lph.tri.nb), lph.hi.low, summary)
# 都道府県別死亡者数データの分析
# データの読み込み
pref.pnt <- readShapePoints("pref_gov.shp")
jpn.pref <- readShapePoly("jpn_pref.shp")
plot(jpn.pref)
points(pref.pnt)
pref_gov <- read.table("pref_gov.txt",sep=",",header=TRUE,row.names=2)
coords <- matrix(0,nrow(pref_gov),2)
coords[,1] <- pref_gov$X
coords[,2] <- pref_gov$Y
# ドロネー三角網図
pref.tri.nb <- tri2nb(coords,row.names=rownames(pref_gov))
plot(pref.tri.nb, coords, add=T)
# Moran's I
moran.test(pref.pnt$diabetes, nb2listw (pref.tri.nb,style="W"))
moran.test(pref.pnt$diabetes, nb2listw (pref.tri.nb,style="B"))
moran.test(pref.pnt$diabetes, nb2listw (pref.tri.nb,style="C"))
moran.test(pref.pnt$diabetes, nb2listw (pref.tri.nb,style="S"))
#
moran.test(pref.pnt$geriatric, nb2listw (pref.tri.nb,style="W"))
moran.test(pref.pnt$malignant, nb2listw (pref.tri.nb,style="W"))
moran.test(pref.pnt$hypertensi, nb2listw (pref.tri.nb,style="W"))
# Geary's C
geary.test(pref.pnt$diabetes, nb2listw(pref.tri.nb, style="W"))
geary.test(pref.pnt$geriatric, nb2listw(pref.tri.nb, style="W"))
geary.test(pref.pnt$malignant, nb2listw(pref.tri.nb, style="W"))
geary.test(pref.pnt$hypertensi, nb2listw(pref.tri.nb, style="W"))
# Join Count
diabates.hi.low <- cut(pref.pnt$diabetes,
breaks=c(0,mean(pref.pnt$diabetes),max(pref.pnt$diabetes)),
labels=c("low","high"))
names(diabates.hi.low) <- pref.pnt$KENCODDE
joincount.test(diabates.hi.low, nb2listw(pref.tri.nb, style="W"))
by(card(pref.tri.nb), diabates.hi.low, summary)
# Local Moran's I
LMI1 <- localmoran(pref.pnt$diabetes, nb2listw(pref.tri.nb, style="W"))
pref.lm <- data.frame(cbind(LMI1[,1],(pref.pnt$diabetes-
mean(pref.pnt$diabetes))/sd(pref.pnt$diabetes)),row.names=pref.pnt$KENCODE)
colnames(pref.lm) <- c("Ii","standardized diabetes")
pref.lm
# 図5.11
plot(pref.lm,xlab="Local Moran's I", ylab="死亡率（標準化済み）", pch=20,
cex=1.3, cex.axis=1.3, cex.lab=1.2)
lines(x=c(-1,2), y=c(0,0), col="grey", lwd=2)
lines(x=c(0,0), y=c(-3,4), col="grey", lwd=2)
# text(pref.lm,rownames(pref.lm), adj=1.2, cex=0.8)
text(1.5, 3.3, "徳島県", cex=1.2)
text(1.16, 1.6, "香川県", cex=1.2)
text(-0.32, 2.1, "大分県", cex=1.2)
text(-0.6, -1.65, "愛知県", cex=1.2)
text(-0.72, 1.7, "富山県", cex=1.2)
text(0.9, -0.57, "東京都", cex=1.2)
text(0.85, -1.91, "神奈川県", cex=1.2)
# 図5.12
library(classInt)
jpn.pref <- readShapePoly("jpn_pref.shp",IDvar="PREF_CODE")
LMI <- as.data.frame(cbind(jpn.pref$PREF_CODE, pref.lm$Ii))
colnames(LMI) <- c("PREF_CODE","Ii")
ID.match <- match(jpn.pref$PREF_CODE, LMI$PREF_CODE)
ID.match
jpn_LMI<- LMI[ID.match,]
summary(jpn_LMI)
jpn_pref_LMI<- spCbind(jpn.pref,jpn_LMI)
pal <- c("white","grey","grey2")
fj_LMI <- classIntervals(round(jpn_pref_LMI$Ii,2), n=5, style="fisher")
plot(fj_LMI, pal=pal)
fj_LMI_Col <- findColours(fj_LMI,pal)
plot(jpn_pref_LMI,col=fj_LMI_Col)
# title("Diabetes mellitus (local Moran's I)")
legend("topleft",fill=attr(fj_LMI_Col,"palette"),
legend=names(attr(fj_LMI_Col,"table")), cex=1.3, bty="n")
# local G
localG.diabates <- localG(pref.pnt$diabetes, nb2listw(pref.tri.nb, style="W"))
localGs.diabates <- localG(pref.pnt$diabetes,
nb2listw(include.self(pref.tri.nb), style="W"))
library(classInt)
jpn.pref <- readShapePoly("jpn_pref.shp",IDvar="PREF_CODE")
localG.mat <- as.data.frame(cbind(jpn.pref$PREF_CODE, localG.diabates))
colnames(localG.mat) <- c("PREF_CODE","localG")
ID.match <- match(jpn.pref$PREF_CODE, localG.mat$PREF_CODE)
ID.match
jpn_localG<- localG.mat[ID.match,]
summary(jpn_localG)
jpn_pref_localG<- spCbind(jpn.pref,jpn_localG)
pal <- c("white", "grey","grey2")
fj_localG <- classIntervals(round(jpn_pref_localG$localG,2), n=5,
style="fisher")
plot(fj_localG, pal=pal)
fj_localG_Col <- findColours(fj_localG,pal)
plot(jpn_pref_localG,col=fj_localG_Col)
# title("Diabetes mellitus (localG)")
legend("topleft",fill=attr(fj_localG_Col,"palette"),
legend=names(attr(fj_localG_Col,"table")), cex=1.3, bty="n")
