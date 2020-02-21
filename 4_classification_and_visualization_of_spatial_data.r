##パッケージ
library(maptools)
library(classInt)

# 地図データの読み込み
jpn_pref <- readShapePoly("jpn_pref.shp",IDvar="PREF_CODE")
plot(jpn_pref, col="grey")
summary(jpn_pref)
# 地図属性テーブルと空間データテーブルとの結合
jpn_COD <- read.table("COD.csv",sep=",",header=TRUE)
summary(jpn_COD)
ID.match <- match(jpn_pref$PREF_CODE, jpn_COD$PREF_CODE)
ID.match
jpn_COD1 <- jpn_COD[ID.match,]
summary(jpn_COD1)
jpn_pref_COD <- spCbind(jpn_pref,jpn_COD1)
names(jpn_pref_COD)
summary(jpn_pref_COD)
# カラーパレットの作成
pal0 <- c("white","grey","grey2")

#### 4.1 等量分類 ####
q_pref <- classIntervals(round(jpn_pref_COD$malignant,2), n=5, style="quantile")
plot(q_pref, pal=pal0, cex.axis=1.3, cex.lab=1.2, lwd=2, main="",
xlab="悪性新生物による死亡者数(人口10万人あたり)")
q_pref_Col <- findColours(q_pref,pal0)
plot(jpn_pref_COD,col=q_pref_Col)
title("悪性新生物による死亡者数(人口10万人あたり) (等量分類)", cex=1.4)
legend("topleft",fill=attr(q_pref_Col,"palette"), cex=1.4,
legend=names(attr(q_pref_Col,"table")),bty="n")

#### 4.2 等間隔分類 ####
eq_pref <- classIntervals(round(jpn_pref_COD$malignant,2), n=5, style="equal")
plot(eq_pref, pal=pal0, cex.axis=1.3, cex.lab=1.2, lwd=2, main="",
xlab="悪性新生物による死亡者数(人口10万人あたり)")
eq_pref_Col <- findColours(eq_pref,pal0)
plot(jpn_pref_COD,col=eq_pref_Col)
title("悪性新生物による死亡者数(人口10万人あたり) (等間隔分類)", cex=1.4)
legend("topleft",fill=attr(eq_pref_Col,"palette"), cex=1.4,
legend=names(attr(eq_pref_Col,"table")),bty="n")

#### 4.3 標準偏差分類 ####
sd_pref <- classIntervals(round(jpn_pref_COD$malignant,2), n=5, style="sd")
plot(sd_pref, pal=pal0, cex.axis=1.3, cex.lab=1.2, lwd=2, main="",
xlab="悪性新生物による死亡者数(人口10万人あたり)")
sd_pref_Col <- findColours(sd_pref,pal0)
plot(jpn_pref_COD,col=sd_pref_Col)
title("悪性新生物による死亡者数(人口10万人あたり) (標準偏差分類)", cex=1.4)
legend("topleft",fill=attr(sd_pref_Col,"palette"), cex=1.4,
legend=names(attr(sd_pref_Col,"table")),bty="n")

#### 4.4 自然階級分類 ####
fj_pref <- classIntervals(round(jpn_pref_COD$malignant,2), style="fisher")
plot(fj_pref, pal=pal0, cex.axis=1.3, cex.lab=1.2, lwd=2, main="",
xlab="悪性新生物による死亡者数(人口10万人あたり)")
fj_pref_Col <- findColours(fj_pref,pal0)
plot(jpn_pref_COD,col=fj_pref_Col)
title("悪性新生物による死亡者数(人口10万人あたり) (Fisher-Jenks法)", cex=1.4)
legend("topleft",fill=attr(fj_pref_Col,"palette"), cex=1.4,
legend=names(attr(fj_pref_Col,"table")),bty="n")

#### 4.5 階級区分を指定した分類 ####
fix_pref <- classIntervals(round(jpn_pref_COD$malignant,2), n=4, style="fixed",
fixedBreaks=c(0, 200, 250, 300, 350))
plot(fix_pref, pal=pal0, cex.axis=1.3, cex.lab=1.2, lwd=2, main="",
xlab="悪性新生物による死亡者数(人口10万人あたり)")
fix_pref_Col <- findColours(fix_pref,pal0)
plot(jpn_pref_COD,col=fix_pref_Col)
title("悪性新生物による死亡者数(人口10万人あたり) (階級区分指定)", cex=1.4)
legend("topleft",fill=attr(fix_pref_Col,"palette"), cex=1.4,
legend=names(attr(fix_pref_Col,"table")),bty="n")

#### 4.6 非階層クラスタリングによる分類 ####
km_pref <- classIntervals(round(jpn_pref_COD$malignant,2), n=5, style="kmeans")
plot(km_pref, pal=pal0, cex.axis=1.3, cex.lab=1.2, lwd=2, main="",
xlab="悪性新生物による死亡者数(人口10万人あたり)")
km_pref_Col <- findColours(km_pref,pal0)
plot(jpn_pref_COD,col=km_pref_Col)
title("悪性新生物による死亡者数(人口10万人あたり) (非階層クラスタリング)",
cex=1.4)
legend("topleft",fill=attr(km_pref_Col,"palette"), cex=1.4,
legend=names(attr(km_pref_Col,"table")),bty="n")

#### 4.7 階層クラスタリングによる分類 ####
hc_pref <- classIntervals(round(jpn_pref_COD$malignant,2), n=5,
style="hclust",  method="complete")
plot(hc_pref, pal=pal0, cex.axis=1.3, cex.lab=1.2, lwd=2, main="",
xlab="悪性新生物による死亡者数(人口10万人あたり)")
hc_pref_Col <- findColours(hc_pref,pal0)
plot(jpn_pref_COD,col=hc_pref_Col)
title("悪性新生物による死亡者数(人口10万人あたり) (階層クラスタリング)",
cex=1.4)
legend("topleft",fill=attr(hc_pref_Col,"palette"), cex=1.4,
legend=names(attr(hc_pref_Col,"table")),bty="n")

#### 階級区分に依存しない分類 ####
pr_pref <- classIntervals(round(jpn_pref_COD$malignant,2), style="pretty")
plot(pr_pref, pal=pal4, main="", cex.axis=1.3, cex.lab=1.2, lwd=2)
pr_pref_Col <- findColours(pr_pref,pal4)
plot(jpn_pref_COD,col=pr_pref_Col)
title("Malignant neoplasms (pretty)", cex=1.4)
legend("topleft",fill=attr(pr_pref_Col,"palette"), cex=1.4,
legend=names(attr(pr_pref_Col,"table")),bty="n")

#### 4.8 ドットマップ ####
# R2.9.0では実行可能
# R2.12.xでは実行できず
# 東北地方を切り出す
#  tohoku_COD <- jpn_pref_COD[jpn_pref_COD$PREF_CODE>=2,]
# tohoku_COD <- tohoku_COD[tohoku_COD$PREF_CODE<=7,]
# データの読み込み
tohoku_COD <- readShapePoly("tohoku_COD.shp")
# 座標の抽出
tohoku_coord <- coordinates(tohoku_COD)
# ドットマップの作成
# 人口5万人に１ドット
tohoku_dots <- dotsInPolys(tohoku_COD, as.integer(tohoku_COD$Pop/50000))
# 図4.8
plot(tohoku_COD, col="grey", border="white", lwd=2)
plot(tohoku_dots, pch=19, cex=1, col="black", add=TRUE)

#### 4.9 比例シンボルマップ ####
# 図4.9
# symbols()関数を使う方法
plot(tohoku_COD, col="grey", border="white", lwd=2)
symbols(x=tohoku_coord[,1], y=tohoku_coord[,2],
circles=tohoku_COD$hypertensi/35, inch=FALSE, bg="black", add=TRUE)
text(x=tohoku_coord[,1], y=tohoku_coord[,2]+0.3, cex=1.3, col="black",
c("青森県","岩手県", "宮城県","秋田県","山形県","福島県"))
# points()のsymbolを使う方法
plot(tohoku_COD, lwd=3)
# points(x=tohoku_coord[,1], y=tohoku_coord[,2],pch=21,
# cex=tohoku_COD$hypertensive/2, bg="grey")
points(x=tohoku_coord[,1], y=tohoku_coord[,2],pch=21,
cex=tohoku_COD$hypertensi/2, bg="black", lwd=0)
text(x=tohoku_coord[,1], y=tohoku_coord[,2]+0.2, cex=1.3,
c("青森県","岩手県", "宮城県","秋田県","山形県","福島県"))

#### 4.10 グラフ表示 ####
# 図4.10
# パッケージTeachingDemosを使用
library(TeachingDemos)
# 棒グラフで表示するデータテーブルを作成
tohoku_COD2 <- cbind(tohoku_COD$Pop_Dens,tohoku_COD$malignant)
# ポリゴンデータを表示
plot(tohoku_COD, lwd=2)
# 棒グラフを表示
for(i in 1:nrow(tohoku_COD)){
 subplot(barplot(tohoku_COD2[i,], yaxt="n", col=c("grey","black")),
 x=tohoku_coord[i,1],y=tohoku_coord[i,2], vadj=0, size=c(0.4,0.6))
}
# 凡例を表示
legend(138, 41.5, c("人口密度", paste("悪性新生物","死亡者数",sep="\n")),
cex=1.3, fill=c("grey","black"), bty="n")


## Tips 
# カラーパレットの作成
pal0 <- c("grey","grey9")
pal1 <- gray.colors(n=5,start=0.9,end=0.3)
pal2 <- rainbow(n=5,start=0.6,end=0.1)
pal3 <- heat.colors(n=5,alpha=1)
pal4 <- topo.colors(n=5,alpha=1)
pal5 <- terrain.colors(n=5,alpha=1)
pal6 <- cm.colors(n=5,alpha=1)
