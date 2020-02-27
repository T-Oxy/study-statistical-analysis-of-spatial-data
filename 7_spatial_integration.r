library(DCluster)

#データの読み込みと作成
data76 <- read.table("data76.csv", sep=",", header=TRUE, row.names=1)
data76_OE <- data.frame(Observed=data76$n.birth)
data76_OE <- cbind(data76_OE,
Expected=data76$pop*sum(data76$n.birth)/sum(data76$pop),
x=data76$Easting, y=data76$Northing)

#### 7.1 ピアソンのχ2検定 ####
achisq.stat(data76_OE, lambda=1)
achisq.test(Observed~offset(log(Expected)), data76_OE, model="poisson", R=100)
data76_achb_pb <- boot(data76_OE, statistic=achisq.pboot, sim="parametric", ran.gen=poisson.sim, R=100)
#図7.1
plot(data76_achb_pb)

#### 7.2 Potthoff & Whittinghillの検定（過分散の有無に関する検定）####
pottwhitt.stat(data76_OE)
pottwhitt.test(Observed~offset(log(Expected)), data76_OE, model="poisson", R=100)
data76_pw_pb <- boot(data76_OE, statistic=pottwhitt.pboot, sim="parametric", ran.gen=poisson.sim, R=100)
plot(data76_pw_pb)

#### 7.3 Stone's test ####
stone.stat(data76_OE, region=which(row.names(data76_OE)=="20"), lambda=1)
stone.test(Observed~offset(log(Expected)), data76_OE, model="poisson", R=100, region=which(row.names(data76_OE)=="20"), lambda=1)
data76_st_pb <- boot(data76_OE, statistic=stone.pboot, sim="parametric", ran.gen=poisson.sim, R=100, region=which(row.names(data76_OE)=="20"))
plot(data76_st_pb)

#### 7.4 Tangoの検定 ####
data76_OE <- cbind(data76_OE, x=data76$Easting, y=data76$Northing)
coords <- as.matrix(data76_OE[,c("x","y")])
dlist <- dnearneigh(coords, 0, Inf)
dlist <- include.self(dlist)
dlist.d <- nbdists(dlist, coords)
col.W.tango <- nb2listw(dlist, glist=lapply(dlist.d, function(x){exp(-x)}), style="C")
tango.stat(data76_OE, col.W.tango, zero.policy=TRUE)
tango.test(Observed~offset(log(Expected)), data76_OE, model="poisson", R=100, list=col.W.tango, zero.policy=TRUE)
data76_tn_pb <- boot(data76_OE, statistic=tango.pboot, sim="parametric", ran.gen=poisson.sim, R=100, listw=col.W.tango, zero.policy=TRUE)
plot(data76_tn_pb)

#### 7.5 Whittermoreの統計量 ####
col.W.whitt <- col.W.tango
whittermore.stat(data76_OE, col.W.whitt, zero.policy=TRUE)
whittermore.test(Observed~offset(log(Expected)), data76_OE, model="poisson", R=100, list=col.W.whitt, zero.policy=TRUE)
data76_wt_pb <- boot(data76_OE, statistic=whittermore.pboot, sim="parametric", ran.gen=poisson.sim, R=100, listw=col.W.whitt, zero.policy=TRUE)
plot(data76_wt_pb)

#### 7.6 Besag-Newell's test ####
data76_bn_perboot <- boot(data76_OE, statistic=besagnewell.boot, R=100, k=20)
plot(data76_bn_perboot)
besagnewell.stat(data76_OE, k=30)

#### 7.7 OpenshawのGAM ####
thegrid <- as(data76_OE, "data.frame")[,c("x","y")]
data76_opg <- opgam(data=as(data76_OE, "data.frame"), thegrid=thegrid,radius=20000, step=1000, alpha=0.05)
data76_opg
#図7.2
plot(data76_OE$x,data76_OE$y, xlab="Easting", ylab="Northing", cex=2,cex.axis=1.2, cex.lab=1.2)
points(data76_opg$x, data76_opg$y, col="red", pch="*", cex=4)

#### 7.8 Kulldorffの空間スキャン統計量 ####
data76 <- read.table("data76.csv", sep=",", header=TRUE, row.names=1)
data76_OE <- data.frame(Observed=data76$n.birth)
data76_OE <- cbind(data76_OE,
Expected=data76$pop*sum(data76$n.birth)/sum(data76$pop),
x=data76$Easting, y=data76$Northing)
dist <- (data76_OE$x-data76_OE$x[1])^2+(data76_OE$y-data76_OE$y[1])^2
index<-order(dist)
kullnagar.stat(data76_OE[index,], fractpop=.5)

## Kulldorff & Nagarwallaの検定 
# Permutation model
data76_kn_perboot <- boot(data76_OE, statistic=kullnagar.boot, R=100, fractpop=.5)
plot(data76_kn_perboot)
# Multinomial model
data76_m_pboot <- boot(data76_OE, statistic=kullnagar.pboot, sim="parametric", ran.gen=multinom.sim, R=100, fractpop=0.5)
plot(data76_kn_mboot)
# Poisson model
data76_kn_pboot <- boot(data76_OE, statistic=kullnagar.pboot, sim="parametric", ran.gen=poisson.sim, R=100, fractpop=0.5)
plot(data76_kn_pboot)
# Poisson-Gamma model
data76_kn_pgboot <- boot(data76_OE, statistic=kullnagar.pboot, sim="parametric",ran.gen=negbin.sim, R=100, fractpop=0.5)
plot(data76_kn_pgboot)
