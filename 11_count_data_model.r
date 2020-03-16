# パッケージ
library(DCluster)
library(pscl)
library(R2jags)
# データ
data84 <- read.table("data84.csv", sep=",", header=TRUE)
summary(data84)
# 図11.1 ヒストグラムの作図
hist(data84$Hospital, 0:200, xlab="一般病院数", ylab="自治体数", cex=1.5, cex.axis=1.5, main="", xlim=c(0,30), ylim=c(0,1000))

## ポアソン・ガンマモデルによる相対リスク
# 相対リスクの計算
r <- sum(data84$TrAcc)/sum(data84$Pop)
data84$TAExpected <- data84$Pop * r
data84$TARR <- data84$TrAcc / data84$TAExpected
# ポアソンモデル
data84$pvalpois <- ppois(data84$TrAcc,data84$TAExpected, lower.tail=FALSE)

## 経験ベイズ推定
# ポアソン・ガンマモデルによる相対リスクの計算
eb <- empbaysmooth(data84$TrAcc,data84$TAExpected)
data84$EBPG <- eb$smthrr

## 階層ベイズ推定
# JAGSコード(poisson_gamma.txt)
# Bivand et al. (2008) "Applied Spatial Data Analysis with R"
model{
for(i in 1:n){
acc.obs[i] ~ dpois(mu[i])
mu[i] <- theta[i] * acc.exp[i]
theta[i] ~ dgamma(nu, alpha)
}
nu~dgamma(0.01, 0.01)
alpha~dgamma(0.01,0.01)
}
# Rコード
n <- nrow(data84)
acc.obs <- data84$TrAcc
acc.exp <- data84$TAExpected
data <- list("n","acc.obs","acc.exp")
inits <- list(nu=1,alpha=1)
parameters<- c("theta","nu","alpha")
model.file <- system.file(package="R2jags", "model", "poisson_gamma.txt")
model84_pg <- jags(data=data, inits=inits, parameters, n.iter=1100,
n.burnin=100, n.chains=1, model.file=model.file)
print(model84_pg, digits=3)
plot(model84_pg)


#### 11.1 ポアソン回帰モデル ####
# 最尤推定法
# glm()関数を使う方法
model84_p <- glm(Hospital~Popd+Pop65r, data=data84, family=poisson)
summary(model84_p)
# ベイズ推定
# MCMCpoisson()関数を使う方法
library(MCMCpack)
model84_p.mcmc <- MCMCpoisson(Hospital~Popd+Pop65r, data=data84, mcmc=3000,
burnin=1000)
summary(model84_p.mcmc)
# ベイズ法
# bayescountパッケージのbayescount()関数を使う方法
# ベイズ推定
# JAGSを使う方法
####
# JAGSコード(poisson_regress.txt)
model{
for(i in 1:n){
y[i] ~ dpois(lambda[i])
lambda[i] <- exp(b0 + b1*x[i,1] + b2*x[i,2])
}
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
b2 ~ dnorm(0.0, 1.0E-6)
}
#
####
#
# Rコード
y <- as.vector(data84$Hospital)
n <- length(y)
x <- as.matrix(cbind(
data84$Popd,
data84$Pop65r))
data <- list("n","y","x")
inits <- list(b0=0, b1=0, b2=0)
parameters<- c("b0","b1","b2")
model.file <- system.file(package="R2jags", "model", "poisson_regress.txt")
model84_p.jags <- jags(data=data, inits=inits, parameters, n.iter=2000,
n.burnin=1000, n.chains=1, model.file=model.file)
print(model84_p.jags, digits=3)
plot(model84_p.jags)
# ポアソン回帰モデル（分散不均一）
####
#
# JAGSコード(poisson_regress_hetero.txt)
model{
for(i in 1:n){
y[i] ~ dpois(lambda[i])
lambda[i] <- exp(b0 + b1*x[i,1]) * u[i]
u[i] ~ dgamma(alpha, alpha)
}
alpha <- 1/ralpha
ralpha ~ dgamma(1, 1)
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
}
#
####
#
# Rコード
y <- as.vector(data84$Hospital)
n <- length(y)
x <- as.matrix(cbind(
data84$Popd,
data84$Pop65r))
data <- list("n","y","x")
inits <- list(b0=0, b1=0)
parameters<- c("b0","b1"
model.file <- system.file(package="R2jags", "model",
"poisson_regress_hetero.txt")
model84_ph.jags <- jags(data=data, inits=inits, parameters, n.iter=1100,
n.burnin=100, n.chains=1, model.file=model.file)
print(model84_ph.jags, digits=3)
plot(model84_ph.jags)

#### 11.2 負の二項分布モデル ####
# 最尤法
# glm.nb()関数を使う方法
model84_nb <- glm.nb(Hospital~Popd+Pop65r, data=data84)
summary(model84_nb)
# glm()関数を使う方法
model84_nb1 <- glm(Hospital~Popd+Pop65r, data=data84,
family=negative.binomial(1))
summary(model84_nb1)
# anova(model84_p, model84_nb1)

# ベイズ法
# bayescountパッケージのbayescount()関数を使う方法

# ベイズ法
# bayesmパッケージのrhierNegbinRw()関数を使う方法

# ベイズ法
# MCMCpackのMCMCmetrop1R()関数を使う方法
# ランダムウォーク・メトロポリス法
negbinfun <- function(theta, y, X){
k <- length(theta)
beta <- theta[1:(k-1)]
alpha <- exp(theta[k])
mu <- exp(X %*% beta)
log.like <- sum(
lgamma(y+alpha) - lfactorial(y) - lgamma(alpha) +
alpha * log(alpha/(alpha+mu)) +
y * log(mu/(alpha+mu))
)
}
x1 <- data84$Popd
x2 <- data84$Pop65r
XX <- cbind(1,x1,x2)
yy <- data84$Hospital
mu <- exp(1.5+x1+2*x2)*rgamma(length(x1),1)
post.samp <- MCMCmetrop1R(negbinfun, theta.init=c(0,0,0,0), y=yy, X=XX,
thin=1, mcmc=10000, burnin=1000,
tune=1.5, verbose=500, logfun=TRUE,
seed=list(NA,1))
raftery.diag(post.samp)
plot(post.samp)
summary(post.samp)

# ベイズ法
# JAGSを使う方法
####
# JAGSコード(negbin.txt)
model{
for(i in 1:n){
y[i] ~ dnegbin(p[i],r)
p[i] <- r/(r+lambda[i])
lambda[i] <- exp(b0 + b1*x[i,1] + b2*x[i,2])
}
r ~ dgamma(0.001, 0.001)
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
b2 ~ dnorm(0.0, 1.0E-6)
}
#
####
#
# Rコード
y <- as.vector(data84$Hospital)
n <- length(y)
x <- as.matrix(cbind(
data84$Popd,
data84$Pop65r))
data <- list("n","y","x")
inits <- list(r=1, b0=0, b1=0, b2=0)
parameters<- c("b0","b1","b2")
model.file <- system.file(package="R2jags", "model", "negbin.txt")
# n.iter=11000のとき、計算時間がかかる
model84_nb.jags <- jags(data=data, inits=inits, parameters, n.iter=2500,
n.burnin=500, n.chains=1, model.file=model.file)
print(model84_nb.jags, digits=3)
plot(model84_nb.jags)

#### 11.3 zero-inflated poisson (ZIP) model ####
# 最尤法
# zeroinfl()関数を使う方法
# 単純なモデル
model84_zip1 <- zeroinfl(Hospital~Popd+Pop65r, data=data84)
summary(model84_zip1)
# model84_zip1 <- zeroinfl(Hospital~Popd+Pop65r | Popd+Pop65r, data=data84)と同じ結果となる。
# 0値の発生確率を定数項で与えるモデル
model84_zip2 <- zeroinfl(Hospital~Popd+Pop65r | 1, data=data84)
summary(model84_zip2)
# カウントデータのモデル(Hospital~Popd+Pop65r+Emp3r)とzero-inflated model(Hospital~Popd+Pop65r)の組み合わせモデル
model84_zip3 <- zeroinfl(Hospital~Popd+Pop65r+Emp3r | Popd+Pop65r, data=data84)
summary(model84_zip3)

# ベイズ法
# zicパッケージのzic()関数を使う方法

# ベイズ法
# bayescountパッケージのbayescount()関数を使う方法

# ベイズ法
# JAGSを使う方法
####
# JAGSコード(zip.txt)
#
model{
for(i in 1:n){
z[i] ~ dpois(mu[i])
mu[i] <- (1-nz[i]) * -(log(1-p[i] + p[i] * exp(-lambda[i]))) +
         nz[i] * -(log(p[i]) -lambda[i] + y[i] * log(lambda[i]) - logfact(y[i]))
logit(p[i]) <- -1 * (a0 +a1*x[i,1] + a2*x[i,2])
lambda[i] <- exp(b0 + b1*x[i,1] + b2*x[i,2])
}
a0 ~ dnorm(0.0, 1.0E-6)
a1 ~ dnorm(0.0, 1.0E-6)
a2 ~ dnorm(0.0, 1.0E-6)
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
b2 ~ dnorm(0.0, 1.0E-6)
}
#
#  model{
#  for(i in 1:n){
#  z[i] ~ dpois(mu[i])
#  mu[i] <- (1-nz[i]) * -(log(p[i] + (1-p[i]) * exp(-lambda[i]))) +
#           nz[i] * -(log(1-p[i]) -lambda[i] + y[i] * log(lambda[i]) - logfact(y[i]))
#  logit(p[i]) <- a0 +a1*x[i,1] + a2*x[i,2]
#  lambda[i] <- exp(b0 + b1*x[i,1] + b2*x[i,2])
#  }
#  a0 ~ dnorm(0.0, 1.0E-6)
#  a1 ~ dnorm(0.0, 1.0E-6)
#  a2 ~ dnorm(0.0, 1.0E-6)
#  b0 ~ dnorm(0.0, 1.0E-6)
#  b1 ~ dnorm(0.0, 1.0E-6)
#  b2 ~ dnorm(0.0, 1.0E-6)
#  }
####
# Rコード
nz <- ifelse(y==0,0,1)
p <- rep(0.5,n)
z <- rep(0,n)
data <- list("n","y","x", "nz","p","z")
inits <- list(a0=0, a1=0, a2=0, b0=0, b1=0, b2=0)
parameters<- c("a0","a1","a2","b0","b1","b2")
model.file <- system.file(package="R2jags", "model", "zip.txt")
model84_zip.jags <- jags(data=data, inits=inits, parameters, n.iter=2500,
n.burnin=500, n.chains=1, model.file=model.file)
print(model84_zip.jags, digits=3)
plot(model84_zip.jags)
#
####

#### 11.4 zeri-inflated negative binomial (ZINB) model ####
# 最尤推定法
# 単純なモデル
model84_zinb1 <- zeroinfl(Hospital~Popd+Pop65r, data=data84, dist="negbin")
summary(model84_zinb1)
# 0値の発生確率を定数項で与えるモデル
model84_zinb2 <- zeroinfl(Hospital~Popd+Pop65r | 1, data=data84, dist="negbin")
summary(model84_zinb2)
# カウントデータのモデル(Hospital~Popd+Pop65r+Emp3r)とzero-inflated model(Hospital~Popd+Emp3r)の組み合わせモデル
model84_zinb3 <- zeroinfl(Hospital~Popd+Pop65r+Emp3r|Popd+Emp3r, data=data84, dist="negbin")
summary(model84_zinb3)

# ベイズ法
# bayescountパッケージのbayescount()関数を使う方法

# ベイズ推定
# JAGSコード(zinb.txt)
model{
for(i in 1:n){
z[i] ~ dpois(mu[i])
logit(p[i]) <- -1 * (a0 +a1*x[i,1] + a2*x[i,2])
mu[i] <- (1-nz[i]) * -(log(1 - p[i]
			   + p[i] * (r[i]/(r[i]+lambda[i]))^r[i] )) +
        nz[i] * -(log(p[i]) + loggam(y[i] + r[i]) - loggam(r[i])
             	   - loggam(y[i] + 1) + r[i] * log(r[i])
    			   + y[i] * log(lambda[i])
			   - (r[i] + y[i]) * log(r[i] + lambda[i]))
lambda[i] <- exp(b0 + b1*x[i,1] + b2*x[i,2])
r[i] ~ dgamma(0.001, 0.001)
}
a0 ~ dnorm(0.0, 1.0E-6)
a1 ~ dnorm(0.0, 1.0E-6)
a2 ~ dnorm(0.0, 1.0E-6)
b0 ~ dnorm(0.0, 1.0E-6)
b1 ~ dnorm(0.0, 1.0E-6)
b2 ~ dnorm(0.0, 1.0E-6)
}
####
# Rコード
nz <- ifelse(y==0,0,1)
p <- rep(0.1,n)
z <- rep(0,n)
r <- rep(1,n)
data <- list("n","y","x", "nz","p","z","r")
inits <- list(a0=0, a1=0, a2=0, b0=0, b1=0, b2=0)
parameters<- c("a0","a1","a2","b0","b1","b2")
#  inits <- list(b0=0, b1=0, b2=0)
#  parameters<- c("b0","b1","b2")
model.file <- system.file(package="R2jags", "model", "zinb.txt")
model84_zinb.jags <- jags(data=data, inits=inits, parameters, n.iter=2500,
n.burnin=500, n.chains=1, model.file=model.file)
print(model84_zinb.jags, digits=3)
plot(model84_zinb.jags)
