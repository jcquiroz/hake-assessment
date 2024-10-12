if(!require('KScorrect')) {
  install.packages('KScorrect')
  library('KScorrect')
}
library(tidyverse)
library(FSA)

# growth Hembras & Machos
nedades = 15
Linf    = (163.72 + 108.89)/2
k       = .9*(0.08 + 0.11)/2
Lo      = 15
t0      = 0
M       = 0.21

cv      = 0.2
tallas  = seq(2,110) #floor(Linf))

mu_edad = sigma_edad = NULL
mu_edad[1] = Lo;

for (i in 2:nedades) {
  mu_edad[i] = mu_edad[i-1]*exp(-k)+Linf*(1-exp(-k))
}
sigma_edad = cv*mu_edad

# mixture
propo1 <- rep(1/nedades,1,nedades)
N = rep(1,1,nedades)
prop <- for (i in 2:nedades) {
  N[i] = N[i-1]*exp(-M)
}
propo2 = N/sum(N)

x <- rmixnorm(n = 5000, mean = mu_edad, sd = sigma_edad, pro = propo2)
hist(x, n=length(tallas), main="random bimodal sample")

# Confirm 'rmixnorm' above produced specified model
require(mclust)
mod <- mclust::Mclust(x)
mod             # Best model (correctly) has two-components with unequal variances
mod$parameters	# and approximately same parameters as specified above
# sd^2            # Note reports var (sigma-squared) instead of sd used above

# Density, distribution, and quantile functions
lenghtll <- seq(0, Linf, .1)
plot(lenghtll, dmixnorm(lenghtll, mean=mu_edad, sd=sigma_edad, pro=propo2),
     type="l", main="Normal mixture density")

plot(lenghtll, pmixnorm(lenghtll, mean=mu_edad, sd=sigma_edad, pro=propo2),
     type="l", main="Normal mixture cumulative")

plot(stats::ppoints(100), qmixnorm(stats::ppoints(100), mean=mu_edad, sd=sigma_edad, pro=propo2),
     type="l", main="Normal mixture quantile")

# 'expand' can be specified to prevent non-calculable quantiles
# Reduce 'expand'. (Values < 0.8 allow correct approximation)
q1 <- qmixnorm(stats::ppoints(100), mean=mu_edad, sd=sigma_edad, pro=propo2, expand=.5)
plot(stats::ppoints(100), q1, type="l", main="Quantile with reduced range")

# Requires functions from the 'mclust' package
# Confirmation that qmixnorm approximates correct solution
mpar   <- mod$param
approx <- qmixnorm(p=ppoints(100), mean = mpar$mean, pro = mpar$pro, sd = sqrt(mpar$variance$sigmasq))
known  <- qnorm(p=ppoints(100), mean=mpar$mean, sd=sqrt(mpar$variance$sigmasq))
cor(approx, known)  # Approximately the same
plot(approx, main="Quantiles for (unimodal) normal")
lines(known)
legend("topleft", legend=c("known", "approximation"), pch=c(NA,1), lty=c(1, NA), bty="n")


## Veamos
dbins = 2
set_tallas <- tibble(length = x) |> filter(length > 0) |> rowwise() |> 
  mutate(age = DLMtool::iVB(t0, k, Linf, length), bins = FSA::lencat(length, w=dbins)) |> 
  drop_na() |> mutate(agec = ceiling(age))

raw    <- xtabs(~bins+agec, data=set_tallas)
WR.key <- prop.table(raw, margin=1)
#View(WR.key)

maxedad = 18
wr.key <- WR.key[,c(1:maxedad)]

FSA::alkPlot(wr.key,"barplot")
FSA::alkPlot(wr.key,"barplot",col="Cork")
FSA::alkPlot(wr.key,"barplot",col=heat.colors(8))
FSA::alkPlot(wr.key,"barplot",showLegend=TRUE)
FSA::alkPlot(wr.key,"area")
FSA::alkPlot(wr.key,"lines")
FSA::alkPlot(wr.key,"splines", ylim = c(0.1,1) )
FSA::alkPlot(wr.key,"splines",span=0.2)
FSA::alkPlot(wr.key,"bubble")
FSA::alkPlot(wr.key,"bubble",col=FSA::col2rgbt("#468e8c",0.5))

len.n <- xtabs(~bins, data = set_tallas)

alkMeanVar(wr.key, length ~ bins + agec, set_tallas, len.n) # Bettoli-Miranda method

alkMeanVar(wr.key, length ~ bins + agec, set_tallas, len.n, method="QuinnDeriso") # Quinn-Deriso method


###########################################################
## Examples of fitting

plot(length ~ agec, data=set_tallas, pch=19)

# Fitting the typical parameterization of the von B function
( vb1 <- vbFuns() )

fit1 <- nls(length ~ vb1(agec, Linf, K, t0), data = set_tallas, 
  start = vbStarts(length ~ agec, data= set_tallas))
summary(fit1,correlation=TRUE)

curve(vb1(x, Linf = coef(fit1)), from=1, to=50, col="red", lwd=3, add=TRUE)


