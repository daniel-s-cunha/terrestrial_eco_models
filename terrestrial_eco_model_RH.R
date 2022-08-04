# This code is used to model terrestrial processes 
# for the purpose of downstream estimation of phenology,
# land cover classification, and carbon stocks
# Annual estimates of harmonic and trend parameters are outputted for this purpose

setwd("Documents/FastCP/")
library(lme4)

create_design <- function(orig_df,h,tot_days,period){
  #assume orig_df has column t for time
  orig_df$year = floor(orig_df$t*tot_days/period)
  for (i in 1:h){
    orig_df = cbind(orig_df, sin(2*pi*i*orig_df$t*tot_days/period),cos(2*pi*i*orig_df$t*tot_days/period))
    names(orig_df)[2+2*i] = paste("sin",i,sep = "")
    names(orig_df)[3+2*i] = paste("cos",i,sep = "")
  }
  return(orig_df)
}

create_random_design <- function(df){
  #Get random effects design matrix
  #Currently for one harmonics
  #model = lmer(y ~ 1 + sin1 + cos1 + sin2 + cos2 + t + (0+sin1|year) + (0+cos1|year), data = df)
  #model = lmer(y ~ 1 + sin1 + cos1 + sin2 + cos2 + t + (0+sin1|year) + (0+cos1|year)+ (0+sin2|year) + (0+cos2|year), data = df)
  model = lmer(y ~ 1 + sin1 + cos1 + sin2 + cos2 + t + (1|year) + (0+sin1|year) + (0+cos1|year), data = df)
  #model = lmer(y ~ 1 + sin1 + cos1 + sin2 + cos2 + t + (0+t|year) + (0+sin1|year) + (0+cos1|year), data = df)
  #model = lmer(y ~ 1 + sin1 + cos1 + sin2 + cos2 + t + (1|year) + (0+t|year) + (0+sin1|year) + (0+cos1|year), data = df)
  Z = as.matrix(getME(model, "Z"))
  return(Z)  
}

maximize_sigma2_ga <- function(Ega2,nyears,nran_eff){
  # Since sigma2_ga is diagonal, return the vector of its diagonal elements
  # "year" is the vector of each random effect's year
  # NEED TO DOUBLE CHECK THAT Egamma IS ORDERED CORRECTLY
  s2g = c()
  for(i in 1:nran_eff){
    s2g   = c(s2g, rep(sum(Ega2[(nyears*(i-1)+1):(nyears*i)])/nyears, nyears)) 
  }
  return(s2g)  
}


period         = 365
h              = 2 #harmonics #currently random design only supports two harmonics
ndvi_ts        = read.csv(file = "data/NDVI_42.3993_-72.0266_14595.csv")
n              = dim(ndvi_ts)[1]
tot_day        = 14595#7303
ndvi_ts        = create_design(ndvi_ts,h,tot_day,period)
Z              = create_random_design(ndvi_ts)
nran_eff       = dim(Z)[2]/max(ndvi_ts$year)

f              = ""
for (i in 1:h){f = paste(f,paste("sin",i,sep = ""),sep=" + ")}
for (i in 1:h){f = paste(f,paste("cos",i,sep = ""),sep=" + ")}
formula        = as.formula(paste("y2 ~ 1 + t",f,sep=""))
maxit = 20
epsilon = 0.1
#Initialize the posterior expectation of random effects
nyears     = length(unique(ndvi_ts$year))
sigma2_ga  = rep(1,nran_eff*nyears)
sigma2     = 2.
Covgamma   = solve(t(Z)%*%Z + sigma2*diag(1/sigma2_ga))
Egamma     = rep(0,nran_eff*nyears)
p          = length(Egamma)
ndvi_ts$y2 = ndvi_ts$y - Z%*%Egamma #working response for XBeta
it = 1

for (it in 1:maxit) {
  # [ M-step ]
  m         = lm(formula, data = ndvi_ts)
  y3        = ndvi_ts$y - m$fitted.values #working response for Zgamma
  sigma2    = (sum(m$residuals^2) + sum(diag(Covgamma%*%t(Z)%*%Z)))/n
  sigma2_ga = maximize_sigma2_ga(diag(Covgamma + Egamma%*%t(Egamma)),nyears,nran_eff) #this should be a vector
  l         = dnorm(ndvi_ts$y, m$fitted + Z%*%Egamma, sqrt(sigma2), log = TRUE) 
  Q_new     = sum(l) - sum(log(diag(Covgamma)))/2 - sum(diag((Covgamma + Egamma%*%t(Egamma))%*%diag(sigma2_ga)))/2
  #
  # [ E-step ]
  Covgamma   = solve(t(Z)%*%Z + sigma2*diag(1/sigma2_ga))
  Egamma     = Covgamma%*%t(Z)%*%y3
  ndvi_ts$y2 = ndvi_ts$y - Z%*%Egamma #working response for XBeta
  #
  # [ Check convergence ]
  message("[", it, "] Q = ", Q_new)
  if (it > 15 && (Q_new - Q) / (Q + .1) < epsilon)
    break
  Q = Q_new
}

#list(m = m, sigma2 = sigma2, Egamma = Egamma)
# range = 400:600
# plot(x = ndvi_ts$t[range],y = (m$fitted + Z%*%Egamma)[range],col='orange')
# points(x = ndvi_ts$t[range],y = ndvi_ts$y[range], col='blue')

npred   = 10000
pred_df = data.frame(y=rnorm(n = (npred+1)),t=seq(0,1,1/npred))
pred_df = create_design(pred_df,h,tot_day,period)
pred_df = pred_df[pred_df$year %in% unique(ndvi_ts$year),]
pred_Z  = create_random_design(pred_df)
dim(pred_Z)
pred_y  = predict.lm(m,pred_df) + pred_Z%*%Egamma

lo = .6; hi = .75;
plot_values = function(t,low,high){(t >= low & t <= high)}
vals0 = plot_values(ndvi_ts$t,lo,hi)
vals1 = plot_values(pred_df$t,lo,hi)
plot( x = ndvi_ts[vals0,]$t,y = ndvi_ts[vals0,]$y, ylim = c(0.,.8),ylab='NDVI',xlab='time',col='blue')
lines(x = pred_df[vals1,]$t,y = pred_y[vals1],     col='orange')


# plot( x = unique(ndvi_ts$year),y = m$coefficients[2] + Egamma[1:nyears], ylim = c(-0.4,0.),ylab='Trend (dNDVI / dYear)',xlab='Year',col='blue')
# lines(x = unique(ndvi_ts$year),y = m$coefficients[2] + Egamma[1:nyears], ylim = c(-0.4,0.),col='blue')
# abline(h=m$coefficients[2],lty=2)


# plot( x = unique(ndvi_ts$year),y = m$coefficients[3] + Egamma[(nyears+1):(2*nyears)], ylim = c(-0.2,0.),ylab='Interseasonal Sine Harmonic',xlab='Year',col='blue')
# lines(x = unique(ndvi_ts$year),y = m$coefficients[3] + Egamma[(nyears+1):(2*nyears)], ylim = c(-0.2,0.),col='blue')
# abline(h=m$coefficients[3],lty=2)

