## Run analysis, write model results


# install.packages("devtools")
# install.packages("FLCore", repo = "http://flr-project.org/R")
# install_github("ices-tools-prod/msy")
# install.packages("ggplotFL", repos="http://flr-project.org/R")
# #install.packages("msy", repos="http://flr-project.org/R")
# install.packages(c("copula","triangle", "coda"))
# install.packages("FLEDA", repos="http://flr-project.org/R")
# install.packages("FLFishery", repos="http://flr-project.org/R")
# 
# install.packages(c("FLa4a", "FLasher", "ggplotFL"), repos="https://flr-project.org/R")
# install.packages(c("ggplot2", "snpar", "foreach", "data.table"))
# devtools::install_github("flr/a4adiags")



library(devtools)
library(FLCore)
library(icesTAF)
library(FLa4a)
library(FLEDA)
library(ggplotFL)
library(gridExtra)
library(icesAdvice)
library(FLFishery)
library(FLasher)
library(msy)
library(a4adiags)
library(devtools)
library(ggplot2)





#  #TO WORK from TAF
# # 
# taf.bootstrap()
# # # # load index and stock
# 
# load(taf.data.path("stock", "meg78_stock.RData"))
# load(taf.data.path("indices", "meg78_indices.RData"))
# 

 #TO WORK from MY PC

setwd("~/GitHub/2024_meg.27.7b-k8abd_assessment")
#
#
load("bootstrap/data/stock/meg78_stock.RData")
load("bootstrap/data/indices/meg78_indices.RData")

mkdir("model")


# specify submodels defined as follows:

#FIT: ONLY SURVEYS

tun.sel2<- tun.sel[c("SP_PORC","CPUE.IRLFRsurvey")]

plot(tun.sel2[[1]]@index[1,])+geom_point()+ggtitle("index age 1")+ylab("Ind./30min haul")
plot(tun.sel2[[2]]@index[1,])+geom_point()+ggtitle("index age 1")+ylab("numbers per 10 hours fished")

#strange outlier, set to NA in IRLFRsurvey 2011

tun.sel2[[2]]@index['1','2011'] <- NA




#we do not really believe that the increase in 1-year-olds in the catch is real so we shouldn't formulate a model that treats it as real. I think the better approach is to 1984:


plot(stock@catch.n['1',])
stock@catch.n['1',as.character(1984:2000)] <- NA

# Also deal with  pretty flat after age 7 so reduce the age where they are 'bound'



qmod <- list(~ I(1 / (1 + exp(-age))),~ I(1 / (1 + exp(-age))))
srmod <- ~ factor(year) 

fmod <- ~ factor(replace(age,age>7,7)) + factor(year)


fit1 <- sca(stock, tun.sel2, fmodel = fmod, qmodel = qmod, srmodel = srmod)

stk1 <- stock + fit1


tun.sel2

################################################


#FIT 1: FIT FINAL

stk1 <- stock + fit1

plot(stk1)



# utilities
SavePlot <- function(plotname, width = 10, height = 7) {
  file <- file.path("report", paste0("meg78_data_", plotname, ".png"))
  dev.print(png, file, width = width, height = height, units = "in", res = 600)
}


save(stk1,stock, tun.sel2, fit1, file = "model/MegFit_FINAL.Rdata")

submodels(fit1)

# fmodel: ~factor(replace(age, age > 7, 7)) + factor(year)
# srmodel: ~factor(year)
# n1model: ~s(age, k = 3)
# qmodel:
#   SP_PORC:          ~I(1/(1 + exp(-age)))
# CPUE.IRLFRsurvey: ~I(1/(1 + exp(-age)))
# vmodel:
#   catch:            ~s(age, k = 3)
# SP_PORC:          ~1
# CPUE.IRLFRsurvey: ~1


name(stk1) <- 'Run1'

AIC(fit1)

BIC(fit1)


idx<-tun.sel2



#retro without smoothers-> 
back <- 5
retro <- split(1:back, 1:back)
retro <- lapply(retro, function(x){
  yr <- range(stock)["maxyear"] - x
  stk <- window(stock, end=yr)
  idx <- window(idx, end=yr)
  stk + sca(stk, idx,fmodel=fmod, qmodel=qmod, srmodel=srmod)
})

#FLStocks(retro)

mkdir("report")

retro$"0" <- stock+fit1

plot(FLStocks(retro),col=1,lwd=1)+ facet_wrap(~qname,scales='free_y')

save(retro, file = "model/MegRetroFit_FINAL.Rdata")

SavePlot('Retro_Fit_FINAL')


####AI: ERROR in Retro value calculation

Retro_F <- data.frame(Y0=c(fbar(retro$`0`))
                      ,Y1=c(fbar(retro$`1`),NA)
                      ,Y2=c(fbar(retro$`2`),NA,NA)
                      ,Y3=c(fbar(retro$`3`),NA,NA,NA)
                      ,Y4=c(fbar(retro$`4`),NA,NA,NA,NA)                
                      ,Y5=c(fbar(retro$`5`),NA,NA,NA,NA,NA)
)


mohn(Retro_F,details=T) 
mohn(Retro_F,plot=T) 

error <- function(years,x,se){
  #  polygon(c(years,rev(years)),c(x-se,rev(c(x+se))),col='#0000FF20',border=NA)
  polygon(c(years,rev(years)),c(x-1.96*se,rev(c(x+1.96*se))),col='#0000FF20',border=NA)
}

pal <- c("#7570B3","#1B9E77","#E6AB02","#D95F02","#E7298A")
years <- stk1@range[4]:stk1@range[5]
nyears <-length(years)
fits <- simulate(fit1, 1000)
stks <- stock + fits
fbar <- fbar(stk1)
fbarse <- apply(fbar(stks),2,sd) # AI: considering 1sd
fbar2sd <- 2*apply(fbar(stks),2,sd) # AI: considering 2sd
ylim <- c(0,max(fbar+1.96*fbarse))
ylim <- c(0,max(fbar+1.96*fbar2sd))
xlim <-range(years)
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='F',main=paste0('Fbar ',stk1@range[6],'-',stk1@range[7]))
error(years,fbar,fbarse)
error(years,fbar,fbar2sd)
for(i in 1:5) lines(years[1:(nyears-i)],fbar(retro[[i]]),lwd=1,col=pal[i])
lines(years,fbar,lwd=2,col=4)
#abline(h=0.28,lty=2)
#abline(h=0.181,lty=3)
#abline(h=0.39,lty=3)
#lines(1986:2005,c(0.3486,0.3176,0.3423,0.3854,0.3835,0.3433,0.271,0.2041,0.2105,0.2664,0.3087,0.3065,0.2512,0.155,0.1261,0.182,0.1981,0.2328,0.2376,0.1622))
#legend('topleft',c('Fbar','Fmsy','Frange'),lty=c(1,2,3),col=c(4,1,1),lwd=c(2,1,1),bty='n',y.intersp=0.8)
SavePlot('Retro_Fit_Fbar')

Retro_SSB <- data.frame(Y0=c(ssb(retro$`0`))
                        ,Y1=c(ssb(retro$`1`),NA)
                        ,Y2=c(ssb(retro$`2`),NA,NA)
                        ,Y3=c(ssb(retro$`3`),NA,NA,NA)
                        ,Y4=c(ssb(retro$`4`),NA,NA,NA,NA)                
                        ,Y5=c(ssb(retro$`5`),NA,NA,NA,NA,NA)
)

ssb <- ssb(stk1)/1000
ssbse <- apply(ssb(stks),2,sd)/1000 # AI: considering 1sd
ssb2sd <- 2*apply(ssb(stks),2,sd)/1000 # AI: considering 2sd
ylim <- c(0,max(ssb+1.96*ssbse))# AI: considering 1sd
ylim <- c(0,max(ssb+ssb2sd))# AI: considering 2sd
xlim <-range(years)
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='F',main=paste0('SSB '))
error(years,ssb,ssbse) #1sd
error(years,ssb,ssb2sd) #2sd
for(i in 1:5) lines(years[1:(nyears-i)],ssb(retro[[i]])/1000,lwd=1,col=pal[i])
lines(years,ssb,lwd=2,col=4)
SavePlot('Retro_Fit_SSB')


mohn(Retro_SSB) 
mohn(Retro_SSB,details=T) 
mohn(Retro_SSB,plot=T) 



recr <- function(x) x@stock.n[1,]
Retro_R <- data.frame(Y0=c(recr(retro$`0`))
                      ,Y1=c(recr(retro$`1`),NA)
                      ,Y2=c(recr(retro$`2`),NA,NA)
                      ,Y3=c(recr(retro$`3`),NA,NA,NA)
                      ,Y4=c(recr(retro$`4`),NA,NA,NA,NA)                
                      ,Y5=c(recr(retro$`5`),NA,NA,NA,NA,NA)
)
mohn(Retro_R) 
mohn(Retro_R,details=T) 
mohn(Retro_R,plot=T) 


AIC(fit1)
BIC(fit1)


###################


###REPORT SECTION

mkdir("report")


# load fit, index and stock

#load("model/MegFit_FINAL.Rdata")

# load retro run
#load("model/MegRetroFit_FINAL.Rdata")


stk1 <- stock + fit1
plot(stk1)
name(stk1) <- 'Run1'

AIC(fit1)
BIC(fit1)



# hacked function from a4aFitresiduals-class
plotr <- function(x, y=missing, auxline="smooth",...){
  args <- list()
  args$data <- as.data.frame(x)
  args$x <- as.formula("data~year|factor(age)*qname")
  args$type=c("p", auxline)
  args$groups <- quote(qname)
  args$cex=0.3
  args$lwd=0
  args$ylab="standardized residuals"
  args$xlab=""
  args$panel=function(x,y,...){
    panel.abline(h=0, col.line="gray80")
    panel.xyplot(x,y,...)
  }
  args$par.settings=list(
    superpose.symbol=list(col=1, pch=19, cex=0.2),
    superpose.line=list(col="gray75", lty=1, lwd=0, col=NA),
    strip.background=list(col="gray90"),
    strip.border=list(col="black"),
    box.rectangle=list(col="gray90"))
  args$main="log residuals of catch and abundance indices by age"
  if(is(latticeExtra::useOuterStrips, "function")) latticeExtra::useOuterStrips(do.call("xyplot", args)) else do.call("xyplot", args)
}

res <- residuals(fit1, stock, tun.sel2)
print(plotr(res, auxline = "none"))
SavePlot('Residuals1',10,6)

res$catch.n[is.na(res$catch.n)] <- 0 #hack
#res$FR_EVHOE[is.na(res$FR_EVHOE)] <- 0 #hack
res$SP_PORC[is.na(res$SP_PORC)] <- 0 #hack
#res$CPUE.Vigo84[is.na(res$CPUE.Vigo84)] <- 0 #hack
#res$CPUE.Vigo99[is.na(res$CPUE.Vigo99)] <- 0 #hack
#res$LPUE.ITBB[is.na(res$LPUE.ITBB)] <- 0 #hack
#res$CPUE.IAM[is.na(res$CPUE.IAM)] <- 0 #hack
res$CPUE.IRLFRsurvey[is.na(res$CPUE.IRLFRsurvey)] <- 0 #hack

bubbles(res)
SavePlot('Residuals2',10,6)

res <- residuals(fit1, stock, tun.sel2)
qqmath(res)
SavePlot('Residuals3')

wireframe(fit1@harvest)
#a <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,0.9),ylab='F',main='Fishing mortality')
a <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,1),ylab='F',main='Fishing mortality')
a
SavePlot('F')

#


fitted <- predict(pars(fit1))

#a1 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchabilty',main='Q FR_EVHOE',ylim=c(0,4e-5))
#a1
a2 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchabilty',main=names(tun.sel2)[1])
a2
#a3 <- xyplot(data~age,groups=year,data=fitted$qmodel[3],type='b',ylab='Catchabilty',main='Q CPUE.Vigo84',ylim=c(0,0.005))
#a3
#a4 <- xyplot(data~age,groups=year,data=fitted$qmodel[4],type='b',ylab='Catchabilty',main='Q CPUE.Vigo99',ylim=c(0,0.0015))
#a4
#a5 <- xyplot(data~age,groups=year,data=fitted$qmodel[5],type='b',ylab='Catchabilty',main='Q LPUE.ITBB',ylim=c(0,0.02))
#a5
#a6 <- xyplot(data~age,groups=year,data=fitted$qmodel[6],type='b',ylab='Catchabilty',main='CPUE.IAM',ylim=c(0,0.03))
#a6
a7 <- xyplot(data~age,groups=year,data=fitted$qmodel[2],type='b',ylab='Catchabilty',main=names(tun.sel2)[2])
a7

grid.arrange(a,a2,a7,ncol=2)
SavePlot('F and Q',8,6)




# simulate

fits <- simulate(fit1, 1000)
flqs <- FLQuants(sim=iterMedians(stock.n(fits)), det=stock.n(fit1))
keylst <- list(points=FALSE, lines=TRUE, space="right")
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l", main="Median simulations VS fit", scales=list(y=list(relation="free")), auto.key=keylst)
#fits@stock.n[,1:2] <- fit1@stock.n[,1:2]
#fits@catch.n[,1:2] <- fit1@catch.n[,1:2]
stks <- stock + fits
plot(stks) + facet_wrap(~qname,scales='free_y')
SavePlot('summary0')


#PROJECTIONS


#summary table for the report (AI: input for the projections accepted WGBIE2020: F NOT SCALED (Average F(2017-2019)), R geometric mean (1984:2017))

years <- stk1@range[4]:stk1@range[5]
nyears <- length(years)
GM <- round(exp(mean(log(c(stk1@stock.n[1,1:(nyears-2)])))),0)/1000 #Geometric mean of recruitment (1984-2017); thousands.
GM

#PROJECTIONS: F SCALED OR UNSCALED

fsq <- fsq <- mean(fbar(stk1)[,nyears-2:0]) #F NOT SCALED (Average F(2020-2022))
#fsq <- fsq <- fbar(stk1)[,nyears]  #F SCALED TO THE LAST YEAR

# AI: 08/06/2022 in the ADG the 2sd is requested and it is included in all xxse variables for example: tsbse <- 2*apply(tsb(stks),2,sd).

#Summary table using coeficient interval of 1sd

tsb <- tsb(stk1)
tsbse <- apply(tsb(stks),2,sd)
ssb <- ssb(stk1)
ssbse <- apply(ssb(stks),2,sd)
ssbInt <- NA #get this from stf
catchobs <- stock@catch
discardsobs <- stock@discards
catch <- stk1@catch
recr <- stk1@stock.n[1,]
recrse <- apply(stks@stock.n[1,],2,sd)
fbar <- fbar(stk1)
fbarse <- apply(fbar(stks),2,sd)

#Summary table using coeficient interval of 2sd

tsb <- tsb(stk1)
tsb2sd <- 2*apply(tsb(stks),2,sd)
ssb <- ssb(stk1)
ssb2sd <- 2*apply(ssb(stks),2,sd)
ssbInt <- NA #get this from stf
catchobs <- stock@catch
discardsobs <- stock@discards
catch <- stk1@catch
recr <- stk1@stock.n[1,]
recr2sd <- 2*apply(stks@stock.n[1,],2,sd)
fbar <- fbar(stk1)
fbar2sd <- 2*apply(fbar(stks),2,sd)



sum1 <-
  data.frame(
    Year = c(years, paste0(max(years) + 1, "*")),
    Lan = c(catchobs - discardsobs, NA),
    Dis = c(discardsobs, NA),
    Cat = c(catchobs, NA),
    CatEst = c(catch, NA),
    Tsb = c(tsb, NA),
    Ssb = c(ssb, ssbInt),
    SsbCv = c(ssbse / ssb, NA),
    Recr = c(recr, GM),
    RecrCv = c(recrse / recr, NA),
    Fbar = c(fbar, fsq),
    FbarCv = c(fbarse / fbar, NA)
  )

knitr::kable(sum1,row.names=F,digits=c(rep(3,12)))
write.csv(sum1,'report/Summary_FINAL.csv',row.names=F)

plot(years,tsb,type='l',ylim=c(0,120))
lines(years,ssb,type='l',col=2)
plot(years,ssb/tsb,type='l',ylim=c(0,1))


#summary for standard graphs template with coeficient interval of 1sd
sum2 <- data.frame(year=c(years,max(years+1))
                   ,recrlo=c(recr-recrse,NA)
                   ,recr=c(recr,GM)
                   ,recrhi=c(recr+recrse,NA)
                   ,tsblo=c(tsb-tsbse,NA)
                   ,tsb=c(tsb,NA)
                   ,tsbhi=c(tsb+tsbse,NA)
                   ,ssblo=c(ssb-ssbse,NA)
                   ,ssb=c(ssb,ssbInt)
                   ,ssbhi=c(ssb+ssbse,NA)
                   ,catch=c(catchobs,NA)
                   ,lan=c(catchobs-discardsobs,NA)
                   ,dis=c(discardsobs,NA)
                   ,ibc=NA
                   ,ur=NA
                   ,yssb=NA
                   ,flo=c(fbar-fbarse,NA)
                   ,f=c(fbar,fsq)
                   ,fhi=c(fbar+fbarse,NA)
)
knitr::kable(sum2,row.names=F,digits=c(rep(0,16),rep(3,3)))
write.csv(sum2,'report/Summary_2_1sd.csv', row.names=F)



#summary for standard graphs template with coeficient interval of 2sd
sum3 <- data.frame(year=c(years,max(years+1))
                   ,recrlo=c(recr-recr2sd,NA)
                   ,recr=c(recr,GM)
                   ,recrhi=c(recr+recr2sd,NA)
                   ,tsblo=c(tsb-tsb2sd,NA)
                   ,tsb=c(tsb,NA)
                   ,tsbhi=c(tsb+tsb2sd,NA)
                   ,ssblo=c(ssb-ssb2sd,NA)
                   ,ssb=c(ssb,ssbInt)
                   ,ssbhi=c(ssb+ssb2sd,NA)
                   ,catch=c(catchobs,NA)
                   ,lan=c(catchobs-discardsobs,NA)
                   ,dis=c(discardsobs,NA)
                   ,ibc=NA
                   ,ur=NA
                   ,yssb=NA
                   ,flo=c(fbar-fbar2sd,NA)
                   ,f=c(fbar,fsq)
                   ,fhi=c(fbar+fbar2sd,NA)
)
knitr::kable(sum3,row.names=F,digits=c(rep(0,16),rep(3,3)))
write.csv(sum3,'report/Summary_3_2sd.csv', row.names=F)


# sensitivity for template - not needed anymore?
ages <- stk1@range[1]:stk1@range[2]
stk3y <- window(stk1, start = max(years) - 2)
p <- apply(landings.n(stk3y) / catch.n(stk3y), 1, mean, na.rm = T)
p <- c(ifelse(is.na(p), 0, p))
sen <- data.frame(
  Age = ages,
  M = c(apply(m(stk3y), 1, mean)),
  Mat = c(apply(mat(stk3y), 1, mean)),
  PF = c(apply(harvest.spwn(stk3y), 1, mean)),
  PM = c(apply(m.spwn(stk3y), 1, mean)),
  Sel = c(apply(harvest(stk3y), 1, mean)) * p,
  WeCa = c(apply(landings.wt(stk3y), 1, mean)),
  Fd = c(apply(harvest(stk3y), 1, mean)) * (1 - p),
  WeCad = c(apply(discards.wt(stk3y), 1, mean)),
  Fi = 0,
  WeCai = 0
)
knitr::kable(sen, row.names = F, digits = c(0, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3))
write.csv(sen, "report/Sen.csv", row.names = F)



plot(fit1, stock)
SavePlot("Fit1")
plot(fit1, tun.sel[1])
SavePlot("Fit2")
plot(fit1, tun.sel[2])
SavePlot("Fit3")
plot(fit1, tun.sel[3])
SavePlot("Fit4")
plot(fit1, tun.sel[4])
SavePlot("Fit5")
plot(fit1, tun.sel[5])
SavePlot("Fit6")
plot(fit1, tun.sel[6])
SavePlot("Fit7")
plot(fit1, tun.sel[7])
SavePlot("Fit8")

fitSumm(fit1)





# Summary plot with reference points

error <- function(years, x, se) {
  #  polygon(c(years,rev(years)),c(x-se,rev(c(x+se))),col='#0000FF20',border=NA)
  polygon(c(years, rev(years)), c(x - 1.96 * se, rev(c(x + 1.96 * se))), col = "#0000FF20", border = NA)
}

pal <- c("#7570B3", "#1B9E77", "#E6AB02", "#D95F02", "#E7298A")


#### KOBE PLOT (AI: updated RP for megrim 78 ICES. 2022a. Benchmark Workshop for selected Megrim Stocks (WKMEGRIM)

kobe <- data.frame(Btrig=40444,Fmsy=0.233,B=rev(c(ssb(stk1))),F=rev(c(fbar(stk1))))
with(kobe[1,],plot(B/Btrig,F/Fmsy,xlim=c(0,3.5),ylim=c(0,3),main='meg78'))
rect(-1,-1,5,5,col='yellow',border=NA)
rect(-1,1,1,5,col='red',border=NA)
rect(1,-1,5,1,col='green',border=NA)
with(kobe[1,],text(B/Btrig,F/Fmsy,max(years),pos=1))
with(kobe[1,],points(B/Btrig,F/Fmsy,pch=16))
with(kobe,lines(B/Btrig,F/Fmsy))
#plot(kobe)
SavePlot('Kobe',6,6)

