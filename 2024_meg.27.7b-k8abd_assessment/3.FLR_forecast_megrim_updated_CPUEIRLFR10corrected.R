# model.R - Run forecast
# 2020_sol.27.4_forecast/model.R

# Copyright Iago MOSQUEIRA (WMR), 2020
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# install.packages("devtools")
# install.packages("FLCore", repo = "http://flr-project.org/R")
# install_github("ices-tools-prod/msy")
# install.packages(c("FLa4a", "FLasher", "ggplotFL"), repos="https://flr-project.org/R")
# install.packages(c("ggplot2", "snpar", "foreach", "data.table"))
# devtools::install_github("flr/a4adiags")



library(FLCore)
library(icesTAF)
library(ggplot2)
library(ggplotFL)
library(FLFishery)
library(FLasher)
library(FLa4a)
library (a4adiags)

getwd()

setwd("~/GitHub/2024_meg.27.7b-k8abd_assessment")



mkdir("model")

# Load assessment

load('Refpoints_WKMEGRIM_megrim78_run3/refPts_CPUEIRLFR10_corrected.RData')#AI: importante cargar primero el archivo de Refpoints y luego el Fit, porque en este ultimo tenemos todos los a?os del assessment.

load('model/MegFit_FINAL.RData')

run <- stk1 # replace XXX with the final FLStock object


#years <- range(stks)["minyear"]:range(stks)["maxyear"]

# --- SETUP

# refpts: Bpa, Fpa, Blim, Flim, BmsyTrig, Fmsy, FmsyLow, FmsyUpp

# MODEL, ADVICE and FORECAST year

dy <- dims(run)$maxyear #last year of assessment input data
ay <- dy + 1
fy <- ay + 1

# TAC & advice current year

#advice <- FLQuant(19184, dimnames=list(age='all', year=ay), units="tonnes")  aqui se pondr?a el consejo para el a?o del grupo de trabajo es decir, en WGBIE 2022 hay que poner el consejo para el 2022. 

#advice <- FLQuant(23596, dimnames=list(age='all', year=ay), units="tonnes") #en WGBIE 2023, advice para el a?o 2023 son 23596 tn.
advice <- FLQuant(23303, dimnames=list(age='all', year=ay), units="tonnes") #en WGBIE 2024, advice para el a?o 2024 son 23303 tn.
tac <- advice



# GEOMEAN but last year
#rec1gm <- exp(mean(log(window(stock.n(run)["1",], end=-1))))

# GEOMEAN for all years minus last 2
rec1gm <- exp(mean(log(window(stock.n(run)["1",], end=-2))))

# GEOMEAN for years 1984:2018. AI: used in WGGBIE 2021 (GEOMEAN for all years minus last 2) # AI: distintas formas de obtener "rec1gm" en el primero por defecto, en este hay que modificar el a?o manualmente.
rec1gm <- exp(mean(log(window(stock.n(run)["1",], start=1984, end=2021))))

#################
# GEOMEAN: para cambiar el reclutamiento del ?ltimo a?o del assessment por geomean si consideramos que tenemos mucha incertidumbre.

#plot(run@stock.n['1',])
#run@stock.n['1','2022'] <- rec1gm #AI: oara la evaluci?n con datos hasta 2022.



# --- SETUP future

# 3 years, 5 years wts/selex, 3 years discards
fut <- stf(run, nyears=3, wts.nyears = 5, fbar.nyears=5, disc.nyears=3)  

# GET F status quo (Fsq)
Fsq <-expand(yearMeans(fbar(fut)[, ac(seq(dy - 2, dy))]), year=2024) #F NOT SCALED (Average F last 3 years)
#Fsq <-FLCore::expand(yearMeans(fbar(fut)[, ac(seq(dy - 2, dy))]), year=2023) # Be careful not to have package (tidyr)instlaed. In case of error use: FLCore::expand
#Fsq <- expand(fbar(fut)[, ac(dy)], year=ay) # #F SCALED TO THE LAST YEAR # LC 2024 use this???

#######################################################################################
# #AI script used in WGBIE 2021 for Fsq
# #PROJECTIONS: F SCALED OR UNSCALED??
# 
# fsq <- fsq <- mean(fbar(stk1)[,nyears-2:0]) #F NOT SCALED (Average F(2017-2019))
# #fsq <- fsq <- fbar(stk1)[,nyears]  #F SCALED TO THE LAST YEAR

########################################################################################

# SET geomean SRR
gmsrr <- predictModel(model=rec~a, params=FLPar(c(rec1gm), units="thousands",
  dimnames=list(params="a", year=seq(ay, length=3), iter=1)))

# > gmsrr   # AI: reemplaza el reclutamiento con rec1gm del a?o intermedio (2023) y los dos siguientes.
# An object of class "FLQuants": EMPTY
# model:  
#   rec ~ a
# 
# params:  
#   An object of class "FLPar"
# year
# params    2022        2023      2024  
# a        221059      221059    221059
# units:  thousands 

# GENERATE targets from refpts
refpts <- FLPar(refPts[1,])
targets <- expand(as(refpts, 'FLQuant'), year=fy)

##############################################################################
#AI: added to include a new variable in the output
# GENERATE targets from refpts   
#need to run the forecast once to see what the SSB at the start of the advice year is: as.numeric(ssb(runs[[1]])[,as.character(fy)])

as.numeric(ssb(runs[[1]])[,as.character(fy)]) 
#83883.42


# Then enter it manually
refPts <- cbind(refPts,  83883) 
dimnames(refPts)[[2]][ncol(refPts)] <- "constSSB"
refpts <- FLPar(refPts[1,])
targets <- expand(as(refpts, 'FLQuant'), year=fy)


#################

#AI: to complete in the advice the assumptions for interim year and forecast.


fbar(runs$Fmsy)[, '2024']
fbar(runs$Fmsy)[, '2025']
ssb(runs$Fmsy)[, '2025']
rec(runs$Fmsy)[, '2025']
catch(runs$Fmsy)[, '2024']
landings(runs$Fmsy)[, '2024']
discards(runs$Fmsy)[, '2024']

#Otro forma de obtener los valores de los a?os intermedios:

# as.numeric(catch(runs[[1]])[,as.character(ay)])
# #[1] 19993.89
# ay
# #[1] 2022
# as.numeric(landings(runs[[1]])[,as.character(ay)])
# #[1] 17123.85
# as.numeric(discards(runs[[1]])[,as.character(ay)])
# #[1] 2870.042


ssb(runs$Fmsy)[, '2024'] # La SSB del a?o intermedio (2024) hay que incluirla en la tabla del advice donde se da la summary table pero este ultimo a?o no hay valor (NA).
#89073

###############################################################################


# --- PROJECT catch options

# Targets 2021
C0 <- FLQuant(0, dimnames=list(age='all', year=fy))
Fiy <- FLQuants(fbar=append(Fsq, C0))
Ciy <- FLQuants(catch=append(tac, C0))


# RUN for Fsq
Fsqrun <- fwd(fut, sr=gmsrr, control=as(Fiy, "fwdControl"))


# TEST if Cay < TACay
if(catch(Fsqrun)[, ac(ay)] <= tac) {
  itarget <- Fiy
} else {
  itarget <- Ciy
  Fsqrun <- fwd(fut, sr=gmsrr, control=as(Ciy, "fwdControl"))
}

# DEFINE catch options
catch_options <- list(

  # FMSY
  Fmsy=FLQuants(fbar=targets["Fmsy",]),

  # lowFMSY
  lFmsy=FLQuants(fbar=targets["Fmsylower_unconstr",]),

  # uppFMSY
  uFmsy=FLQuants(fbar=targets["Fmsyupper",]),

  # F0
  F0=FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=fy))),

  # Fpa
  Fpa=FLQuants(fbar=targets["Fpa",]),

  # Flim
  Flim=FLQuants(fbar=targets["Flim",]),

  # Bpa
  Bpa=FLQuants(ssb_flash=targets["Bpa",]),

  # Blim
  Blim=FLQuants(ssb_flash=targets["Blim",]),

  # MSYBtrigger
  MSYBtrigger=FLQuants(ssb_flash=targets["MSYBtrigger",]),

  # F sq
  Fsq=FLQuants(fbar=expand(fbar(Fsqrun)[, ac(ay)], year=fy)),
  
  # SB sq
  SBsq=FLQuants(ssb_flash=expand(ssb(Fsqrun)[, ac(ay+1)], year=fy)),
  
  # TAC
  rotac=FLQuants(catch=expand(tac, year=fy))
)

# C0
F0 <- FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=fy + 1)))

# CONVERT to fwdControl

fctls <- lapply(catch_options, function(x) {
  as(FLQuants(c(window(itarget, end=ay), x, F0)), "fwdControl")
}
)

# RUN!

runs <- FLStocks(lapply(fctls, function(x) fwd(fut, sr=gmsrr, control=x)))

# SSB 2025




# COMPARE

Map(compare, runs, fctls)


# --- PROJECT F levels

flevels <- seq(0, 0.50, 0.01)

control <- as(as(c(lapply(window(itarget, end=ay), propagate, length(flevels)),
  FLQuants(fbar=FLQuant(flevels,  dimnames=list(year=fy, iter=seq(length(flevels))))),
  lapply(F0, propagate, length(flevels))), "FLQuants"), "fwdControl")

f_runs <- divide(fwd(fut, sr=gmsrr, control=control), names=flevels)

# SAVE

save(runs, f_runs, rec1gm, tac, advice, file="model/runs_CPUEIRLFR10corrected.RData")





