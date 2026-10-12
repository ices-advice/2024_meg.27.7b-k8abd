#############################################################################################
# catch_table.R - DESC
# /catch_table.R

# Copyright Iago MOSQUEIRA (WMR), 2021
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

### Clean slate
rm(list=ls())
#path <- ("D:/ICES benchmarks/2022 WKMEGRIM/RefPts/Meg89/")

path <- ("~/GitHub/2024_meg.27.7b-k8abd_assessment/")
setwd(path)

library(icesAdvice)
library(flextable)
library(data.table)
library(writexl)

load(paste0(path,"model/runs_CPUEIRLFR10corrected.RData"))

forecastYr <- range(runs[[1]])["maxyear"]-1

discardFages <- 1:2  #is this correct?
landFages <- range(runs[[1]])["minfbar"]:range(runs[[1]])["maxfbar"] #It is 3-6.

Funwanted <- function(x, ages=discardFages) { 
  quantMeans((discards.n(x)[ac(ages),] / catch.n(x)[ac(ages),]) *
               harvest(x)[ac(ages)])
}

Fwanted <- function(x, ages=landFages) {  
  quantMeans((landings.n(x)[ac(ages),] / catch.n(x)[ac(ages),]) *
               harvest(x)[ac(ages)])
}

metrics(window(runs[[1]], start=forecastYr, end=forecastYr),
        list(catch=catch, wanted=landings, unwanted=discards,
             F=fbar, Fwanted=Fwanted, Funwanted=Funwanted))

unlist(lapply(runs, is, "FLStock"))

# TABLE. Advice sheet annual catch options

tab_options <- lapply(runs, function(x) {
  model.frame(
    FLQuants(c(metrics(window(x, start=forecastYr, end=forecastYr),
                       list(catch=catch, wanted=landings, unwanted=discards,
                            F=fbar, Fwanted=Fwanted, Funwanted=Funwanted)),
               metrics(window(x, start=forecastYr+1, end=forecastYr+1),
                       list(SSB=ssb)))), drop=TRUE)})


tab_options <- rbindlist(tab_options, idcol="basis")

# ADD SSB change, 100 - (old / new) * 100
ssbiy <- c(ssb(runs[[1]])[, as.character(forecastYr)])
tab_options[, ssbchange:=(SSB - ssbiy) / ssbiy * 100]

# ADD TAC change
tab_options[, tacchange:=(catch - c(tac)) / c(tac) * 100]

# ADD advice change
tab_options[, advicechange:=(catch - c(advice)) / c(advice) * 100]

# FIX zeroes

tab_options[, (2:ncol(tab_options)) := lapply(.SD, function(x) ifelse(abs(x) < 0.0001, 0, x)), .SDcols=2:ncol(tab_options)]

# SET first column
# tab_options$basis <- c("MSY approach: F[MSY]", "F=MAP F[MSY lower]",
#                        "F=MAP F[MSY upper]", "F=0", "F[pa]", "F[lim]", 
#                        paste0("SSB (",forecastYr+1,")=B[pa]"),
#                        paste0("SSB(",forecastYr+1,")=B[lim]"), 
#                        paste0("SSB(",forecastYr+1,")=MSY B[trigger]"),
#                        paste0("SSB(",forecastYr+1,")=SSB(",forecastYr,")"),
#                        "F[2020]", "Roll-over TAC")

tab_options$basis <- c("MSY approach: F[MSY]", "F=MAP F[MSY lower]",
                       "F=MAP F[MSY upper]", "F=0", "F[pa]", "F[lim]", 
                       paste0("SSB (",forecastYr+1,")=B[pa]"),
                       paste0("SSB(",forecastYr+1,")=B[lim]"), 
                       paste0("SSB(",forecastYr+1,")=MSY B[trigger]"),
                       "F[2023]",                                             #"F[2023]" hay que cambiarlo manualmente, en el WGBIE 2023 ponemos F(2023).
                       paste0("SSB(",forecastYr+1,")=SSB(",forecastYr,")"),
                       "Roll-over TAC")                     

# CALL round and icesRound
tab_options[, c(2,3,4,8) := lapply(.SD, round, digits=0), .SDcols = c(2,3,4,8)]
tab_options[, c(5,6,7) := lapply(.SD, round, digits=3), .SDcols = c(5,6,7)]
tab_options[, c(9,10,11) := lapply(.SD, icesRound), .SDcols = c(9,10,11)]

## ICES Rounding for F
x <- as.logical(tab_options[,"F"]>=0.2)
tab_options[x,"F"] <- round(tab_options[x,"F"],2)
x <- as.logical(tab_options[,"Fwanted"]>=0.2)
tab_options[x,"Fwanted"] <- round(tab_options[x,"Fwanted"],2)
x <- as.logical(tab_options[,"Funwanted"]>=0.2)
tab_options[x,"Funwanted"] <- round(tab_options[x,"Funwanted"],2)

# CREATE table
#ft <- flextable(tab_options[c(1,2,3,1,11,4,5,6,7,8,9,10,12),])
ft <- flextable(tab_options)

# SET colnames

ftabops <- set_header_labels(ft, 
                             basis="Basis", 
                             catch=paste0("Total catch (",forecastYr,")"), 
                             wanted=paste0("Wanted catch (",forecastYr,")"),
                             unwanted=paste0("Unwanted catch (",forecastYr,")"), 
                             F=paste0("F[total] (ages ",min(landFages),"-",max(landFages),") (",forecastYr,")"),
                             Fwanted=paste0("F[wanted] (ages ",min(landFages),"-",max(landFages),") (",forecastYr,")"), 
                             Funwanted= paste0("F[unwanted] (ages ",min(discardFages),"-",max(discardFages),") (",forecastYr,")"), 
                             SSB=paste0("SSB (",forecastYr+1,")"), 
                             ssbchange="% SSB change", 
                             tacchange="% TAC change",
                             advicechange="% Advice change")

dimnames(tab_options)[[2]] <- c("Basis", 
                               paste0("Total catch (",forecastYr,")"), 
                               paste0("Wanted catch (",forecastYr,")"),
                               paste0("Unwanted catch (",forecastYr,")"), 
                               paste0("F[total] (ages ",min(landFages),"-",max(landFages),") (",forecastYr,")"),
                               paste0("F[wanted] (ages ",min(landFages),"-",max(landFages),") (",forecastYr,")"), 
                               paste0("F[unwanted] (ages ",min(discardFages),"-",max(discardFages),") (",forecastYr,")"), 
                               paste0("SSB (",forecastYr+1,")"), 
                               "% SSB change", "% TAC change", "% Advice change")


writexl::write_xlsx(tab_options,
                    path=paste0(path,"report/catch_options_CPUEIRLFR10corrected.xlsx"))


