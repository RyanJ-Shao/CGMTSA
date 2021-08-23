# The MIT License (MIT)
# Copyright (c) 2020 UCAS

#' Calculate the metrics of CGM data.
#' @param inputdir Path of input directory containing CGM files.
#' @param outputdir Path of output directory where metrics will be stored.
#' @param useig Logical. If TRUE the imputed data will be used.
#' @param diffnum The number of differencing.
#' @param threshold The threshold of sd used in mage calculation.
#' @param bthreshold The minimum glucose of tir range.
#' @param athreshold The maximum glucose of tir range.
#' @param interval The interval of CGM data.
#' @export

cgmmetrics <- function(inputdir, outputdir, useig = FALSE, diffnum = 1, threshold =1, bthreshold = 3.9, athreshold = 10, interval = 15){
  fileNames = list.files(inputdir)
  freq = 1440/interval
  fnamevec <- c()
  sdvec <- c()
  meanvec <- c()
  cvvec <- c()
  gmisvec <- c()
  lbgivec <- c()
  hbgivec <- c()
  magevec <- c()
  tirvec <- c()
  moddvec <- c()
  acf_1vec <- c()
  acf_2vec <- c()
  acf_3vec <- c()
  acf_4vec <- c()
  acf_5vec <- c()
  pacf_1vec <- c()
  pacf_2vec <- c()
  pacf_3vec <- c()
  pacf_4vec <- c()
  pacf_5vec <- c()
  for(f in fileNames){
    fname <- unlist(strsplit(f, split = "\\."))[1]
    fnamevec <- append(fnamevec, fname)
    print(paste("processing file:", f))
    cgmtsall <- read.csv(paste(inputdir, "/", f, sep = ''),stringsAsFactors= FALSE)
    cgmtsall <- cgmtsall[order(lubridate::ymd_hm(cgmtsall$timestamp)),]
    vectimestamp <- as.vector(cgmtsall$timestamp)
    vectimestamp <- unlist(strsplit(vectimestamp,split=" "))
    maxtimestamp <- matrix(vectimestamp,ncol=2,byrow=T)[,1]
    cgmtsall <- dplyr::mutate(cgmtsall, timedate = maxtimestamp)

    sd <- sdall(cgmtsall = cgmtsall, useig = useig)
    sdvec <- append(sdvec, sd)
    mean <- meanall(cgmtsall = cgmtsall,useig = useig)
    meanvec <- append(meanvec, mean)
    cv <- cvall(cgmtsall = cgmtsall, useig = useig)
    cvvec <- append(cvvec, cv)
    gmis <-gmi(cgmtsall = cgmtsall, useig = useig)
    gmisvec <- append(gmisvec, gmis)
    lhbgi <-lhbgiall(cgmtsall = cgmtsall, useig = useig)
    lbgi <- lhbgi[1]
    lbgivec <- append(lbgivec, lbgi)
    hbgi <- lhbgi[2]
    hbgivec <- append(hbgivec, hbgi)
    mage <- mageall(cgmtsall = cgmtsall, useig = useig, threshold = threshold)
    magevec <- append(magevec, mage)
    tir <- tirall(cgmtsall = cgmtsall, useig = useig, athreshold = athreshold, bthreshold = bthreshold)
    tirvec <- append(tirvec, tir)
    modd <- round(mean(modd(cgmtsall = cgmtsall, useig = useig)),2)
    moddvec <- append(moddvec, modd)
    gacf <- acfcoff(cgmtsall, useig = useig, diffnum = diffnum, interval = interval)
    gpacf <- pacfcoff(cgmtsall, useig = useig, diffnum = diffnum, interval = interval)
    acf_1vec <- append(acf_1vec, gacf[1])
    acf_2vec <- append(acf_2vec, gacf[2])
    acf_3vec <- append(acf_3vec, gacf[3])
    acf_4vec <- append(acf_4vec, gacf[4])
    acf_5vec <- append(acf_5vec, gacf[5])
    pacf_1vec <- append(pacf_1vec, gpacf[1])
    pacf_2vec <- append(pacf_2vec, gpacf[2])
    pacf_3vec <- append(pacf_3vec, gpacf[3])
    pacf_4vec <- append(pacf_4vec, gpacf[4])
    pacf_5vec <- append(pacf_5vec, gpacf[5])
    mtcdf <- mtcgrpday(cgmtsall,useig = useig, threshold = threshold,  bthreshold = bthreshold, athreshold = athreshold, freq = freq)
    write.csv(mtcdf, paste(outputdir,fname, "_metricsByDay.csv",sep = ""),row.names = FALSE)
  }
  summetrics <- data.frame(ID = fnamevec, SD = sdvec, Mean = meanvec, CV = cvvec, GMI = gmisvec,
                           LBGI = lbgivec, HBGI = hbgivec, MAGE = magevec, TIR = tirvec,
                           MODD = moddvec,acf_1 = acf_1vec, acf_2 = acf_2vec, acf_3 = acf_3vec,
                           acf_4 = acf_4vec, acf_5 = acf_5vec, pacf_1 = pacf_1vec, pacf_2 = pacf_2vec,
                           pacf_3 = pacf_3vec, pacf_4 = pacf_4vec, pacf_5 = pacf_5vec)
  write.csv(summetrics, paste(outputdir,"metricsSummary.csv",sep = ""),row.names = FALSE)
}


mtcgrpday <- function(cgmtsall,useig = FALSE, threshold = 1,  bthreshold = 3.9, athreshold = 10, freq = 96){
  coldate <- unique(cgmtsall$timedate)
  sdcol <- sdgrpbyday(cgmtsall = cgmtsall, useig = useig)
  meancol <- meangrpbyday(cgmtsall = cgmtsall, useig = useig)
  cvcol <- cvgrpbyday(cgmtsall = cgmtsall, useig = useig)
  lhbgilist <- lhbgi(cgmtsall = cgmtsall, useig = useig)
  lbgicol <- lhbgilist$lbgivec
  hbgicol <- lhbgilist$hbgivec
  magecol <- mage(cgmtsall = cgmtsall, useig = useig, threshold = threshold)
  tircol <- tir(cgmtsall = cgmtsall, useig = useig,  bthreshold = bthreshold, athreshold = athreshold, freq = freq)
  moddcol <- c(NA)
  moddcol <- append(moddcol, modd(cgmtsall = cgmtsall, useig = useig))
  metricsdf <- data.frame(Day = coldate, SD = sdcol, Mean = meancol, CV = cvcol,
                           LBGI = lbgicol, HBGI = hbgicol, MAGE = magecol, TIR = tircol,
                           MODD = moddcol)
  return(metricsdf)
}


sdgrpbyday <- function(cgmtsall, useig = FALSE){
  coldate <- unique(cgmtsall$timedate)
  sdvec <- c()
  for(d in coldate){
    cgmts <- dplyr::filter(cgmtsall, timedate == d)
    if(useig){
      sdvec <- append(sdvec, round(sd(cgmts$imglucose), 2))
    }else{
      sdvec <- append(sdvec, round(sd(cgmts$sglucose), 2))
    }
  }
  return(sdvec)
}

sdall <- function(cgmtsall, useig = FALSE){
  sdvec = 0
  if(useig){
    sdvec = round(sd(cgmtsall$imglucose), 2)
  }else{
    sdvec = round(sd(cgmtsall$sglucose), 2)
  }
  return(sdvec)
}

meangrpbyday <- function(cgmtsall, useig = FALSE){
  coldate <- unique(cgmtsall$timedate)
  meanvec <- c()
  for(d in coldate){
    cgmts <- dplyr::filter(cgmtsall, timedate == d)
    if(useig){
      meanvec <- append(meanvec, round(mean(cgmts$imglucose), 2))
    }else{
      meanvec <- append(meanvec, round(mean(cgmts$sglucose), 2))
    }
  }
  return(meanvec)
}

meanall <- function(cgmtsall, useig = FALSE){
  meanall = 0
  if(useig){
    meanall = round(mean(cgmtsall$imglucose), 2)
  }else{
    meanall = round(mean(cgmtsall$sglucose), 2)
  }
  return(meanall)
}

cvgrpbyday <- function(cgmtsall, useig = FALSE){
  coldate <- unique(cgmtsall$timedate)
  cvvec <- c()
  for(d in coldate){
    cgmts <- dplyr::filter(cgmtsall, timedate == d)
    if(useig){
      cvvec <- append(cvvec, round(sd(cgmts$imglucose)/mean(cgmts$imglucose), 2))
    }else{
      cvvec <- append(cvvec, round(sd(cgmts$sglucose)/mean(cgmts$sglucose), 2))
    }
  }
  return(cvvec)
}


cvall <- function(cgmtsall, useig = FALSE){
  cvall = 0
  if(useig){
    cvall = round(sd(cgmtsall$imglucose)/mean(cgmtsall$imglucose), 2)
  }else{
    cvall = round(sd(cgmtsall$sglucose)/mean(cgmtsall$sglucose), 2)
  }
  return(cvall)
}

#https://doi.org/10.2337/dc18-1581
gmi <- function(cgmtsall, useig = FALSE){
  mgmi = 0
  if(useig){
    mgmi <- 3.31 + 0.02392 * mean(cgmtsall$imglucose*18)
  }else{
    mgmi <- 3.31 + 0.02392 * mean(cgmtsall$sglucose*18)
  }
  return(round(mgmi,2))
}

#Assessment of Risk for Severe Hypoglycemia Among Adults With IDDM
#Symmetrization of the blood glucose measurement scale and its applications
lhbgi <- function(cgmtsall, useig = FALSE){
  coldate <- unique(cgmtsall$timedate)
  lbgivec <- c()
  hbgivec <- c()
  for(d in coldate){
    cgmts <- dplyr::filter(cgmtsall, timedate == d)
    if(useig){
      fbg <- 1.794 * (log(cgmts$imglucose) ** 1.026 - 1.861)
    }else{
      fbg <- 1.794 * (log(cgmts$sglucose) ** 1.026 - 1.861)
    }
    rlbg <- fbg[fbg < 0]
    rhbg <- fbg[fbg > 0]
    lbgivec <- append(lbgivec, round(mean(10 * (rlbg ** 2)), 2))
    hbgivec <- append(hbgivec, round(mean(10 * (rhbg ** 2)), 2))
  }
  out <- list(lbgivec=lbgivec, hbgivec=hbgivec)
  return(out)
}


lhbgiall <- function(cgmtsall, useig = FALSE){
  if(useig){
    fbg <- 1.794 * (log(cgmtsall$imglucose) ** 1.026 - 1.861)
  }else{
    fbg <- 1.794 * (log(cgmtsall$sglucose) ** 1.026 - 1.861)
  }
  rlbg <- fbg[fbg < 0]
  rhbg <- fbg[fbg > 0]
  lbgivec <- round(mean(10 * (rlbg ** 2)), 2)
  hbgivec <- round(mean(10 * (rhbg ** 2)), 2)
  out <- c(lbgivec, hbgivec)
  return(out)
}

#Mean Amplitude of Glycemic Excursions, a Measure of Diabetic Instability
mage <- function(cgmtsall, useig = FALSE, threshold = 1){
  coldate <- unique(cgmtsall$timedate)
  magevec <- c()
  for(d in coldate){
    cgmts <- dplyr::filter(cgmtsall, timedate == d)

    daysd <- 0
    dayglucose <- c()
    if(useig){
      daysd <- sd(cgmts$imglucose)
      dayglucose <- cgmts$imglucose
    }else{
      daysd <- sd(cgmts$sglucose)
      dayglucose <- cgmts$sglucose
    }
    if(daysd ==0){
      magevec <- append(magevec, NA)
      next
    }
    daysd <- threshold * daysd
    diffglucose <- diff(dayglucose)

    inorder <- TRUE
    sflag <- TRUE
    turnpoints <- c()
    for(i in seq_along(diffglucose)){
      if(diffglucose[i] ==0){
        next
      }
      if(sflag){
        if(diffglucose[i] < 0){
          inorder <- TRUE
        }else{
          inorder <- FALSE
        }
        sflag <- FALSE
        next
      }
      if(inorder){
        if(diffglucose[i] < 0){
          next
        }else{
          turnpoints <- append(turnpoints, dayglucose[i])
          inorder <- FALSE
        }

      }else{
        if(diffglucose[i] < 0){
          turnpoints <- append(turnpoints, dayglucose[i])
          inorder <- TRUE
        }else{
          next
        }
      }
    }
    if(length(turnpoints) <=2){
      magevec <- append(magevec, NA)
      next
    }
    tpdiff <- diff(turnpoints)
    tpdiff <- tpdiff[abs(tpdiff) > daysd]
    if(length(tpdiff) ==0){
      magevec <- append(magevec, NA)
      next
    }
    if(tpdiff[1] > 0){
      tpdiff <- tpdiff[tpdiff > 0]
    }else if(tpdiff[1] < 0){
      tpdiff <- tpdiff[tpdiff < 0]
    }
    magevec <- append(magevec, abs(round(mean(tpdiff),2)))
  }

  return(magevec)

}
mageall <- function(cgmtsall, useig = FALSE, threshold = 1){
  mage <- 0
  daysd <- 0
  glucosevec <- c()
  if(useig){
    daysd <- sd(cgmtsall$imglucose)
    glucosevec <- cgmtsall$imglucose
  }else{
    daysd <- sd(cgmtsall$sglucose)
    glucosevec <- cgmtsall$sglucose
  }
  daysd <- threshold * daysd
  diffglucose <- diff(glucosevec)
  inorder <- TRUE
  sflag <- TRUE
  turnpoints <- c()
  for(i in seq_along(diffglucose)){
    if(diffglucose[i] ==0){
      next
    }
    if(sflag){
      if(diffglucose[i] < 0){
        inorder <- TRUE
      }else{
        inorder <- FALSE
      }
      sflag <- FALSE
      next
    }
    if(inorder){
      if(diffglucose[i] < 0){
        next
      }else{
        turnpoints <- append(turnpoints, glucosevec[i])
        inorder <- FALSE
      }
    }else{
        if(diffglucose[i] < 0){
          turnpoints <- append(turnpoints, glucosevec[i])
          inorder <- TRUE
        }else{
          next
        }
      }
    }
  tpdiff <- diff(turnpoints)
  tpdiff <- tpdiff[abs(tpdiff) > daysd ]
  if(tpdiff[1] > 0){
    tpdiff <- tpdiff[tpdiff > 0]
  }else if(tpdiff[1] < 0){
    tpdiff <- tpdiff[tpdiff < 0]
  }
  mage <- abs(round(mean(tpdiff),2))
  return(mage)
}


tir <- function(cgmtsall, useig = FALSE, bthreshold = 3.9, athreshold = 10, freq){
  coldate <- unique(cgmtsall$timedate)
  tirvec <- c()
  for(d in coldate){
    cgmts <- dplyr::filter(cgmtsall, timedate == d)
    fcgmts <- ""
    if(useig){
      fcgmts <- dplyr::filter(cgmts, imglucose > bthreshold & imglucose < athreshold)
    }else{
      fcgmts <- dplyr::filter(cgmts, sglucose > bthreshold & sglucose < athreshold)
    }
    tirvec <- append(tirvec, round(nrow(fcgmts)/freq,2))
  }
  return(tirvec)
}

tirall <- function(cgmtsall, useig = FALSE, bthreshold = 3.9, athreshold = 10){
  tirvec <- 0
  fcgmts <- ""
  if(useig){
    fcgmts <- dplyr::filter(cgmtsall, imglucose > bthreshold & imglucose < athreshold)
  }else{
    fcgmts <- dplyr::filter(cgmtsall, sglucose > bthreshold & sglucose < athreshold)
  }
  tirvec <- round(nrow(fcgmts)/nrow(cgmtsall),2)
  return(tirvec)
}

modd <- function(cgmtsall, useig = FALSE){
  moddvec <- c()
  coldate <- unique(cgmtsall$timedate)
  for(i in seq_along(coldate)){
    if(i == 1){
      next
    }
    lastday <- dplyr::filter(cgmtsall,timedate == coldate[i-1])
    cuday <- dplyr::filter(cgmtsall,timedate == coldate[i])
    if(useig){
      md <- abs(round(mean(cuday$imglucose - lastday$imglucose),2))
    }else{
      md <- abs(round(mean(cuday$sglucose - lastday$sglucose),2))
    }
    moddvec <- append(moddvec, md)
  }
  return(moddvec)
}


acfcoff <- function(cgmtsall, useig = TRUE, diffnum = 1, interval = 15){
  glucosets = NULL
  if(useig){
    glucosets <- cgmtsall$imglucose
  }else{
    glucosets <- cgmtsall$sglucose
  }
  glucosets <- ts(glucosets, frequency = 1440/interval)
  gacf <- acf(diff(glucosets, differences = diffnum), plot = FALSE)
  return(gacf$acf[2:6])
}

pacfcoff <- function(cgmtsall, useig = TRUE, diffnum = 1, interval = 15){
  glucosets = NULL
  if(useig){
    glucosets <- cgmtsall$imglucose
  }else{
    glucosets <- cgmtsall$sglucose
  }
  glucosets <- ts(glucosets, frequency = 1440/interval)
  gpacf <- pacf(diff(glucosets, differences = diffnum), plot = FALSE)
  return(gpacf$acf[2:6])
}




