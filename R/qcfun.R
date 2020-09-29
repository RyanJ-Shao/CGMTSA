# The MIT License (MIT)
# Copyright (c) 2020 University of Chinese Academy of Sciences
#
#' Quality Control for CGM data, include imputation and outlier detection

qcfun<- function(cgmts, outlierdet = TRUE, interval = 15, imputation = FALSE, immethod = "linear",
                 maxgap = 60, compeleteday = TRUE, removeday = FALSE, transunits = FALSE, removeflday = TRUE){
  cgmts <- cgmts[order(lubridate::ymd_hms(cgmts$timestamp)),]
  vectimestamp <- as.vector(cgmts$timestamp)
  vectimestamp <- unlist(strsplit(vectimestamp,split=" "))
  maxtimestamp <- matrix(vectimestamp,ncol=2,byrow=T)[,1]
  cgmts <-  dplyr::mutate(cgmts, timedate = maxtimestamp)
  coldate <- unique(cgmts$timedate)
  freq = 1440/interval
  #remove first day and last day
  fday <- coldate[1]
  lday <- coldate[length(coldate)]
  if(removeflday){
    cgmts[cgmts$timedate ==fday,]$timedate <- NA
    cgmts <- cgmts[!is.na(cgmts$timedate),]
    cgmts[cgmts$timedate ==lday,]$timedate <- NA
    cgmts <- cgmts[!is.na(cgmts$timedate),]
  }

  if(transunits){
    cgmts <- dplyr::mutate(cgmts, sglucose = round(sglucose/18,2))
  }

  coldate <- unique(cgmts$timedate)
  if(compeleteday){
    for (d in coldate){
      if (length(cgmts[cgmts$timedate == d,]$timedate) < freq){
        cgmts[cgmts$timedate ==d,]$timedate <- NA
        cgmts <- cgmts[!is.na(cgmts$timedate),]
      }
    }
    cgmts <- cgmts[!is.na(cgmts$timedate),]
  }else if(imputation){
    if(immethod == "linear"){
      print("linear imputation")
      cgmts <- dplyr::mutate(cgmts, imglucose = imputeTS::na.interpolation(cgmts$sglucose,maxgap = as.integer(maxgap/interval)))
    }else if(immethod == "seadec"){
      print("SEADEC imputation")
      tstype <- ts(cgmts$sglucose,frequency = freq)
      tstype <- imputeTS::na.seadec(tstype,algorithm = "interpolation",maxgap = as.integer(maxgap/interval))
      cgmts <- dplyr::mutate(cgmts, imglucose = tstype)
    }else if(immethod == "arima"){
      print("ARIMA imputation")
      cgmts <- dplyr::mutate(cgmts, imglucose = imputeTS::na.kalman(cgmts$sglucose,model = "auto.arima",maxgap = as.integer(maxgap/interval)))

    }
    if(removeday == TRUE){
      is.na.rle <- rle(is.na(cgmts$imglucose))
      is.na.rle$values <- is.na.rle$values & is.na.rle$lengths > as.integer(maxgap/interval)
      reov = cgmts[inverse.rle(is.na.rle), ]
      reovdate = unique(reov$timedate)
      for (d in reovdate){
        cgmts[cgmts$timedate ==d,]$timedate <- NA
        cgmts <- cgmts[!is.na(cgmts$timedate),]
      }

    }
  }else{
    if(removeday == TRUE){
      is.na.rle <- rle(is.na(cgmts$sglucose))
      is.na.rle$values <- is.na.rle$values & is.na.rle$lengths > as.integer(maxgap/interval)
      reov = cgmts[inverse.rle(is.na.rle), ]
      reovdate = dplyr::unique(reov$timedate)
      for (d in reovdate){
        cgmts[cgmts$timedate ==d,]$timedate <- NA
      }
      cgmts <- cgmts[!is.na(cgmts$timedate),]
    }
  }
  if(outlierdet == TRUE){
    cgmts <- dplyr::mutate(cgmts, outliers = NA)
    unidate = unique(cgmts$timedate)
    for (d in unidate){
      udcgm <- cgmts[cgmts$timedate ==d,]
      if(!any(is.na(udcgm$sglucose))){
        udcgmglucose <- udcgm$sglucose
        udcgmts <- ts(udcgmglucose, frequency = freq)
        tsmodel <- forecast::auto.arima(udcgmts)
        ao <- TSA::detectAO(tsmodel)
        io <- TSA::detectIO(tsmodel)
        aodict <- c(ao$lambda2)
        names(aodict) <- ao$ind

        iodict <- c(io$lambda1)
        names(iodict) <- io$ind
        indinc <- intersect(ao$ind, io$ind)
        #remove ao or io index according to theri lambda
        for(ic in indinc){
          if(aodict[[as.character(ic)]] < iodict[[as.character(ic)]]){
            removeind = which(c(ic) %in% names(aodict))
            aodict <- aodict[-removeind]
          }else{
            removeind <- which(c(ic) %in% names(iodict) )
            iodict <- iodict[-removeind]
          }
        }

        for(i in seq_along(aodict)){
          ind <- as.numeric(names(aodict)[i])
          udcgm[c(ind),]$outliers <- "AO"
        }
        for(i in seq_along(iodict)){
          ind <- as.numeric(names(iodict)[i])
          udcgm[c(ind),]$outliers <- "IO"
        }
        cgmts[cgmts$timedate ==d,]$outliers <- udcgm$outliers
      }
    }

  }
  return(cgmts)
}
