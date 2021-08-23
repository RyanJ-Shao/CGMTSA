# The MIT License (MIT)
# Copyright (c) 2020 UCAS

#' Preprocess CGM data.
#' @param inputdir Path of input directory containing CGM files.
#' @param outputdir Path of output directory where preprocessed CGM data will be stored.
#' @param outlierdet Logical. If TRUE the outliers will be detected.
#' @param interval The interval of CGM data.
#' @param imputation Logical. If TRUE the missing data will be imputed.
#' @param imethod The method that will be used to impute missing data.
#' @param maxgap If the missing gap is greater than max gap, the missing gap will be ignore when imputed.
#' @param compeleteday Logical. If TRUE the day with missing data will be filtered.
#' @param removeday Logical. If TRUE the day with missing gap greater than maxgap will be filtered.
#' @param removeday Logical. If TRUE the day with missing gap greater than maxgap will be filtered.
#' @param device Device type: 0 (manual format); 1 (Abbott libre freestyle); 2 (Medtronic ipro2); 3 (Dexcom G6), default 0.
#' @param transunits Logical. If TURE the glucose values will be divided by 18.
#' @param removeflday Logical. If TRUE the data of first and last day will be filter.
#' @export
prepro <- function(inputdir="", outputdir="", outlierdet = TRUE, interval = 15, imputation = FALSE,
                   immethod = "linear", maxgap = 60, compeleteday = TRUE, removeday = FALSE, device = 0, transunits = FALSE, removeflday = TRUE){
	fileNames = list.files(inputdir)
	for(f in fileNames){
	  print(paste("processing file:", f))
	  cgmts <- fformat(fpath = paste(inputdir, f, sep = ''), device = device)
	  colnm <- colnames(cgmts)
	  if(colnm[1] != "timestamp" || colnm[2] != "sglucose" || colnm[3] != "bglucose"){
	    stop(paste("The format fo file '",f ,"' is incorrect and cannot be read.",sep = ""))
	  }
	  cgmts  = qcfun(cgmts, outlierdet, interval, imputation,immethod, maxgap, compeleteday,removeday, transunits, removeflday)
	  write.csv(cgmts, paste(outputdir,f,sep=''),row.names = FALSE)
	  }
}

fformat <- function(fpath, device = 0){
  cgmts = ""
  if(device == 0){
    cgmts <- read.csv(fpath)
    return(cgmts)
  }else if(device == 1){
    cgmts <- read.table(fpath, sep = "\t", skip = 3, encoding = "UTF-8")
    if(length(names(cgmts)) >4 ){
      cgmts <- dplyr::select(cgmts, 2,4,5)
      cgmts[is.na(cgmts$V4),]$V4 <- cgmts[is.na(cgmts$V4),]$V5
      cgmts <- cgmts[!is.na(cgmts$V4),]
      cgmts <- dplyr::select(cgmts, 1,2)
    }else{
      cgmts <- dplyr::select(cgmts, 2,4)
    }
    names(cgmts) <- c("timestamp", "sglucose")
    cgmts <-  dplyr::mutate(cgmts, timestamp = gsub("/","-", cgmts$timestamp))
    cgmts <-  dplyr::mutate(cgmts, bglucose = NA)
  }else if(device == 2){
    cgmts <- read.table(fpath, sep = "\t", skip = 11,header = TRUE)
    cgmts <- dplyr::select(cgmts, 4,10)
    names(cgmts) <- c("timestamp", "sglucose")
    x <- as.character(lubridate::dmy_hms(cgmts$timestamp))
    cgmts$timestamp <- x
    cgmts <-  dplyr::mutate(cgmts, bglucose = NA)
  }else if(device == 3){
    cgmts <- read.table(fpath, sep = ",", skip = 11)
    cgmts <- dplyr::select(cgmts, 2, 8)
    names(cgmts) <- c("timestamp", "sglucose")
    cgmts <-  dplyr::mutate(cgmts, timestamp = gsub("T"," ", cgmts$timestamp))
    x <- as.character(lubridate::ymd_hms(cgmts$timestamp))
    cgmts$timestamp <- x
    cgmts <-  dplyr::mutate(cgmts, bglucose = NA)
  }
  return(cgmts)
}

