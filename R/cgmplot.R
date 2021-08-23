# The MIT License (MIT)
# Copyright (c) 2020 UCAS

#' Generate plots of CGM data.
#' @param inputdir Path of input directory containing CGM files.
#' @param outputdir Path of output directory where plots will be stored.
#' @param useig Logical. If TRUE the imputed data will be used.
#' @param markoutliers Logical. If TRUE the outliers will be labeled in different colors.
#' @param interval The interval of CGM data.
#' @param diffnum The number of differencing.
#' @param html Logical.If TRUE the interactive plot will be exported, else the static pdf plot will exported.
#' @export

cgmplot <- function(inputdir, outputdir, useig= FALSE, markoutliers= TRUE, interval = 15,
                    diffnum = 1, html = TRUE){
  fileNames = list.files(inputdir)
  for(f in fileNames){
    fname = unlist(strsplit(f, split = "\\."))[1]
    print(paste("processing file:", f))
    cgmtsall = read.csv(paste(inputdir, "/", f, sep = ''),stringsAsFactors= FALSE)
    vectimestamp <- as.vector(cgmtsall$timestamp)
    vectimestamp <- unlist(strsplit(vectimestamp,split=" "))
    maxtimestamp <- matrix(vectimestamp,ncol=2,byrow=T)[,1]
    cgmtsall <- dplyr::mutate(cgmtsall, timedate = maxtimestamp)
    #print(head(cgmtsall))
    print("plottinig ACF")
    cgmACF(cgmtsall, fname, outputdir, useig, diffnum, interval)
    print("plottinig PACF")
    cgmPACF(cgmtsall, fname, outputdir, useig, diffnum, interval)
    print("plottinig 3d")
    cgm3d(cgmtsall, fname, outputdir, useig, interval)
    print("plottinig decom")
    cgmdecom(cgmtsall, fname, outputdir, useig,interval, html = html)
    print("plottinig trace")
    cgmtrace(cgmtsall, fname, outputdir,useig, markoutliers, html = html)
  }

}

cgmACF <- function(cgmtsall, fname, outputdir, useig = TRUE, diffnum = 1, interval = 15){
  glucosets = NULL
  if(useig){
    glucosets <- cgmtsall$imglucose
  }else{
    glucosets <- cgmtsall$sglucose
  }
  glucosets <- ts(glucosets, frequency = 1440/interval)
  pdf(paste(outputdir, fname,"_","acf", ".pdf",sep = ""))
  adftest <- tseries::adf.test(diff(glucosets, differences = diffnum))
  pv <- adftest$p.value
  title = paste("ACF, ADF p_value:", as.character(pv), sep = "")
  acf(diff(glucosets, differences = diffnum), main = title)
  dev.off()
}


cgmPACF <- function(cgmtsall, fname, outputdir, useig = TRUE, diffnum = 1, interval = 15){
  glucosets = NULL
  if(useig){
    glucosets = cgmtsall$imglucose
  }else{
    glucosets = cgmtsall$sglucose
  }
  glucosets = ts(glucosets, frequency = 1440/interval)
  pdf(paste(outputdir, fname,"_","pacf", ".pdf",sep = ""))
  adftest <- tseries::adf.test(diff(glucosets, differences = diffnum))
  pv <- adftest$p.value
  title = paste("PACF, ADF p_value:", as.character(pv), sep = "")
  pacf(diff(glucosets, differences = diffnum), main= title)
  dev.off()
}

cgm3d <- function(cgmtsall, fname, outputdir, useig = TRUE,interval){
  freq = 1440/interval
  hm <- merge(0:23, seq(0, 60 - interval, by = interval))
  x <- data.frame('INTERVAL' = chron::chron(time = paste(hm$x, ':', hm$y, ':', 0)))
  x <- data.frame('INTERVAL' = x[order(x$INTERVAL), ])
  x <- as.character(x$INTERVAL)

  if(useig){
    z = matrix(cgmtsall$imglucose,ncol=freq,byrow=T)
  }else{
    z = matrix(cgmtsall$sglucose,ncol=freq,byrow=T)
  }
  fig <- plotly::plot_ly(x = x, z = z)
  fig <- plotly::add_surface(fig)
  fig <- plotly::layout(fig,
      scene = list(
        xaxis = list(
        title = "Time",
        dtick = 10,
        tick0 = 0
        #tickmode = "array",
        #type = "date",
        #tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Days"
      ),
      zaxis = list(title = "Glucose Values")
    )
    )
  htmlwidgets::saveWidget(fig, paste(outputdir, fname,"_","3d", ".html",sep = ""))
}


cgmdecom <- function(cgmtsall, fname, outputdir, useig = TRUE,interval = 15, html = FALSE){
  freq = 1440/interval
  uniday = unique(cgmtsall$timedate)
  #remove uncompelete day
  for(d in uniday){
    if(useig){
      if(any(is.na(cgmtsall[cgmtsall$timedate ==d,]$imglucose))){
        cgmtsall[cgmtsall$timedate ==d,]$timedate <- NA
      }
    }else{
      if(any(is.na(cgmtsall[cgmtsall$timedate ==d,]$sglucose))){
        cgmtsall[cgmtsall$timedate ==d,]$timedate <- NA
      }
    }
  }
  cgmtsall <- cgmtsall[!is.na(cgmtsall$timedate),]
  uniday = unique(cgmtsall$timedate)
  if(useig){
    gts <- ts(cgmtsall$imglucose, frequency = freq)
  }else{
    gts <- ts(cgmtsall$sglucose, frequency = freq)
  }
  stlgts <- stl(gts,s.window = "periodic")
  seasonal <- stlgts$time.series[,1]
  trend <- stlgts$time.series[,2]
  remainder <- stlgts$time.series[,3]
  #plot seasonal
  seafig <- plotly::plot_ly(
    type = "scatter",
    x = c(cgmtsall$timestamp),
    y = c(seasonal),
    mode = "lines"
  )
  seafig <- plotly::layout(seafig,
      xaxis = list(
        title = "Time",
        dtick = 10,
        tick0 = 0,
        tickmode = "array",
        type = "date",
        tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Glucose seansonal component"
      )
    )
  #plot trend
  trfig <- plotly::plot_ly(
    type = "scatter",
    x = c(cgmtsall$timestamp),
    y = c(trend),
    mode = "lines"
  )
  trfig <- plotly::layout(trfig,
      xaxis = list(
        title = "Time",
        dtick = 10,
        tick0 = 0,
        tickmode = "array",
        type = "date",
        tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Glucose trend component"
      )
    )

  #plot trend
  refig <- plotly::plot_ly(
    type = "scatter",
    x = c(cgmtsall$timestamp),
    y = c(remainder),
    mode = "lines"
  )
  refig <- plotly::layout(refig,
      xaxis = list(
        title = "Time",
        dtick = 10,
        tick0 = 0,
        tickmode = "array",
        type = "date",
        tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Glucose remainder component"
      )
    )

  if(html){
    htmlwidgets::saveWidget(seafig, paste(outputdir,fname,"_","seasonal", ".html",sep = ""))
    htmlwidgets::saveWidget(trfig, paste(outputdir, fname,"_","trend", ".html",sep = ""))
    htmlwidgets::saveWidget(refig, paste(outputdir, fname,"_","remainder", ".html",sep = ""))
  }else{
    oldworkdir = getwd()
    setwd(outputdir)
    plotly::orca(seafig, paste(fname,"_","seasonal", ".pdf",sep = ""))
    plotly::orca(trfig, paste( fname,"_","trend", ".pdf",sep = ""))
    plotly::orca(refig, paste(fname,"_","remainder", ".pdf",sep = ""))
    setwd(outputdir)
  }

  #orca(seafig, paste(outputdir,"/", fname,"_","seasonal", ".pdf",sep = ""))
  #orca(trfig, paste(outputdir,"/", fname,"_","trend", ".pdf",sep = ""))
  #orca(refig, paste(outputdir,"/", fname,"_","remainder", ".pdf",sep = ""))

  #return(seafig)
}




#fname = ryan
#cgmtrace(cgmts, "ryan", "Desktop/cgm_software/CGMTS/plotoutput")
cgmtrace <- function(cgmtsall, fname, outputdir,useig = TRUE, markoutliers = TRUE, html = FALSE){

  cgmdate <- unique(cgmtsall$timedate)
  for (d in cgmdate){
    cgmts <- dplyr::filter(cgmtsall, timedate == d)
    vectimestamp <- as.vector(cgmts$timestamp)
    vectimestamp <- unlist(strsplit(vectimestamp,split=" "))
    hms <- matrix(vectimestamp,ncol=2,byrow=T)[,2]
    gts <- c()
    gy <-hms
    if(useig){
      gts <- cgmts$imglucose
    }else{
      gts <- cgmts$sglucose
    }

    fig <- plotly::plot_ly(
      type = "scatter",
      x = c(gy),
      y = c(gts),
      name = 'Glucose Trace',
      mode = "lines"
      )

    if(markoutliers){
      out <- dplyr::filter(cgmts, outliers == "IO" |  outliers == "AO" )
      for(i in 1:nrow(out)){
        y=0
        x = ""
        if(useig){
          y = out[i,]$imglucose
        }else{
          y = out[i,]$sglucose
        }
        fig <- plotly::add_trace(fig,
            type = "scatter",
            x = unlist(strsplit(out[i,]$timestamp, split = " "))[2],
            y = y,
            name = out[i,]$outliers,
            mode = "markers"
          )
      }
    }


    fig <- plotly::layout(fig,
        xaxis = list(
          title = "Time",
          dtick = 10,
          tick0 = 0
        ),
        yaxis = list(
          title = "Glucose value(mmol/L)"
        )
        )
    if(html){
      htmlwidgets::saveWidget(fig, paste(outputdir, fname,"_", d , ".html",sep = ""))
    }else{
      oldworkdir = getwd()
      setwd(outputdir)
      plotly::orca(fig, paste(fname, "_", d , ".pdf",sep = ""))
      setwd(oldworkdir)
    }

  }
}
