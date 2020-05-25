# Continuous Glucose Monitoring Time Series Data Analysis (CGMTSA) User Guide
## 1.	Packages which CGMTSA depends on.
```
dplyr
lubridate
plotly
imputeTS
forecast
TSA
```

## 2.	Install “orca”
CGMTSA package needs orca to export pdf file of CGM plot. Orca is an Electron app that generates images and reports of Plotly things like plotly.js graphs, dash apps, dashboards from the command line. The installation guide of orca can be seen in https://github.com/plotly/orca.
#### Installation guide of orca in MacOS:
- Unzip the mac-release.zip file.
- Double-click on the orca-X.Y.Z.dmg file. This will open an installation window.
- Drag the orca icon into the Applications folder.
- Open finder and navigate to the Applications/ folder.
- Right-click on the orca icon and select Open from the context menu.
- A password dialog will appear asking for permission to add orca to your system PATH.
- Enter your password and click OK.
- This should open an Installation Succeeded window.
- Open a new terminal and verify that the orca executable is available on your PATH.
#### Installation guide of orca in Windows:
- Extract the windows-release.zip file.
- In the release folder, double-click on orca Setup X.Y.Z, this will create an orca icon on your Desktop.
- Right-click on the orca icon and select Properties from the context menu.
- From the Shortcut tab, copy the directory in the Start in field.
- Add this Start in directory to your system PATH (see below).
- Open a new Command Prompt and verify that the orca executable is available on your PATH.
## 3.	License
This project is licensed under The MIT License.
## 4. Install CGMTS Package in R
The package can be installed from GitHub directly:
```
library(“devtools”)
install_github(“RyanJ-Shao/CGMTS”)
```
## 5. Data Format
CGMTSA package accept CSV file as input. The input file includes three columns: "timestamp","sglucose","bglucose", every column are separated by comma. For example:
timestamp |	sglucose | bglucose
-------|--------|--------
2019-12-06 19:50 |	10.9	| NA
2019-12-06 20:05 |	11.9 |	NA
2019-12-06 20:20 |	12.5 |	NA
2019-12-06 20:35 |	13.1 |	NA
2019-12-06 20:50 |	14 |	NA
2019-12-06 21:05 |	14.8 |	NA

The “sglucose” column is glucose from CGM sensor, “bglucose” is glucose from calibration glucose or SMBG. If sglucose is missing, the timestamp of this point still needs to be recorded. It also can directly read original data from Abbott FreeStyle Libre, Dexcom G6 and Medtronic Ipro2.

## 6. Running CGMTS
There are three main functions in CGMTSA: prepro, cgmmetrics and cgmplot. The prepro() function preprocess CGM data, its main function including: read all CSV files from the input directory, detect outliers from CGM data and impute the missing data. The cgmmetrics() function calculate common CGM metrics, such as standard deviation (SD), coefficient of variation (CV), mean glucose, etc. The cgmplot() function generates many different plots of CGM data, such as glucose trace, autocorrelation function plot, partial autocorrelation function plot, etc.
#### (1) Preprocess CGM data
- The prepro() function takes a directory which contains CGM files as input, and export them to an output directory. The prepro() function can directly read original data from Abbott FreeStyle Libre, Dexcom G6 and Medtronic Ipro2 or manually format. If parameters “device” is set to 0, it indicates the file is manually format, 1 indicates the file is FreeStyle Libre format, 2 indicates the file is Dexcom G6, 3 indicates the file is Ipro2 format. In this function we will remove the first and end day of CGM data, because the data quality of these two days is bad in general. The prepro() function includes two main functions: outliers detection and imputation.
- The prepro() function can detect the innovational outlier (IO) and additive outlier (AO). These two types of outliers represent a large portion of outliers likely to be found in practice. The outliers may cause serious bias in estimating autocorrelations, partial autocorrelations and autoregressive moving average parameters. The AO only affect the time when it happened, the IO represents an extraordinary shock in time T influencing the sequence after time T. We use Chang’s method () to detect AO and IO, and mark them in CGM trace plot.
- The prepeo() function will remove incomplete days if the parameters “completeday” is set to TRUE. If “compeleteday” is FALSE and “imputation” is TRUE, then prepro() will use the method in “immethod” to impute the missing gap that is less than “maxgap”. The “immethod” contains “linear”, “seadec” and “arima” to impute missing gap. If “removeday” is TRUE, then the days which missing gap is greater than “maxgap” will removed.
#### (2)	Calculate CGM metrics
- The cgmmetrics () function takes a directory which contains preprocessed CGM files as input, and export the metrics to an output directory. The cgmmetrics() function can calculate the blood glucose metrics that is common and recommended by the consensus, including standard deviation (SD), mean glucose, coefficient of variation (CV), glucose management index(GMI), low blood glucose index (LBGI), high blood glucose index (HBGI), mean amplitude of glycemic excursions (MAGE), time in range (TIR) and mean of daily differences (MODD). The cgmmetrics() can calculate everyday blood glucose metrics and the overall blood glucose metrics. If the “useig” parameters is set to TRUE, this function will use imputed data.
#### (3)	Visualize CGM data
- The cgmplot () function takes a directory which contains preprocessed CGM files as input, and export the plots to an output directory. The cgmplot() function can calculate types of plot, including CGM trace plot, CGM 3 dimensional plot, CGM decomposition plot and ACF, PACF plot. The cgmplot() function will generate HTML file as default except the ACF and PACF plot. If “html” parameters is set to FALSE, it will generate PDF file instead. If the “useig” parameters is set to TRUE, this function will use imputed data.
## 7.	Examples
The CGMTSA package contains an example CGM file, we can load it and test the functions of CGMTSA package on it. First, use prepro() function to detect outliers and impute missing in CGM data.
1. Manual format
```
library(CGMTS)
datadir <- system.file("extdata", package = "CGMTS")
prepro(inputdir=paste(datadir, "/manualFormat/", sep = ""), outputdir, outlierdet = TRUE, interval = 15, imputation = TRUE, immethod = "linear", maxgap = 60, compeleteday = FALSE, removeday = FALSE, device = 0, transunits = FALSE, removeflday = TRUE)
cgmmetrics(inputdir, outputdir ,useig = TRUE, threshold =1, bthreshold = 3.9, athreshold = 10, interval = 15)
cgmplot(inputdir, outputdir, useig= TRUE, markoutliers= TRUE, interval = 15, diffnum = 1, html = TRUE)
```
2. Freestyle Libre format
```
prepro(inputdir=paste(datadir, "/FreeStyleLibre/", sep = ""), outputdir, outlierdet = TRUE, interval = 15, imputation = TRUE, immethod = "linear", maxgap = 60, compeleteday = FALSE, removeday = FALSE, device = 1, transunits = FALSE, removeflday = TRUE)
cgmmetrics(inputdir, outputdir ,useig = TRUE, threshold =1, bthreshold = 3.9, athreshold = 10, interval = 15)
cgmplot(inputdir, outputdir, useig= TRUE, markoutliers= TRUE, interval = 15, diffnum = 1, html = TRUE)
```
3. Medtronic ipro2 format
```
prepro(inputdir=paste(datadir, "/Medtronicipro2/", sep = ""), outputdir, outlierdet = TRUE, interval = 5, imputation = TRUE, immethod = "linear", maxgap = 60, compeleteday = FALSE, removeday = FALSE, device = 2, transunits = TRUE, removeflday = TRUE)
cgmmetrics(inputdir, outputdir ,useig = TRUE, threshold =1, bthreshold = 3.9, athreshold = 10, interval = 5)
cgmplot(inputdir, outputdir, useig= TRUE, markoutliers= TRUE, interval = 5, diffnum = 1, html = TRUE)
```
4. Dexcom G6 format
```
prepro(inputdir=paste(datadir, "/DexcomG6/", sep = ""), outputdir, outlierdet = TRUE, interval = 5, imputation = TRUE, immethod = "linear", maxgap = 60, compeleteday = FALSE, removeday = FALSE, device = 3, transunits = TRUE, removeflday = TRUE)
cgmmetrics(inputdir, outputdir ,useig = TRUE, threshold =1, bthreshold = 3.9, athreshold = 10, interval = 5)
cgmplot(inputdir, outputdir, useig= TRUE, markoutliers= TRUE, interval = 5, diffnum = 1, html = TRUE)
```
