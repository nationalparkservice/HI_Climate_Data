# HI_Climate_Data
Scripts to access and create climate futures from HI climate data sets

# Climate Data sets
Climate data were provided by Dr. Abby Frazier in 2022. Descriptions of the data sets are provided below.

*Zhang (dynamical)*
* ncdf difficult to work with and don't behave like normal ncdf - can send monthly rasters
* name of years is wrong in historical (should be 1990-2009 but left as 2080-2099)
* raw nc
* processed
* For each island present (1990-2009), rcp45 and 85 (2080-2099)

# Data prep scripts
## Zhang extraction
_Zhang_extraction.R_
Extracts Zhang data and puts into RDS files. Requires external HD with Zhang files to run. Output goes to local park file.

## Summarize rainfall observations
_Mean_Rainfall_Obs.R_
Takes data from Rainfall.Atlas and summarizes by year and month for bias correction process. Output is geoTIFF files to be read in to BC script

## Bias correction
_Bias_correction.R_
Reads in OBS data and converts to RDS files. Obs data stored on external HD and output goes to local park file.
  * *Monthly temperature:* Observation-based present day means use the period 1957-1981. For more information, please see the Climate of Hawaii website; Giambelluca et al. 2014; http://climate.geography.hawaii.edu/. Data are in degC

  * *Monthly precipitation:* The observation-based present day mean maps use the period 1990-2009. These are calculated from the month-year rainfall maps, Frazier et al. 2016; available at http://rainfall.geography.hawaii.edu. Present-day means are calculated for each season: Annual, Wet Season (Nov-Apr), and Dry Season (May-Oct). Means are converted to different units: Rainfall is converted from mm to inches (divide mm by 25.4).
  
Variables created
  * ANNPrecip_Obs (output are in mm)
  * ANNPrecipDelta_rcp45/85
  * ANNPrecipBC_rcp45/85...
    * JanPrecip 
    * FebPrecip
    * MarPrecip...
    * WetPrecip
    * DryPrecip
  * ANNTmean... (Output are in Fahrenheit)

Bias correction is Future delta (from Zhang extraction; downscaled future - downscaled obs) + Obs data

Wet/Dry season BC for plots
  * Wet is Nov-Apr, Dry is May-Oct per Zhang README from Frazier
    * Precip is mean delta of each month

## Variable creation
_Variable_creation.R_
Reads in bias-corrected data, calculates ts files for variables to be read into plotting scripts
  * Output: Zone-Monthly-data.csv (Output are in Inches and Fahreneit)

## SPI
_....R_
Only SPI can be calculated b/c evapotranspiration calculations require more sophisticated methods than what can be done with T and P (Thornthwaite), thus ET cannot be calculated from available projection data. 

  * Calculates SPI on stars objects for Obs and each CF.
  * Save mean() stars objects for SPI
  * Reduce ts for each climate zone and whole park
  * Calculate characteristics for each climate zone
  * Plot Maps w/ SPI ts for Obs and each CF
  * Standard panel plot for each zone

# Plotting scripts

## Delta maps for whole park with ts avg of park
_Maps_ts_plots.R_
Plots map of whole park for each CF with ts plot that is average of the park.
Input data are Mean RDS for whole park and annual timeseries
  * Annual Tmean (ANNTmean)
  * Annual Precip (ANNPrecip)
  * Wet/Dry Tmean/Precip (WetTmean) (DryPrecip)

## Wet/Dry seaonal maps and dotplots by zone
_Maps_seasonal_dotplots.R_
  * Tmean and Precip 
  
  





