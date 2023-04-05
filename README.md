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
Extracts Zhang data and puts into RDS files

Abs and delta files.

## Bias correction
Reads in OBS data and converts to RDS files.
  * *Monthly temperature:* Observation-based present day means use the period 1957-1981. For more information, please see the Climate of Hawaii website; Giambelluca et al. 2014; http://climate.geography.hawaii.edu/. Data are in degC

  * *Monthly precipitation:* The observation-based present day mean maps use the period 1990-2009. These are calculated from the month-year rainfall maps, Frazier et al. 2016; available at http://rainfall.geography.hawaii.edu. Present-day means are calculated for each season: Annual, Wet Season (Nov-Apr), and Dry Season (May-Oct). Means are converted to different units: Rainfall is converted from mm to inches (divide mm by 25.4).
  
Variables created
  * ANNPrecip_Obs
  * ANNPrecipDelta_rcp45/85
  * ANNPrecipBC_rcp45/85...
    * JanPrecip
    * FebPrecip
    * MarPrecip...
    * WetPrecip
    * DryPrecip

Bias correction is Future delta (from Zhang extraction; downscaled future - downscaled obs) + Obs data

## Variable creation
Reads in bias-corrected data, calculates ts files for variables to be read into plotting scripts

## SPI
Only SPI can be calculated b/c evapotranspiration calculations require more sophisticated methods than what can be done with T and P (Thornthwaite), thus ET cannot be calculated from available projection data. 

  * Calculates SPI on stars objects for Obs and each CF.
  * Save mean() stars objects for SPI
  * Reduce ts for each climate zone and whole park
  * Calculate characteristics for each climate zone
  * Plot Maps w/ SPI ts for Obs and each CF
  * Standard panel plot for each zone

# Plotting scripts

## Maps_ts_plots.R
Plots map of whole park for each CF with ts plot that is average of the park.
Input data are Mean RDS for whole park and annual timeseries
- Would it make more sense to do RDS that's annual means for whole park and in script can extract the ts?
  * Annual Tmean
  * Annual Precip





