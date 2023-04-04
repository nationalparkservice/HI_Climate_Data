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
  
# *Xue seasonal attached (dynamical)*
# * daily available from NCAR
# * 10-year period in historical and future (pseudo-global warming method)
# * Only 1 Cf
# * Rain and T2
# 
# - Right now Xue doesn't add anything from Zhang. It's average over whole time periods. Try getting daily files from NCAR to look at extreme precip and temps.

# Data prep scripts
## Zhang extraction
Extracts Zhang data and puts into RDS files

Abs and delta files.

## Bias correction
Reads in OBS data and converts to RDS files.
  * *Monthly temperature:* Observation-based present day means use the period 1957-1981. For more information, please see the Climate of Hawaii website; Giambelluca et al. 2014; http://climate.geography.hawaii.edu/. Data are in degC

  * *Monthly precipitation:* The observation-based present day mean maps use the period 1990-2009. These are calculated from the month-year rainfall maps, Frazier et al. 2016; available at http://rainfall.geography.hawaii.edu. Present-day means are calculated for each season: Annual, Wet Season (Nov-Apr), and Dry Season (May-Oct). Means are converted to different units: Rainfall is converted from mm to inches (divide mm by 25.4).

  * Bias correction is Future delta (from Zhang extraction) (downscaled future - downscaled obs) + Obs data

## Variable creation
Reads in bias-corrected data, calculates RDS files for each variable and ts files for variables to be read into plotting scripts

## SPEI
Calculates SPEI on base data and calculates SPEI on each grid cell


# Plotting scripts

## Maps_ts_plots.R
Plots map of whole park for each CF with ts plot that is average of the park.
Input data are Mean RDS for whole park and annual timeseries
- Would it make more sense to do RDS that's annual means for whole park and in script can extract the ts?
  * Annual Tmean
  * Annual Precip





