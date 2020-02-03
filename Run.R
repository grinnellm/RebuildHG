##### Header #####
# 
# Author:       Matthew H. Grinnell
# Affiliation:  Pacific Biological Station, Fisheries and Oceans Canada (DFO) 
# Group:        Quantitative Assessment Methods Section
# Address:      3190 Hammond Bay Road, Nanaimo, BC, Canada, V9T 6N7
# Contact:      e-mail: Matthew.Grinnell@dfo-mpo.gc.ca | tel: (250) 756.7055
# Project:      Herring
# Code name:    Run.R
# Version:      1.0
# Date started: Mar 9, 2018
# Date edited:  Mar 8, 2019
# 
# Overview: 
# Source the 'SpatialAnalysis.R' script for each requested region.
# 
# Requirements: 
# Output from the 'Summary.R' script.
# 
# Notes: 
# 
#
# References:
# 

##### Housekeeping #####

# General options
rm( list=ls( ) )      # Clear the workspace
sTime <- Sys.time( )  # Start the timer
graphics.off( )       # Turn graphics off

# Install missing packages and load required packages (if required)
UsePackages <- function( pkgs, locn="https://cran.rstudio.com/" ) {
  # Reverse the list 
  rPkgs <- rev( pkgs )
  # Identify missing (i.e., not yet installed) packages
  newPkgs <- rPkgs[!(rPkgs %in% installed.packages( )[, "Package"])]
  # Install missing packages if required
  if( length(newPkgs) )  install.packages( newPkgs, repos=locn )
  # Loop over all packages
  for( i in 1:length(rPkgs) ) {
    # Load required packages using 'library'
    eval( parse(text=paste("suppressPackageStartupMessages(library(", rPkgs[i], 
                "))", sep="")) )
  }  # End i loop over package names
}  # End UsePackages function

# Make packages available
UsePackages( pkgs=c("tidyverse", "sp", "scales", "ggforce", "lubridate", 
        "cowplot", "GGally", "magick", "ggrepel", "readxl", "xtable", 
        "viridis", "zoo") )

##### Controls #####

# Select region(s): major (HG, PRD, CC, SoG, WCVI); or minor (A27, A2W)
spRegions <- c( "A2W" )

# File name for dive transect XY
diveFN <- file.path( "Data", "dive_transects_with_lat_long_June2_2017.xlsx" )

# File name for q parameters
qFN <- file.path( "Data", "qPars.csv" )

# # Model: AM1 and/or AM2 (Note: This has not been tested with > 1 model)
# mNames <- "AM2"

# Generate GIF---this can take a long time---use 64 bit R.
makeGIF <- FALSE

# Reference years
refYrsAll <- read_csv( file=
        "SAR, Start, End
        HG, 1951, 2018 
        PRD, 1951, 2018 
        CC, 1951, 2018 
        SoG, 1951, 2018 
        WCVI, 1990, 1999
        A27, 1951, 2018
        A2W, 1951, 2018", 
    col_types=cols("c", "i", "i") )

##### Parameters #####

# Spawn index threshold (tonnes; NA for none)
siThreshold <- NA  # 15000

# Minimum number of consecutive years
nYrsConsec <- 3

# Buffer distance (m; to include locations that are outside the region polygon)
maxBuff <- 10000

# Intended harvest rate
intendU <- 0.2

# First year of intended harvest rate
intendUYrs <- 1983

# Plot quality (dots per inch)
pDPI <- 600

##### Functions #####

# Load helper functions
source( file=file.path("..", "HerringFunctions", "Functions.R") )

# Latex bold (e.g., for table column names)
boldLatex <- function( x )  paste( '\\textbf{', x, '}', sep ='' )

# Latex math (e.g., for table contents)
mathLatex <- function( x )  paste( '$', x, '$', sep ='' )

# Latex subscript in math mode (e.g., change q1 to q_1)
subLatex <- function( x )  
  gsub( pattern="([a-zA-Z]+)([0-9]+)", replacement='\\\\mli{\\1}_{\\2}', x=x )

##### Data #####

# Load q parameters (from the latest assessment)
qPars <- read_csv( file=qFN, col_types=cols() )

# Region names
allRegionNames <- list( 
    major=c("Haida Gwaii (HG)", "Prince Rupert District (PRD)", 
        "Central Coast (CC)", "Strait of Georgia (SoG)", 
        "West Coast of Vancouver Island (WCVI)"), 
    minor=c("Area 27 (A27)", "Area 2 West (A2W)") )

# Cross-walk table for SAR to region and region name
regions <- read_csv(file=
        "SAR, Region, RegionName, Major
        1, HG, Haida Gwaii, TRUE
        2, PRD, Prince Rupert District, TRUE
        3, CC, Central Coast, TRUE
        4, SoG, Strait of Georgia, TRUE
        5, WCVI, West Coast of Vancouver Island, TRUE
        6, A27, Area 27, FALSE
        7, A2W, Area 2 West, FALSE",
    col_types=cols("i", "c", "c", "l") )

# Determine if the region is minor or major
isMajor <- regions[regions$Region ==spRegions,]$Major

##### Main ##### 

# Message
cat( "Investigate", length(spRegions), "region(s):", PasteNicely(spRegions), 
    "\n" )

# Start a loop over region(s)
for( reg in 1:length(spRegions) ) {
  # Get the ith region
  region <- spRegions[reg]
  # Spatial unit: Region, StatArea, Section, or Group
  if( region == "HG" )    spUnitName <- "Group"
  if( region == "PRD" )   spUnitName <- "StatArea"
  if( region == "CC" )    spUnitName <- "StatArea"
  if( region == "SoG" )   spUnitName <- "Region"
  if( region == "WCVI" )  spUnitName <- "StatArea"
  if( region == "A27" )   spUnitName <- "StatArea"
  if( region == "A2W" )   spUnitName <- "Region"
  # Extract reference years for required region
  refYrs <- refYrsAll %>%
      filter( SAR==region )
  # Error if no reference years present
  if( nrow(refYrs) == 0 ) 
    stop( "Specify reference years for biomass threshold (refYrs): ", region,
        call.=FALSE )
  # Message re spatial info
  cat( "\nInvestigate", region, "by", spUnitName, "\n" )
  # Run the spatial analysis
  source( file="SpatialAnalysis.R" )
}  # End reg loop over region(s)

##### Tables #####

##### Output #####

# Save the workspace image 
save.image( file="Image.RData" ) 

##### End ##### 

# Print end of file message and elapsed time
cat( "\nEnd of file Run.R: ", sep="" ) ;  print( Sys.time( ) - sTime )
