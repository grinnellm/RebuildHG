##### Header #####
# 
# Author:       Matthew H. Grinnell
# Affiliation:  Pacific Biological Station, Fisheries and Oceans Canada (DFO) 
# Group:        Quantitative Assessment Methods Section
# Address:      3190 Hammond Bay Road, Nanaimo, BC, Canada, V9T 6N7
# Contact:      e-mail: Matthew.Grinnell@dfo-mpo.gc.ca | tel: (250) 756.7198
# Project:      Herring
# Code name:    Rebuild.R
# Version:      1.0
# Date started: Feb 3, 2020
# Date edited:  Feb 4, 2020
# 
# Overview: 
# Plot catch, biosample, and spawn data on a finer spatial scale (i.e.,
# statistical area, section, location code, or group).
# 
# Requirements: 
# The saved image from 'Summary.R', and the saved image from the stock
# assessment report 'herringsr'.
# 
# Notes: 
# Averages and other statistics are for the locations that have values. For
# example, 'LayersMean' in 'siYr' is the average number of egg layers for
# all the locations that have spawn reported. That is to say, it excludes
# locations with no spawn (which have zero layers).
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

# Select region: major (HG, PRD, CC, SoG, WCVI); or minor (A27, A2W)
region <- c( "A2W" )

# Spatial unit: Region, StatArea, Section, or Group
spUnitName <- "Region"

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

##### Files #####

# File name for dive transect XY
diveFN <- file.path( "Data", "dive_transects_with_lat_long_June2_2017.xlsx" )

# File name for q parameters
qFN <- file.path( "Data", "qPars.csv" )

# File name for reference years
refYrsAll <- refFN <- file.path("Data", "RefYrs.csv")

##### Functions #####

# Load helper functions
source( file=file.path("..", "HerringFunctions", "Functions.R") )

##### Data #####

# Load q parameters (from the latest assessment)
qPars <- read_csv( file=qFN, col_types=cols() ) %>%
  filter( Region == region ) %>%
  rename( qLower=Lower, qMedian=Median, qUpper=Upper ) %>%
  select( qLower, qMedian, qUpper, Survey )

# Load reference years
refYrs <- read_csv( file=refFN, col_types=cols("c", "i", "i") ) %>%
  filter( Region == region )

# Error if no reference years present
if( nrow(refYrs) == 0 ) 
  stop( "Specify reference years for biomass threshold (refYrs): ", region,
        call.=FALSE )

# If old directory exists
if( region %in% list.files() ) {
  # Remove the old directory
  unlink( region, recursive=TRUE )
  # Warning: remove previous summary output
  warning( "Removed existing directory '", region, "'", call.=FALSE )
  # Create the main directory for output
  dir.create( region )
} else {  # End if directory exists, otherwise
  # Create the main directory for output
  dir.create( region )
}  # End if directory doesn't exists

# Get data from the data summary
LoadSavedObjects <- function( loc ) {
  # Load the saved object
  load( file=file.path("..", "DataSummaries", loc,
                       paste("Image.", loc, ".RData", sep="")) )
  # load( file=file.path("Data", paste("Image.", loc, ".RData", sep="")) )
  # Confirm the region (from the saved data -- should match!)
  if( regName != region )
    stop( "Region mismatch -- check saved image", call.=FALSE )
  # Return objects to the main environment
  newSurvYr <<- newSurvYr  # Year in which survey type changed
  yrRange <<- yrRange  # Range of years to consider
  inCRS <<- inCRS  # Coordinate reference system (input)
  outCRS <<- outCRS  # Coordinate reference system (output)
  shapes <<- shapes  # Shapefiles etc
  figWidth <<- figWidth*1.2  # Figure size
  catch <<- catch  # Catch data
  UpdateCatchData <<- UpdateCatchData  # Function to update catch data
  CalcWeightAtAge <<- CalcWeightAtAge  # Function to calculate weight-at-age
  CalcLengthAtAge <<- CalcLengthAtAge  # Function to calculate length-at-age
  CalcBiomassSOK <<- CalcBiomassSOK  #  Function to calculate biomass from SOK
  inclTestCatch <<- inclTestCatch  # Include catch from test fishery
  geoProj <<- geoProj  # Text for maps: projection
  areas <<- areas  # Spatial info
  bio <<- bio  # Biological data
  spawnRaw <<- spawnRaw  # Spawn data
  catchRaw <<- catchRaw  # Raw catch data (for SOK harvest)
  ciLevel <<- ciLevel  # Confidence interval
  ageRange <<- ageRange  # Ages
  qYrs <<- qYrs  # Survey years (q1 and q2)
  smLine <<- smLine  # Smoothing function for line
  nRoll <<- nRoll  # Number of years for rolling mean
  ageShow <<- ageShow  # Age to highlight on the x-at-age plots
  yrBreaks <<- yrBreaks  # Years to show in x-axes
  convFac <<- convFac  # Unit conversion factors
  parsProd <<- parsProd  # Productivity parameters
  ECF <<- ECF  # Egg conversion factor
}  # End LoadSavedObjects function

# Load saved data (directly!)
LoadSavedObjects( loc=region )

# Get small area table
# TODO: Are both of these required? aSmall and areasSM
aSmall <- areas %>%
  select( Region, RegionName, StatArea, Group, Section ) %>%
  distinct( )

# Calculate SOK harvest
harvest <- catchRaw %>%
  filter( DisposalCode == 2, Source=="SOK" ) %>%
  mutate( BiomassSOK=CalcBiomassSOK(SOK=Catch*convFac$lb2kg, 
                                    eggKelpProp=parsProd$eggKelpProp, 
                                    eggBrineProp=parsProd$eggBrineProp, 
                                    eggWt=parsProd$eggWt, ECF=ECF),
          HarvSOK=Catch*convFac$lb2kg/1000 ) %>%
  left_join( y=aSmall, by=c("Region", "StatArea", "Section") ) %>%
  select( Year, Region, StatArea, Group, Section, BiomassSOK, HarvSOK )

# Function to load transect spatial info
LoadTransectXY <- function( loc ) {
  # Load the data and wrangle
  dat <- read_excel( path=loc, sheet=1 ) %>%
    rename( LocationCode=LOC_CODE ) %>%
    mutate( LocationCode=as.integer(LocationCode) ) %>%
    group_by( LocationCode ) %>%
    summarise( 
      Longitude=MeanNA(c(StartLong, MidLong, EndLong)),
      Latitude=MeanNA(c(StartLat, MidLat, EndLat)) ) %>%
    ungroup( ) %>%
    filter( !is.na(Longitude), !is.na(Latitude), LocationCode!=0 )
  # Grab the spatial info (X and Y)
  locSP <- dat %>%
    transmute( X=Longitude, Y=Latitude )
  # Put X and Y into a spatial points object
  locPts <- SpatialPointsDataFrame( coords=locSP, 
                                    data=data.frame(LocationCode=dat$LocationCode), proj4string=CRS(inCRS) )
  # Convert X and Y from WGS to Albers
  locPtsAlb <- spTransform( x=locPts, CRSobj=CRS(outCRS) )
  # Extract spatial info
  dfAlb <- as_tibble( locPtsAlb ) %>%
    rename( Eastings=X, Northings=Y )
  # Return the data
  return( dfAlb )
}  # End LoadTransectXY function

# Load 'auxiliary' dive transect spatial data
transectXY <- LoadTransectXY( loc=diveFN )

# Get a short list of areas: sections and groups
areasSm <- areas %>%
  select( Region, StatArea, Group, Section ) %>%
  distinct( )

# Get spawn index data
GetSI <- function( allSI, loc, XY ) { 
  # Get a subset and wrangle
  raw <- allSI %>%
    replace_na( list(Group="Other") ) %>%
    mutate( Decade=GetDecade(Year), Area=Length*Width, 
            Survey=ifelse(Year < newSurvYr, "Surface", "Dive"),
            YrsSurv=ifelse(Year < newSurvYr, length(yrRange[yrRange<newSurvYr]),
                           length(yrRange[yrRange>=newSurvYr])) ) %>%
    rowwise( ) %>%
    mutate( LyrsMean=MeanNA(c(SurfLyrs, MacroLyrs, UnderLyrs)),
            SITotal=SumNA(c(SurfSI, MacroSI, UnderSI)) ) %>%
    ungroup( ) %>%
    group_by( Decade ) %>%
    mutate( YrsDecade=length(unique(Year)) ) %>%
    ungroup( )
  # Clip the extent
  df <- ClipExtent( dat=raw, spObj=shapes$regSPDF, bufDist=maxBuff, 
                    silent=TRUE )
  # Subset data with 'good' X and Y
  dfNotNA <- df %>%
    filter( !is.na(Eastings) & !is.na(Northings) )
  # Subset data with 'bad' X or Y, and try to fill in using transect X and Y
  dfNA <- df %>%
    filter( is.na(Eastings) | is.na(Northings) ) %>%
    select( -Eastings, -Northings ) %>%
    left_join( y=XY, by="LocationCode" )
  # Re-combine the two subsets
  df2 <- bind_rows( dfNotNA, dfNA )
  # Clip the extent (again)
  res <- ClipExtent( dat=df2, spObj=shapes$regSPDF, bufDist=maxBuff )
  # Stop if we're missing rows
  if( nrow(raw) != nrow(res) )  stop( "Missing rows!", call.=FALSE )
  # Return the spatial data
  return( res )
}  # End GetSI function

# Get spawn index
siAll <- GetSI( allSI=spawnRaw, loc=region, XY=transectXY )

##### Privacy #####

# Load catch and harvest privacy info
LoadPrivacy <- function( sp, sc, fn ) {
  # Path to privacy data
  privPath <- file.path( "Data", "Privacy" )
  # Privacy name 
  privName <- paste( "Priv", fn, sep="" )
  # Generate catch/harvest privacy file name
  privFN <- paste( fn, sp, sc, ".csv", sep="" )
  # If there is a privacy file
  if( privFN %in% list.files(privPath) ) {
    # Load the privacy data
    privDat <- read_csv( file=file.path(privPath, privFN), col_types=cols() ) %>%
      mutate( Private=TRUE ) %>%
      rename( !!privName:=Private )  # Weird but works..
    # Print a message
    cat( "Loading", fn, "privacy info for", sc, "by", sp, "\n" )
  } else {  # End if there is privacy data, otherwise
    # No data
    privDat <- NULL
  }  # End if there is no catch privacy data
  # Return the data
  return( privDat )
}  # End LoadPrivacy function

# Load catch privacy data (if any)
privCatchDat <- LoadPrivacy( sp=spUnitName, sc=region, fn="Catch" )

# Load harvest privacy data (if any)
privHarvDat <- LoadPrivacy( sp=spUnitName, sc=region, fn="Harvest" )

# Deal with privacy issues for catch data (if any)
if( !is.null(privCatchDat) ) {
  # Identify which catch values are private (TRUE)
  catch <- catch %>%
    full_join( y=privCatchDat ) %>%
    replace_na( replace=list(PrivCatch=FALSE) )
} else {  # End if privacy issues, otherwise
  # Add a dummy column
  catch <- catch %>%
    mutate( PrivCatch=FALSE )
}  # End if no privacy issues (catch)

# Deal with privacy issues for harvest data (if any)
if( !is.null(privHarvDat) ) {
  # Identify which harvest values are private (TRUE)
  harvest <- harvest %>%
    full_join( y=privHarvDat ) %>%
    replace_na( replace=list(PrivHarvest=FALSE) )
} else {  # End if privacy issues, otherwise
  # Add a dummy column
  harvest <- harvest %>%
    mutate( PrivHarvest=FALSE )
}  # End if no privacy issues (harvest)

##### Spatial ##### 

# Update shapefiles: allowed spatial units
if( spUnitName %in% c("Region", "StatArea", "Group", "Section") ) {
  # If spatial unit is region
  if( spUnitName == "Region" ) {
    # Make new data frames
    shapes$spUnitDF <- shapes$regDF %>% rename_( SpUnit=spUnitName ) 
    shapes$spUnitCentDF <- shapes$regCentDF %>% 
      rename_( SpUnit=spUnitName ) %>%
      filter( SpUnit == region )
  }  # End if statistical areas
  # If spatial unit is statistical areas
  if( spUnitName == "StatArea" ) {
    # Make new data frames
    shapes$spUnitDF <- shapes$saDF %>% rename_( SpUnit=spUnitName ) 
    shapes$spUnitCentDF <- shapes$saCentDF %>% rename_( SpUnit=spUnitName )
  }  # End if statistical areas
  # If spatial unit is groups
  if( spUnitName == "Group" ) {
    # Make new data frames
    shapes$spUnitDF <- shapes$grpDF %>% rename_( SpUnit=spUnitName ) 
    shapes$spUnitCentDF <- shapes$grpCentDF %>% rename_( SpUnit=spUnitName )
  }  # End if groups
  # If spatial unit is sections
  if( spUnitName == "Section" ) {
    # Make new data frames
    shapes$spUnitDF <- shapes$secDF %>% rename_( SpUnit=spUnitName ) 
    shapes$spUnitCentDF <- shapes$secCentDF %>% rename_( SpUnit=spUnitName )
  }  # End if groups
} else {  # End if spatial units are defined, otherwise
  # Error
  stop( "Variable spUnitName = ", spUnitName, " not defined", call.=FALSE )
}  # End if spatial units are not defined

# Rename variables to set spatial resolution
areas <- areas %>%
  mutate( Section=formatC(Section, width=3, flag="0"),
          StatArea=formatC(StatArea, width=2, flag="0") ) %>%
  rename_( SpUnit=spUnitName )
bio <- bio %>%
  mutate( Section=formatC(Section, width=3, flag="0"),
          StatArea=formatC(StatArea, width=2, flag="0") ) %>%
  rename_( SpUnit=spUnitName )
catch <- catch %>%
  mutate( Section=formatC(Section, width=3, flag="0"),
          StatArea=formatC(StatArea, width=2, flag="0") ) %>%
  rename_( SpUnit=spUnitName )
harvest <- harvest %>%
  mutate( Section=formatC(Section, width=3, flag="0"),
          StatArea=formatC(StatArea, width=2, flag="0") ) %>%
  rename_( SpUnit=spUnitName )
siAll <- siAll %>%
  mutate( Section=formatC(Section, width=3, flag="0"),
          StatArea=formatC(StatArea, width=2, flag="0") ) %>%
  rename_(SpUnit=spUnitName )

##### Main #####

##### Tables #####

##### Output #####

##### End ##### 

# Print end of file message and elapsed time
cat( "\nEnd of file Rebuild.R: ", sep="" ) ;  print( Sys.time( ) - sTime )