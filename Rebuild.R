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
rm(list = ls()) # Clear the workspace
sTime <- Sys.time() # Start the timer
graphics.off() # Turn graphics off

# Install missing packages and load required packages (if required)
UsePackages <- function(pkgs, locn = "https://cran.rstudio.com/") {
  # Reverse the list
  rPkgs <- rev(pkgs)
  # Identify missing (i.e., not yet installed) packages
  newPkgs <- rPkgs[!(rPkgs %in% installed.packages()[, "Package"])]
  # Install missing packages if required
  if (length(newPkgs)) install.packages(newPkgs, repos = locn)
  # Loop over all packages
  for (i in 1:length(rPkgs)) {
    # Load required packages using 'library'
    eval(parse(text = paste("suppressPackageStartupMessages(library(", rPkgs[i],
                            "))",
                            sep = ""
    )))
  } # End i loop over package names
} # End UsePackages function

# Make packages available
UsePackages(pkgs = c(
  "tidyverse", "sp", "scales", "ggforce", "lubridate",
  "cowplot", "GGally", "magick", "ggrepel", "readxl", "xtable",
  "viridis", "zoo", "SpawnIndex"
))

##### Controls #####

# Select region: major (HG, PRD, CC, SoG, WCVI); or minor (A27, A2W)
region <- c("All")

# Spatial unit: Region, StatArea, Section, or Group
spUnitName <- "Region"

##### Parameters #####

# Apply privacy restrictions (i.e., rule of 3)
doPrivacy <- TRUE

# Spawn index threshold (tonnes; NA for none)
siThreshold <- NA # 15000

# Minimum number of consecutive years
nYrsConsec <- 3

# Buffer distance (m; to include locations that are outside the region polygon)
maxBuff <- 10000

# Intended harvest rate
intendU <- 0.2

# First year of intended harvest rate
intendUYrs <- 1983

# Plot quality (dots per inch)
pDPI <- 320

##### Files #####

# File name for dive transect XY
diveFN <- file.path("Data", "dive_transects_with_lat_long_June2_2017.xlsx")

# File name for q parameters
qFN <- file.path("Data", "qPars.csv")

# File name for reference years
refYrsAll <- refFN <- file.path("Data", "RefYrs.csv")

##### Functions #####

# Load helper functions
source(file = file.path("..", "HerringFunctions", "Functions.R"))

##### Data #####

# Load required data

# Message
cat("Investigate", region, "by", spUnitName, "\n")

# Load saved data from the stock assessment model
LoadAssessment <- function(loc) {
  # Load the saved object
  load(file = file.path(
    "..", "herringsr", "models", loc, "AM2",
    paste("aaa_gfiscam", ".RData", sep = "")
  ))
  # Grab quantiles for q
  qPars <- t(model$mcmccalcs$q.quants) %>%
    as_tibble() %>%
    rename(qLower = `5%`, qMedian = `50%`, qUpper = `95%`) %>%
    mutate(Survey = c("Surface", "Dive"))
  # Make a list of things to return
  res <- list(qPars = qPars)
  # Return the list
  return(res)
} # End LoadAssessment function

# If there is an assessment
if(region != "All") {
  # Load assessment values (directly!)
  assessOutput <- LoadAssessment(loc = region)
  # Load q from the assessment
  qPars <- assessOutput$qPars
} else { # End if loading assesment output, otherwise
  # Fill q = 1
  qPars <- tibble(
    qLower=c(1, 1), qMedian=c(1, 1), qUpper=c(1, 1),
    Survey=c("Surface", "Dive")
  )
}

# Load reference years
refYrs <- read_csv(file = refFN, col_types = cols("c", "i", "i")) %>%
  filter(Region == region)

# Error if no reference years present
if (nrow(refYrs) == 0) {
  stop("Specify reference years for biomass threshold (refYrs): ", region,
       call. = FALSE
  )
}

# If old directory exists
if (region %in% list.files()) {
  # Remove the old directory
  unlink(region, recursive = TRUE)
  # Create the main directory for output
  dir.create(region)
} else { # End if directory exists, otherwise
  # Create the main directory for output
  dir.create(region)
} # End if directory doesn't exists

# Get data from the data summary
LoadSavedObjects <- function(loc) {
  # Load the saved object
  load(file = file.path(
    "..", "DataSummaries", loc,
    paste("Image.", loc, ".RData", sep = "")
  ))
  # load( file=file.path("Data", paste("Image.", loc, ".RData", sep="")) )
  # Confirm the region (from the saved data -- should match!)
  if (regName != region) {
    stop("Region mismatch -- check saved image", call. = FALSE)
  }
  # Return objects to the main environment
  newSurvYr <<- newSurvYr # Year in which survey type changed
  yrRange <<- yrRange # Range of years to consider
  inCRS <<- inCRS # Coordinate reference system (input)
  outCRS <<- outCRS # Coordinate reference system (output)
  shapes <<- shapes # Shapefiles etc
  figWidth <<- figWidth * 1.2 # Figure size
  catch <<- catch # Catch data
  UpdateCatchData <<- UpdateCatchData # Function to update catch data
  CalcWeightAtAge <<- CalcWeightAtAge # Function to calculate weight-at-age
  CalcLengthAtAge <<- CalcLengthAtAge # Function to calculate length-at-age
  inclTestCatch <<- inclTestCatch # Include catch from test fishery
  geoProj <<- geoProj # Text for maps: projection
  areas <<- areas # Spatial info
  bio <<- bio # Biological data
  spawnRaw <<- spawnRaw # Spawn data
  catchRaw <<- catchRaw # Raw catch data (for SOK harvest)
  ciLevel <<- ciLevel # Confidence interval
  ageRange <<- ageRange # Ages
  qYrs <<- qYrs # Survey years (q1 and q2)
  smLine <<- smLine # Smoothing function for line
  nRoll <<- nRoll # Number of years for rolling mean
  ageShow <<- ageShow # Age to highlight on the x-at-age plots
  yrBreaks <<- yrBreaks # Years to show in x-axes
  convFac <<- convFac # Unit conversion factors
  ECF <<- ECF # Egg conversion factor
} # End LoadSavedObjects function

# Load saved data (directly!)
LoadSavedObjects(loc = region)

# Get small area table
# TODO: Are both of these required? aSmall and areasSM
aSmall <- areas %>%
  select(Region, RegionName, StatArea, Group, Section) %>%
  distinct()

# Calculate SOK harvest
harvest <- catchRaw %>%
  filter(DisposalCode == 2, Source == "SOK") %>%
  mutate(
    BiomassSOK = CalcBiomassSOK(SOK = Catch * convFac$lb2kg),
    HarvSOK = Catch * convFac$lb2kg / 1000
  ) %>%
  left_join(y = aSmall, by = c("Region", "StatArea", "Section")) %>%
  select(Year, Region, StatArea, Group, Section, BiomassSOK, HarvSOK)

# Function to load transect spatial info
LoadTransectXY <- function(loc) {
  # Load the data and wrangle
  dat <- read_excel(path = loc, sheet = 1) %>%
    rename(LocationCode = LOC_CODE) %>%
    mutate(LocationCode = as.integer(LocationCode)) %>%
    group_by(LocationCode) %>%
    summarise(
      Longitude = MeanNA(c(StartLong, MidLong, EndLong)),
      Latitude = MeanNA(c(StartLat, MidLat, EndLat))
    ) %>%
    ungroup() %>%
    filter(!is.na(Longitude), !is.na(Latitude), LocationCode != 0)
  # Grab the spatial info (X and Y)
  locSP <- dat %>%
    transmute(X = Longitude, Y = Latitude)
  # Put X and Y into a spatial points object
  locPts <- SpatialPointsDataFrame(
    coords = locSP,
    data = data.frame(LocationCode = dat$LocationCode), proj4string = CRS(inCRS)
  )
  # Convert X and Y from WGS to Albers
  locPtsAlb <- spTransform(x = locPts, CRSobj = CRS(outCRS))
  # Extract spatial info
  dfAlb <- as_tibble(locPtsAlb) %>%
    rename(Eastings = X, Northings = Y)
  # Return the data
  return(dfAlb)
} # End LoadTransectXY function

# Load 'auxiliary' dive transect spatial data
transectXY <- LoadTransectXY(loc = diveFN)

# Get a short list of areas: sections and groups
areasSm <- areas %>%
  select(Region, StatArea, Group, Section) %>%
  distinct()

# Get spawn index data
GetSI <- function(allSI, loc, XY) {
  # Get a subset and wrangle
  raw <- allSI %>%
    replace_na(list(Group = "Other")) %>%
    mutate(
      Duration = as.numeric(End - Start + 1),
      Decade = GetDecade(Year), Area = Length * Width, Week = week(Start),
      Survey = ifelse(Year < newSurvYr, "Surface", "Dive"),
      YrsSurv = ifelse(Year < newSurvYr, length(yrRange[yrRange < newSurvYr]),
                       length(yrRange[yrRange >= newSurvYr])
      )
    ) %>%
    rowwise() %>%
    mutate(
      LyrsMean = MeanNA(c(SurfLyrs, MacroLyrs, UnderLyrs)),
      SITotal = SumNA(c(SurfSI, MacroSI, UnderSI))
    ) %>%
    ungroup() %>%
    group_by(Decade) %>%
    mutate(YrsDecade = length(unique(Year))) %>%
    ungroup()
  # Clip the extent
  df <- ClipExtent(
    dat = raw, spObj = shapes$regSPDF, bufDist = maxBuff,
    silent = TRUE
  )
  # Subset data with 'good' X and Y
  dfNotNA <- df %>%
    filter(!is.na(Eastings) & !is.na(Northings))
  # Subset data with 'bad' X or Y, and try to fill in using transect X and Y
  dfNA <- df %>%
    filter(is.na(Eastings) | is.na(Northings)) %>%
    select(-Eastings, -Northings) %>%
    left_join(y = XY, by = "LocationCode")
  # Re-combine the two subsets
  df2 <- bind_rows(dfNotNA, dfNA)
  # Clip the extent (again)
  res <- ClipExtent(dat = df2, spObj = shapes$regSPDF, bufDist = maxBuff)
  # Stop if we're missing rows
  if (nrow(raw) != nrow(res)) stop("Missing rows!", call. = FALSE)
  # Return the spatial data
  return(res)
} # End GetSI function

# Get spawn index
siAll <- GetSI(allSI = spawnRaw, loc = region, XY = transectXY)

# Check for weird durations
if (any(siAll$Duration > 20) | any(siAll$Duration < 0)) {
  # Count how many
  oddDuration <- c(which(siAll$Duration > 20), which(siAll$Duration < 0))
  # Warning
  cat("Set", length(oddDuration), "duration(s) to NA\n")
  # Set to NA
  siAll <- siAll %>%
    mutate(Duration = ifelse(Duration %in% 0:20, Duration, NA))
} # End checking duration

##### Privacy #####

# Load catch and harvest privacy info
LoadPrivacy <- function(sp, sc, fn) {
  # Path to privacy data
  privPath <- file.path("Data", "Privacy")
  # Privacy name
  privName <- paste("Priv", fn, sep = "")
  # Generate catch/harvest privacy file name
  privFN <- paste(fn, sp, sc, ".csv", sep = "")
  # If there is a privacy file
  if (privFN %in% list.files(privPath)) {
    # Load the privacy data
    privDat <- read_csv(file = file.path(privPath, privFN), col_types = cols()) %>%
      mutate(Private = TRUE) %>%
      rename(!!privName := Private) # Weird but works..
    # Print a message
    cat("Loading", fn, "privacy data\n")
  } else { # End if there is privacy data, otherwise
    # No data
    privDat <- NULL
    # Print a message
    cat("No", fn, "privacy data found\n")
  } # End if there is no catch privacy data
  # Return the data
  return(privDat)
} # End LoadPrivacy function

# Load catch privacy data (if any)
privCatchDat <- LoadPrivacy(sp = spUnitName, sc = region, fn = "Catch")

# Load harvest privacy data (if any)
privHarvDat <- LoadPrivacy(sp = spUnitName, sc = region, fn = "Harvest")

# Names for privacy data (in join below)
if (spUnitName == "Region") {
  namesPriv <- "Year"
} else {
  namesPriv <- c("Year", spUnitName)
}

# Deal with privacy issues for catch data (if any)
if (!is.null(privCatchDat) & doPrivacy) {
  # Identify which catch values are private (TRUE)
  catch <- catch %>%
    full_join(y = privCatchDat, by = namesPriv) %>%
    replace_na(replace = list(PrivCatch = FALSE))
} else { # End if privacy issues, otherwise
  # Add a dummy column
  catch <- catch %>%
    mutate(PrivCatch = FALSE)
} # End if no privacy issues (catch)

# Deal with privacy issues for harvest data (if any)
if (!is.null(privHarvDat) & doPrivacy) {
  # Identify which harvest values are private (TRUE)
  harvest <- harvest %>%
    full_join(y = privHarvDat, by = namesPriv) %>%
    replace_na(replace = list(PrivHarvest = FALSE))
} else { # End if privacy issues, otherwise
  # Add a dummy column
  harvest <- harvest %>%
    mutate(PrivHarvest = FALSE)
} # End if no privacy issues (harvest)

##### Spatial #####

# Update shapefiles: allowed spatial units
if (spUnitName %in% c("Region", "StatArea", "Group", "Section")) {
  # If spatial unit is region
  if (spUnitName == "Region") {
    # Make new data frames
    shapes$spUnitDF <- shapes$regDF %>% rename(SpUnit = spUnitName)
    shapes$spUnitCentDF <- shapes$regCentDF %>%
      rename(SpUnit = spUnitName) %>%
      filter(SpUnit == region)
  } # End if statistical areas
  # If spatial unit is statistical areas
  if (spUnitName == "StatArea") {
    # Make new data frames
    shapes$spUnitDF <- shapes$saDF %>% rename(SpUnit = spUnitName)
    shapes$spUnitCentDF <- shapes$saCentDF %>% rename(SpUnit = spUnitName)
  } # End if statistical areas
  # If spatial unit is groups
  if (spUnitName == "Group") {
    # Make new data frames
    shapes$spUnitDF <- shapes$grpDF %>% rename(SpUnit = spUnitName)
    shapes$spUnitCentDF <- shapes$grpCentDF %>% rename(SpUnit = spUnitName)
  } # End if groups
  # If spatial unit is sections
  if (spUnitName == "Section") {
    # Make new data frames
    shapes$spUnitDF <- shapes$secDF %>% rename(SpUnit = spUnitName)
    shapes$spUnitCentDF <- shapes$secCentDF %>% rename(SpUnit = spUnitName)
  } # End if groups
} else { # End if spatial units are defined, otherwise
  # Error
  stop("Variable spUnitName = ", spUnitName, " not defined", call. = FALSE)
} # End if spatial units are not defined

# Rename variables to set spatial resolution
areas <- areas %>%
  mutate(
    Section = formatC(Section, width = 3, flag = "0"),
    StatArea = formatC(StatArea, width = 2, flag = "0")
  ) %>%
  rename(SpUnit = spUnitName)
bio <- bio %>%
  mutate(
    Section = formatC(Section, width = 3, flag = "0"),
    StatArea = formatC(StatArea, width = 2, flag = "0")
  ) %>%
  rename(SpUnit = spUnitName)
catch <- catch %>%
  mutate(
    Section = formatC(Section, width = 3, flag = "0"),
    StatArea = formatC(StatArea, width = 2, flag = "0")
  ) %>%
  rename(SpUnit = spUnitName)
harvest <- harvest %>%
  mutate(
    Section = formatC(Section, width = 3, flag = "0"),
    StatArea = formatC(StatArea, width = 2, flag = "0")
  ) %>%
  rename(SpUnit = spUnitName)
siAll <- siAll %>%
  mutate(
    Section = formatC(Section, width = 3, flag = "0"),
    StatArea = formatC(StatArea, width = 2, flag = "0")
  ) %>%
  rename(SpUnit = spUnitName)

##### Main #####

# Aggregate catch by year and spatial unit
catchYrSp <- catch %>%
  group_by(Year, SpUnit) %>%
  summarise(Catch = SumNA(Catch), PrivCatch = unique(PrivCatch)) %>%
  ungroup() %>%
  arrange(Year, SpUnit) %>%
  mutate(CatchShow = ifelse(PrivCatch, 0, Catch))

# Aggregate harvest by year and spatial unit
harvYrSp <- harvest %>%
  group_by(Year, SpUnit) %>%
  summarise(HarvSOK = SumNA(HarvSOK), PrivHarvest = unique(PrivHarvest)) %>%
  ungroup() %>%
  arrange(Year, SpUnit) %>%
  mutate(HarvSOKShow = ifelse(PrivHarvest, 0, HarvSOK))

# Add harvest data to catch table
chYrSp <- full_join(catchYrSp, harvYrSp, by = c("Year", "SpUnit"))

# Update years (use the same year to compare timing among years)
year(siAll$Start) <- 0000
year(siAll$End) <- 0000

# Make it long (for plots)
siAllLong <- siAll %>%
  gather("Start", "End", key = "Timing", value = "Date") %>%
  mutate(
    Survey = factor(Survey, levels = c("Surface", "Dive")),
    Timing = factor(Timing, levels = c("Start", "End"))
  )

# # Spawn duration
# siDuration <- siAll %>%
#   # Add 1 so that that spawns that start and end same day are 1 day long
#   mutate(Duration = End - Start + 1) %>%
#   group_by(Year, Survey, SpUnit) %>%
#   summarise(
#     DurationMean = MeanNA(Duration),
#     DurationSD = sd(Duration, na.rm = TRUE),
#     DurationNum = n()
#   ) %>%
#   ungroup() %>%
#   mutate(
#     Lower = DurationMean - DurationSD,
#     Upper = DurationMean + DurationSD,
#     Survey = factor(Survey, levels = c("Surface", "Dive"))
#   )

# Aggregate spawn index by year and spatial unit
siYrSp <- siAll %>%
  group_by(Year, SpUnit) %>%
  summarise(
    NumLocs = n_distinct(LocationCode),
    LengthTotal = SumNA(Length),
    WidthMean = MeanNA(Width),
    AreaTotal = SumNA(Area),
    LayersMean = MeanNA(LyrsMean),
    SITotal = SumNA(SITotal),
    DateFirst = MinNA(Start, End),
    DateLast = MaxNA(Start, End),
    DateDiff = DateLast - DateFirst,
    Survey = unique(Survey)
  ) %>%
  group_by(SpUnit) %>%
  mutate(NConsec = CountConsecutive(Year)) %>%
  ungroup() %>%
  arrange(Year, SpUnit) %>%
  left_join(y = qPars, by = "Survey") %>%
  mutate(
    BiomassLower = SITotal / qUpper, BiomassMedian = SITotal / qMedian,
    BiomassUpper = SITotal / qLower,
    Survey = factor(Survey, levels = c("Surface", "Dive"))
  )

# Combine catch with spawn index by year and spatial unit
allYrSp <- full_join(x = chYrSp, y = siYrSp, by = c("Year", "SpUnit")) %>%
  arrange(Year, SpUnit) %>%
  mutate(Survey = ifelse(Year < newSurvYr, "Surface", "Dive")) %>%
  replace_na(replace = list(NConsec = -1)) %>%
  mutate(
    Survey = factor(Survey, levels = c("Surface", "Dive")),
    HarvestLower = Catch / (Catch + BiomassUpper),
    HarvestMedian = Catch / (Catch + BiomassMedian),
    HarvestUpper = Catch / (Catch + BiomassLower)
  ) %>%
  filter(!is.na(SpUnit))

# Determine ratio of max SOK harvest to max spawn index
rSOK <- max(allYrSp$HarvSOK, na.rm = TRUE) / max(allYrSp$SITotal, na.rm = TRUE)

# Count the number of fish aged by year (and as a proportion) by seine gear:
# use the 'SampWt' column to fix unrepresentative sampling if identified
npAgedYear <- bio %>%
  filter(GearCode == 29) %>%
  select(Year, SpUnit, Age, SampWt) %>%
  na.omit() %>%
  group_by(Year, SpUnit, Age) %>%
  summarise(Number = SumNA(SampWt)) %>%
  mutate(Proportion = Number / SumNA(Number)) %>%
  ungroup() %>%
  arrange(Year, SpUnit, Age)

# Calculate weight-at-age by year and area
weightAge <- bio %>%
  group_by(SpUnit) %>%
  do(CalcWeightAtAge(.)) %>%
  ungroup() %>%
  select(Year, SpUnit, Age, Weight) %>%
  arrange(Year, SpUnit, Age)

# Calculate running mean weight-at-age by year (if data exist)
if (exists("weightAge")) {
  muWeightAge <- weightAge %>%
    group_by(SpUnit, Age) %>%
    mutate(muWeight = rollmean(x = Weight, k = nRoll, align = "right", na.pad = TRUE)) %>%
    ungroup() %>%
    mutate(Age = factor(Age))
}

# Calculate length-at-age by year and area
lengthAge <- bio %>%
  group_by(SpUnit) %>%
  do(CalcLengthAtAge(.)) %>%
  ungroup() %>%
  select(Year, SpUnit, Age, Length) %>%
  arrange(Year, SpUnit, Age)

# Calculate running mean length-at-age by year (if data exist)
if (exists("lengthAge")) {
  muLengthAge <- lengthAge %>%
    group_by(SpUnit, Age) %>%
    mutate(muLength = rollmean(x = Length, k = nRoll, align = "right", na.pad = TRUE)) %>%
    ungroup() %>%
    mutate(Age = factor(Age))
}

# Calculate proporiton of total spawn by week and spatial unit
propWeekSI <- siAll %>%
  filter(!is.na(Week)) %>%
  group_by(SpUnit, Week) %>%
  summarise(SITotal = SumNA(SITotal)) %>%
  group_by(SpUnit) %>%
  mutate(SIProp = SITotal / SumNA(SITotal)) %>%
  ungroup()

##### Figures #####

# Spawn index time series
siPlot <- ggplot(
  data = allYrSp,
  aes(x = Year, group = Survey)
) +
  geom_vline(xintercept = newSurvYr - 0.5, linetype = "dashed", size = 0.25) +
  scale_x_continuous(breaks = seq(from = 1000, to = 3000, by = 10)) +
  expand_limits(x = yrRange) +
  myTheme +
  facet_grid(SpUnit ~ ., scales = "free_y") +
  theme(legend.position = "top")

# Spawn index plot: with catch
siPlotCatch <- siPlot +
  geom_line(aes(y = SITotal), na.rm = TRUE) +
  geom_point(aes(y = SITotal, shape = Survey), na.rm = TRUE) +
  labs(y = "Spawn index and catch (t)") +
  scale_y_continuous(labels = comma) +
  geom_col(data = filter(allYrSp, !PrivCatch), aes(y = Catch), alpha = 0.5) +
  geom_point(
    data = filter(allYrSp, PrivCatch), aes(y = CatchShow), shape = 8,
    na.rm = TRUE
  ) +
  ggsave(
    filename = file.path(region, "SpawnIndexCatch.png"),
    height = min(8.75, n_distinct(allYrSp$SpUnit) * 1.9 + 1),
    width = figWidth
  )

# Spawn index plot with SOK harvest (harvest is in lbs -- need to scale)
siPlotHarv <- siPlot +
  geom_line(aes(y = SITotal), na.rm = TRUE) +
  geom_point(aes(y = SITotal, shape = Survey), na.rm = TRUE) +
  labs(y = "Spawn index (t)") +
  scale_y_continuous(
    labels = comma,
    sec.axis = sec_axis(~ . * rSOK, labels = comma, name = "SOK harvest (t)")
  ) +
  geom_col(data = filter(allYrSp, !PrivHarvest), aes(y = HarvSOK / rSOK), alpha = 0.5) +
  geom_point(
    data = filter(allYrSp, PrivHarvest), aes(y = HarvSOKShow), shape = 8,
    na.rm = TRUE
  ) +
  ggsave(
    filename = file.path(region, "SpawnIndexHarv.png"),
    height = min(8.75, n_distinct(allYrSp$SpUnit) * 1.9 + 1),
    width = figWidth
  )

# Spawn timing by year and spatial unit
timingPlot <- ggplot(data = filter(siAllLong, !is.na(Survey)), aes(x = Year)) +
  geom_point(aes(y = Date, shape = Survey),
             alpha = 0.5,
             na.rm = TRUE
  ) +
  geom_vline(xintercept = newSurvYr - 0.5, linetype = "dashed", size = 0.25) +
  scale_x_continuous(breaks = seq(from = 1000, to = 3000, by = 10)) +
  expand_limits(x = yrRange) +
  labs(y = "Date") +
  facet_grid(SpUnit ~ Timing) +
  myTheme +
  theme(legend.position = "top") +
  ggsave(
    filename = file.path(region, "SpawnTiming.png"),
    height = min(8.75, n_distinct(siAll$SpUnit) * 1.9 + 1), width = figWidth
  )

# Spawn duration by year and spatial unit
durationPlot <- ggplot(data = siAll, mapping = aes(y = Duration, x = Year)) +
  geom_point(mapping = aes(shape = Survey), alpha = 0.5, na.rm = TRUE) +
  # geom_boxplot(mapping = aes(group=Year), na.rm = TRUE) +
  geom_vline(xintercept = newSurvYr - 0.5, linetype = "dashed", size = 0.25) +
  scale_x_continuous(breaks = seq(from = 1000, to = 3000, by = 10)) +
  expand_limits(x = yrRange, y = 0) +
  labs(y = "Duration (days)") +
  facet_grid(SpUnit ~ .) +
  myTheme +
  theme(legend.position = "top") +
  ggsave(
    filename = file.path(region, "SpawnDuration.png"),
    height = min(8.75, n_distinct(siAll$SpUnit) * 1.9 + 1), width = figWidth
  )

# Number of spawns by year and spatial unit
numSpawnPlot <- ggplot(data = siAll, mapping = aes(x = Year)) +
  geom_bar() +
  geom_vline(xintercept = newSurvYr - 0.5, linetype = "dashed", size = 0.25) +
  scale_x_continuous(breaks = seq(from = 1000, to = 3000, by = 10)) +
  expand_limits(x = yrRange, y = 0) +
  labs(y = "Number of spawns") +
  facet_grid(SpUnit ~ .) +
  myTheme +
  theme(legend.position = "top") +
  ggsave(
    filename = file.path(region, "SpawnNumber.png"),
    height = min(8.75, n_distinct(siAll$SpUnit) * 1.9 + 1), width = figWidth
  )

# Spawn index for spawn number by year and spatial unit
siNumPlot <- ggplot(data = siAll, mapping = aes(y = SITotal, x = Year)) +
  geom_point(mapping = aes(shape = Survey), alpha = 0.5, na.rm = TRUE) +
  geom_vline(xintercept = newSurvYr - 0.5, linetype = "dashed", size = 0.25) +
  scale_x_continuous(breaks = seq(from = 1000, to = 3000, by = 10)) +
  scale_y_continuous(labels = comma) +
  # scale_y_log10(labels = comma) +
  expand_limits(x = yrRange, y = 0) +
  labs(y = "Spawn index (t)") +
  facet_grid(SpUnit ~ .) +
  myTheme +
  theme(legend.position = "top") +
  ggsave(
    filename = file.path(region, "SpawnIndexNum.png"),
    height = min(8.75, n_distinct(siAll$SpUnit) * 1.9 + 1), width = figWidth
  )

# Plot weight-at-age by year (if data exist)
if (exists("muWeightAge")) {
  weightAgePlot <- ggplot(data = muWeightAge) +
    geom_line(aes(x = Year, y = muWeight, group = Age, colour = Age),
              size = 1, na.rm = TRUE
    ) +
    scale_colour_viridis(guide = guide_legend(nrow = 1), discrete = TRUE) +
    labs(y = "Weight-at-age (g)") +
    scale_x_continuous(breaks = yrBreaks) +
    #    coord_cartesian( ylim=wtRange ) +
    expand_limits(x = yrRange) +
    facet_grid(SpUnit ~ .) +
    myTheme +
    theme(legend.position = "top") +
    ggsave(
      filename = file.path(region, "WeightAge.png"), width = figWidth,
      height = min(9, n_distinct(npAgedYear$SpUnit) * 2 + 1)
    )
}

# Plot length-at-age by year
if (exists("muLengthAge")) {
  lengthAgePlot <- ggplot(data = muLengthAge) +
    geom_line(aes(x = Year, y = muLength, group = Age, colour = Age),
              size = 1, na.rm = TRUE
    ) +
    scale_colour_viridis(guide = guide_legend(nrow = 1), discrete = TRUE) +
    labs(y = "Length-at-age (mm)") +
    scale_x_continuous(breaks = yrBreaks) +
    #    coord_cartesian( ylim=lenRange ) +
    expand_limits(x = yrRange) +
    facet_grid(SpUnit ~ .) +
    myTheme +
    theme(legend.position = "top") +
    ggsave(
      filename = file.path(region, "LengthAge.png"), width = figWidth,
      height = min(9, n_distinct(npAgedYear$SpUnit) * 2 + 1)
    )
}

# Plot proportion of spawn by week
propWeekSIPlot <- ggplot(
  data = propWeekSI, mapping = aes(x = Week, y = SIProp, group = Week)
) +
  geom_bar(stat = "identity") +
  labs(x = "Week of the year", y = "Proportion of spawn index") +
  scale_x_continuous(breaks = seq(from = 0, to = 50, by = 5)) +
  facet_grid(SpUnit ~ .) +
  myTheme +
  ggsave(
    filename = file.path(region, "PropWeekSI.png"), width = figWidth,
    height = min(9, n_distinct(npAgedYear$SpUnit) * 2 + 1)
  )

##### Tables #####

##### Output #####

##### End #####

# Print end of file message and elapsed time
cat("\nEnd of file Rebuild.R: ", sep = "")
print(Sys.time() - sTime)
