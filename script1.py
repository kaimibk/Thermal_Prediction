#!/usr/bin/env python 3
import os
import random
import numpy as np
from get_data import get_data

###############################################
# Script 1
###############################################

print("Script1 successully started!")

###############################################
# Variables
###############################################

# The wavelengths in micrometers we are
# measuring at.
# There are two at 3.959 to reflect the
# fact that they have different dynamic
# ranges
# (not like I do anything with it though)
wavelengths = [3.959,3.959,1.64,11.03,12.02]
radiances = [[],[],[],[],[]]

# A list to store the times
UNIXtimes = []

# A list to store the NTI
NTIs = []

# Constants for the Spectral Radiance
# formula; c1 has units W * micrometer^4
# / m^2 whereas c2 has units micrometer
# * Kelvin
c1 = 3.74151e8
c2 = 1.43879e4

# Absolute zero in Celsius
absolute_zero = -273.15

# Current working directory
cwd = os.getcwd() + "/"

###############################################
# Data Collection
###############################################

# Gather time and NTI data with pandas
df = get_data(2019,10,10000,-156,-154,19,20).sort_values('UNIX_Time').reset_index()
df.rename(index=str, columns={"Mo": "Month", "Dy": "Day", "Hr": "Hour", "Mn" : "Minute"}, inplace=True)
print(f"...{df.shape[0]} measurements retrieved...")

# Convert these to easy-to-use lists
radiances[0] = list(df["B21"])
radiances[1] = list(df["B22"])
radiances[2] = list(df["B6"])
radiances[3] = list(df["B31"])
radiances[4] = list(df["B32"])
UNIXtimes = list(df["UNIX_Time"])
NTIs = list(df["Ratio"])
excesses = list(df["Excess"])
temps = list(df["Temp"])
errors = list(df["Err"])

# N is how many data entries we have
N = len(UNIXtimes)

###############################################
# Models
###############################################

# The temperature D (Celsius) of a lava flow
# over time t (hours) be modeled by the
# following equation:
#
# D = a log(t) + b
#
# where a and b are empirical constants
# For surface termperatures:
#     a = -140
#     b = 303
#
# This seems to be accurate for any time
# after 0.05 hours
#
# Hon 1994

# TimeCooled - hours
# SurfaceTemp - Kelvin
def getSurfaceTemp(time_cooled):
    return -140 * np.log10(time_cooled) \
           + 303 - absolute_zero

# The inverse of the above function

# TimeCooled - hours
# SurfaceTemp - Kelvin
def getTimeCooled(surface_temperature):
    return np.exp((surface_temperature \
            + absolute_zero - 303.0)\
            /-140.0)


# Apply Planck's Blackboy Radiation
# Law to get the theoretical
# spectral radiance at some
# temperature and wavelength
#
# Wright 2016

# Temperature - Kelvin
# Wavelength - micrometer
# SpectralRadiance - W / (micrometer * m^2)
def getSpectralRadiance(temperature,wavelength):

    candidate_radiance =\
            c1 / (np.pi * (wavelength**5) \
            * (np.exp(c2 / (wavelength * \
                      temperature)) - 1.0))

#   If the radiance is too large, the value
#   overfills and returns -10.00
#   if (candidate_radiance > 99.99):
#       return -10.00
#   else:
#       return candidate_radiance
    return candidate_radiance

# The inverse of the above function

# Temperature - Kelvin
# Wavelength - micrometer
# SpectralRadiance - W / (micrometer * m^2)
def getTemperature(spectral_radiance,wavelength):

    candidate_temperature =\
            ((wavelength/c2)*np.log(1.0\
            + c1 / (spectral_radiance \
            * np.pi * wavelength**5)))**(-1)

    return candidate_temperature

# Radiant flux can also be approximated
# by spectral radiances from the
# 4 micrometer band from the pixel (ON)
# and a pixel next to it (OFF)
#
# Wright 2016

# RadianceON - W / (micrometer * m^2)
# RadianceOFF - W / (micrometer * m^2)
# RadiantFlux - W / (micrometer * m^2)
def getRadiantFlux(radianceON,radianceOFF):

    return (1.89e7)*\
            (radianceON-radianceOFF)

# Calculate the NTI given the spectral
# radiances at two wavelengths
# (4 micrometer and 12 micrometer)
# as specified in MODVOLC

# Radiance4 - W / (micrometer * m^2)
# Radiance12 - W / (micrometer * m^2)
# NTI - unitless
def getNTI(radiance4,radiance12):

    return (radiance4 - radiance12) / \
           (radiance4 + radiance12)

# The inverse of the above function

# Radiance4 - W / (micrometer * m^2)
# Radiance12 - W / (micrometer * m^2)
# NTI - unitless
def getRadiance4(NTI,radiance12):

    return radiance12*(1.0+NTI)/(1.0-NTI)

# Use geospatial trajectory modeling
# to calculate how often flyovers occur
# for a particular location
#
# This is too difficult at this moment
# so the minimum "twice per day" is used

# Latitude - degrees
# Longitude - degees
# FlyoverPeriod - hours
def getFlyoverPeriod(latitude,longitude):
    return 1.5 * 24.0 / 2

###############################################
# Data Analysis
###############################################



###############################################
# Experiment 1
###############################################
#
# The goal is to confirm how the MODVOLC
# researchers computed the excess radiance
# (flux) from the spectral radiances and
# the average surface temperature
#
###############################################
#
# Hypothesis:
#     Follow exactly the procedure the
#     paper mentioned (Wright 2016)
#
###############################################

Ntrials = 100

# Pick a random number of these measurements
indexesRand = random.sample(range(N),Ntrials)

radiance21sRand = [radiances[0][i] \
        for i in indexesRand]
radiance22sRand = [radiances[1][i] \
        for i in indexesRand]

tempsRand = [temps[i] for i in indexesRand]
errorsRand = [errors[i] for i in indexesRand]

# These are the empirical results
excessesRand = [excesses[i] \
        for i in indexesRand]

# This will be the theoretical results
A = []
for i in range(Ntrials):
    if (radiance22sRand[i] < 0):
        A.append(getRadiantFlux(\
                 radiance21sRand[i],\
                 getSpectralRadiance(\
                 tempsRand[i],\
                 wavelengths[0])))
    else:
        A.append(getRadiantFlux(\
                 radiance22sRand[i],\
                 getSpectralRadiance(\
                 tempsRand[i],\
                 wavelengths[1])))

# Now let's see how they differ
exp001file = open(cwd+"exp001.dat","w")
exp001file.write("#Flux          ")
exp001file.write("# Emp     Theo  ")
for i in range(Ntrials):
    exp001file.write((\
            "{:7.3f} {:7.1f} | "\
            ).format(\
            excessesRand[i],A[i]))
exp001file.close()

###############################################
# Experiment 2
###############################################
#
# The goal is to analyze how the cooling
# rate bounds how quickly the temperature
# of a pixel may decrease and consequently
# how many false negatives may be missing
#
###############################################
#
# Hypothesis:
#     Assume the flyover period is constant;
#     increment the START time for the
#     volcano cooling (because the cooling
#     rate slows over time)
#
#     Use this temperature drop to
#     calculate a respective radiance drop
#     to use as a threshold for identifying
#     false negatives
#
###############################################

# Vary the start times (in hours)
#
# Because of the band saturation, we get
# overfills for temperatures above 500 K
# so I set the minimum start time to be
# the time associated with 500 K
#
# There may be also be some more
# assumptions we can make here but I
# am not making any yet

minSTARTtime = getTimeCooled(500.0)

STARTtimes = minSTARTtime + [\
              0.00,\
              4.00,\
              8.00,\
             16.00\
              ]

Ntrials = len(STARTtimes)

# Recorded when we think data is
# missing
missingData = [[] \
        for i in range(Ntrials)]

# Record when we think data has
# an anomolous drop
anomolousData = [[] \
        for i in range(Ntrials)]

for j in range(Ntrials):

    # The blackbody radiance law is
    # increasing in terms of temperature
    # so we can calculate the maximum
    # and minimum radiance from the
    # maximum and minimum temperatures
    deltaRadiance_max = \
        getSpectralRadiance(\
        getSurfaceTemp(STARTtimes[j]),\
        wavelengths[0])\
      - getSpectralRadiance(\
        getSurfaceTemp(STARTtimes[j]\
          + getFlyoverPeriod(0.0,0.0)),\
        wavelengths[0])

    print("")
    print("deltaRadiance_max:",\
            deltaRadiance_max)
    print("T1:",\
        getSurfaceTemp(STARTtimes[j]))
    print("T2:",\
        getSurfaceTemp(STARTtimes[j]\
          + getFlyoverPeriod(0.0,0.0)))
    print("radiance1:",\
        getSpectralRadiance(\
        getSurfaceTemp(STARTtimes[j]),\
        wavelengths[0]))
    print("radiance2:",\
        getSpectralRadiance(\
        getSurfaceTemp(STARTtimes[j]\
          + getFlyoverPeriod(0.0,0.0)),\
        wavelengths[0]))

    # We aren't really sure about this
    # but, again, we'll use the
    # minimum flyover period of twice
    # per day
    deltaTime_max = \
            getFlyoverPeriod(0.0,0.0)

    # Records how many unique flyovers
    # occured
    Nunique = 0

    # Record the maximum radiance for
    # the 4 micrometer for some
    # UNIXtime
    maxRadianceData = []
    
    # Record how many hotspot pixels
    # were recorded at that time
    NhotspotData = []
    
    # Record the UNIX time for each
    # unique data point
    UNIXtimeData = []

    # Initialize
    maxRadiance = 0
    Nhotspot = 0
    missingData[j].append(1)
    anomolousData[j].append(1)
    maxRadianceData.append(0)
    NhotspotData.append(0)
    UNIXtimeData.append(0)

    for i in range(N-1):

        # Get the radiance (use band
        # 22 preferentially)
        if (radiances[1][i] < 0):
            currentRadiance = \
                radiances[0][i]
        else:
            currentRadiance = \
                radiances[1][i]

        if (currentRadiance\
            > maxRadiance):
            maxRadiance =\
                currentRadiance
            otherRadiance =\
                    radiances[4][i]

        # We look to see if there is
        # an overlap or gap in data
        deltaTime = (UNIXtimes[i+1]\
                    - UNIXtimes[i])\
                    / (60.0*60.0)

        # If there is an overlap,
        # continue reading more
        # data
        if (deltaTime == 0):
            Nhotspot = \
                Nhotspot + 1
            continue

        # Otherwise, calculate the
        # maximum radiance and
        # record this on our time
        # series
        else: 
            maxRadianceData.append(
                maxRadiance)
            NhotspotData.append(
                Nhotspot + 1)
            UNIXtimeData.append(
                UNIXtimes[i])

        # If there is a gap, we know
        # that the NTI is below
        # threshold.
        # Assume the NTI is exactly at
        # the threshold and assume the
        # radiance of band 32 has not
        # changed: calculate the
        # radiance of band 21/22
        if (deltaTime > deltaTime_max):
            nextRadiance = \
                getRadiance4(-0.8,\
                    otherRadiance)

            # If there is missing data
            # we also need to check to
            # see if it is anomolous
            deltaRadiance = \
                maxRadiance-nextRadiance

            if (deltaRadiance > \
                deltaRadiance_max):
                anomolousData[j].append(1)
            else:
                anomolousData[j].append(0)

            missingData[j].append(1)

        # Otherwise, do nothing
        else:
            missingData[j].append(0)

        # This was one unique flyover
        Nunique = Nunique + 1

        # Reset the values
        maxRadiance = 0
        Nhotspot = 0

        # If there was a gap in data
        # previously, no need to check
        # the previous time slot for
        # an anomoly
        if (missingData[j][Nunique-1] == 1):
            continue

        # If there is no gap in data
        # previously, we can now check
        # if there was an anomoly from
        # the previous time slot
        deltaRadiance = \
            maxRadianceData[Nunique-1]\
            - maxRadianceData[Nunique]

        if (deltaRadiance > \
            deltaRadiance_max):
            anomolousData[j].append(1)
        else:
            anomolousData[j].append(0)

###############################################
# Data Post-Processing
###############################################

print("")
print(len(UNIXtimeData),\
      len(maxRadianceData),\
      len(NhotspotData),\
      len(missingData[0]),\
      len(anomolousData[0]),\
      Nunique)

# Output our data onto some file
exp002file = open(cwd+'exp002.dat',"w")

# Limit ourselves to a manageable
# number of data points
Nhandicap = 200
Nunique = min(Nunique,Nhandicap)

longstring = "{:10d}   {:7.3f} {:3d}"
for j in range(Ntrials):
    longstring = longstring + \
            " {:1d} {:1d}"

longarguments = []
for i in range(Nunique):
    longarguments.append([])
    longarguments[i].append(\
            UNIXtimeData[i+1])
    longarguments[i].append(\
            maxRadianceData[i+1])
    longarguments[i].append(\
            NhotspotData[i+1])
    for j in range(Ntrials):
        longarguments[i].append(\
            missingData[j][i+1])
        longarguments[i].append(\
            anomolousData[j][i+1])

for i in range(Nunique):
    exp002file.write((longstring+\
                      "\n").format(\
                          *longarguments[i]))
#                      UNIXtimeData[i+1],\
#                      maxRadianceData[i+1],\
#                      NhotspotData[i+1],\
#                      missingData[0][i+1],\
#                      anomolousData[0][i+1],\
#                      missingData[1][i+1],\
#                      anomolousData[1][i+1],\
#                      missingData[2][i+1],\
#                      anomolousData[2][i+1]))

exp002file.close()

###############################################
# Data Visualization
###############################################

# Make the gnuplot script that will visualize our data
gnuplotfile = open(cwd+"gnuplotfile","w")
def gw(aline):
    gnuplotfile.write(aline)
    return

gw('set term pngcairo size 2400,3600\n')
gw('set output "'+'exp002.png"\n')
gw('set tmargin 0\n')
gw('set bmargin 0\n')
gw('set lmargin 1\n')
gw('set rmargin 1\n')
gw(f'set multiplot layout {Ntrials},1 ' + \
    'columnsfirst margins 0.1,0.95,.1,.9 ' + \
    'spacing 0.1,0 title ' + \
    '"Maximum Radiances Detected with Varying Threshold" ' + \
    'font ",36" offset 0,3\n')

gw('unset xlabel\n')
gw('unset xtics\n')
gw('set grid x\n')
gw('set ylabel "Maximum Radiance (W/micrometer*m^2)"' + \
        'font ",24"\n')
gw('set yrange [-1:]\n')

for i in range(Ntrials):
    if (i == Ntrials - 1):
        gw('set xlabel "UNIX Time (s)" ' + \
                'font ",24"\n')
        gw('set xtics\n')
    gw('plot "'+cwd+'exp002.dat" u 1:2 w l t "",\\\n')
    gw('     "'+cwd+'exp002.dat" u ' + \
            f'1:(${2*i+4}==1?$2:1/0) ' + \
            'w p lc "blue" pt 7 ps 2 t "Missing",\\\n')
    gw('     "'+cwd+'exp002.dat" ' + \
            f'u 1:(${2*i+5}==1?$2:1/0) ' + \
            'w p lc "red" pt 7 ps 2 t "Anomolous",\\\n')
    gw('     "'+cwd+'exp002.dat" ' + \
            f'u 1:(${2*i+4}==1&&${2*i+5}==1?$2:1/0) ' + \
            'w p lc "green" pt 7 ps 2 t "Both"\n')

gnuplotfile.close()

os.system(("gnuplot < "+ cwd + "gnuplotfile"))

print("Script 1 successully exited!")
