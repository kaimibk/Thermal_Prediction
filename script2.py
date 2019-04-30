#!/usr/bin/env python 3
import os
import random
import time
import numpy as np
from get_data import get_data

###############################################
# Script 2
###############################################

print("Script2 successully started!")

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

# The center of whatever area we are observing
center_long = -155.0
center_lat = 19.5

# The range of this area of observation
dcenter_long = 1.0
dcenter_lat = 0.5

# Functions to return the length in latitude
# or longitude to kilometers, and vice versa
def lat_to_km(latitude,longitude,dlatitude):
    return dlatitude * 110.574

def long_to_km(latitude,longitude,dlongitude):
    return dlongitude * 111.320 * \
            np.cos(latitude * np.pi / 180.0)

def km_to_lat(latitude,longitude,dkm):
    return dkm / 110.574

def km_to_long(latitude,longitude,dkm):
    return dkm / (111.320 * \
            np.cos(latitude * np.pi / 180.0))

# The resolution of our heatmap
resolution = 0.25
dlat = km_to_lat(center_lat,center_long,\
        resolution)
dlong = km_to_long(center_lat,center_long,\
        resolution)
Nlat = int(dcenter_lat * 2.0 / dlat)
Nlong = int(dcenter_long * 2.0 / dlong)

# Gather time and NTI data with pandas
df = get_data(\
        2019,10,10000,\
        center_long-dcenter_long,center_long+dcenter_long,\
        center_lat-dcenter_lat,center_lat+dcenter_lat\
        ).sort_values('UNIX_Time').reset_index()
df.rename(index=str, columns={"Mo": "Month", "Dy": "Day", "Hr": "Hour", "Mn" : "Minute"}, inplace=True)
print(f"...{df.shape[0]} measurements retrieved...")

# Convert these to easy-to-use lists
radiances[0] = list(df["B21"])
radiances[1] = list(df["B22"])
radiances[2] = list(df["B6"])
radiances[3] = list(df["B31"])
radiances[4] = list(df["B32"])
latitudes = list(df["Latitude"])
longitudes = list(df["Longitude"])
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
    return -140.0 * np.log10(time_cooled) \
           + 303.0 - absolute_zero

# The inverse of the above function

# TimeCooled - hours
# FractionCovered - between 0 and 1
# BackgroundTemp - Kelvin
# SurfaceTemp - Kelvin
#def getTimeCooled(surface_temperature):
#    return 10.0**((surface_temperature \
#            + absolute_zero - 303.0)\
#            /-140.0)
def getTimeCooled(surface_temperature,\
        background_temp,fraction_covered):
#   lava_temperature = \
#           (surface_temperature - \
#           background_temp) / \
#           fraction_covered + 1.0
#   return 10.0**((lava_temperature \
    return 10.0**((surface_temperature \
            - 303.0) /-140.0)


# Apply Planck's Blackboy Radiation
# Law to get the theoretical
# spectral radiance at some
# temperature and wavelength
#
# Wright 2016

# Temperature - Kelvin
# Wavelength - micrometer
# SpectralRadiance - W / (micrometer * m^2)
def getSpectralRadiance1(temperature,\
        wavelength):

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

# Same as above but account for the
# fact that some of the surface area
# is one temperature, and some is
# another, so the total radiance
# (which is integrated over the area)
# is the weighted sum of the radiance
# of the radiance at the two
# temperatures

# Temperature1 - Kelvin (lava)
# Temperature2 - Kelvin (background)
# FractionCovered - between 0 and 1
# Wavelength - micrometer
# SpectralRadiance - W / (micrometer * m^2)
def getSpectralRadiance2(temperature1,\
        temperature2,fraction_covered,\
        wavelength):

    candidate_radiance =\
            c1 / (np.pi * (wavelength**5)) \
            * (fraction_covered / \
            (np.exp(c2 / (wavelength * \
            temperature1)) - 1.0) + \
            (1.0 - fraction_covered) / \
            (np.exp(c2 / (wavelength * \
            temperature2)) - 1.0))

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
    return 1.1 * 24.0 / 2

###############################################
# Experiment 3
###############################################
#
# The goal is to track the radiances pixel
# by pixel both on a graph and on a visual
# heatmap; we will also make a gif to see
# its evolution over time
#
###############################################
#
# Hypothesis:
#
###############################################

time_per_frame = \
        getFlyoverPeriod(0.0,0.0) / 2
frames_per_second = 4

minSTARTtime = getTimeCooled(500.0,27.0,0.90)
#minSTARTtime = getTimeCooled(1400.0,27.0,0.90)

STARTtime = 0.00

STARTtime = STARTtime + minSTARTtime

gif_start_time = UNIXtimes[0]
gif_end_time = gif_start_time + \
        100 * 60*60*getFlyoverPeriod(0.0,0.0)

os.system(("rm -r "+cwd+"png/"))
os.system(("mkdir "+cwd+"png/"))

def makeImages(local_observable,\
        start_index,total_indexes):

###############################################
# Data Processing
###############################################
    
#   print("")
#   print(len(UNIXtimeData),\
#         len(maxRadianceData),\
#         len(NhotspotData),\
#         len(missingData),\
#         len(anomalousData),\
#         len(suspiciousData),\
#         Nunique)

    print("New frame!")
    print("")

    min_radiance = 0.0
    max_radiance = 15.0
    
    # Output our data onto some file
    exp003file = open(cwd+\
            'exp003heatmap.dat',"w")

    for j in range(Nlat):
        for k in range(Nlong):
            longitude = center_long -\
                    dcenter_long + \
                    dlong * (0.5 + k)
            latitude = center_lat -\
                    dcenter_lat + \
                    dlat * (0.5 + j)
            radiance = \
                    local_observable[j][k]

            radiance = min(radiance,
                    max_radiance)
            radiance = max(radiance,
                    min_radiance)
            
#           min_radiance = min(\
#                   min_radiance,\
#                   radiance)
#           max_radiance = max(\
#                   max_radiance,\
#                   radiance)

            exp003file.write((\
#                   '{:12.6f} ' + \
#                   '{:12.6f} ' + \
                    '{:6d} ' + \
                    '{:6d} ' + \
                    '{:12.6f}\n').format(\
#                   latitude,longitude,\
                    j,k,\
#                   longitude,latitude,\
                    radiance))

#           exp003file.write(\
#                   f'{latitude} '+ \
#                   f'{longitude} '+ \
#                   f'{radiance}\n')
        exp003file.write('\n')
    
    exp003file.close()
 
    # Output our data onto some file
    exp003file = open(cwd+'exp003.dat',"w")

    for j in range(Nunique):
        UNIXtime = UNIXtimeData[j+1]
        radiance = maxRadianceData[j+1]
        exp003file.write(\
                f'{UNIXtime} ' + \
                f'{radiance}\n')
    
    exp003file.close()
    
###############################################
# Data Visualization
###############################################
     
    for n in range(total_indexes):
        nframe = start_index + n + 1

        # Make the gnuplot script that will
        # visualize our data
        gnuplotfile = open(cwd+"gnuplotfile","w")
        def gw(aline):
            gnuplotfile.write(aline)
            return

        gw('set term pngcairo size 2400,1200\n')
        gw(('set output "' + cwd + \
                'png/{:04d}.png"\n').format(nframe))
        gw(f'set multiplot\n')
        gw('set title ' + \
            '"Radiances over Kileaua" ' + \
            'font ",36" offset 0,3\n')
        gw('set size 0.5, 1.0\n')
        gw('set origin 0.0, 0.0\n')

        gw('set lmargin at screen 0.05\n')
        gw('set rmargin at screen 0.45\n')
        gw('set bmargin at screen 0.10\n')
        gw('set tmargin at screen 0.95\n')

#       gw('set pm3d map\n')
        gw('set view map\n')
        gw('unset key\n')
        
        gw(f'xmin = {center_long-dcenter_long}\n')
        gw(f'xmax = {center_long+dcenter_long}\n')
        gw(f'ymin = {center_lat-dcenter_lat}\n')
        gw(f'ymax = {center_lat+dcenter_lat}\n')
        gw(f'cbmin = {min_radiance}\n')
        gw(f'cbmax = {max_radiance}\n')
#       gw('set xrange [xmin:xmax]\n')
#       gw('set yrange [ymin:ymax]\n')
        gw(f'set yrange [0:{Nlong-1}]\n')
        gw(f'set xrange [0:{Nlat-1}]\n')
        gw('set zrange [cbmin:cbmax]\n')
        gw('set cbrange [cbmin:cbmax]\n')
        gw('set palette defined (' + \
                '0 "white", ' + \
                '0.2 "yellow", ' + \
                '1 "red"' + \
                ')\n')
        gw('set xlabel "Longitude" font ",24"\n')
        gw('set ylabel "Latitude" font ",24"\n')
        gw('set cblabel "Radiance" font ",24"\n')
    #   gw('set format x ""\n')
    #   gw('set grid xtics lt 0 lw 2\n')

        gw('unset ytics\n')
        gw('set ytics (' + \
                f'"{center_long-dcenter_long}" 0, ' + \
                f'"{center_long+dcenter_long}" {Nlong-1}' + \
                ')\n')
        gw('unset xtics\n')
        gw('set xtics (' + \
                f'"{center_lat-dcenter_lat}" 0, ' + \
                f'"{center_lat+dcenter_lat}" {Nlat-1}' + \
                ')\n')
    
#       gw('plot "'+cwd+'exp003heatmap.dat" ' + \
#               'matrix u 1:2:3 w image\n')
        gw('plot "'+cwd+'exp003heatmap.dat" ' + \
                'u 1:2:3 w image\n')
    
        gw('unset lmargin\n')
        gw('unset rmargin\n')
        gw('unset bmargin\n')
        gw('unset tmargin\n')
        gw('set size 0.5, 0.5\n')
        gw('set origin 0.5, 0.0\n')
        gw(f'xmin = {gif_start_time}\n')
        gw(f'xmax = {gif_end_time}\n')
        gw('set xrange [xmin:xmax]\n')
        gw('set yrange [cbmin:cbmax]\n')
        gw('set xlabel "UNIXtime" font ",24"\n')
        gw('set ylabel "Radiance" font ",24"\n')
        gw('unset ytics\n')
        gw('set ytics\n')
        gw('unset xtics\n')
        gw('set xtics\n')
        gw('plot "'+cwd+'exp003.dat" ' + \
                'u 1:2 w l')
        
        gnuplotfile.close()
        
        os.system(("gnuplot < "+ cwd + "gnuplotfile"))

###############################################
# Data Analysis
###############################################

# We aren't really sure about this
# but, again, we'll use the
# minimum flyover period of twice
# per day
deltaTime_max = \
        getFlyoverPeriod(0.0,0.0)

# Records how many unique flyovers
# occured
Nunique = 0
Nframes_total = 0

# Record the maximum radiance for
# the 4 micrometer for some
# UNIXtime
maxRadianceData = [0]

# Record how many hotspot pixels
# were recorded at that time
NhotspotData = [0]

# Record the UNIX time for each
# unique data point
UNIXtimeData = [0]

# Initialize
background_4radiance = \
        getSpectralRadiance2(\
        10.0,temps[0],0.0,\
        wavelengths[0])
local_4radiances = [\
        [background_4radiance\
        for i in range(Nlong)]\
        for j in range(Nlat)]

maxRadiance = 0
Nhotspot = 0
missingData = [1]
#anomalousData = [1]
anomalousData = []
suspiciousData = [1]

latitude1 = -155.6
longitude1 = 19.2

for i in range(N-1):

    if (UNIXtimes[i] > gif_end_time): break

    # Get the radiance (use band
    # 22 preferentially)
    if (radiances[1][i] < 0):
        currentRadiance = \
            radiances[0][i]
    else:
        currentRadiance = \
            radiances[1][i]

    nlat_min = int(np.floor(\
            (latitudes[i] - 5*dlat -\
            (center_lat - \
            dcenter_lat)) / dlat))
    nlat_max = int(np.floor(\
            (latitudes[i] + 5*dlat -\
            (center_lat - \
            dcenter_lat)) / dlat))
    if (nlat_min < 0): \
            nlat_min = 0
    if (nlat_max < 0): \
            nlat_max = 0
    if (nlat_min >= Nlat) : \
            nlat_min = Nlat - 1
    if (nlat_max >= Nlat) : \
            nlat_max = Nlat - 1

    nlong_min = int(np.floor(\
            (longitudes[i] - 5*dlat -\
            (center_long - \
            dcenter_long)) / dlong))
    nlong_max = int(np.floor(\
            (longitudes[i] + 5*dlat -\
            (center_long - \
            dcenter_long)) / dlong))
    if (nlong_min < 0): \
            nlong_min = 0
    if (nlong_max < 0): \
            nlong_max = 0
    if (nlong_min >= Nlong) : \
            nlong_min = Nlong - 1
    if (nlong_max >= Nlong) : \
            nlong_max = Nlong - 1

    print(Nlat)
    print(Nlat*dlat)
    print(nlat_min)
    print(nlat_max)
    print(Nlong)
    print(Nlong*dlong)
    print(nlong_min)
    print(nlong_max)
    print(currentRadiance)
    print("")
    
    for j in range(nlat_max-nlat_min+1):
        for k in range(nlong_max-nlong_min+1):
            nlat = nlat_min + j
            nlong = nlong_min + k
            local_4radiances[nlat][nlong] =\
                    currentRadiance

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
    maxRadianceData.append(
        maxRadiance)
    NhotspotData.append(
        Nhotspot + 1)
    UNIXtimeData.append(
        UNIXtimes[i])

    # This was one unique flyover
    Nunique = Nunique + 1

    Nframes = np.floor(deltaTime/\
            time_per_frame)

    makeImages(local_4radiances,\
            int(Nframes_total),\
            int(Nframes))

    Nframes_total = \
            Nframes_total + Nframes

    if (i < N-1): \
    background_4radiance = \
            getSpectralRadiance2(\
            10.0,temps[i+1],0.0,\
            wavelengths[0])
    local_4radiances = [\
            [background_4radiance\
            for i in range(Nlong)]\
            for j in range(Nlat)]

    # There is a corresponding
    # maximum drop in radiance
    # for this period of time
    deltaRadiance_max = \
        getSpectralRadiance2(\
        getSurfaceTemp(STARTtime),\
        temps[i],0.90,\
        wavelengths[0])\
      - getSpectralRadiance2(\
        getSurfaceTemp(STARTtime\
          + deltaTime),\
        temps[i],0.90,\
        wavelengths[0])

    # And if we have dropped this
    # amount, record this as an
    # anomaly
    deltaRadiance = \
        maxRadianceData[Nunique-1]\
        - maxRadianceData[Nunique]

    if (deltaRadiance > \
        deltaRadiance_max):
        anomalousData.append(1)
    else:
        anomalousData.append(0)

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

        # This missing data is
        # suspicious if the previous
        # radiance was very high
        deltaRadiance = \
            maxRadiance-nextRadiance

        if (deltaRadiance > \
            deltaRadiance_max):
            suspiciousData.append(1)
        else:
            suspiciousData.append(0)

        missingData.append(1)

    # Otherwise, do nothing
    else:
        missingData.append(0)
        suspiciousData.append(0)

    # Reset the values
    maxRadiance = 0
    Nhotspot = 0

os.system(f'ffmpeg -r {frames_per_second} ' + \
        '-f image2 -s 1920x1080 ' + \
        '-i ' + cwd +'png/%04d.png ' + \
        '-vcodec libx264 ' + \
        '-crf 25  -pix_fmt yuv420p ' + \
        cwd + 'exp003.mp4')

print("Script 2 successully exited!")
