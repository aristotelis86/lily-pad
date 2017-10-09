#!/usr/bin/env python3
# Use this script to produce graphics and extract information for a specific
# sheet and/or control point. (Needs modification to call the functions for the 
# required analysis, example calls are given)

import sheetprocessing as myLib

sheetN = 0; # sheet (numbering of sheets from 0)
pointN = 16; # control point from previous sheet (numbering of cpoints from 0)

myLib.read_sheet_info( sheetN ) # just print available information on the sheet

myLib.energy_plot( sheetN ) # create+save energy evolution figure

myLib.length_plot( sheetN ) # create+save sheet's length evolution figure

myLib.point_plots( sheetN, pointN ) # create+save plots for position, velocity,
                                    # force and phase space in x,y directions 
                                    # for one control point.

myLib.frequency_plot( sheetN, pointN=-1 ) # create+save plots of the FFT analysis
                               # of the time series from both x and y 
                               # position of one control point.

myLib.normalMode_frequency_plot( sheetN, pointN=-1, mode=1, stretchRatio=0.1, align="y" )
# create+save the plot of the FFT analysis of the vibration (pinned-pinned)
# tracked by one control point. 
                                       
myLib.impulse_frequency_analysis( sheetN, newL=44.03, pointN=-1, g=10, align="y" )
# create+save the plot of the FFT analysis of the vibration (pinned-free)
# tracked by one control point. 
                                       
