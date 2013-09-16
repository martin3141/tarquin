#! /usr/bin/python

from math import *

# 1.5T
#ft = 63.684

# 3T
ft = 127.788

# 7T
#ft = 300

# 9.4T
#ft = 400

sig = ( "Lip13a", "Lip13b", "Lip09", "MM09", "Lip20", "Lip20", "Lip20", "MM20", "MM20", "MM20", "MM20", "MM12", "MM14", "MM17")
FWHM = ( 0.15, 0.089, 0.14, 0.14, 0.15, 0.15, 0.2, 0.15, 0.2, 0.15, 0.2, 0.15, 0.17, 0.15)

print "LCM MM/Lip Basis"

for n in range(0, len(FWHM)):
    print sig[n],"\t", FWHM[n],"\t", -pow((FWHM[n]*ft)*pi/2,2)/log(0.5)

#sig = ( "MM09a", "MM09b", "MM12", "MM14", "MM17", "MM20a", "MM20b", "MM20c", "MM20d", "MM20e", "MM20f","MM30", "MM32", "MM38a", "MM38b", "MM38c", "MM43" )
#FWHM = ( 0.1, 0.04, 0.1, 0.08, 0.07, 0.12, 0.08, 0.08, 0.09, 0.2, 0.05, 0.07, 0.10, 0.26, 0.15, 0.15, 0.16 )

sig = ( "MM09a", "MM09b", "MM12", "MM14", "MM17", "MM20a", "MM20b", "MM20c", "MM20d", "MM20e", "MM20f","MM30", "MM32", "MM38a", "MM38b", "MM38c" )
FWHM = ( 0.1, 0.04, 0.1, 0.08, 0.07, 0.12, 0.08, 0.08, 0.09, 0.2, 0.05, 0.07, 0.10, 0.26, 0.15, 0.15 )

print ""
print "High res MM Basis"

for n in range(0, len(FWHM)):
    print sig[n],"\t", FWHM[n],"\t", -pow((FWHM[n]*ft)*pi/2,2)/log(0.5)
