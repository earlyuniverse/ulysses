#!/usr/bin/env python3
#
#----------------------------------------------------------------------
# $Id: example.py 1302 2021-12-06 10:09:55Z salam $
#
# Copyright (c) 2018-, Frederic A. Dreyer, Keith Hamilton, Alexander Karlberg,
# Gavin P. Salam, Ludovic Scyboz, Gregory Soyez, Rob Verheyen
#
#----------------------------------------------------------------------
# This file is part of FastJet contrib.
#
# It is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# It is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this code. If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------
#
# Load a sample file and plot it.
#
# Usage:
#   python3 example.py [--file filename] [--bkg file_bkg]
#                      [--njet njet]  [--npxl npixels]
#

import read_lund_json as lund
#from  import LundImage
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Plot lund images')
parser.add_argument('--file', action = 'store', default = 'jets.json')
parser.add_argument('--njet', type   = int,     default = 2,  help='Maximum number of jets to analyse')
parser.add_argument('--npxl', type   = int,     default = 25, help="Number of pixels in each dimension of the image")

args = parser.parse_args()

# set up the reader and get array from file
xval   = [0.0, 3.0]
yval   = [-3.0, 5.0]

# start by creating a reader for the json file produced by example.cc
# (one json entry per line, correspond to one jet per json entry)
reader = lund.Reader(args.file, args.njet)

# Then examine the jets it contains
print ("Contents of the file", args.file)
for jet in reader:
    # jet is an array of declusterings.
    # The jet's pt can be obtained by looking at the first declustering (jet[0])
    # and extracting the subjet-pair pt ("p_pt")
    print("  Jet with pt = {:6.1f} GeV with {:3d} primary Lund-plane declusterings".format(jet[0]["p_pt"], len(jet)))
print()

# Reset the reader to the start and use it with a helper
# class to extract an image for each jet
reader.reset()
image_generator = lund.LundImage(reader, args.njet, args.npxl, xval, yval)
images = image_generator.values()

# Get the average of the images
print("Now creating average lund image from the {} jets".format(len(images)))
avg_img = np.average(images,axis=0)

# Plot the result
fig=plt.figure(figsize=(6, 4.5))
plt.title('Averaged Lund image')
plt.xlabel('$\ln(R / \Delta)$')
plt.ylabel('$\ln(k_t / \mathrm{GeV})$')
plt.imshow(avg_img.transpose(), origin='lower', aspect='auto',
           extent=xval+yval, cmap=plt.get_cmap('BuPu'))
plt.colorbar()

print("Close the viewer window to exit")
plt.show()
