#! /usr/bin/env python

'''
ABOUT:
This program calculates image statistics for a rectangular annulus.

The path to the image to be anaylzed must be the first argument for command
line execution.  The path to a file that contains X and Y coordinates for the 
two rectangles must be the second argument for command line execution.  
For example:

python imstat_rect.py image.fits coords.dat

where coords.dat contains:

#box1x1 box1x2 box1y1 box1y2 box2x1 box2x2 box2y1 box2y2
1.0 3.0 2.0 4.0 0.5 3.5 1.0 5.0
12.5 13.5 6.3 9.6 0 10 0 10
etc.

Note that, if box 2 coordinates are 0 0 0 0, the program will return statistics
for all pixels within box 1. Also note that, if there are two boxes, box 1 must 
be the inner box, and box2 must be the outer box.

The program will save the statistics to a file named <image>_<coords>.dat.
'''

import argparse
from astropy.io import ascii
from astropy.io import fits as pyfits
import numpy as np
import os
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------

class ImStatBox():
    '''
    Parent class.
    '''

    # -------------------------------------------------------------------------

    def __init__(self, image, coord_list):
        '''
        Assigns argument variables to class instances.
        '''

        self.image = image
        self.coord_list = coord_list

    # -------------------------------------------------------------------------

    def determine_region(self):
        '''
        Determines whether the statistics region is a box or an annlulus based
        on the coordinates of 'box2'
        '''

        # Determine if box2 coordinates are all 0s. If they are, set
        # region_flag to True
        self.region_flag_list = []
        box2_list = [self.box2x1, self.box2x2, self.box2y1, self.box2y2]
        for box2 in box2_list:
            if all(element == 0 for element in box2) == True:
                self.region_flag_list.append(True)
            else:
                self.region_flag_list.append(False)
        
        # If there are two boxes, stats region is annulus. If there is one
        # box, stats region is a box.
        if all(self.region_flag_list) == True:
            self.region = 'box'
        else:
            self.region = 'annulus'

    # -------------------------------------------------------------------------

    def get_coordinates(self):
        '''
        Reads in the data from the coordinates file.
        '''

        self.data = ascii.read(self.coord_list, data_start=0, names=['box1x1', 
            'box1x2', 'box1y1', 'box1y2', 'box2x1', 'box2x2', 'box2y1',
            'box2y2'])
        self.box1x1 = self.data['box1x1']
        self.box1x2 = self.data['box1x2']
        self.box1y1 = self.data['box1y1']
        self.box1y2 = self.data['box1y2']
        self.box2x1 = self.data['box2x1']
        self.box2x2 = self.data['box2x2']
        self.box2y1 = self.data['box2y1']
        self.box2y2 = self.data['box2y2']

    # -------------------------------------------------------------------------

    def get_image(self):
        '''
        Uses pyfits to read in the zeroth extension image data.
        '''

        self.frame = pyfits.open(self.image)[0].data

    # -------------------------------------------------------------------------

    def perform_statistics(self):
        '''
        Calculates basic statistics of the region.
        '''

        self.npix_list = [region.size for region in self.region_list]
        self.mean_list = [np.mean(region) for region in self.region_list]
        self.midpt_list = [np.median(region) for region in self.region_list]
        self.stdev_list = [np.std(region) for region in self.region_list]
        self.min_list = [np.min(region) for region in self.region_list]
        self.max_list = [np.max(region) for region in self.region_list]

    # -------------------------------------------------------------------------

    def write_statistics(self):
        '''
        Writes statistics to output file.
        '''

        num_regions = len(self.mean_list)
        filename = '{}_{}.dat'.format(self.image.split('.')[0], 
            self.coord_list.split('.')[0])
        with open(filename, 'w') as f:
            for i in range(num_regions):
                f.write(
                    'Region {} npix: {}\n'.format(i, self.npix_list[i]) +
                    'Region {} mean: {}\n'.format(i, self.mean_list[i]) +
                    'Region {} midpt: {}\n'.format(i, self.midpt_list[i]) +
                    'Region {} stdev: {}\n'.format(i, self.stdev_list[i]) +
                    'Region {} min: {}\n'.format(i, self.min_list[i]) +
                    'Region {} max: {}\n'.format(i, self.max_list[i])
                    )

    # -------------------------------------------------------------------------
    # The main controller
    # -------------------------------------------------------------------------

    def imstat_box_main(self):
        '''
        The main controller.
        '''

        self.get_image()
        self.get_coordinates()
        self.determine_region()

        if self.region == 'annulus':

            # Set small boxes to all 0s
            for box1y1, box1y2, box1x1, box1x2 in zip(self.box1y1, self.box1y2,
                self.box1x1, self.box1x2):
                    self.frame[box1y1:box1y2,box1x1:box1x2] = 0.0

            # Build large box
            self.large_box_list = [
                self.frame[b2y1:b2y2,b2x1:b2x2] for b2y1,b2y2,b2x1,b2x2 in 
                zip(self.box2y1,self.box2y2,self.box2x1,self.box2x2)]

            # Find indices whose values are not 0.0:
            self.region_list = [large_box[large_box != 0.0] for large_box in 
                self.large_box_list]

        elif self.region == 'box':
            self.region_list = [
                self.frame[b1y1:b1y2,b1x1:b1x2] for b1y1,b1y2,b1x1,b1x2 in
                zip(self.box1y1,self.box1y2,self.box1x1,self.box1x2)]

        self.perform_statistics()
        self.write_statistics()

# -----------------------------------------------------------------------------
# For command line execution
# -----------------------------------------------------------------------------

def parse_args():
    '''
    Parse command line arguments, returns args object.
    '''

    # Create help strings
    image_help = 'Path to image to be analyzed.'
    coord_list_help = 'Path to file containing rectangle coordinates.'

    # Add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('image', type=str, help=image_help)
    parser.add_argument('coord_list', type=str, help=coord_list_help)

    # Parse args
    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------

def test_args(args):
    '''
    Ensure valid command line arguments.
    '''

    # Assert image and coords list exists.
    assert os.path.exists(args.image) == True, \
        'File {} does not exist'.format(args.image)
    assert os.path.exists(args.coord_list) == True, \
        'File {} does not exist'.format(args.coord_list)

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    args = parse_args()
    test_args(args)

    imstat_box = ImStatBox(args.image, args.coord_list)
    imstat_box.imstat_box_main()