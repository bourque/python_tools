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
12.5 13.5 6.3 9.6 0 0 0 0
etc.

Note that, if box 2 coordinates are 0 0 0 0, the program will return statistics
for all pixels within box 1. Also note that, if there are two boxes, box 1 must 
be the inner box, and box2 must be the outer box.
'''

import argparse
from astropy.io import ascii
from astropy.io import fits as pyfits
import numpy as np
#np.set_printoptions(threshold=np.nan)
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

    def build_boxes(self):
        '''
        Constructs the large box and small box (2D arrays) based on the
        coordinates.
        '''

        self.large_box_list = [
            self.frame[b2x1:b2x2,b2y1:b2y2] for b2x1,b2x2,b2y1,b2y2 in 
            zip(self.box2x1,self.box2x2,self.box2y1,self.box2y2)]
        self.small_box_list = [
            self.frame[b1x1:b1x2,b1y1:b1y2] for b1x1,b1x2,b1y1,b1y2 in 
            zip(self.box1x1,self.box1x2,self.box1y1,self.box1y2)]

    # -------------------------------------------------------------------------

    def calc_avg(self):
        '''
        Calculates the average row and column value.
        '''

        self.avg_row_list = [ext_data[x1:x2,y1:y2].mean(axis=1) for ext_data, 
            x1, x2, y1, y2 in zip(self.ext_data_list, self.x1s, self.x2s, 
            self.y1s, self.y2s)]
        self.avg_col_list = [ext_data[x1:x2,y1:y2].mean(axis=0) for ext_data, 
            x1, x2, y1, y2 in zip(self.ext_data_list, self.x1s, self.x2s, 
            self.y1s, self.y2s)]

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
    # The main controller
    # -------------------------------------------------------------------------

    def imstat_box_main(self):
        '''
        The main controller.
        '''

        self.get_image()
        self.get_coordinates()

        # Set small boxes to all 0s
        for box1x1, box1x2, box1y1, box1y2 in zip(self.box1x1, self.box1x2, 
            self.box1y1, self.box1y2):
                self.frame[box1x1:box1x2,box1y1:box1y2] = 0.0

        # Build large box
        self.large_box_list = [
            self.frame[b2x1:b2x2,b2y1:b2y2] for b2x1,b2x2,b2y1,b2y2 in 
            zip(self.box2x1,self.box2x2,self.box2y1,self.box2y2)]

        # Find indices whose values are not 0.0:
        annulus_list = [large_box[large_box != 0.0] for large_box in 
            self.large_box_list]

        # Perform statistics
        mean_list = [np.mean(annulus) for annulus in annulus_list]
        midpt_list = [np.median(annulus) for annulus in annulus_list]
        stdev_list = [np.std(annulus) for annulus in annulus_list]
        min_list = [np.min(annulus) for annulus in annulus_list]
        max_list = [np.max(annulus) for annulus in annulus_list]

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

    # Add time arguments
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