#! /usr/bin/env python

'''
ABOUT:
This program creates average row and/or average column plots on FITS images.
Plots can be created for individual images or a list of images placed in a
text file, one image per line.  One can specify to plot average row and/or
average column with the -p or --plot_type argument.  Specifying the
-a or --all_switch will create one plot per image (if "off") or one plot
for all images in the list (if "on").  The save desitation of the plots
can be specified with the -s or --save_dst argument.  By default, both
average columns and average row plots will be created, one per image, and
saved in the current working directory.

The image of interest of the list of images of interest must the the first
argument for command line execution and must contain the extension and
indices to plot over.  For example:

python avg_row_col.py abcdefgh_flt.fits[1][100:300,550:650]
'''

import argparse
import numpy as np
import os
import pyfits
import matplotlib.pyplot as plt

class AvgRowCol():
    '''
    Parent class.
    '''

    # -------------------------------------------------------------------------

    def __init__(self, images, plot_type, all_switch, save_dst):
        '''
        Assigns argument variables to class instances.
        '''

        self.images = images
        self.plot_type = plot_type
        self.all_switch = all_switch
        self.save_dst = save_dst

    # -------------------------------------------------------------------------

    def get_image_list(self):
        '''
        Reads in one image or a list of images and returns a python list of the
        image(s).  Also ensures each image string contains FITS extension and
        indices.
        '''

        # For single image:
        if '.fits' in self.images:
            self.images = [self.images]

        # For a list of images:
        else:
            with open(self.images, 'r') as image_file:
                self.image_list = image_file.readlines()
            self.image_list = [line.strip() for line in self.image_list]

        # Ensure images have extension and indices by checking "[" and "]".
        for image in self.image_list:
            assert len([i for i in image if i == '[']) == 2, 'Missing or ' + \
                'Invalid extension or indices.'
            assert len([i for i in image if i == ']']) == 2, 'Missing or ' + \
                'Invalid extension or indices.'

    # -------------------------------------------------------------------------

    def plot_all_data(self, values):
        '''
        Creates the average col or row plot, plotting all images on the 
        same plot.
        '''

        fig = plt.figure()
        plt.minorticks_on()
        plt.grid()
        plt.xlabel(self.descrip + ' (pixels)', labelpad=10)
        plt.title('Average of ' + self.anti_descrip + 's')
        lengths = [len(value) for value in values]
        plt.xlim([0, max(lengths) - 1])
        
        for frame, ext, value in zip(self.frames, self.exts, values):
            plt.plot(value, markersize=4, markerfacecolor='none', 
                     label=frame + '[' + ext + ']')
        plt.legend()
        filename = os.path.join(self.save_dst, 
            'avg_' + self.anti_descrip.lower() + '_ext' + ext + '.png')
        plt.savefig(filename)
        print 'Saved figure to ' + filename

    # -------------------------------------------------------------------------

    def plot_single_data(self, values):
        '''
        Creates the average col or row plot for each image.
        '''

        for frame, ext, value in zip(self.frames, self.exts, values):
            fig = plt.figure()
            plt.minorticks_on()
            plt.grid()
            plt.xlabel(self.descrip + ' (pixels)', labelpad=10)
            plt.title('Average of ' + self.anti_descrip + 's')
            plt.xlim([0, len(values[0])])
            plt.plot(value, 'k', markersize=4, markerfacecolor='none')
            filename = os.path.join(save_dst, frame.split('.')[0] + '_avg_' + \
                                    self.anti_descrip.lower() + '_ext' + ext + \
                                    '.png')
            plt.savefig(filename)
            print 'Saved figure to ' + filename

    # -------------------------------------------------------------------------
    # The main controller
    # -------------------------------------------------------------------------

    def avg_row_col_main(self):
        '''
        The main controller.
        '''

        # Construct list of images to be examined
        self.get_image_list()

        # Parse filename, ext, and indices.
        self.frames = [image.split('[')[0] for image in self.image_list]
        self.exts = [image.split('[')[1][0] for image in self.image_list]
        y_indices = [image.split('[')[-1].split(',')[0] for image in 
                     self.image_list]
        x_indices = [image.split('[')[-1].split(',')[1].strip(']') 
                     for image in self.image_list]
        x1s = [int(x.split(':')[0]) for x in x_indices]
        x2s = [int(x.split(':')[1]) for x in x_indices]
        y1s = [int(y.split(':')[0]) for y in y_indices]
        y2s = [int(y.split(':')[1]) for y in y_indices]

        # Read in data
        ext_data_list = []
        for frame,ext in zip(self.frames, self.exts):
            open_frame = pyfits.open(frame)
            ext_data_list.append(open_frame[int(ext)].data)
            open_frame.close()

        # Calculate average row or column
        avg_row_list = [ext_data[x1:x2,y1:y2].mean(axis=1) for ext_data, 
                        x1, x2, y1, y2 in zip(ext_data_list, x1s, x2s, 
                        y1s, y2s)]
        avg_col_list = [ext_data[x1:x2,y1:y2].mean(axis=0) for ext_data, 
                        x1, x2, y1, y2 in zip(ext_data_list, x1s, x2s, 
                        y1s, y2s)]

        # Set plotting parameters
        plt.rcParams['legend.fontsize'] = 10
        plt.rcParams['font.family'] = 'Helvetica'
        plt.minorticks_on()

        # Plot the data
        if self.plot_type == 'row' or self.plot_type == 'both':
            self.descrip = 'Row'
            self.anti_descrip = 'Column'
            if self.all_switch == 'off':
                self.plot_single_data(avg_row_list)
            elif self.all_switch == 'on':
                self.plot_all_data(avg_row_list)
        if self.plot_type == 'col' or self.plot_type == 'both':
            self.descrip = 'Column'
            self.anti_descrip = 'Row'
            if self.all_switch == 'off':
                self.plot_single_data(avg_col_list)
            elif self.all_switch == 'on':
                self.plot_all_data(avg_col_list)


# -----------------------------------------------------------------------------
# For command line execution
# -----------------------------------------------------------------------------

def parse_args():
    '''
    Parse command line arguments, returns args object.
    '''

    # Create help string
    image_help = 'The image(s) (with extension and indices) to be ' + \
                 'examined.  This could be a single image or a text file' + \
                 'containing multiple images.'
    plot_type_help = 'The type of plot to be produced. This can be "row"' + \
                     ' for an average row plot, "col" for an average ' + \
                     'column plot, or "both" to produce both types.'
    all_switch_help = 'If "on", will plot all images on the same plot. ' + \
                      'Conversely, if "off", will create a separate plot ' + \
                      'for each image.'
    save_dst_help = 'The path to where the plots will be saved.  If no ' + \
                    'value is given, the plots will be displayed on the ' + \
                    'screen.'

    # Add time arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('images', type=str, help=image_help)
    parser.add_argument('-p', '--plot_type', dest='plot_type', 
                        action='store', type=str, required=False, 
                        default='both', help=plot_type_help)
    parser.add_argument('-a', '--all_switch', dest = 'all_switch', 
                        action='store', type=str, required=False, 
                        default='off', help=all_switch_help)
    parser.add_argument('-s', '--save_dst', dest='save_dst', 
                        action='store', type=str, required=False, 
                        default=os.getcwd() + '/', help=save_dst_help)

    # Parse args
    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------

def test_args(args):
    '''
    Ensure valid command line arguments.
    '''

    # Assert image or image list exists.
    assert os.path.exists(args.images) == True, 'File ' + images + 'does' + \
        'not exist.'

    # Assert plot_type is "row", "col", or "both".
    valid_plot_types = ['row', 'col', 'both']
    assert args.plot_type in valid_plot_types, 'Invalid plot_type. ' + \
        'plot_type can be "row", "col", or "both".'

    # Assert all_switch is "on" or "off".
    valid_all_switches = ['on', 'off']
    assert args.all_switch in valid_all_switches, 'Invalid all_switch. ' + \
        'all_switch can be "on" or "off".'

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    args = parse_args()
    test_args(args)

    avg_row_col = AvgRowCol(args.images, args.plot_type, args.all_switch,
                            args.save_dst)
    avg_row_col.avg_row_col_main()