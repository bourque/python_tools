# ! /usr/bin/env python
'''
ABOUT:
This program performs basic aperture photometry though PyRAF.  Photometry
is performed on a list of images within image_list.txt.  image_list.txt
should contain (5) columns: filename, exposure time, zeropoint, fwhm
and sigma.

DEPENDS:
Python 2.7.1
Numpy 1.5.1

AUTHOR:
Matthew Bourque
Space Telescope Science Institute
bourque@stsci.edu

LAST UPDATED:
05/14/12 (Bourque)
'''

import os
import sys
import glob
import pyraf
from pyraf import iraf
from iraf import noao
from iraf import digiphot
from iraf import daophot
import numpy as np
from numpy import *
import pyfits
import matplotlib.pyplot as mpl


class Photometry():
    '''
    Parent class
    '''
    
    # ------------------------------------------------------------------------
    
    def __init__(self, root, data_list, threshold, tolerance, apertures):
        '''
        Loads in information about images to be used for photometry.  Reads in
        data_list, which should be a file containing (5) columns: filename,
        exposure time, photometric zeropoint, fwhm and standard deviation of
        the background.
        '''
        
        self.root = root
        self.data_list = data_list
        self.threshold = threshold
        self.tolerance = tolerance
        self.apertures = apertures
        
        self.imlist = loadtxt(self.data_list, dtype='str', unpack=True, \
                              usecols=[0])
        self.explist, self.zptlist, self.fwhmlist, self.sigmalist \
                    = loadtxt(self.data_list, unpack=True, usecols=[1,2,3,4])
        
    # ------------------------------------------------------------------------
    
    def run_all(self):
        '''
        A wrapper function that calls other functions in the appropriate
        order.
        '''
        
        self.remove_existing_files()
        self.unit_converter()
        self.source_find()
        self.source_match()
        self.do_photometry()
        self.textdump()
        self.mag_calc()
        self.apcorr_plot()
        self.write_results()
        
    # ------------------------------------------------------------------------
        
    def remove_existing_files(self):
        '''
        Removes specific photometry files, if they exist, in order to avoid
        PyRAF crashing
        '''
        
        file_list = glob.glob('*')
        remove_list = ['_counts.fits', '_sources.dat', '_match.dat']
        remove_list += ['_phot.dat']
        
        for file in file_list:
            if file[9:] in remove_list:
                if os.path.exists(file):
                    os.remove(file)
                    
    # ------------------------------------------------------------------------
    
    def unit_converter(self):
        '''
        Converts DRZ images from units of counts/sec to counts
        '''
        
        for image, exptime in zip(self.imlist, self.explist):
            iraf.imarith(operand1=image + '[1]', op='*', operand2=exptime, \
                         result=image[:9] + '_counts.fits')
                         
    # ------------------------------------------------------------------------
    
    def source_find(self):
        '''
        Uses iraf.daofind to find sources.  datapars and findpars parameters
        are hardcoded here.
        '''
        
        cntlist = glob.glob(self.root + '/*counts.fits')
        
        for image, fwhm, sigma in zip(cntlist, self.fwhmlist, self.sigmalist):
            iraf.datapars.fwhmpsf = fwhm
            iraf.datapars.sigma = sigma
            iraf.datapars.gain = 'ccdgain'
            iraf.datapars.datamax = 60000
            iraf.findpars.threshold = self.threshold
            
            iraf.daofind(image=image, output=image[:-12] + '_sources.dat', \
                         verify='no', verbose='no')
                         
    # ------------------------------------------------------------------------
        
    def source_match(self):
        '''
        Uses iraf.xyxymatch to match sources between images.
        '''
        
        countlist = glob.glob(self.root + '/*counts.fits')
        sourcelist = glob.glob(self.root + '/*sources.dat')
        refimage = sourcelist[0]
        
        for image, source in zip(countlist, sourcelist):
            iraf.xyxymatch(input=source, reference=refimage, \
                           output=image[:-12] + '_match.dat', \
                           tolerance=self.tolerance, verbose='no')
                           
    # ------------------------------------------------------------------------
    
    def do_photometry(self):
        '''
        Uses iraf.daophot to perform aperture photometry.  centerpars,
        fitskypars and photpars parameters are hard coded here.
        '''
        
        coordinates = glob.glob(self.root + '/*match.dat')
        countlist = glob.glob(self.root + '/*counts.fits')
        
        for image, coords, zpt in zip(countlist, coordinates, self.zptlist):
            iraf.photpars.apertures = self.apertures
            iraf.photpars.zmag = zpt
            iraf.daophot.phot(image=image, coords=coords, \
                              output=image[:-12] + '_phot.dat', verify='no', \
                              verbose='no')
                              
    # ------------------------------------------------------------------------
    
    def textdump(self):
        '''
        Prepares photometry data for science
        '''
        
        countlist = glob.glob(root + '/*counts.fits')
        photlist = glob.glob(root + '/*phot.dat')
        
        for image, phot_file in zip(countlist, photlist):
            saveout = sys.stdout
            fsock = open(image[:-12] + '_txtdump_tmp.dat', 'w')
            sys.stdout = fsock
            iraf.txdump(textfiles=phot_file, \
                        fields='xcenter,ycenter,id,msky,area,flux,mag,merr', \
                        expr='yes', headers='yes', parameters='no')
            sys.stdout = saveout
            fsock.close()
            sfile = open(image[:-12] + '_txtdump_tmp.dat')
            rfile = open(image[:-12] + '_txtdump.dat', 'w')
            for s in sfile:
                rfile.write(s.replace('INDEF', ' 99.99 '))
            sfile.close()
            rfile.close()
            os.remove(image[:-12] + '_txtdump_tmp.dat')
            
    # ------------------------------------------------------------------------
    
    def mag_calc(self):
        '''
        Calculates the instrumental magnitude and difference in magnitude 
        between apertures.  Calculates the zeropoint offset.  Calculates the
        VEGAmag magnitudes.  
        '''
        
        txtdumplist = glob.glob(root + '/*txtdump.dat')
    
        for image, exptme, zpt in zip(txtdumplist, self.explist, self.zptlist):
            counts1, counts2 = loadtxt(image[:-12] + '_txtdump.dat', \
                                       unpack=True, usecols=[6,7])
            
            insmag1 = -2.5*np.log10(counts1/exptme)
            insmag2 = -2.5*np.log10(counts2/exptme)
            delta_mag = (insmag1 - insmag2)  
            
            sum = 0.0
            tot = 0
            for width in insmag1:
                if -6.6 < width < -3.29:
                    for height in delta_mag:
                        if 0.035 < height < 0.087:
                            sum = sum + height
                            tot = tot + 1
            zpt_off =  (sum / tot)
            vegamag = zpt - zpt_off + insmag1
            
            return counts1, counts2, insmag1, insmag2, delta_mag, vegamag
            
    # ------------------------------------------------------------------------
    
    def apcorr_plot(self):
        '''
        Creates a plot of delta_mag versus instru_mag to determine in which
        region(s) to compute zpt_off
        '''
        
        counts1, counts2, insmag1, insmag2, delta_mag, vegamag = \
            self.mag_calc()
        
        for image in self.imlist:
            mpl.rcParams['font.family'] = 'Times New Roman'
            mpl.rcParams['font.size'] = 12.0
            mpl.rcParams['xtick.major.size'] = 10.0
            mpl.rcParams['xtick.minor.size'] = 5.0
            mpl.rcParams['ytick.major.size'] = 10.0
            mpl.rcParams['ytick.minor.size'] = 5.0
            mpl.minorticks_on()
            mpl.ion()
            mpl.title(image)
            mpl.ylabel('$\Delta$ mag (r=' + self.apertures[0] + ', ' + \
                       self.apertures[-1] + ')')
            mpl.xlabel('-2.5 log(flux)')
            mpl.scatter(insmag1, delta_mag, s=1, c='k')
            mpl.savefig(image[:-9] + '_apcorr.png')
            mpl.clf()
            
    # ------------------------------------------------------------------------
    
    def write_results(self):
        '''
        Writes desired results to a file.
        '''
        
        counts1, counts2, insmag1, insmag2, delta_mag, vegamag = \
            self.mag_calc()
            
        for image in self.imlist:
            np.savetxt(image[:9] + '_results.dat', zip(counts1, counts2, \
                       insmag1, insmag2, delta_mag, vegamag), \
                       fmt='%f %f %f %f %f %f')
            file = open(image[:9] + '_results.dat', 'r+')
            body = file.read()
            file.seek(0)
            header = '# COUNTS(r=' + self.apertures[0] + ')'
            header += ' COUNTS(r=' + self.apertures[-1] + ')'
            header += ' MAG(r=' + self.apertures[0] + ')'
            header += ' MAG(r=' + self.apertures[-1] + ')'
            header += ' DELTAMAG  VEGAMAG\n'
            header += body
            file.write(header)
            file.close()
        
    # ------------------------------------------------------------------------

if __name__ == '__main__':
    
    root = os.getcwd()
    data_list = raw_input('Data list: ')
    threshold = raw_input('Source finding threshold: ')
    tolerance = raw_input('Matching tolerance: ')
    apertures = raw_input('Aperture radii (Separate w/ comma): ')
    
    p = Photometry(root, data_list, threshold, tolerance, apertures)
    p.run_all()