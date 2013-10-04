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
04/20/12 (Bourque)
'''

import os
import sys
import glob
import pyraf
from pyraf import iraf
from iraf import noao, digiphot, daophot
import numpy as np
from numpy import *
import pyfits
import matplotlib.pyplot as mpl


def run_all(root, data_list, threshold, tolerance, apertures):
    '''
    The parent function that runs all other functions.
    '''
    
    imlist, explist, zptlist, fwhmlist, sigmalist = load_data(data_list)
    
    remove_existing_files()
    unit_converter(imlist, explist, pam)
    source_find(root, threshold, fwhmlist, sigmalist)
    source_match(root, imlist, tolerance)
    do_photometry(root, apertures, zptlist)
    textdump(root)
    insmag_calc(root, explist, zptlist)
    determine_apcorr(root, explist, zptlist, imlist, apertures)
    write_results(root, explist, zptlist, imlist, apertures)
 
 
def load_data(data_list):
    '''
    Loads in information about images to be used for photometry.  Reads in
    data_list, which should be a file containing (5) columns: filename,
    exposure time, photometric zeropoint, fwhm and standard deviation of
    the background.
    '''
 
    imlist = loadtxt(data_list, dtype='str', unpack=True, usecols=[0])
    explist, zptlist, fwhmlist, sigmalist \
           = loadtxt(data_list, unpack=True, usecols=[1,2,3,4])
    return imlist, explist, zptlist, fwhmlist, sigmalist
           
           
def remove_existing_files():
    '''
    Removes sourcelists and photometry files, if they exist.
    '''
    file_list = glob.glob('*')
    remove_list = ['_counts.fits', '_sources.dat', '_match.dat', '_phot.dat']
    
    for file in file_list:
        if file[9:] in remove_list: 
            if os.path.exists(file):
                os.remove(file)

    
def unit_converter(imlist, explist, pam):
    '''
    Converts DRZ images from units of counts/sec to counts
    '''
    
    for image, exptime in zip(imlist, explist):
        hdu = pyfits.open(image)
        hdr = hdu[0].header
        detector = hdr['DETECTOR']
        hdu.close()
        iraf.imarith(operand1=image + '[1]', op='*', operand2=exptime, \
                     result=image[:9] + '_counts.fits')
    
        
def source_find(root, threshold, fwhmlist, sigmalist):
    '''
    Uses iraf.daofind to find sources.  datapars and findpars parameters
    are hard coded here.
    '''
    countlist = glob.glob(root + '/*counts.fits')
    
    for image, fwhm, sigma in zip(countlist, fwhmlist, sigmalist):
        iraf.datapars.fwhmpsf = fwhm
        iraf.datapars.sigma = sigma
        iraf.datapars.gain = 'ccdgain'
        iraf.datapars.datamax = 60000
        iraf.findpars.threshold = threshold
    
        iraf.daofind(image=image, output=image[:-12] + '_sources.dat', \
                     verify='no', verbose='no')

    
def source_match(root, imlist, tolerance):
    '''
    Uses iraf.xyxymatch to match sources between images.
    '''
    
    countlist = glob.glob(root + '/*counts.fits')
        
    sourcelist = glob.glob(root + '/*sources.dat')
    refimage = sourcelist[0]
    for image, sources in zip(countlist, sourcelist):
        iraf.xyxymatch(input=sources, reference=refimage, \
                       output=image[:-9] + '_match.dat', tolerance=tolerance, \
                       verbose = 'no')
    
    
def do_photometry(root, apertures, zptlist):                
    '''
    Uses iraf.daophot to perform aperture photometry.  centerpars,
    fitskypars and photpars parameters are hard coded here.
    '''
    
    coordinates = glob.glob(root + '/*match.dat')
    countlist = glob.glob(root + '/*counts.fits')
    
    for image, coords, zpt in zip(countlist, coordinates, zptlist):
        iraf.photpars.apertures = apertures
        iraf.photpars.zmag = zpt
        iraf.daophot.phot(image=image, coords=coords, \
                          output=image[:-12] + '_phot.dat', verify='no', \
                          verbose = 'no')
    
    
def textdump(root):
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
    
    
def insmag_calc(root, explist, zptlist):
    '''
    Calculates the instrumental magnitude and difference in magnitude 
    between apertures.  Calculates the zeropoint offset.  Calculates the
    VEGAmag magnitudes.  
    '''
    
    txtdumplist = glob.glob(root + '/*txtdump.dat')
    
    for image, exptime, zpt in zip(txtdumplist, explist, zptlist):
        counts1, counts2 = loadtxt(image[:-12] + '_txtdump.dat', unpack=True, \
                                   usecols=[6,7])
        
        instru_mag1 = -2.5*np.log10(counts1/exptime)
        instru_mag2 = -2.5*np.log10(counts2/exptime)
        delta_mag = (instru_mag1 - instru_mag2)  
        
        return counts1, counts2, instru_mag1, instru_mag2, delta_mag
        
        
def determine_apcorr(root, explist, zptlist, imlist, apertures):# Make more interactive
    '''
    Creates a plot of delta_mag versus instru_mag to determine in which
    region(s) to compute zpt_off
    '''
    
    counts1, counts2, instru_mag1, instru_mag2, delta_mag, = \
        insmag_calc(root, explist, zptlist)
    
    for image, zpt in zip(imlist, zptlist):
        mpl.rcParams['font.family'] = 'Times New Roman'
        mpl.rcParams['font.size'] = 12.0
        mpl.rcParams['xtick.major.size'] = 10.0
        mpl.rcParams['xtick.minor.size'] = 5.0
        mpl.rcParams['ytick.major.size'] = 10.0
        mpl.rcParams['ytick.minor.size'] = 5.0
        mpl.minorticks_on()
        mpl.ion()
        mpl.title(image)
        mpl.ylabel('$\Delta$ mag (r=' + apertures[0] + ', ' + \
                   apertures[-1] + ')')
        mpl.xlabel('-2.5 log(flux)')
        mpl.scatter(instru_mag1, delta_mag, s=1, c='k')
        mpl.savefig(image[:-9] + '_apcorr.png')
        mpl.clf()
        # Show plot on screen
        
        zpt_off_calc(instru_mag1, delta_mag, left, right, top bottom, zpt)
        

def zpt_off_calc(instru_mag1, delta_mag, left, right, top, bottom, zpt):

    sum = 0.0
    tot = 0
    for width in instru_mag1:
        if left < width < right:
            for height in delta_mag:
                if top < height < bottom:
                    sum = sum + height
                    tot = tot + 1
    zpt_off =  (sum / tot)
    vegamag = zpt - zpt_off + instru_mag1


def write_results(root, explist, zptlist, imlist, apertures):
    '''
    Writes desired results to a file
    '''
    
    counts1, counts2, instru_mag1, instru_mag2, delta_mag = \
        mag_calc(root, explist, zptlist)
    vegamag = 
        
    for image in imlist:
        np.savetxt(image[:9] + '_results.dat', zip(counts1, counts2, \
                   instru_mag1, instru_mag2, delta_mag, vegamag), \
                   fmt='%f %f %f %f %f %f')
        file = open(image[:9] + '_results.dat', 'r+')
        body = file.read()
        file.seek(0)
        header = '# COUNTS(r=' + apertures[0] + ')'
        header += ' COUNTS(r=' + apertures[-1] + ')'
        header += ' MAG(r=' + apertures[0] + ')'
        header += ' MAG(r=' + apertures[-1] + ')'
        header += ' DELTAMAG  VEGAMAG\n'
        header += body
        file.write(header)
        file.close()
    

if __name__ == '__main__':

    root = os.getcwd()
    data_list = raw_input('Data list: ')
    threshold = raw_input('Source finding threshold: ')
    tolerance = raw_input('Matching tolerance: ')
    apertures = raw_input('Aperture radii (Separate w/ comma): ')
    
    run_all(root, data_list, threshold, tolerance, apertures)