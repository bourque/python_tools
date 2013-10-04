import argparse
import numpy as np
import os
import pyfits
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------

def plot_all_data(frames, exts, values, save_dst, descrip):
	'''
	Creates the average col or row plot, plotting all images on the same plot.
	'''

	fig = plt.figure()
	plt.minorticks_on()
	plt.grid()
	plt.xlabel(descrip + ' (pixels)', labelpad = 10)
	plt.title('Average of ' + descrip + 's')
	lengths = [len(value) for value in values]
	plt.xlim([0, max(lengths) - 1])
	
	for frame, ext, value in zip(frames, exts, values):
		plt.plot(value, markersize = 4, markerfacecolor = 'none', 
				 label = frame + '[' + ext + ']')
	plt.legend()
	filename = save_dst + 'avg_' + descrip.lower() + '_ext' + ext + '.png'
	plt.savefig(filename)
	print 'Saved figure to ' + filename

# -----------------------------------------------------------------------------

def plot_single_data(frames, exts, values, save_dst, descrip):
	'''
	Creates the average col or row plot for each image.
	'''

	for frame, ext, value in zip(frames, exts, values):
		fig = plt.figure()
		plt.minorticks_on()
		plt.grid()
		plt.xlabel(descrip + ' (pixels)', labelpad = 10)
		plt.title('Average of ' + descrip + 's')
		plt.xlim([0, len(values[0])])
		plt.plot(value, 'k', markersize = 4, markerfacecolor = 'none')
		filename = save_dst + frame.split('.')[0] + '_avg_' + \
				   descrip.lower() + '_ext' + ext + '.png'
		plt.savefig(filename)
		print 'Saved figure to ' + filename

# -----------------------------------------------------------------------------
# The main controller
# -----------------------------------------------------------------------------

def avg_row_col_main(image, plot_type, all_switch, save_dst):
	'''
	The main controller.
	'''

	# Construct list of images to be examined
	if '.fits' in image:
		images = [image]
	else:
		with open(image, 'r') as image_file:
			images = image_file.readlines()
		images = [line.strip() for line in images]

	# Parse filename, ext, and indices.
	frames = [image.split('[')[0] for image in images]
	exts = [image.split('[')[1][0] for image in images]
	x_indices = [image.split('[')[-1].split(',')[0] for image in images]
	y_indices = [image.split('[')[-1].split(',')[1].strip(']') 
				 for image in images]
	x1s = [int(x.split(':')[0]) for x in x_indices]
	x2s = [int(x.split(':')[1]) for x in x_indices]
	y1s = [int(y.split(':')[0]) for y in y_indices]
	y2s = [int(y.split(':')[1]) for y in y_indices]

	# Read in data
	ext_data_list = []
	for frame,ext in zip(frames, exts):
		open_frame = pyfits.open(frame)
		ext_data_list.append(open_frame[int(ext)].data)
		open_frame.close()

	# Calculate average row or column
	avg_row_list = [ext_data[y1:y2,x1:x2].mean(axis = 1) for ext_data, y1, y2, 
				    x1, x2 in zip(ext_data_list, y1s, y2s, x1s, x2s)]
	avg_col_list = [ext_data[y1:y2,x1:x2].mean(axis = 0) for ext_data, y1, y2, 
					x1, x2 in zip(ext_data_list, y1s, y2s, x1s, x2s)]

	# Set plotting parameters
	plt.rcParams['legend.fontsize'] = 10
	plt.rcParams['font.family'] = 'Helvetica'
	plt.minorticks_on()

	# Plot the data
	if plot_type == 'row' or plot_type == 'both':
		descrip = 'Row'
		if all_switch == 'off':
			plot_single_data(frames, exts, avg_row_list, save_dst, descrip)
		elif all_switch == 'on':
			plot_all_data(frames, exts, avg_row_list, save_dst, descrip)
	if plot_type == 'col' or plot_type == 'both':
		descrip = 'Column'
		if all_switch == 'off':
			plot_single_data(frames, exts, avg_col_list, save_dst, descrip)
		elif all_switch == 'on':
			plot_all_data(frames, exts, avg_col_list, save_dst, descrip)


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
    parser.add_argument('-i', dest = 'images', action = 'store', 
    					type = str, required = True, help = image_help)
    parser.add_argument('-p', dest = 'plot_type', action = 'store',
                        type = str, required = False, default = 'both',
                        help = plot_type_help)
    parser.add_argument('-a', dest = 'all_switch', action = 'store',
    					type = str, required = False, default = 'off',
    					help = all_switch_help)
    parser.add_argument('-s', dest = 'save_dst', action = 'store',
    					type = str, required = False, 
    					default = os.getcwd() + '/',
    					help = save_dst_help)

    # Parse args
    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------

def test_args(args):
	'''
	Ensure valid command line arguments.
	'''

	# Assert image or image list exists.
	assert os.path.exists(args.images) == True, 'File ' + images + 'Does' + \
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

	avg_row_col_main(args.images, args.plot_type, args.all_switch, 
		args.save_dst)