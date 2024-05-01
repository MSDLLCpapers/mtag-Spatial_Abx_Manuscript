# bug_config.py
#
# The user may freely change the variables in this file.
#     These variables are used by all three python scripts.

# To process only one study, specify its name.  Tp process all studies, set to null ('').
PROCESS_ONE_STUDY = '3454_COLON'

# If we have a list of study names to input, which I usually don't because I run the "runall.py" code April 10
READ_STUDY_FILENAMES_FROM_FILE = 0

# To run the plot_heatmap script on all files using run_all_files.py and to toggle between 
# running plot_heatmap for each individual file or for the study as a summary 
WHICH_HEATMAP = "python3 plot_image_heatmap_3.py" # python3 plot_study_heatmap_3.py or python3 plot_image_heatmap_3.py

# To write the neighborhood list files, set to 1.  (These are not presently used.)
WRITE_NHOOD_LISTS = 1  # write neighborhood bug lists for each refBug?

# To write info on all 12 replicate and fov files for each study, set to 1.
WRITE_ALL_FILES = 1    # write all replicate and fov file info?

RADIUS_IN_MICRONS = 5  # radius of neighborhood in microns (5 microns = 5000 nm)
PIXEL_SIZE_IN_NM = 70  # diameter of an image pixel in nm

DEBUG = 0       # if 1, write trace output (no output if 0)
CHECK_DATA = 0  # if 1, print info on the data as it is read or processed.
USE_PANDAS = 0  # if 1, write a pandas dataframe (1) or write a bare csv (0)

# set a limit on max bugID (array size) then check it's not exceeded
MAX_BUG_ID = 500  # Bug IDs are 1 to 500

# input dir
INPUT_PATH = '../data/Spatial_and_Abundance'

# output dir for results
# RESULTS_PATH = '~/../results'
RESULTS_PATH = '../Results'

# Do we want a subfolder in results to keep output datasets distinct?
#     Specify a name here and the results folders will be 'results/round02/STUDYPREFIX'
#     rather than results/STUDYPREFIX
SUB_FOLDER = 'BugCountsR5_April102023'  # use either null ('') or a name (like 'round02')

# column number from Kanvas CSV for bug ID, X, Y
IDCOL = 5
XCOL = 2
YCOL = 3
