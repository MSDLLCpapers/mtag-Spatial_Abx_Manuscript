"""
search_nhoods_2.py:

    This is the second step in a 3 script process.  It creates neighborhood count data for
    entire study, as well as for each file in a 12 file study.  So 13 files are written for
    each study (unless the global variable PROCESS_ALL_FILES is set to 0.  Then only a
    single count file for the study will be written).

To process only one study (rather than all available studies in a loop), set the
    PROCESS_ONE_STUDY constant to the name of the ONE study to process, like 3402_CECUM.
    If you set PROCESS_ONE_STUDY to '', it will process all studies in the input folder, taking
    more than 24 hours to finish.

This code processes all studies that exist in the input folder -- the 12 files in each study,
    and then combines them into one -- so it takes 13 passes through the study.
For each pass, there are two variables written:
- a 2D count of nhood bugs that surround each reference bug
  For one entire study and for each image, this produces 13 files.
- For every referenceBug in image, a file is written containing a list of all that bug's
  neighborhod bug IDs that are found within X microns.
  So if there are 40 bugs in a study, we write 40 files at study scale (summing all 12 images).
  Then we write ~40 files for each image file, or another ~480 more files.
  These are written to a subfolder called 'STUDYNAME/radius_5/nhood_lists', if the radius is 5 microns.

This code writes out two variables:

1) for each study's 12 files, write one file that contains a 500x500 matrix of counts (integers)
    for each possible reference bug (bugA), how many other candidateBugs appear
        within RADIUS_IN_PIXELS of bugA (are in its neighborhood)?
    This 500x500 value matrix is like a heat map of bug counts.
    The output filename will be [input filename prefix] / 'nhood_counts_radius_X_500.csv'
        where X is RADIUS_IN_NM.
    An example filename:
results/round02/3402_CECUM/nhood_counts_radius_5_500.csv  (if the nhood radius is 5 nm)

2) Let's say there are 50 unique bugs in a 12 file study, so I write 50 neighborhood list
    files for each study -- one file for each reference bug ID.
    This file contains all the neighborhood bugs found surrounding each instance of the
    reference bug.  If this bug (bug ID 7) occurs 45 times in the study, there will be 45
    lines in the output file (nhood_lists/bug_7.csv).  Each line contains the following:

    7, reference#, fov#, then the list of nearby bug IDs in that neighborhood

    (This file is for future use, none has been identified yet)

    Each file describes the neighborhoods of one referenceBug (refBug)
    Each row in the file lists the memberBugs found in each of the refBug's
        neighborhoods.  If refBug is found 75 times in the study, there will be 75
        rows in the file, each listing the IDs of all memberBugs found within
        RADIUS_IN_NM of each instance of refBug.
    So if there are 50 different bugs in a study, there will be 50 'nhood' files.
    Each of these files is named:
        results/round02/STUDYPREFIX/radius_R/bugID_X_nhood_lists.csv"
        The R is the RADIUS_IN_NM and the X is the reference bugID.
    The first 3 columns in each row are: refBugID, replicate number, FOV number
        the remaining columns are the list of all memberBugIDs in that row's neighborhood.

An example STUDYNAME is: 3402_CECUM_tissue_replicate_1_fov_1_cell_information_counts.csv
An example STUDYPREFIX is: 3402_CECUM

read in all 12 Kanvas input files and fill in a 500x500 matrix with
    the count of each bugID found
process all the Kanvas csv files in a folder and write two kinds of files per study:
1) one NxN neighborhood csv of counts, where each row shows, for bugA,  the count of each
    other bug in bugA's neighborhoods
2) N files, as many as there are different bugIDs, one file per bug ID,
    listing all the memberBugIDs found in each of the refBugID's neighborhoods.

The algorithm:
For each Kanvas study:
    create an empty square NxN matrix that has the same column and row labels
        (include all N possible bug IDs)
    Read in one input CSV file from Kanvas into a pandas data frame.
    Sort the frame's rows on centroid_X (column 3), to order them left to right
        along the X axis in image.
Select the first bug in the sequence as the first refBug -- that's the bug on row 1.

Search the neighborhood surrounding each refBug to find all nearby memberBugs:
    Search backward from refBug until you find a row whose Xcentroid is larger than
        RADIUS_IN_PIXELS to the left (along X axis) -- this is startRow in the neighborhood search.
    Now search forward through sorted list of bugs until you find row with Xcentroid
        that's over RADIUS_IN_PIXELS to right (on X axis) -- this is endRow.
    As you search forward from startRow, see if the location of each possible memberBug
        (candBug) is within RADIUS_IN_PIXELS of refBug (using a pythagorean L2 distance).
    If candBug IS within RADIUS_IN_PIXELS distance of refBug, then increment the value at
        counts[refBugID,candBugID], and add this candBug to a list of bugIDs in the
        current neighborhood surrounding refBug.
    Then we proceed to the next row down the matrix until we reach the right end of the image
        (the end of the sorted file).
    Now write the list of candBugIDs to the file:
        STUDYPREFIX/nhood_lists_radius_5/bug_X_nhood_lists_rep_Y_fov_Z.csv in a single row
        (one neighborhood), where X is refBugID's ID number.

Proceed to the next row (the next refBug) in the input file.  Repeat the preceding
    steps to find the memberBugs within RADIUS_IN_PIXELS of this refBug.

When finished with all bugAs in the input file, process the next CSV file in this study
    Once all the study's replicates and fovs are processed and added to the 500x500
    count matrix, write the matrix of counts to a csv file.  An example filename is:
    results/STUDYPREFIX/radius_5/nhood_counts_radius_5_500.csv
The each row is written to the nhoods files after the search around each refBug is
    complete.  These files are appended to each time a new refBug is processed.

ChangeLog:
April 7
- changed fileCount to imageCount everywhere
- changed writerows(list(zip(*array))) to writerows(array) for nhood counts

"""

# import matplotlib.pyplot as plt
import pandas
import math
import os
import numpy
import csv
import glob
import shutil


# load constants from a file...
exec(open('./bug_config.py').read(), globals())

# PROCESS_ONE_STUDY, WRITE_ALL_FILES, WRITE_NHOOD_LISTS
# DEBUG, CHECK_DATA, USE_PANDAS
# RADIUS_IN_MICRONS
# PIXEL_SIZE_IN_NM
# MAX_BUG_ID
# INPUT_PATH, RESULTS_PATH, SUB_FOLDER
# IDCOL, XCOL, YCOL

# -------------- functions --------------

def write_study_nhood_count_matrix(lastStudyPrefix):
    # write the study-level NxN matrix of counts for all neighborhood bugs on each reference bug's row
    if DEBUG:
        print('\nWriting 500x500 study count matrix for ', lastStudyPrefix)

    studyPath = outPath + os.sep + lastStudyPrefix  # for study
    radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)
    outName = radPath + os.sep + 'nhood_counts_500.csv'
    # example filename:
    #    results / [round02 /] 3454_COLON / radius_5 / nhood_counts_500.csv (without spaces)

    # if results/round02/STUDYPREFIX folder does not exist, create it
    if not os.path.exists(studyPath):
        print('Error: ', studyPath, ' missing in write_study_nhood_count_matrix')
        quit()

    if USE_PANDAS:
        outdf = pandas.DataFrame(countMatrix)  # no column or row labels, already in matrix
        outdf.to_csv(outName)
    else:    # write results to csv, not dataframe
        with open(outName, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # writer.writerows(list(zip(*studyCountMatrix))) # changed April 7
            writer.writerows(studyCountMatrix)

    # if we finished the only study we want to process, quit
    if PROCESS_ONE_STUDY != '':
        print('Processed only a single study ', PROCESS_ONE_STUDY, ', Quitting')
        quit()


def write_image_nhood_count_matrix(lastStudyPrefix, replicate, fov):
    # write the file-level NxN matrix of counts for all neighborhood bugs on each reference bug's row
    if DEBUG:
        print('\nWriting 500x500 image count matrix for ', lastStudyPrefix, ' rep=', replicate, ' fov=', fov)

    studyPath = outPath + os.sep + lastStudyPrefix  # for study
    radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)
    outName =  radPath + os.sep + 'nhood_counts_rep_' + replicate + '_fov_' + fov + '_500.csv'
    # example filename:
    #    results / [round02 /] 3454_COLON / radius_5 / nhood_counts_rep_3_fov_2_500.csv (without spaces)

    # if results/round02/STUDYPREFIX folder does not exist, create it
    if not os.path.exists(nlPath):
        print('Error: ', nlPath, ' missing in write_image_nhood_count_matrix')
        quit()

    if USE_PANDAS:
        outdf = pandas.DataFrame(imageCountMatrix)  # no column or row labels, already in matrix
        outdf.to_csv(outName)
    else:    # write results to csv, not dataframe
        with open(outName, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # writer.writerows(list(zip(*imageCountMatrix))) # changed April 7
            writer.writerows(imageCountMatrix)

def append_to_study_nhood_list_file(refID, studyNhoodBugIDlist):
    # write the list of bugIDs in this refID bug's neighborhood to its study nhood list file.
    # append a comma delimited list of memberBugIDs that make up all of refBugID's nhoods (1 nhood list on one row)
    # results / [round02 /] 3454_COLON / radius_5 / nhood_lists (without spaces)

    if DEBUG:
        print('\nAppending to study nhood list file for ', str(refID) + ' - ', lastStudyPrefix)

    studyPath = outPath + os.sep + studyPrefix  # for study
    radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)   # for each bug's neighborhood lists
    nlPath = radPath + os.sep + 'nhood_lists'  # for each bug's neighborhood lists

    if DEBUG == 2:
        print('Writing ', nlPath, ' nhood bug list, len=', len(studyNhoodBugIDlist))

    # create list starting with refBug ID, replicate, fov
    outlist = [refID, int(replicate), int(fov)]
    for id in studyNhoodBugIDlist:
        outlist.append(id)

    if DEBUG == 2:
        print( outlist)

    nhoodName = nlPath + os.sep + 'bug_' + str(refID) + '.csv'
    # e.g. results / [round02 /] 3454_COLON / radius_5 / nhood_counts_500.csv (without spaces)
    with open(nhoodName, 'a+', newline='') as csvfile:
        csvWriter = csv.writer(csvfile)
        csvWriter.writerow(outlist) # write a 2xN array (2 rows)

    # instead should this be...
    # for row in outList:
    #     csvWriter.writerow(outList)

def append_to_file_nhood_list_file(refID, fileNhoodBugIDlist):
    # write the list of bugIDs in this file's neighborhoods to refbug's nhood list file
    # append a comma delimited list of memberBugIDs that make up all of refBugID's nhoods (1 nhood list on one row)
    # results / [round02 /] 3454_COLON / radius_5 / nhood_lists (without spaces)

    if DEBUG:
        print('\nAppending to nhood list file for ', str(refID), ' - ', lastStudyPrefix)

    studyPath = outPath + os.sep + studyPrefix  # for study
    radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)
    nlPath = radPath + os.sep + 'nhood_lists'  # for each bug's neighborhood lists

    if DEBUG == 2:
        print('Writing ', nlpath, ' nhood bug list, len=', len(fileNhoodBugIDlist))

    # start the list with refBugID, replicate, fov, then the list of bugs in the nhood
    outlist = [refID, int(replicate), int(fov)]
    for id in fileNhoodBugIDlist:
        outlist.append(id)

    if CHECK_DATA:
        print( outlist)

    nhoodName = nlPath + os.sep + 'bug_' + str(refID) + '_rep_' + replicate + '_fov_' + fov + '.csv'
    # results / [round02 /] 3454_COLON / radius_5 / nhood_counts_rep_3_fov_2_500.csv (without spaces)
    with open(nhoodName, 'a+', newline='') as csvfile:
        csvWriter = csv.writer(csvfile)
        csvWriter.writerow(outlist)  # write a 2xN array (two rows)

    # instead should this be...
    # for row in outList:
    #     csvWriter.writerow(outList)

def delete_nhood_lists(nlPath):
    # remove all nhood list CSVs for this study and current radius (all file and study level files)
    #studyPath = outPath + os.sep + studyPrefix  # for study
    #radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)
    #nlPath = radPath + os.sep + 'nhood_lists'  # for each bug's neighborhood lists
    # results / [round02 /] 3454_COLON / radius_5 / nhood_lists (without spaces)

    # if folder does not exist, stop
    if os.path.exists(nlPath) == False:
        print('NOT deleting nhood list files for ', nlPath)
        return

    if 1 or DEBUG:
        print('Deleting nhood list files for ', nlPath)

    # remove study 'radius_X' folder and all contents
    shutil.rmtree(nlPath)


#--------------- end functions --------------

# if __name__ == '__main__':

# so the number of pixels in the radius is: 5000 nm / (70 nm/pixels), or 71.4286 pixels
RADIUS_IN_PIXELS = RADIUS_IN_MICRONS * 1000 / PIXEL_SIZE_IN_NM

# if an output SUB_FOLDER 'round02' is needed (to keep timepoints or datasets apart)
if SUB_FOLDER != '':
    outPath = RESULTS_PATH + os.sep + SUB_FOLDER  # optional folder for timepoints or rounds
else:
    outPath = RESULTS_PATH

if PROCESS_ONE_STUDY != '':
    print('Processing one study: ', PROCESS_ONE_STUDY)
else:
    print( 'Processing all studies')

if DEBUG:
    print('Making output folders')

# round2 input CSV file format is:
# Cell_Label, Cell_Size, Centroid_X, Centroid_Y, Intensity, ASV_Cluster_ID, ASV_Cluster_Scientific_Name
# int[1-4883], int[1-1079], float[0-1998.2], float[0-1998.6], float [0-1], int[1-439],string

# round1 format is:
# Cell_Label,Cell_Size,Centroid_X,Centroid_Y,Intensity,NCBI_Taxonomy_ID,Species_Scientific_Name

# read input file names from folder or from file?
if READ_STUDY_FILENAMES_FROM_FILE:
    # Get the list of Kanvas input files from a file, not via a directory of the folder.
    # The default input file is 'study_name_list.txt' in the current directory.
    #
    # Each line should include the file's full path, like:
    # ~/../data-2022-nov07-round02/csv/3415_COLON_tissue_replicate_1_fov_4_cell_information.csv
    #
    # STUDY_FILE_LIST_FILENAME shoud be present in the current folder.
    # The file will be created for you if you set READ_STUDY_FILENAMES_FROM_FILE to 1 then run this code again.

    if STUDY_FILE_LIST_FILENAME == '':
        print('Warning: STUDY_FILE_LIST_FILENAME is null.  I assume it is study_filename_list.txt and will proceed.')
        STUDY_FILE_LIST_FILENAME = 'study_filename_list.txt'

    if os.path.exists(STUDY_FILE_LIST_FILENAME) == False:
        print( '\nWarning: READ_STUDY_FILENAMES_FROM_FILE is 1, but ' + STUDY_FILE_LIST_FILENAME + ' not found.')
        print( 'So I am creating it.')

        # create 'study_file_list.txt' in the current directory
        ls_command = 'ls -1 ' + INPUT_PATH + os.sep + 'MER*.csv > ' + STUDY_FILE_LIST_FILENAME
        os.system(ls_command)

        print( '\nRun this code again to use this list file for the input filenames.')
        print( 'You can disable use of list by setting READ_STUDY_FILENAMES_FROM_FILE to 0 in bug_config.py.')
        print( 'Specify the name of the list file STUDY_FILE_LIST_FILENAME in bug_config.txt.')
        #print( '\nTo generate another list file, change STUDY_FILE_LIST_FILENAME in bugs_config.py,')
        #print( '    and re-run this code to create that file.')
        quit()

    with open(STUDY_FILE_LIST_FILENAME) as file:
        origNameList = [line.rstrip() for line in file]
else:
    # match all files that end in csv
    origNameList = glob.glob(INPUT_PATH + os.sep + '*.csv')

numFilesFound = len(origNameList)
if numFilesFound < 5:
    print('ERROR: Only ', numFilesFound, ' study input files found with wildcard')
    quit()

# sort names in nameList so each study is processed in full before next begun
fileList = sorted(origNameList)

if DEBUG == 2:
    for nm in fileList:
        print(nm)

firstName = fileList[0]
path,studyName = os.path.split(firstName)
lastStudyPrefix = studyName[0:16] # 3402_COLON
studyActive = False

# studyPath = outPath + os.sep + lastStudyPrefix  # for study
# radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)
# nlPath = radPath + os.sep + 'nhood_lists'  # for each bug's neighborhood lists

# create imageCountMatrix and studyCountMatrix
if DEBUG:
    print('creating countMatrixes')

# create count matrixes for first file in list
if USE_PANDAS:
    studyCountArray = numpy.zeros((MAX_BUG_ID+1, MAX_BUG_ID+1), dtype=int)
    imageCountArray = numpy.zeros((MAX_BUG_ID+1, MAX_BUG_ID+1), dtype=int)

    studyCountMatrix = pandas.DataFrame(studyCountArray, columns=range(0,MAX_BUG_ID+1))
    imageCountMatrix = pandas.DataFrame(imageCountArray, columns=range(0,MAX_BUG_ID+1))
else:
    studyCountMatrix = numpy.zeros((MAX_BUG_ID+1, MAX_BUG_ID+1), dtype=int)  # was [] no dtype
    imageCountMatrix = numpy.zeros((MAX_BUG_ID+1, MAX_BUG_ID+1), dtype=int)  # was [] no dtype

# Store the possible BUG IDs into the 0th row and 0th column of the count matrix
studyCountMatrix[0,:] = range(0,MAX_BUG_ID+1)
studyCountMatrix[:,0] = range(0,MAX_BUG_ID+1)
imageCountMatrix[0,:]  = range(0,MAX_BUG_ID+1)
imageCountMatrix[:,0]  = range(0,MAX_BUG_ID+1)

# iterate through a sorted list of Kanvas input files, usually 12 files per study
for fullFilePath in fileList:
    # split full name into path and filename
    path,fileName = os.path.split(fullFilePath)

    # get prefix in filename (to find other files in that study)
    # typical name: '3402_CECUM_tissue_replicate_1_fov_1_cell_information_counts.csv'
    studyPrefix = fileName[0:16] # i.e. 3402_CECUM

    if DEBUG:
        print('\nLooking at ', fileName)
        print('studyPrefix ', studyPrefix)

    # extract replicate # and fov # from filename
    replicate = fileName[34]
    fov       = fileName[40]

    if DEBUG:
        print( '\nstart studyPrefix=', studyPrefix, ' replicate=', replicate, ' fov=', fov)
        print('    lastStudyPrefix=', lastStudyPrefix)
        print('    studyActive=', studyActive)

    # did the name of the study just change?
    #     If so, and we're NOT processing the first file in studyNameList, write results.
    if studyPrefix != lastStudyPrefix and studyActive == True:
        write_study_nhood_count_matrix(lastStudyPrefix)
        lastStudyPrefix = studyPrefix
        studyActive = False

    # if not actively processing a study and name does not match, skip this file
    if PROCESS_ONE_STUDY != '' and studyPrefix != PROCESS_ONE_STUDY:
        print('Skipping ', studyPrefix, ' fov=', fov, ' replicate=', replicate)
        lastStudyPrefix = studyPrefix
        continue

    # if starting a new study, create output folders and init count vars and nhood lists
    if studyActive == False:
        # create output folders for study/onefile
        studyPath = outPath + os.sep + studyPrefix  # for study
        print('Creating dir ', studyPath)
        if os.path.exists(studyPath) == False:
            os.mkdir(studyPath)

        radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)
        print('Creating dir ', radPath)
        if os.path.exists(radPath) == False:
            os.mkdir(radPath)

        # delete then create the nhood list folder even if WRITE_NHOOD_LISTS turned off
        #   I don't want old nhood lists from a prior run to remain.
        nlPath = radPath + os.sep + 'nhood_lists'
        print( 'About to delete ', nlPath)

        # for a new study, clear all nhood lists
        delete_nhood_lists(nlPath)

        print('Creating dir ', nlPath)
        # recreate empty dir
        os.mkdir(nlPath)

        # zero out countMatrix for next study
        if DEBUG:
            print('zeroing countMatrix')

        # zero out studyCountMatrix for new study
        if USE_PANDAS:
            studyCountMatrixArray[:] = 0
            studyCountMatrix = pandas.DataFrame(studyCountMatrixArray, columns=range(0,MAX_BUG_ID+1))
        else:
            studyCountMatrix[:] = 0
            # put the possible BUG IDs into the 0th row and 0th column of the count matrix
            studyCountMatrix[0,:] = range(0,MAX_BUG_ID+1)  # store bugIDs 0-500 in 501 columns
            studyCountMatrix[:,0] = range(0,MAX_BUG_ID+1)

    if DEBUG:
        print('\nProcessing matching study ', studyPrefix, ' fov=', fov, ' replicate=', replicate)

    # we are starting to process a new file, so zero out file's countMatrix
    imageCountMatrix[:] = 0
    imageCountMatrix[0,:] = range(0,MAX_BUG_ID+1)  # store bugIDs 0-500 in 501 columns
    imageCountMatrix[:,0] = range(0,MAX_BUG_ID+1)
    lastStudyPrefix = studyPrefix
    studyActive = True

    # first row in matrix is the bugID, as is the first column.
    # must add these to non-pandas 2D matrix.
    # outmatrix is always 501x501 to hold all possible bug IDs (in columns 1 to 500)

    # print('count matrix initted')
    # print('shape of matrix = ', countMatrix.shape)

    # read next Kanvas input file into a pandas dataframe
    if DEBUG:
        print( 'Reading Kanvas CSV ', fullFilePath)
    kanvasDF = pandas.read_csv(fullFilePath)

    # sort the dataframe on Centroid_X (col 3, py index 2) so rows are ordered from low to high along X coord
    #     this simplifies a neighborhood search, since I can exclude all bugs w/ X coord too low or high # TODO:
    #     fall within RADIUS_IN_PIXELS of the referenceBug at the center of the neighborhood.
    sortedDF = kanvasDF.sort_values(by=['Centroid_X'], ascending=True, inplace=False)

    numRows, numCols = sortedDF.shape  # 7 columns, 1013 rows are typical

    endstr = str(numRows)  # how many rows of input to process?

    # step through each row in the Kanvas dataframe, from low to high X value, (left to right in image)
    #    get refBugID from the refRow we're processing
    #    search the nhood around refBug's X value, first to the left of refBug along X, then to the right
    for refRow in range(0, numRows):   # skips the first label row
        if refRow % 1000 == 1 and DEBUG < 2:  # print this only once per input file
            print('Processing ', str(refRow), ' of ', endstr, ' ', studyPrefix, ' fov=', fov, ' replicate=', replicate)

        refID = sortedDF.iloc[refRow, IDCOL]
        if refID == -1:  # skip invalid bug IDs
            continue

        refX = sortedDF.iloc[refRow, XCOL]
        refY = sortedDF.iloc[refRow, YCOL]

        # print first 40 rows of kanvas file to see what bugs should be in nhood
        if DEBUG > 1:
            guideDF = sortedDF.iloc[0:40, :]
            print(guideDF)

        if DEBUG > 1:
            print('refX=', refX, ' refY=', refY, ' refID=', refID)

        if False:
            r = 0
            refBugID = sortedDF.iloc[refRow, IDCOL]  # reference bug ID
            refX     = sortedDF.iloc[refRow, XCOL]  # X coord of refBug
            refY     = sortedDF.iloc[refRow, YCOL]  # Y coord of refBug

            print('refBugID=', refBugID, ' refX=', refX, ' refY=', refY)

            r = 1
            refBugID = sortedDF.iloc[refRow, IDCOL]  # reference bug ID
            refX     = sortedDF.iloc[refRow, XCOL]  # X coord of refBug
            refY     = sortedDF.iloc[refRow, YCOL]  # Y coord of refBug

            print('refBugID=', refBugID, ' refX=', refX, ' refY=', refY)

            r = 2
            refBugID = sortedDF.iloc[refRow, IDCOL]  # reference bug ID
            refX     = sortedDF.iloc[refRow, XCOL]  # X coord of refBug
            refY     = sortedDF.iloc[refRow, YCOL]  # Y coord of refBug

            print('refBugID=', refBugID, ' refX=', refX, ' refY=', refY)

            r = 3
            refBugID = sortedDF.iloc[r, IDCOL]  # reference bug ID
            refX     = sortedDF.iloc[r, XCOL]      # X coord of refBug
            refY     = sortedDF.iloc[r, YCOL]      # Y coord of refBug

            print('refBugID=', refBugID, ' refX=', refX, ' refY=', refY)
            quit()

        # list of bugIDs in this refBug's nhood (used for both the file and study nhood lists)
        nhoodBugIDlist = []

        # Do NOT add this refBug to its own nhoodBugIDlist, but do add other bugs with the same ID
        #     to the list when they are nearby
        # AND increment refBug in the 2D matrix -- at matrix[refID, refID] -- even for a nhood of 1

        # ---- PROCESS THIS refBugID bug's NHOOD ----
        # search BACKWARD from refBug's row through rows until cand.Centroid_X + RADIUS_IN_PIXELS < ref.Centroid_X
        # calc the cand bug's L2 distance from ref bug
        # if distance is less than RADIUS_IN_PIXELS, increment value in output matrix[refBugID, candBugID]
        #     and append candID to nhoodList
        row = refRow - 1
        if row < 0:
            row = 0

        if DEBUG == 2:
            print('Entering decr loop, row=', row)

        while row >= 0 and abs( refX - sortedDF.iloc[row, XCOL]) <= RADIUS_IN_PIXELS:
            diffx = refX - sortedDF.iloc[row, XCOL]
            diffy = refY - sortedDF.iloc[row, YCOL]

            if DEBUG == 2:
                print( 'decr row=', row, ' candX=', sortedDF.iloc[row, XCOL], 'candY=', sortedDF.iloc[row, YCOL])

            L2radius = math.sqrt( diffx * diffx + diffy * diffy)
            if L2radius <= RADIUS_IN_PIXELS:
                candID = sortedDF.iloc[row, IDCOL]
                if candID > 0:
                    if DEBUG == 2:
                        print('added to decr list candID=', candID, ' rad=', L2radius, ' diffx=', diffx, ' diffy=', diffy)
                    imageCountMatrix[refID, candID]  += 1
                    studyCountMatrix[refID, candID] += 1

                    nhoodBugIDlist.append(candID)  # add candID to list of bugs in this nhood

            row -= 1  # set the next candID row to the PRIOR row

        if DEBUG:
            print( 'decr nhood bug list for ', refID)
            for bug in fileNhoodBugIDlist:
                print( bug, end=',')
            print('\n')

        # search FORWARD from refBug's row through the rows until cand.Centroid_X - RADIUS_IN_PIXELS > ref.Centroid_X
        # calc the candBug's L2 distance from ref bug
        # if distance is less than RADIUS_IN_PIXELS, increment value in output matrix at row=refBugID col=candBugID
        #     and write candID to nhoodList
        if DEBUG == 2:
            print('Entering incr loop, row=', refRow + 1)

        row = refRow + 1
        while row < numRows and abs( refX - sortedDF.iloc[row, XCOL]) <= RADIUS_IN_PIXELS:
            diffx = refX - sortedDF.iloc[row, XCOL]
            diffy = refY - sortedDF.iloc[row, YCOL]

            if DEBUG == 2:
                print( 'incr row=', row, ' candX=', sortedDF.iloc[row, XCOL], 'candY=', sortedDF.iloc[row, YCOL])

            L2radius = math.sqrt( diffx * diffx + diffy * diffy)
            if L2radius <= RADIUS_IN_PIXELS:
                # if refBugID row does not exist in outmatrix, create it
                candID = sortedDF.iloc[row, IDCOL]
                if candID > 0:
                    if DEBUG == 2:
                        print('added to incr list candID=', candID, ' rad=', L2radius, ' diffx=', diffx, ' diffy=', diffy)
                    imageCountMatrix[refID, candID]  += 1
                    studyCountMatrix[refID, candID] += 1
                    nhoodBugIDlist.append(candID)  # add candID to list of bugs in this file nhood

            row += 1  # set the next candID row to the NEXT row

        if DEBUG:
            print( refRow, ' final bug nhood list for refBug', refID)
            for bug in nhoodBugIDlist:
                print( bug, end=',')
            print('\n')

        if DEBUG == 3:
            print('ENDED NHOOD LIST LOOPS, RADIUS_IN_PIXELS=', RADIUS_IN_PIXELS, '\n')

        # if nhood empty, add refID as only member, and increment its count on diagonal of matrix
        if len(nhoodBugIDlist) == 0:
            if DEBUG:
                print( 'Adding solitary refBug ', refID, ' to nhood list')

            # DO NOT add a solo bug nhood to the master count, but do add it to the nhood list file
            # countMatrix[refID, refID] += 1
            nhoodBugIDlist.append(refID)   # add refID to the list of bugs in this nhood

        if WRITE_NHOOD_LISTS:
            append_to_study_nhood_list_file(refID, nhoodBugIDlist)
            #if WRITE_ALL_FILES:  # MANUALLY_DISABLED
                # append_to_file_nhood_list_file(refID, nhoodBugIDlist)
                # if uncommented, this writes nhood_list files for all 12 of the image files
                # up to 480 files (assuming 40 bugs per image and 12 images)

        # now clear the nhood ID list and move on to the next refID row (refRow += 1)
        # nhoodBugIDlist = []

    # now that one file has been read, write out its count matrix
    if WRITE_ALL_FILES:
        write_image_nhood_count_matrix(lastStudyPrefix, replicate, fov)

# do I want to divide all bug count values by 2 before writing it out? since each A:B and B:A
#     relation is counted twice?
# NO, because I incremented count in different cells in matrix for A:B vs B:A.
# The counts are correct (and the same) in the upper or lower triangle of count matrix.

lastStudyPrefix = studyPrefix

# write results from the final study (since there was no subsequent study to triggger the write)
printf( 'Writing final count matrix for ', lastStudyPrefix)
write_study_nhood_count_matrix(lastStudyPrefix)

print('Done')
