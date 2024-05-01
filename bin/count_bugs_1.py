"""
count_bugs_1.py
    The first step in a 3 script process.  It writes the count of each bug for
    each file in a 12 file study as well as a combined file for th eentire study.
    The output is a 500x2 CSV file of bugIDs and counts.

This script can process an entire 12 files study, as well as each of the 12 files in that study,
    producing one file that sums all 12 files (bug_counts_500.csv) and 12 files in
    the format 'bug_counts_rep_2_fov_3_500.csv' (for rep 2 and fov 3).  It also produces
    a version of each that contains only the bugs present in that image (no 0 count bugs).
    All output files appear in the main STUDYNAME folder.

If PROCESS_ONE_STUDY is set to null (''), it processes all studies in the input folder.
    This variable is set in bug_config.py.

The sequence of steps:

1) Count the number of each bug ID in all 12 fov+replicate files, then sum the ID across
    all 12 files in the entire study.
2) Write the summed count for each bug into one CSV per study, named for that study,
    and write the count for each replicate+fov file.
    Thus the summed file for the study (bug_counts_500.csv) contains 500 rows (0-499), one row
        per bug, even if the bug count is 0.
    I also create a version of each of these 13 files that excludes the bugs with a count of 0.
        These are names the same, but instead of '_500', each filename ends in '_no0'.
    Thus each line of the 500 line output file is: bugID, count.

The output file's name is results/STUDY_NAME/bug_counts_500.csv.  If PROCESS_ONE_STUDY
is set to null (''), this script will process all the studies in the Kanvas data directory.

There is no downside to running this script repeatedly.  Existing output files will
    be overwritten.
"""

import os
import numpy
import pandas
import csv
import glob


# load constants from a file
exec(open('./bug_config.py').read(), globals())

# PROCESS_ONE_STUDY, WRITE_ALL_FILES, WRITE_NHOOD_LISTS
# DEBUG, CHECK_DATA, USE_PANDAS
# RADIUS_IN_MICRONS
# PIXEL_SIZE_IN_NM
# MAX_BUG_ID
# INPUT_PATH, RESULTS_PATH, SUB_FOLDER
# IDCOL, XCOL, YCOL


def process_one_file():
    if DEBUG:
        print('Processing... ', studyName)
        # print("studyName=", studyName)

    # print( 'reading image ', fullName)

    # read in CSV file data
    csvDataframe = pandas.read_csv(fullName)

    # read through each csv, for each row, add 1 to that bugID
    for idx, row in csvDataframe.iterrows():
        bugID = row["ASV_Cluster_ID"]

        if DEBUG == 2:
            print("Bug ID is ", bugID)

        # MAX_BUG_ID is actually MAX_BUG_ID-1 -- count array is 0 to MAX_BUG_ID
        if bugID >= MAX_BUG_ID:
            print('ERROR: bugID ', bugID, ' exceeds ', MAX_BUG_ID-1)
            quit()

        if bugID > 0:  # unrecognized bugs have ID = -1, so skip it
            countArray[bugID] += 1
            totCountArray[bugID] += 1


def write_one_result():
    # write countArray for one file in the 12 file study.

    if DEBUG:
        print('\nWriting one file ', studyName)
        # a typical studyName is: 3402_CECUM_tissue_replicate_1_fov_1_cell_information.csv

    repStr = studyName[34]
    fovStr = studyName[40]

    # create dataframe from four arrays
    idxs = numpy.array(range(0,len(countArray)))

    numBugs = sum(countArray)
    print( 'Num bugs ', numBugs, ' rep ', repStr, ' fov ', fovStr)  # same as number of rows

    outArray = numpy.column_stack((idxs, countArray))

    # now array is a 500x2 matrix of bug IDs and counts

    if USE_PANDAS:
        odf = pandas.DataFrame(outArray, columns=['bugID', 'count'], dtype=int)

    studyPath = outPath + os.sep + studyPrefix
    if os.path.exists(studyPath) == False:
        os.mkdir(studyPath)

    outName = studyPath + os.sep + 'bug_counts_rep_' + repStr + '_fov_' + fovStr + '_500.csv'
    # like 3454_COLON/bug_counts_rep_2_fov_3_500.csv

    if USE_PANDAS:
        odf.to_csv(outName)
    else:    # write results to csv, not dataframe
        with open(outName, 'w', newline='') as csvfile:
            csvWriter = csv.writer(csvfile)
            for idx in range(len(countArray)):
                csvWriter.writerow([idx, int(countArray[idx]) ])

    # # now write countArray without the zero columns and rows (only the bugs that are present)
    # outName = studyPath + os.sep + 'bug_counts_rep_' + repStr + '_fov_' + fovStr + '_no0.csv'
    #
    # if USE_PANDAS:
    #     odf.to_csv(outName)
    # else:    # write results to csv, not dataframe
    #     with open(outName, 'w', newline='') as csvfile:
    #         writer = csv.writer(csvfile)
    #         for idx in range(len(countArray)):
    #             if countArray[idx] >= 1:
    #                 writer.writerow([idx, int(countArray[idx]) ])
    #
    # # for nonzero counts, print bugID and count
    # if DEBUG:
    #     for idx in range(len(countArray)):
    #         if countArray[idx] >= 1:
    #             print( idx, int(countArray[idx]))
    #     print('\n')


def write_study_results():
    # write totCountArray for entire study (all 12 files)
    print('Writing study ', lastStudyPrefix) # , ' num files processed ', numFilesProcessed)

    # create dataframe from four arrays
    idxs = numpy.array(range(0,len(totCountArray)))  # 0 to 500

    numBugs = sum(totCountArray)
    print( 'Num bugs in study ', int(numBugs), '\n')  # same as number of rows

    # insert a row 0 to include all bug IDs from 0 to 500
    outArray = numpy.column_stack((idxs, totCountArray))

    # now array is a 500x2 matrix of bug counts

    if USE_PANDAS:
        odf = pandas.DataFrame(outArray, columns=['bugID', 'count'], dtype=int)

    studyPath = outPath + os.sep + lastStudyPrefix
    if os.path.exists(studyPath) == False:
        os.mkdir(studyPath)

    # write 500x2 count Array
    outName = studyPath + os.sep + 'bug_counts_500.csv' #  like 3454_COLON.csv

    header_names = ['bug ID', 'count']  # I don't use these

    if USE_PANDAS:
        odf.to_csv(outName)
    else:    # write results to csv, not dataframe
        with open(outName, 'w', newline='') as csvfile:
            csvWriter = csv.writer(csvfile)
            # writer.writerow(header_names)  # 0th row is header names
            for idx in range(len(totCountArray)):
                csvWriter.writerow([idx, int(totCountArray[idx]) ])

    # for nonzero counts, print bugID and count
    if CHECK_DATA:
        for idx in range(len(totCountArray)):
            if totCountArray[idx] >= 1:
                print( idx, int(totCountArray[idx]))
        print('\n')

    # # now write totCountArray without the zero columns and rows (only the bugs that are present)
    # outName = studyPath + os.sep + 'bug_counts_no0.csv' #  like 3454_COLON.csv
    #
    # if USE_PANDAS:
    #     odf.to_csv(outName)
    # else:    # write results to csv, not dataframe
    #     with open(outName, 'w', newline='') as csvfile:
    #         writer = csv.writer(csvfile)
    #         #writer.writerow(header_names)  # 1st row is header names
    #         for idx in range(len(totCountArray)):
    #             if totCountArray[idx] >= 1:
    #                 writer.writerow([idx, int(totCountArray[idx]) ])

    # if we finished the only study we want to process, quit
    if PROCESS_ONE_STUDY != '':
        print( 'Processed only a single study: ', lastStudyPrefix, ', Quitting')
        quit()



#------------ no more functions ------------

# if __name__ == '__main__':

if SUB_FOLDER != '':
    outPath = RESULTS_PATH + os.sep + SUB_FOLDER  # name of sub-results folder
else:
    outPath = RESULTS_PATH

# each study's output file will be named: results/round02/3402_CECUM_bug_counts.csv

# if output folder does not exist, create it
if os.path.exists(RESULTS_PATH) == False:
    os.mkdir(RESULTS_PATH)

# round02
if os.path.exists(outPath) == False:
    os.mkdir(outPath)

# round2 input CSV file format is:
# Cell_Label, Cell_Size, Centroid_X, Centroid_Y, Intensity, ASV_Cluster_ID, ASV_Cluster_Scientific_Name
# int[1-4883], int[1-1079], float[0-1998.2], float[0-1998.6], float [0-1], int[1-439],string

# round1 format is:
# Cell_Label,Cell_Size,Centroid_X,Centroid_Y,Intensity,NCBI_Taxonomy_ID,Species_Scientific_Name

# get min and max of each column
# maxvals = csvData.max()
# print(maxvals)

# match all files that end in csv
nameList = glob.glob(INPUT_PATH + os.sep + '*.csv')

numFiles = len(nameList)
if numFiles < 5:
    print('ERROR: Only ', numFiles, ' files found')
    quit()

# sort names in nameList
nameList = sorted(nameList)
# print(nameList)

# create 1D array for count of each bugID (up to 500 bugIDs, so array is 0 to 501)
totCountArray = numpy.zeros([MAX_BUG_ID])  # total bug counts for all 12 files in study
countArray    = numpy.zeros([MAX_BUG_ID])  # bug counts for 1 file in study

# set name to first name in list
firstName = nameList[0]
path,studyName = os.path.split(firstName)
lastStudyPrefix = studyName[0:16] # e.g. 3402_CECUM
studyActive = False

for fullName in nameList:
    # split full name into path and filename
    path,studyName = os.path.split(fullName)
    # typical studyName is 3402_CECUM_tissue_replicate_1_fov_1_cell_information.csv

    # get prefix in filename (to find other files in that study)
    # typical name: '3402_CECUM_tissue_replicate_1_fov_1_cell_information_counts.csv'
    studyPrefix = studyName[0:16]  # i.e. 3402_CECUM

    if studyPrefix != lastStudyPrefix and studyActive == True:
        write_study_results()  # all replicate + fov files combined

        # zero out totCountArray to process next study
        totCountArray[:] = 0
        lastStudyPrefix = studyPrefix
        studyActive = False

    # if single file is specified, and name does not match, then skip file
    if PROCESS_ONE_STUDY != '' and studyPrefix != PROCESS_ONE_STUDY:
        print( 'Skipping ', studyName)
        continue

    process_one_file()  # one replicate + fov file

    lastStudyPrefix = studyPrefix

    write_one_result()  # one replicate + fov file

    studyActive = True
    countArray[:] = 0

# write results from the final study (since there no subsequent study file
#     was found in list to triggger the write).
write_study_results()  # all replicate + fov files combined

print('Done')
