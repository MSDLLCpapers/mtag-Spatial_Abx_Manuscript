"""
plot_study_heatmap_3.py
    The third step in a 3 script process.

This script generates a heatmap ONLY for an entire study, not for each of the 12 image files in a study.
Use plot_image_heatmap_3.py to get heatmaps for all 12 image files.

THIS SCRIPT DOES NOT REPLACE 'plot_heatmap_3.py'.  It cannot generate heatmap for an entire study,
    only for each of the 12 files in a study.
Run plot_heatmap_3.py to get a heatmap for the entire study.

1) read in 2D nhood count matrix from script step 2 into count_matrix
2) read in the bug counts file from script step 1 into bug_count list (2xN dimension, not really a list)
3) remove all zero columns and rows in the nhood bug count matrix (those bug IDs that do not exist)
   name those much smaller matrixes as VARIABLE_No0
4) calc the bug affinity ratios (actual count / expected count)
5) create heatmap, write it to disk

totalNumBugs = sum up the total number of bugs in the file from bug_count array (sum the 2nd row)

for one row in the 2D matrix (refBug), process each col (menmberBug)
totalNumRefBugNhoodBugs    = in 2D matrix from script #2, sum all memberBug counts in the refBug's row
totalNumOfMemberBugInStudy = in study_counts file from script #1, get the count from the bug's row

and write into expected_count matrix.

Then calc the 'affinity ratio' for the bug and store in affinity matrix:
    affinity = actual count of the bug / expected count of the bug

Do this for all bugs.

Create a 50x50 heatmap from 50x50 matrix using Seaborn.  (50 is estimate.  It usu ranges from 30 to 50.)
Display the heatmap, and write it to a PDF and CSV.

"""

import os
import csv
import pandas
import math
import numpy
import seaborn as sns
import matplotlib.pyplot as plt
from copy import copy, deepcopy


# load constants from a file...
exec(open('./bug_config.py').read(), globals())

# PROCESS_ONE_STUDY, WRITE_ALL_FILES, WRITE_NHOOD_LISTS
# DEBUG, CHECK_DATA, USE_PANDAS
# RADIUS_IN_MICRONS
# PIXEL_SIZE_IN_NM
# MAX_BUG_ID
# INPUT_PATH, RESULTS_PATH, SUB_FOLDER
# IDCOL, XCOL, YCOL


def remove_zeros_from_matrix(matrix):
    # remove all columns from matrix that contain only zeros, then all rows that contain only zeros
    if DEBUG:
        print('Enter remove_zeros_from_matrix')

    # remove zero rows and zero cols from NxN matrix
    nr, nc = matrix.shape

    if DEBUG:
        print( 'remove_zeros_from_matrix: numRows=', nr, ', numcols=', nc)

    # initialize a matrix with all row data from the 1st column (the row bug IDs)
    noZeroCols = matrix[:,0]
    noZeroCols = noZeroCols.reshape(nr, 1)  # transpose row into column

    if CHECK_DATA:
        print( 'shape of 1st col in noZeroCols: ', noZeroCols.shape)

    # step through matrix, look for zero columns, append col if not all zeros
    colCount = 1
    for c in range(1,nc):
        sumOfCol = numpy.sum(matrix[1:nr, c])  # sum of all rows except 1st

        # if column is NOT all zeros, append column to growing 2D matrix (grow wider)
        if sumOfCol > 0:
            acol = matrix[:,c]  # extract the column, but becomes a row (a list).
            acol = acol.reshape(nr,1)
            noZeroCols = numpy.hstack((noZeroCols, acol))  # append the transposed row (column)
            if CHECK_DATA:
                print( sumOfCol, ' shape of moving noZeroCols matrix: ', noZeroCols.shape)
            colCount += 1

    if CHECK_DATA:
        print( 'saved ', colCount, ' nonzero columns')
        print( 'cut   ', nc - colCount + 1, ' zero columns')
        print( 'shape of noZeroCol matrix: ', noZeroCols.shape)

    numCols = colCount  # get a new value for numCols

    # initialize a matrix with all col data from the NEW 1st row (the col bug IDs)
    noZeroRows = noZeroCols[0, 0:colCount]  # this is a 1D list of bug IDs

    if DEBUG:
        print( 'shape of 1st row in noZeroRows0: ', noZeroRows.shape)

    noZeroRows = noZeroRows.reshape(1, colCount)

    if DEBUG:
        print( 'shape of 1st row in noZeroRows1: ', noZeroRows.shape)

    # step through matrix, look for zero rows, append row if not all zeros
    rowCount = 1
    for r in range(1,nr):
        sumOfRow = numpy.sum(noZeroCols[r, 1:nc])  # sum of all cols except 1st

        # if row is NOT all zeros, append row to growing 2D matrix (grow taller)
        if sumOfRow > 0:
            arow = noZeroCols[r,:]
            if DEBUG:
                print( 'shape of arow: ', arow.shape)
            noZeroRows = numpy.vstack((noZeroRows, arow))  # append the row
            if DEBUG:
                print( sumOfRow, ' shape of noZeroRows matrix: ', noZeroRows.shape)

    if DEBUG:
        print( 'saved ', rowCount, ' nonzero rows')
        print( 'shape of final noZeroRows matrix: ', noZeroRows.shape)

    return noZeroRows


def create_heatmap(affinityRatios, radPath, studyPrefix):
    # show heatmap of (actual) count / predicted count for each nhood bug in a ref bug's row
    # red is higher then expected > 1, and blue is lower than expected < 1
    nr, nc = affinityRatios.shape

    if DEBUG:
        print( 'affinityRatios NumRows: ', nr, ' NumCols: ', nc)  # includes 0 row and col

    sns.set_theme()

    # if number of bugs > 20 need to make scale smaller
    sns.set(font_scale=0.4)  # 1 is default
    # Do not set a global font or the number of rows and cols cannot auto-adjust.
    # Some labels will disappear from the plot

    # find max ratio (exclude row 0 and col 0)
    maxratio = numpy.amax(affinityRatios[1:nr,1:nc])
    if DEBUG:
        print( "Max value in affinityRatios is ", maxratio)

    colNameFloats = affinityRatios[0, 1:nc]  # bug IDs are in row 0 of affinity matrix
    colNameInts = [int(x) for x in colNameFloats]
    colNames = ["%3d" % x for x in colNameInts]
    df = pandas.DataFrame(affinityRatios[1:nr, 1:nc], index=colNames, columns = colNames)

    if CHECK_DATA:
        for c in range(0, len(colNames)):
            print(c+1, colNames[c])
        print('\n')

    # seaborn.color_palette("coolwarm", as_cmap=True)  # blue to white to red
    ax = sns.heatmap(df, vmin=0, vmax=maxratio, square=True, cmap='coolwarm')

    # ax.set(xlabel='Centroid ID', ylabel='Neighbor ID', title=studyPrefix)
    heatmapTitle = studyPrefix + '_radius_' + str(RADIUS_IN_MICRONS)
    ax.set_title(heatmapTitle, fontsize=12)
    ax.set_xlabel('Neighbor Bug ID', fontsize=10)
    ax.set_ylabel('Center Bug ID', fontsize=10)
    ax.tick_params(axis='x', rotation= 90)
    ax.tick_params(axis='y', rotation=0)

    # write the heatmaps to the main folder: RESULTS_DIR / STUDYNAME
    outName = radPath + os.sep + 'heatmap'

    # write plot to file, png and pdf
    plt.savefig(outName + '.png')
    plt.savefig(outName + '.pdf')


def find_nans(matrix):  # not used
    print('looking for nans in 2D matrix')
    for r in range(0,matrix.shape[0]):
        for c in range(0,matrix.shape[1]):
            if numpy.isnan(matrix[r,c]):
                print( 'nan found at row:', r, ' col:', c)


def find_infs(matrix):  # not used
    print('looking for infs in 2D matrix')
    for r in range(0,matrix.shape[0]):
        for c in range(0,matrix.shape[1]):
            if numpy.isinf(matrix[r,c]):
                print( 'inf found at row:', r, ' col:', c)


#--------------- end functions --------------

# if __name__ == '__main__':

if PROCESS_ONE_STUDY == '':
    print('ERROR: You need to set PROCESS_ONE_STUDY in bug_config.py.')
    quit()

# this code processes only one study at a time - should specify the study name in bug_config.py
studyPrefix = PROCESS_ONE_STUDY  # specified in bug_config.py

RADIUS_IN_PIXELS = RADIUS_IN_MICRONS * 1000 / PIXEL_SIZE_IN_NM

if SUB_FOLDER != '':
    RESULTS_PATH = RESULTS_PATH + os.sep + SUB_FOLDER #MT replaced 'round02' with SUBFOLDER

# what paths do we read files from?
studyPath = RESULTS_PATH + os.sep + studyPrefix  # for study
radPath = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS)

# create heatmaps and write the matrix no0 CSVs for each (with zero rows/cols removed)
# the path to one file's neighborhood count matrix for each bug
countMatrixFilename = radPath + os.sep + 'nhood_counts_500.csv'
# example filename:
#    results / [round02 /] 
# 3454_COLON / radius_5_nhood_counts_500.csv (without spaces)

if not os.path.exists(countMatrixFilename):
    print( 'Nhood matrix file is missing ', countMatrixFilename)
    quit()

print( 'reading countMatrixFilename=', countMatrixFilename)

# read in 500x500 nhoodCounts matrix
#reader = csv.reader(open(countMatrixFilename, 'rb'), delimiter=',')
#x = list(reader)
#origNhoodCountMatrix = numpy.array(x).astype('float')

# This 500x500 nhood count matrix includes the counts for all input files
#     in the study, usually 12.  That summing across those files was done in
#     nhood_search_2.py.
nhoodCountMatrix500 = numpy.loadtxt(open(countMatrixFilename, 'rb'), delimiter=',')

#print( '    len: ', len(origNhoodCountMatrix))
if CHECK_DATA:
    print('    shape= ', origNhoodCountMatrix.shape)

    # make sure matrix col 0 and row 0 contain the bug IDs
    print('nhoodCountMatrix500:')
    print( '    cell 0,0 = ', nhoodCountMatrix500[0, 0])
    print( '    cell 0,1 = ', nhoodCountMatrix500[0, 1])
    print( '    cell 1,0 = ', nhoodCountMatrix500[1, 0])
    print( '    cell 1,1 = ', nhoodCountMatrix500[1, 1])
    print( '    cell 2,0 = ', nhoodCountMatrix500[2, 0])
    print( '    cell 2,1 = ', nhoodCountMatrix500[2, 1])

    print( 'reading bugCountFilename=', bugCountFilename)

# read in number of each bug found in study (500 IDs in list)
#reader = csv.reader(open(bugCountFilename, 'rb'), delimiter=',')
#x = list(reader)
#bugCount500 = numpy.array(x).astype('float')

# name of bug count file for one file
bugCountFilename = studyPath + os.sep + 'bug_counts_500.csv'

if not os.path.exists(bugCountFilename):
    print( 'Bug count file is missing ', bugCountFilename)
    quit()

if DEBUG:
    print('reading bugCountFilename= ', bugCountFilename)

# a 500x500 element matrix with many ZERO rows and cols
bugCounts500 = numpy.loadtxt(open(bugCountFilename, 'rb'), delimiter=',')

if DEBUG == 2:
    print('    shape= ', bugCounts500.shape)  # 500,2

    print( 'bugCounts500:')
    print( '    row 0 = ', bugCounts500[0,0], bugCounts500[0,1])
    print( '    row 1 = ', bugCounts500[1,0], bugCounts500[1,1])
    print( '    row 2 = ', bugCounts500[2,0], bugCounts500[2,1])
    print( '    row 3 = ', bugCounts500[3,0], bugCounts500[3,1])
    print( '    row 4 = ', bugCounts500[4,0], bugCounts500[4,1])
    print( '    row 5 = ', bugCounts500[5,0], bugCounts500[5,1])
    print( '    row 6 = ', bugCounts500[6,0], bugCounts500[6,1])
    print( '    row 7 = ', bugCounts500[7,0], bugCounts500[7,1])
    print( '    row 8 = ', bugCounts500[8,0], bugCounts500[8,1])
    print( '    row 9 = ', bugCounts500[9,0], bugCounts500[9,1])
    print( '    row 10 = ', bugCounts500[10,0], bugCounts500[10,1])

# remove all rows and columns in 2D nhood matrix that contain only zeros
nhoodCountMatrixNo0 = remove_zeros_from_matrix(nhoodCountMatrix500)

# what is size of trimmed matrix?  row 0 and col 0 are the bugIDs,
#    so matrix has 1 extra row and 1 extra column
# numRows and numCols set only once, done here
numRowsNo0, numColsNo0 = nhoodCountMatrixNo0.shape

# this is the number of bugs in the study (with a count of at least 1)
numBugs = numRowsNo0 - 1

if DEBUG:  # matrix's typical size: 43, 43
    print('Writing nhoodCountMatrixNo0, rows=', numRowsNo0, 'cols=', numColsNo0)

# write file to RESULTS_PATH / STUDYNAME (not into RESULTS_PATH / STUDYNAME / radius_5)
outName = radPath + os.sep + 'nhood_counts_no0.csv'
with open(outName, 'w', newline='') as csvfile:
    csvWriter = csv.writer(csvfile, delimiter=',')
    csvWriter.writerows(nhoodCountMatrixNo0)

if CHECK_DATA:
    print( 'shape of bugCounts500: ', bugCounts500.shape)  # 500,2

bugCountsNo0 = remove_zeros_from_matrix(bugCounts500)  # return a Nx2 matrix of ID,count

if CHECK_DATA:
    nr, nc = bugCountsNo0.shape
    print( 'trimmed list bugCountsNo0: numRows=', nr, ' numCols=', nc)  # 51,2

    # Now check number of bugs in the input files.  Number (and bug IDs) must match.
    nr1, nc1 = bugCountsNo0.shape

    print('bugCountsNo0 IDs: ', nr1)
    for i in range(1, nr1):
        print(i, bugCountsNo0[i,0], bugCountsNo0[i,1])
    print('\n')

    nr2, nc2 = nhoodCountMatrixNo0.shape

    print('nhoodCountMatrixNo0 IDs: ', nr2)
    for i in range(1, nr2):
        print(i, nhoodCountMatrixNo0[0, i])
    print('\n')

    # check to see if bugs in bug_counts_500.csv match the bugs in radius_5_nhood_counts_500.csv
    if nr1 != nr2:
        print('Error: number of bugs in these files do not match:')
    print( '    bug_counts_500.csv (', nr1, ') and radius_5_nhood_counts_500.csv (', nr2, ')')

# write out this file of bug_counts with 0 count bugs removed
# write file to RESULTS_PATH / STUDYNAME (not into RESULTS_PATH / STUDYNAME / radius_5)
if DEBUG:
    print('Writing bugCountsNo0')

outName = studyPath + os.sep + 'bug_counts_no0.csv'
with open(outName, 'w', newline='') as csvfile:
    csvWriter = csv.writer(csvfile, delimiter=',')
    #writer.writerows(list(zip(*bugCountsNo0)))
    csvWriter.writerows(bugCountsNo0)

# this is the total number of bugs in the file (bug counts are on the second row (1))
totalNumOfBugs = sum(bugCountsNo0[1:numRowsNo0,1])  # 5327

if CHECK_DATA:
    for i in range(1,numRowsNo0):
        print('Bug ', bugCountsNo0[i,0], ' count ', bugCountsNo0[i,1])
    print('\nbugCounts shape: ', bugCountsNo0.shape)
    print('total Number Of Bugs: ', totalNumOfBugs)

# convert actual count matrix into dataframe for seaborn
# countMatrixNo0df = pandas.DataFrame(trimCountMatrix)
if DEBUG:
    print( 'Creating empty arrays/lists')

# create empty 2D matrix of predicted counts (same memberBug IDs)
expectedNhoodCountMatrix = numpy.zeros([numRowsNo0, numRowsNo0])

# convert predicted count matrix into dataframe for seaborn
# expectedNhoodCountDF = pandas.DataFrame(predictedCounts)

# create empty 2D matrix of actual/predicted (same memberBug IDs in col 0 and row 0)
affinityRatios = numpy.zeros([numRowsNo0, numRowsNo0])

# calc the expected fraction of each nhood bug in the study based on
#     that bug's count / total num bugs.
#     Summing all these fractions for all bugs will equal 1.0
nhoodBugFraction = bugCountsNo0[:,1] / totalNumOfBugs

if CHECK_DATA:
    print( 'shape of nhoodBugFraction=', nhoodBugFraction.shape)  # 51, 1
    for n in range(1,nhoodBugFraction.shape[0]):
        if nhoodBugFraction[n] < 0.001:
            print( 'nhoodBugFraction is zero at: ', n)

# calc 2D matrix of the expected counts for each nhood bug in each ref bug's row
# first sum up the total num nhood bugs for each refbug
#     (sum of counts in a refBug's row minus the refBug's entry for itself)
refBugNhoodCountTotal = numpy.zeros(numRowsNo0)

if DEBUG:
    print( 'Filling refBugNhoodCountTotal, numrows=', numRowsNo0)

# sum up the count of all neighbor bugs (columns) for each refBug (row)
# subtract neighborhood of 1 (self) from count
for b in range(1,numRowsNo0):
    refBugNhoodCountTotal[b] = sum(nhoodCountMatrixNo0[b, 1:numRowsNo0])
        # - actualNhoodCountMatrix[b,b]
        # DO NOT subtract the refID's self, since a solo refID's was
        #   never added to the nhoods in nhood_search.py (line 481)

if CHECK_DATA:
    print( 'refBugNhoodCountTotal:')
    for n in range(1,numRowsNo0):
        print('bug ', n, ' countTotal ', refBugNhoodCountTotal[n])

# check there are no zeros in count matrix.  All refBug rows should contain some nhood bugs.
for n in range(1,refBugNhoodCountTotal.shape[0]):
    if refBugNhoodCountTotal[n] < 0.001:
        print( 'Error: refBugNhoodCount is zero for refBugID ', \
        nhoodCountMatrixNo0[n,0], ', row: ', n)

if DEBUG or CHECK_DATA:
    print( 'Shape of refBugNhoodCountTotal ', refBugNhoodCountTotal.shape)

# now create the expected nhood count matrix using the nhood bug bug fractions
#    and the total num of bugs in each refbug's nhoods
# this tells how many of each nhood bug *should* be in all of each refbug's nhoods
if DEBUG:
    print( 'Filling expectedNhoodCountMatrix ', studyPrefix)

# transfer bug IDs into row 0 and col 0 of expected count matrix
expectedNhoodCountMatrix[0,1:numColsNo0] = nhoodCountMatrixNo0[0,1:numColsNo0]
expectedNhoodCountMatrix[1:numColsNo0,0] = nhoodCountMatrixNo0[1:numColsNo0,0]

zeros = 0
for r in range(1,expectedNhoodCountMatrix.shape[0]):      # refbug row
    for n in range(1,expectedNhoodCountMatrix.shape[1]):  # nhood bug column
        expectedNhoodCountMatrix[r,n] = nhoodBugFraction[n] * refBugNhoodCountTotal[r]
        if expectedNhoodCountMatrix[r,n] == 0:
            zeros += 1
            print( zeros, ' Zero in expectedCount Matrix at: ', r, n)

if DEBUG:
    print( 'shape of expectedNhoodCountMatrix=', expectedNhoodCountMatrix.shape)

# write expected count array
outName = radPath + os.sep + 'expected_nhood_count_matrix_no0.csv'
with open(outName, 'w', newline='') as csvfile:
    csvWriter = csv.writer(csvfile)
    # csvWriter.writerows(list(zip(*expectedNhoodCountMatrix)))
    csvWriter.writerows(expectedNhoodCountMatrix)

# now we can create the affinity matrix showing how the actual counts compare to the expected

# copy bug IDs into row 0 and col 0 of the affinityRatio matrix
if DEBUG:
    print( 'Initting affinityRatios col 0 and row 0')

for i in range(1,numRowsNo0):
    affinityRatios[0,i] = bugCountsNo0[i,0]
    affinityRatios[i,0] = bugCountsNo0[i,0]

if DEBUG:
    print( 'Filling affinityRatios')

if DEBUG:
    print('shape of affinityRatios=', affinityRatios.shape) # 51, 51
    print('shape of nhoodCountMatrixNo0=', nhoodCountMatrixNo0.shape)  # 43, 43
    print('shape of expectedNhoodCountMatrix=', expectedNhoodCountMatrix.shape) # 51, 51

for r in range(1,affinityRatios.shape[0]):
    for c in range(1,affinityRatios.shape[1]):  # this is actual / expected
        affinityRatios[r,c] = nhoodCountMatrixNo0[r,c] / expectedNhoodCountMatrix[r,c]

# write affinityMatrix to a CSV
# studyPath = outPath + os.sep + studyPrefix  # path to study results folder

if DEBUG:
    print('Writing affinity matrix for ', studyPrefix)

outName = radPath + os.sep + 'affinity_ratio_matrix_no0.csv'
if USE_PANDAS:
    outdf = pandas.DataFrame(affinityRatios)  # no column or row labels, already in matrix
    outdf.to_csv(outName)
else:    # write results to csv, not dataframe
    with open(outName, 'w', newline='') as csvfile:
        csvWriter = csv.writer(csvfile, delimiter=',')
        # writer.writerows(list(zip(*affinityRatios)))
        csvWriter.writerows(affinityRatios)

# suppress large values in heatmap for visualization
#     create temp matrix, convert all large values to a max of 25
#     display truncated array but save original
maxRatio = 2
mask = affinityRatios > maxRatio
affinityRatiosMasked = deepcopy(affinityRatios)
affinityRatiosMasked[mask] = maxRatio

# assign bug IDs again to row 0 and col 0
affinityRatiosMasked[0,0:numRowsNo0] = affinityRatios[0,0:numRowsNo0]
affinityRatiosMasked[0:numRowsNo0,0] = affinityRatios[0:numRowsNo0,0]

# display the heatmap and write the image to a png and PDF
if DEBUG:
    print( 'Creating heatmap')

create_heatmap(affinityRatiosMasked, radPath, studyPrefix)

# write CSV of affinity ratios
# outName = studyPath + os.sep + 'radius_' + str(RADIUS_IN_MICRONS) + '_heatmap.csv'
# example filename:
#    results / [round02 /] 3454_COLON / radius_5_heatmap.csv (without spaces)

if CHECK_DATA:
    # look for nan in matrix
    nanmask = numpy.isnan(affinityRatios)
    infmask = numpy.isinf(affinityRatios)

    find_infs(affinityRatios)
    find_nans(affinityRatios)

    print( 'max in affinityRatios: ', numpy.nanmax(affinityRatios))

    # sort affinityRatios to see how many are > 100
    affList = affinityRatios[1:numRows, 1:numRows].flatten()
    saffList = sorted(affList, reverse=True)

    for c in range(0,100):
        print(saffList[c])

# sort affinityRatios to see top values, where they are
#     debug why those values are so high
# why aren't the 2D matrices symmetric?
#     row is ref bug, col is nhood bug -- not the same

print('Done')

# for each refBug (row) and memberBug (col) in (actual) nhoodCountMatrixNo0 2D matrix:
#    what is the predicted num of each memberBug in each cell?
#
# compute predictedNumOf MemberBug and write into predictedMatrix
#
# compute this one scalar val once for entire study:
# totalNumOfBugsInStudy = in _bug_counts_500 file from script #1, sum up all bug counts
#
# compute this val once for each refBug, write to 500 element vector
# totalNumOfRefBugNhoodBugs = in 2D matrix from script #2, sum all memberBug counts in the refBug's row
# write into totalNumOfRefBugNhoodBugs(1-500)
#
# compute this val once for each memberBug, write to 500 element vector
# totalNumOfBug = in _bug_counts_500 file from script #1, get count from the bug's row
# write into totalNumOfBug(1-500)
