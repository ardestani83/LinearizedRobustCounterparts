# Section 6
Instructions on how to produce the figures and tables from Section 6.

# Figure 1:
- Run performFig1.m

# Figures 2 and 3:
- Run the commands in performFigure23_all.m
(Files used: performFigure23_Tests.m, performFigure23_Analysis.m)

# Figure 4:
- Run the commands in performFigure4_all.m
(File used: performFigure4_Tests.m, performFigure4_Analysis.m)

# Table 9:
- Run the commands in performTable9_All.m
- Import the data from the ./results_table9/NewsVendorTable.csv file in Excel file NewsVendorTable.xlsx (in "Raw" sheet) to generate the table from the paper
(File used: performTable9Analysis.m, performTable9Tests.m)

# Details about folders:
./results_figure23: contains results generated to produce figures 2 and 3
./results_figure4: contains results generated to produce figure 4./results_table9: contains results generated to produce Table 9# Other files:
perform_allTests.m : performs all tests included in this folder

testNewsvendor5.m : Applies the different solutions scheme to multi-item newsvendor problem with budgeted uncertainty set
testNewsvendorEll3.m : Applies the different solutions scheme to multi-item newsvendor problem with uncertainty set that intersects box with ellipsoidvineBeta.m : Generates a random correlation matrix with random correlation with controlled spreadmergeMatResults.m : used in analysing results to merge the result files from the different runsplotPerformGivenX.m : used to plot figures 2, 3, and 4