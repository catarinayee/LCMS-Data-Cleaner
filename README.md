The LCMS Data Cleaner takes the CSV output of MS-DIAL, and cleans + normalizes + filters the data according to workflow specifications.

1) Drops all unknown features from the dataset
2) Cleans up the metabolite names so that they are parseable by MetaboAnalyst
3) Calculates the CV value and then drops all features that exceed the threshold
4) Takes away all columns that aren't the sample concentrations
5) Figures out how many samples there are, and then gets that slice of the dataset so that we can normalize it.
6) Subtracts the blanks' average from the sample concentration.
7) Drops all features that have too many missing values or didn't appear in over 50% of samples
8) Applies a half-minimum value to all of the remaining missing values so they normalize correctly
9) Normalizes the data by Log2
10) Re-attaches the Metabolite name to the corresponding row of the feature
11) This gets rid of duplicates; putting this at the end ensures no false negatives from the first instance of the feature being an outlier.
12) This scrapes the Sample Group from the end of the sample name after the last underscore, ie if your sample is named 345_L_V, it would be assigned the V group. It then adds it as a new row under the column headers.
13) Outputs a new CSV file ready to be put into the Biomarker Analysis tool
