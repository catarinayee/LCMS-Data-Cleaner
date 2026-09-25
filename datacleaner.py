import sys
import os
import traceback
import numpy as np
import pandas as pd

infile = pd.read_excel("Untargeted metabolomics_CSF.xlsx", skiprows=4)

notunknown = infile[infile["Metabolite name"] != "Unknown"]
notunknown["CV"] = notunknown["QC.1"] / notunknown["QC"]

nobadcv = notunknown[notunknown["CV"] <= 0.3]
noxtracols = nobadcv.drop(columns=["Alignment ID", "Average Rt(min)", "Adduct type", "Post curation result", "Fill %", "MS/MS assigned", "Reference RT", "Reference m/z", "Formula", "Ontology", "INCHIKEY", "SMILES", "Annotation tag (VS1.0)", "RT matched", "m/z matched", "MS/MS matched", "Comment", "Manually modified for quantification", "Manually modified for annotation", "Isotope tracking parent ID", "Isotope tracking weight number", "RT similarity", "m/z similarity", "Simple dot product", "Weighted dot product", "Reverse dot product", "Matched peaks count", "Matched peaks percentage", "Total score", "S/N average", "Spectrum reference file name", "MS1 isotopic spectrum", "MS/MS spectrum"])
notdupes = noxtracols.drop_duplicates(subset=["Metabolite name"])

colnames = notdupes.columns.to_list()
startcol = 0
stopcol = 0

for col in colnames:
    if "Pool-QC_4" in col:
        targetcol = colnames.index(col) + 1

    elif col == "Blank":
        stopcol = colnames.index(col) + 1

blanknorm = notdupes.iloc[:, targetcol:stopcol].sub(notdupes["Blank"], axis=0)

#lognorm = np.log2(blanknorm)

blanknorm.to_csv("test.csv", index=False)