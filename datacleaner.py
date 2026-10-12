import sys
import os
import traceback
import numpy as np
import pandas as pd

infile = pd.read_excel("Untargeted metabolomics_CSF.xlsx", skiprows=4)

notunknown = infile[infile["Metabolite name"] != "Unknown"]

notunknown["Metabolite name"] = notunknown["Metabolite name"].str.removeprefix("low score: ")
notunknown["Metabolite name"] = notunknown["Metabolite name"].str.removeprefix("no MS2: ")

notunknown["CV"] = notunknown["QC.1"] / notunknown["QC"]
nobadcv = notunknown[notunknown["CV"] <= 0.3]

noxtracols = nobadcv.drop(columns=["Alignment ID", "Average Rt(min)", "Adduct type", "Post curation result", "Fill %", "MS/MS assigned", "Reference RT", "Reference m/z", "Formula", "Ontology", "INCHIKEY", "SMILES", "Annotation tag (VS1.0)", "RT matched", "m/z matched", "MS/MS matched", "Comment", "Manually modified for quantification", "Manually modified for annotation", "Isotope tracking parent ID", "Isotope tracking weight number", "RT similarity", "m/z similarity", "Simple dot product", "Weighted dot product", "Reverse dot product", "Matched peaks count", "Matched peaks percentage", "Total score", "S/N average", "Spectrum reference file name", "MS1 isotopic spectrum", "MS/MS spectrum"])

colnames = noxtracols.columns.to_list()
startcol = 0
stopcol = 0


for col in colnames:
    if "Pool-QC_4" in col:
        targetcol = colnames.index(col) + 1

    elif col == "Blank":
        stopcol = colnames.index(col)

blanknorm = noxtracols.iloc[:, targetcol:stopcol].sub(noxtracols["Blank"], axis=0)

maskednorm = blanknorm.mask(blanknorm <= 0)

filteredmask = maskednorm.dropna(axis=0, thresh=len(maskednorm.columns)*.5)

hmv = filteredmask.min(axis=1) / 2 

minvalmask = filteredmask.T.fillna(hmv).T

lognorm = np.log2(minvalmask)

lognorm.insert(0, "Metabolite name", noxtracols.loc[lognorm.index, "Metabolite name"].values)
#lognorm.insert(1, "m.z", noxtracols.loc[lognorm.index, "Average Mz"].values) 

nodupes = lognorm.drop_duplicates(subset=["Metabolite name"])

groupslist =[" "]
nameslist = nodupes.columns.to_list()

for i in nameslist[1:]:
    groupslist.append(i.split('_')[-1])

newframe = pd.Series(groupslist, index=nodupes.columns).to_frame().T

finaldf = pd.concat([newframe, nodupes.iloc[0:, :]])

finaldf.to_csv("test.csv", index=False)

