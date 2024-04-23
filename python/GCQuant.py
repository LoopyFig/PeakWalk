#!/usr/bin/python3

# os operation libraries
import sys
import os
import getopt

from multiprocessing import Process, Queue, Lock

# data libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from sklearn.linear_model import LinearRegression

# get input options
try:
  opts, args = getopt.getopt(sys.argv[1:], "i:d:s:q:")
except getopt.GetoptError as err:
  print(err)
  sys.exit(2)
intensities = None
dilutions = None
stdlib = None
quantifications = None
for o, a in opts:
  if o == "-i":
    print(a)
    intensities = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-d":
    dilutions = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-s":
    stdlib = a
  if o == "-q":
    quantifications = os.path.abspath(a)


# set up transposed dataframe for looping regression
subjects = dilutions["sample"][(dilutions.type == "subject")].sort_values().tolist()
dilutions = dilutions[(dilutions["type"] == stdlib)].sort_values(["dilu", "sample"])
dilutionIntensities = intensities[dilutions["sample"].tolist()].transpose()
concentrations = dilutions.dilu.tolist()
sampleIntensities = intensities[subjects].transpose()
sampleIntensities.columns = dilutionIntensities.columns

# added to handle lists where concentrations are handled
# on a target by target basis
concentrationModifiers = intensities["concentration"].tolist()

dsamp = list()
dconc = list()
beta = list()
rsquared = list()
ii = 0
for i in dilutionIntensities:
  fragment = dilutionIntensities[[i]]
  fragment["concentration"] = concentrations
  fragment = fragment[(fragment[i] > 0)]
  dsamp.append(len(fragment.index))
  dconc.append(len(fragment.concentration.unique()))
  if len(fragment.index) > 0:
    x = np.array(fragment[i]).reshape(-1, 1)
    y = np.array(fragment.concentration).reshape(-1, 1)
    curve = LinearRegression(fit_intercept = False).fit(x, y)
    curveBeta = curve.coef_.item(0)
    beta.append(curveBeta)
    rsquared.append(curve.score(x, y))

    for j in sampleIntensities.index:
      sampleIntensities.at[j, i] = curveBeta * sampleIntensities.at[j, i] * concentrationModifiers[ii]

  else:
    beta.append(np.nan)
    rsquared.append(np.nan)

    for j in sampleIntensities.index:
      sampleIntensities.at[j, i] = np.nan
  ii = ii + 1

sampleIntensities = sampleIntensities.transpose()
quant = intensities.iloc[:,1:14]
quant["dsamp"] = dsamp
quant["dconc"] = dconc
quant["beta"] = beta
quant["rsquared"] = rsquared

quant = pd.concat([quant, sampleIntensities], axis = 1)
quant.to_csv(quantifications)


