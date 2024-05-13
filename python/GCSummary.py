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

# get input options
try:
  opts, args = getopt.getopt(sys.argv[1:], "i:r:m:s:")
except getopt.GetoptError as err:
  print(err)
  sys.exit(2)
intensities = None
rtimes = None
masscharges = None
for o, a in opts:
  if o == "-i":
    print(a)
    intensities = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-r":
    rtimes = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-m":
    masscharges = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-s":
    summaryname = os.path.abspath(a)

# create initial fragment summaries
summaries = intensities[["tmz", "id", "subid", "name", "monoisotopic", "cas", "formula"]]

# summarize rtimes
rtimes = rtimes.iloc[:,14:len(rtimes.columns)]
rtimes[rtimes == 0] = np.nan
summaries["rtMed"] = rtimes.median(axis=1)
summaries["rtMin"] = rtimes.min(axis=1)
summaries["rtMax"] = rtimes.max(axis=1)
summaries["rtRange"] = summaries["rtMax"].sub(summaries["rtMin"])

masscharges = masscharges.iloc[:,14:len(masscharges.columns)]
masscharges[masscharges == 0] = np.nan
summaries["mzMed"] = masscharges.median(axis=1)
summaries["mzMin"] = masscharges.min(axis=1)
summaries["mzMax"] = masscharges.max(axis=1)
summaries["mzRange"] = summaries["mzMax"].sub(summaries["mzMin"])

# summarize intensities
intensities = intensities.iloc[:,14:len(intensities.columns)]
intensities[intensities == 0] = np.nan
summaries["iMean"] = intensities.mean(axis=1)
summaries["iStd"] = intensities.std(axis=1)
summaries["iCV"] = summaries["iStd"].div(summaries["iMean"])
summaries["detected"] = intensities.count(axis=1)
summaries["detectedfraction"] = summaries["detected"].div(len(intensities.columns))

# identification quality
ids = summaries["id"].drop_duplicates()
summaries["fragCnt"] = 0
summaries["qualityCnt"] = 0
for i in ids:
  iiloc = np.where(summaries["id"] == i)[0].tolist()
  summaries.loc[summaries.index[iiloc], "fragCnt"] = len(iiloc)
  simulCnt = intensities.iloc[iiloc].count(axis=0)
  summaries.loc[summaries.index[iiloc], "qualityCnt"] = (simulCnt >= 3).sum()

# rearrange columns, add back intensities to summary
summaries = summaries[["mzMed", "rtMed", "tmz", "id", "subid", "name", "monoisotopic", "formula", "mzMin", "mzMax", "mzRange", "rtMin", "rtMax", "rtRange", "iMean", "iStd", "iCV", "detected", "detectedfraction", "fragCnt", "qualityCnt"]]
summaries = pd.concat([summaries, intensities], axis=1)

# filter rows
summaries = summaries[summaries["detected"] > 0]

# write summary file
summaries.to_csv(summaryname)

