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
  opts, args = getopt.getopt(sys.argv[1:], "i:r:s:")
except getopt.GetoptError as err:
  print(err)
  sys.exit(2)
intensities = None
rtimes = None
for o, a in opts:
  if o == "-i":
    print(a)
    intensities = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-r":
    rtimes = pd.read_csv(os.path.abspath(a), sep=",")
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

# summarize intensities
intensities = intensities.iloc[:,14:len(intensities.columns)]
intensities[intensities == 0] = np.nan
summaries["iMean"] = intensities.mean(axis=1)
summaries["iStd"] = intensities.std(axis=1)
summaries["iCV"] = summaries["iStd"].div(summaries["iMean"])
summaries["detected"] = intensities.count(axis=1)
summaries["detectedfraction"] = summaries["detected"].div(len(intensities.columns))

# rearrange columns, add back intensities to summary
summaries = summaries[["tmz", "rtMed", "id", "subid", "name", "monoisotopic", "formula", "rtMin", "rtMax", "rtRange", "iMean", "iStd", "iCV", "detected", "detectedfraction"]]
summaries = pd.concat([summaries, intensities], axis=1)

# filter rows
summaries = summaries[summaries["detected"] > 0]

# write summary file
summaries.to_csv(summaryname)

