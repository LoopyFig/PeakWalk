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
  opts, args = getopt.getopt(sys.argv[1:], "i:r:m:q:n:b:l:s:xh")
except getopt.GetoptError as err:
  print(err)
  sys.exit(2)
# required parameters
intensities = None
rtimes = None
masscharges = None
quants = None
summaryname = None
# additional parameters
batchinfo = None # batch file
labels = ["subject"] # sample type filter
stdlib = None # standards library
hasrep = False # toggles replicate summarization
for o, a in opts:
  if o == "-h":
    print("PeakWalk: Automated GC Identification and Quantification")
    print("-i intensity file, required")
    print("-r retention time file, required")
    print("-m masscharge file, required")
    print("-q quantification file, optional, overrides intensity with concentration, requires -b")
    print("-n summary file name, required")
    print("-b batch file, optional")
    print("-l filter labels, optional, comma-separated, requires -b")
    print("-s standards library, optional, function not implemented")
    print("-x toggle replicate summarization, optional, requires -b")
    print("-h help, optional")
    sys.exit()
  if o == "-i":
    print(a)
    intensities = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-r":
    rtimes = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-m":
    masscharges = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-q":
    quants = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-n":
    summaryname = os.path.abspath(a)
  if o == "-b":
    batchinfo = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-l":
    labels = a.split(",")
  if o == "-s":
    stdlib = a
  if o == "-x":
    hasrep = True

# get list of samples based on input
samples = []
sources = []
libsamples = []
libsources = []
if batchinfo is not None and labels:
  # get samples and sources for selected sample types
  for label in labels:
    samples = samples + batchinfo[batchinfo["type"] == label]["sample"].tolist()
    if hasrep:
      sources = sources + batchinfo[batchinfo["type"] == label]["id"].unique().tolist()
  # get samples and sources for standard library
  if stdlib:
    libsamples = batchinfo[batchinfo["type"] == stdlib]["sample"].tolist()
    if hasrep:
      libsources = batchinfo[batchinfo["type"] == stdlib]["id"].unique().tolist()
else:
  samples = intensities.columns[14:].tolist()
  samples.sort()

# create initial fragment summaries
summaries = masscharges[["tmz", "id", "subid", "name", "monoisotopic", "cas", "formula"]]

# summarize rtimes
rtimes = rtimes[samples]
rtimes[rtimes == 0] = np.nan
if hasrep:
  # get medians for replicates
  replicates = pd.DataFrame()
  for source in sources:
    replicates[source] = rtimes[batchinfo[batchinfo["id"] == source]["sample"].tolist()].median(axis=1)
  rtimes = replicates
summaries["rtMed"] = rtimes.median(axis=1)
summaries["rtMin"] = rtimes.min(axis=1)
summaries["rtMax"] = rtimes.max(axis=1)
summaries["rtRange"] = summaries["rtMax"].sub(summaries["rtMin"])

# summarize masses
masscharges = masscharges[samples]
masscharges[masscharges == 0] = np.nan
if hasrep:
  # get medians for replicates
  replicates = pd.DataFrame()
  for source in sources:
    replicates[source] = masscharges[batchinfo[batchinfo["id"] == source]["sample"]].median(axis=1)
  masscharges = replicates
summaries["mzMed"] = masscharges.median(axis=1)
summaries["mzMin"] = masscharges.min(axis=1)
summaries["mzMax"] = masscharges.max(axis=1)
summaries["mzRange"] = summaries["mzMax"].sub(summaries["mzMin"])

# summarize standards if library provided
if libsamples:
  libintensities = intensities[libsamples]
  libintensities[libintensities == 0] = np.nan
  if hasrep:
    replicates = pd.DataFrame()
    for source in libsources:
      replicates[source] = libintensities[batchinfo[batchinfo["id"] == source]["sample"]].median(axis=1)
    libintensities = replicates
  # basic info on detection in the standard
  summaries["libDetectCnt"] = libintensities.count(axis=1)
  summaries["libDetectFrac"] = summaries["libDetectCnt"].div(len(libintensities.columns))
  # more complex info about quality in the standard
  ids = summaries[(summaries["libDetectCnt"] > 0)]["id"].drop_duplicates()
  summaries["libFragCnt"] = 0
  summaries["libQualityCnt"] = 0
  for i in ids:
    # get the number of detected fragments per id in std
    # get number of quality detections per id in std
    #   a quality detection is defined as 3 detected fragments in one sample
    frags = summaries.index[(summaries["id"] == i) & (summaries["libDetectCnt"] > 0)].tolist()
    summaries.loc[frags, "libFragCnt"] = len(frags)
    simulCnt = libintensities.loc[frags].count(axis=0)
    summaries.loc[frags, "libQualityCnt"] = (simulCnt >= 3).sum()
  summaries["libQualityFrac"] = summaries["libQualityCnt"].div(len(libintensities.columns))

# summarize intensities
intensities = intensities[samples]
intensities[intensities == 0] = np.nan
if hasrep:
  # get medians for replicates
  replicates = pd.DataFrame()
  for source in sources:
    replicates[source] = intensities[batchinfo[batchinfo["id"] == source]["sample"]].median(axis=1)
  intensities = replicates
summaries["iMean"] = intensities.mean(axis=1)
summaries["iStd"] = intensities.std(axis=1)
summaries["iCV"] = summaries["iStd"].div(summaries["iMean"])
summaries["detectCnt"] = intensities.count(axis=1)
summaries["detectFrac"] = summaries["detectCnt"].div(len(intensities.columns))

# identification quality and best fragment identification
ids = summaries[(summaries["detectCnt"] > 0)]["id"].drop_duplicates()
summaries["fragCnt"] = 0
summaries["qualityCnt"] = 0
summaries["bestFlag"] = 0
for i in ids:
  # get the number of detected fragments per id
  # get number of quality detections per id
  #   a quality detection is defined as 3 detected fragments in one sample
  frags = summaries.index[(summaries["id"] == i) & (summaries["detectCnt"] > 0)].tolist()
  summaries.loc[frags, "fragCnt"] = len(frags)
  simulCnt = intensities.loc[frags].count(axis=0)
  summaries.loc[frags, "qualityCnt"] = (simulCnt >= 3).sum()
  # identify the best fragment for analysis via detected count and subid
  #   in the future alternative metrics like CV might be considered
  detectCntMax = summaries.loc[frags]["detectCnt"].max()
  subidMin = summaries.loc[(summaries["id"] == i) & (summaries["detectCnt"] == detectCntMax), "subid"].min()
  summaries.loc[(summaries["id"] == i) & (summaries["detectCnt"] == detectCntMax) & (summaries["subid"] == subidMin), "bestFlag"] = 1
summaries["qualityFrac"] = summaries["qualityCnt"].div(len(intensities.columns))

# if concentration is provided add associated details
if quants is not None:
  summaries["beta"] = quants["beta"]
  quants = quants[samples]
  quants[quants == 0] = np.nan
  if hasrep:
    # get medians for replicates
    replicates = pd.DataFrame()
    for source in sources:
      replicates[source] = quants[batchinfo[batchinfo["id"] == source]["sample"]].median(axis=1)
    quants = replicates

# select and rearrange columns based on selected inputs
summaryFields = ["mzMed", "rtMed", "tmz", "id", "subid", "name", "monoisotopic", "formula", "cas", "mzMin", "mzMax", "mzRange", "rtMin", "rtMax", "rtRange", "iMean", "iStd", "iCV"]
if libsamples:
  summaryFields = summaryFields + ["libFragCnt", "libQualityFrac"]
summaryFields = summaryFields + ["detectCnt", "detectFrac", "fragCnt", "qualityCnt", "qualityFrac", "bestFlag"]
if quants is not None:
  summaryFields = summaryFields + ["beta"]
summaries = summaries[summaryFields]

# add intensities or concentrations to summary
if quants is None:
  summaries = pd.concat([summaries, intensities], axis=1)
else:
  summaries = pd.concat([summaries, quants], axis=1)

# filter rows
summaries = summaries[summaries["detectCnt"] > 0]

# write summary file
summaries.to_csv(summaryname)

