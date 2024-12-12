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
  opts, args = getopt.getopt(sys.argv[1:], "f:i:r:m:q:n:b:l:s:c:xh", ["quant"])
except getopt.GetoptError as err:
  print(err)
  sys.exit(2)
# required parameters
featuresLoc = ""
intensitiesLoc = "feature.sample.i.csv"
rtimesLoc = "feature.sample.rt.csv"
masschargesLoc = "feature.sample.mz.csv"
quantsLoc = "feature.sample.quant.csv"
summaryname = [
  "feature.sample.summary.csv","feature.subject.summary.csv",
  "feature.sample.qsummary.csv","feature.subject.qsummary.csv"
  ] # summary file naming options based on input
defaultnameToggle = True
# additional parameters
batchinfo = None # batch file
labels = ["subject"] # sample type filter
stdlib = None # standards library
repToggle = False # toggles replicate summarization
quantToggle = False # toggles quantification options
corrToggle = False
detectMin = sys.float_info.min
for o, a in opts:
  if o == "-h":
    print("PeakWalk: Automated GC Identification and Quantification")
    print("-f feature files location, optional")
    print("  default feature files location is working directory")
    print("  sets location for input and output")
    print("-i intensity file, optional")
    print("-r retention time file, optional")
    print("-m masscharge file, optional")
    print("  -irm are optional, but feature files are required, see -f")
    print("  default file names are feature.sample.{i, mz, rt}.csv")
    print("-q quantification file, optional, requires -b")
    print("  overrides intensity with concentration")
    print("-n summary file name, optional")
    print("-b batch file, optional")
    print("-l filter labels, optional, comma-separated, requires -b")
    print("-s standards library, optional, function not implemented")
    print("-x toggle replicate summarization, optional, requires -b")
    print("-c sets correlation filter, optional")
    print("-h help, optional")
    sys.exit()
  if o == "-f":
    featuresLoc = a + "/"
  if o == "-i":
    print(a)
    intensitiesLoc = a
  if o == "-r":
    rtimesLoc = a
  if o == "-m":
    masschargesLoc = a
  if o == "-q":
    quantToggle = True
    quantsLoc = a
  if o == "--quant":
    quantToggle = True
  if o == "-n":
    defaultnameToggle = False
    summaryname = a
  if o == "-b":
    batchinfo = pd.read_csv(os.path.abspath(a), sep=",")
  if o == "-l":
    labels = a.split(",")
  if o == "-s":
    stdlib = a
  if o == "-x":
    repToggle = True
  if o == "-c":
    corrToggle = True
    corrMin = float(a)

# read in data based on input
intensities = pd.read_csv(os.path.abspath(featuresLoc + intensitiesLoc), sep=",")
rtimes = pd.read_csv(os.path.abspath(featuresLoc + rtimesLoc), sep=",")
masscharges = pd.read_csv(os.path.abspath(featuresLoc + masschargesLoc), sep=",")
quants = None
if quantToggle:
  quants = pd.read_csv(os.path.abspath(featuresLoc + quantsLoc), sep=",")
if defaultnameToggle:
  if quantToggle:
    if repToggle:
      summaryname = summaryname[3]
    else:
      summaryname = summaryname[2]
  else:
    if repToggle:
      summaryname = summaryname[1]
    else:
      summaryname = summaryname[0]
summaryname = featuresLoc + summaryname

# get list of samples based on input
samples = []
sources = []
libsamples = []
libsources = []
if batchinfo is not None and labels:
  # get samples and sources for selected sample types
  for label in labels:
    samples = samples + batchinfo[batchinfo["type"] == label]["sample"].tolist()
    if repToggle:
      sources = sources + batchinfo[batchinfo["type"] == label]["id"].unique().tolist()
  # get samples and sources for standard library
  if stdlib:
    libsamples = batchinfo[batchinfo["type"] == stdlib]["sample"].tolist()
    if repToggle:
      libsources = batchinfo[batchinfo["type"] == stdlib]["id"].unique().tolist()
else:
  samples = intensities.columns[14:].tolist()
  samples.sort()

# create initial fragment summaries
summaries = masscharges[["tmz", "trt", "id", "subid", "name", "monoisotopic", "cas", "formula"]]

# summarize rtimes
rtimes = rtimes[samples]
rtimes[rtimes == 0] = np.nan
if repToggle:
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
if repToggle:
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
  if repToggle:
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
if repToggle:
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

# quality-based filtering
summaries["badFlag"] = 0
summaries.loc[(summaries["detectFrac"] < detectMin), "badFlag"] = 1
# optional correlation-based filtering
if corrToggle:
  summaries.loc[(summaries["detectCnt"] < 3), "badFlag"] = 1
  ids = summaries[(summaries["badFlag"] == 0)]["id"].drop_duplicates()
  summaries["corrMean"] = 0
  for i in ids:
    frags = summaries.index[(summaries["id"] == i) & (summaries["badFlag"] == 0)].to_numpy()
    fragCnt = len(frags)
    print(str(i) + " " + str(frags) + " " + str(fragCnt))
    if fragCnt > 1:
      fragCorr = intensities.loc[frags].T.corr().values
      np.fill_diagonal(fragCorr, np.nan)
      noOverlap = np.isnan(fragCorr).all(axis=1)
      summaries.loc[frags[noOverlap], "badFlag"] = 1
      frags = frags[~noOverlap]
      fragCorr = fragCorr[~noOverlap, :][:, ~noOverlap]
      print(str(fragCorr))
      fragCnt = fragCnt - np.count_nonzero(noOverlap)
      while fragCnt > 1:
        meanCorr = np.nanmean(fragCorr, axis=0)
        worstFrag = np.argmin(meanCorr)
        print(worstFrag)
        print(meanCorr)
        if meanCorr[worstFrag] < corrMin:
          summaries.loc[frags[worstFrag], "badFlag"] = 1
          frags = np.delete(frags, worstFrag)
          fragCorr = np.delete(fragCorr, worstFrag, 0)
          fragCorr = np.delete(fragCorr, worstFrag, 1)
          fragCnt -= 1
        else:
          for j in range(0, len(frags)):
            summaries.loc[frags[j], "corrMean"] = meanCorr[j]
          break
    if fragCnt <= 1:
      summaries.loc[frags, "badFlag"] = 1

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
  summaries.loc[(summaries["id"] == i) & (summaries["subid"] == subidMin), "bestFlag"] = 1
summaries["qualityFrac"] = summaries["qualityCnt"].div(len(intensities.columns))

# if concentration is provided add associated details
if quantToggle:
  summaries["beta"] = quants["beta"]
  quants = quants[samples]
  quants[quants == 0] = np.nan
  if repToggle:
    # get medians for replicates
    replicates = pd.DataFrame()
    for source in sources:
      replicates[source] = quants[batchinfo[batchinfo["id"] == source]["sample"]].median(axis=1)
    quants = replicates
  # identify the best fragment for quantification via detected count and subid and presence of standards
  ids = summaries[(summaries["detectCnt"] > 0) & (summaries["beta"] > 0)]["id"].drop_duplicates()
  summaries["bestQuantFlag"] = 0
  for i in ids:
    frags = summaries.index[(summaries["id"] == i) & (summaries["detectCnt"] > 0) & (summaries["beta"] > 0)].tolist()
    detectCntMax = summaries.loc[frags]["detectCnt"].max()
    subidMin = summaries.loc[(summaries["id"] == i) & (summaries["detectCnt"] == detectCntMax) & (summaries["beta"] > 0), "subid"].min()
    summaries.loc[(summaries["id"] == i) & (summaries["subid"] == subidMin), "bestQuantFlag"] = 1

# select and rearrange columns based on selected inputs
summaryFields = ["mzMed", "rtMed", "tmz", "trt", "id", "subid", "name", "monoisotopic", "formula", "cas", "mzMin", "mzMax", "mzRange", "rtMin", "rtMax", "rtRange", "iMean", "iStd", "iCV"]
if libsamples:
  summaryFields = summaryFields + ["libFragCnt", "libQualityFrac"]
summaryFields = summaryFields + ["detectCnt", "detectFrac", "fragCnt", "qualityCnt", "qualityFrac", "bestFlag", "badFlag"]
if corrToggle:
  summaryFields = summaryFields + ["corrMean"]
if quantToggle:
  summaryFields = summaryFields + ["beta", "bestQuantFlag"]
summaries = summaries[summaryFields]

# add intensities or concentrations to summary
if not quantToggle:
  summaries = pd.concat([summaries, intensities], axis=1)
else:
  summaries = pd.concat([summaries, quants], axis=1)
  summaries = summaries[summaries["beta"] > 0]

# filter rows
summaries = summaries[summaries["badFlag"] == 0]
summaries.drop('badFlag', axis=1, inplace=True)

# write summary file
summaries.to_csv(summaryname)

