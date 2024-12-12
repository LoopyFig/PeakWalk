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

# import dtw
from fastdtw import fastdtw

# set global variables
drtMax = 0.02

# DTW similarity score
def DTW(trtdelta, rtdelta):
  _, path = fastdtw(rtdelta, trtdelta)
  length = max(len(rtdelta), len(trtdelta))
  similarity = 1 - len(path) / length
  return similarity

# Search for targets on a single sample
# using DWT criteria
def matchDTW(matches, idSubids, mzRts):
  # Get Relative distances for targets and matches
  trtdeltas = abs(idSubids.trt.values[:, None] - idSubids.trt.values).tolist()
  rtdeltas = abs(mzRts.rt.values[:, None] - mzRts.rt.values).tolist()

  # Get DWT scores for all relevant matches
  for idSubid in idSubids.idSubid:
    submatches = matches[(matches.idSubid == idSubid)].index
    trtdelta = trtdeltas[idSubid]
    for match in submatches:
      mzRt = matches.mzRt.at[match]
      rtdelta = rtdeltas[mzRt]
      matches.dtw.loc[match] = abs(DTW(trtdelta, rtdelta))

  # Perform a queued match based on DTW score
  idSubidQueue = idSubids.idSubid.tolist()
  idSubidQueue.reverse()
  while idSubidQueue:
    idSubid = idSubidQueue.pop()
    submatches = matches[(matches.idSubid == idSubid)].index
    dtwBest = np.inf
    mzRtBest = -1
    matchBest = None
    for match in submatches:
      dtw = matches.dtw.at[match]
      mzRt = matches.mzRt.at[match]
      if dtw < dtwBest:
        idSubidContest = mzRts.idSubid.iat[mzRt]
        if idSubidContest == -1:
          dtwBest = dtw
          mzRtBest = mzRt
          matchBest = match
        # FOR LATER: consider validity of contest, especially for split peaks
        elif dtw < matches.dtw.loc[(matches.idSubid == idSubidContest) & (matches.mzRt == mzRt)].iat[0]:
          dtwBest = dtw
          mzRtBest = mzRt
          idSubids.mzRt.iat[idSubidContest] = -1
          mzRts.idSubid.iat[mzRt] = -1
          idSubidQueue.append(idSubidContest)
    if mzRtBest > -1:
      idSubids.mzRt.iat[idSubid] = mzRtBest
      mzRts.idSubid.iat[mzRtBest] = idSubid
      matches.flag.at[match] = 1

  # return modified dataframes with match selections
  return (matches, idSubids, mzRts)

# Calculate inferred rt (irt) and delta rt (drt)
# If drop, removes previous irt and drt
# Modifies dataframes and returns a reasonable deltart bound (drtBound)
def irtFilter(ids, idSubids, mzRts):
  # If drop, removes the previous irt and drt
  if "irt" in ids.columns:
    ids.drop(columns="irt",inplace=True)
    idSubids.drop(columns=["irt","drt","mz","rt"],inplace=True)

  # Calculate inferred retention time per id (irt) as low median of the subid rts
  # Calculate delta rt (drt) from abs(irt-rt)
  # NOTE: right now the apply uses a "low median" method, but it's not clear
  # that this is always best practice. a slightly more sophisticated method would
  # look for a "better median" rather than a low one, so we'll try that later
  irts = pd.merge(idSubids[idSubids.mzRt != -1], mzRts[mzRts.idSubid != -1],
    on=["idSubid","mzRt"]).groupby("id")["rt"].apply(lambda x: x.sort_values().iat[int(len(x)/2 + len(x)%2)-1])
  irts.name = "irt"
  ids = pd.merge(ids, irts, on="id", how="left")
  idSubids = pd.merge(idSubids, ids, on="id")
  idSubids = pd.merge(idSubids, mzRts[["mzRt","mz","rt"]], on="mzRt", how="left")
  idSubids["drt"] = abs(idSubids.irt - idSubids.rt)
  idSubids = idSubids.sort_values(by="idSubid")

  # mzRts should match idSubids
  # print("mzRts: "  + str(len(mzRts[mzRts.idSubid != -1])))
  # print("idSubids: "  + str(len(idSubids[idSubids.mzRt != -1])))

  # Clean-up matches outside a reasonable rt bound
  # Bound is approximated via double the 3rd quantile of deltart
  drtBound = 2*idSubids[~np.isnan(idSubids.drt)].drt.quantile(0.75)
  if drtBound > drtMax:
    drtBound = drtMax
  return (ids, idSubids, drtBound)

# Perform first shot matching based on rt
# Get inferred rt and delta rt bound
# Filter out matches outside delta rt bound
def firstShot(matches, ids, idSubids, mzRts, drtBoundLimit = 0):
  # Get best first-shot match based on rt
  # Scan through unique ids to get all candidate matches
  # Scan through unique subids to get all subcandidate matches
  # Assign first free mzRt to idSubid in order of rt
  # Link mzRt to idSubid for second shot matching
  for id in ids.id:
    #print("loop")
    candidates = matches[(matches.id == id) & (matches.flag == 0)]
    #print(candidates)
    subids = idSubids[(idSubids.id == id) & (idSubids.mzRt == -1)].subid.unique()
    #print(subids)
    for i in range(0, len(subids)):
      subid = subids[i]
      subcandidates = candidates[candidates.subid == subid]
      for j in subcandidates.index:
        mzRt = subcandidates.mzRt[j]
        idSubid = subcandidates.idSubid[j]
        if mzRts.idSubid.iloc[mzRt] == -1:
          mzRts.idSubid.iat[mzRt] = idSubid
          idSubids.mzRt.iat[idSubid] = mzRt
          matches.flag.at[j] = 2
          break

  # Calculate inferred retention time per id (irt) as median of the subid rts
  # Calculate delta rt (drt) from abs(irt-rt)
  # Get deltart bound (drtBound)
  ids, idSubids, drtBound = irtFilter(ids, idSubids, mzRts)
  # print("drtBound: " + str(drtBound))
  if drtBoundLimit > drtBound:
    drtBound = drtBoundLimit

  # Clean-up matches outside a reasonable rt range
  # Range is approximated via double the 3rd quantile of deltart
  for i in idSubids[idSubids.drt > drtBound].index:
    mzRt = idSubids.mzRt[i]
    idSubid = idSubids.idSubid[i]
    mzRts.idSubid.iat[mzRt] = -1
    idSubids.mzRt.at[i] = -1
    # NOTE: prior version only flags filtered results
    # now we just just flag all prior assignments above
    # this prevents getting "stuck" on single bad solution when running one-shot
    # at the cost of potentially missing some matches that arise through corrections
    # explore this later
    # matches.flag.loc[(matches.mzRt == mzRt) & (matches.idSubid == idSubid)] = 1

  # Drop outdated irt and drt and update
  # We keep the same bound, as we want to avoid a filtering loop
  # A mzRt candidate is a idSubid match within the drtBound
  ids, idSubids, drtHypothesis = irtFilter(ids, idSubids, mzRts)
  # print("drtHypothesis One: " + str(drtHypothesis))

  # Return updated dataframes and bound
  return (matches, ids, idSubids, mzRts, drtBound, drtHypothesis)

# Perform second-shot correction of matches
def secondShot(matches, ids, idSubids, mzRts, drtBound):
  # Selection occurs in reverse rt order
  # Method:
  # 1. Check if idSubid slot is empty
  # True:
  #   2a. Check if there are any unassigned mzRt candidates
  #   True:
  #     3aa. Choose best (lowest drt) candidate to assign
  #   False:
  #     3ab. Check if there are any assigned mzRt candidates
  #     4ab. Compare if better drt can be achieved
  #     True:
  #       5aba: Choose best (lowest drt) candidate to steal
  # False:
  #   2b. Check if there are any unassigned mzRt candidates with better drt
  #   3b. Compare if better drt can be achieved
  #   True:
  #     4ba. Choose best (lowest drt) candidate to replace
  for id in ids.id[~np.isnan(ids.irt)].iloc[::-1]:
    candidates = matches[(matches.id == id)]
    for i in idSubids.loc[idSubids.id == id].iloc[::-1].index:
      subcandidates = candidates[candidates.subid == idSubids.subid[i]]
      take = [-1, 0]
      steal = [-1, 0]
      if idSubids.mzRt[i] == -1:
        for j in subcandidates.iloc[::-1].index:
          mzRt = subcandidates.mzRt[j]
          drt = abs(idSubids.irt[i] - mzRts.rt.iloc[mzRt])
          if drt > drtBound:
            continue
          idSubid = mzRts.idSubid.iloc[mzRt]
          if idSubid == -1 and (take[0] == -1 or take[1] > drt):
            take = [mzRt, drt]
          elif take[0] == -1 and idSubids.drt.iloc[idSubid] > drt and (steal[0] == -1 or steal[1] > drt):
            steal = [mzRt, drt]
        if take[0] != -1:
          mzRts.idSubid.iloc[take[0]] = idSubids.idSubid[i]
          idSubids.mzRt.loc[i] = take[0]
          idSubids.drt.loc[i] = take[1]
        elif steal[0] != -1:
          mzRts.idSubid.iloc[take[0]] = idSubids.idSubid[i]
          idSubid = mzRts.idSubid.iloc[steal[0]]
          idSubids.mzRt.iloc[idSubid] = -1
          idSubids.irt.iloc[idSubid] = -1
          idSubids.drt.iloc[idSubid] = -1
          idSubids.mzRt.loc[i] = steal[0]
          idSubids.drt.loc[i] = steal[1]
      else:
        for j in subcandidates.iloc[::-1].index:
          mzRt = subcandidates.mzRt[j]
          drt = abs(idSubids.irt[i] - mzRts.rt.iloc[mzRt])
          if drt > drtBound:
            continue
          idSubid = mzRts.idSubid.iloc[mzRt]
          if idSubid == -1 and (take[0] == -1 or take[1] > drt):
            take = [mzRt, drt]
        if take[0] != -1 and idSubids.drt[i] > take[1]:
          mzRts.idSubid.iloc[take[0]] = idSubids.idSubid[i]
          idSubids.mzRt.loc[i] = take[0]
          idSubids.drt.loc[i] = take[1]

  # Update irt, get drtHypothesis
  ids, idSubids, drtHypothesis = irtFilter(ids, idSubids, mzRts)

  # Return updated dataframes and drtHypothesis
  return (matches, ids, idSubids, mzRts, drtHypothesis)

# Search for targets on a single sample
def targetSearch(targets, sample, boundTestLimit = np.inf, dtw = True):
  # Search original dataframe
  # NOTE: in this implementation, Order and Fragment shares the same window as DTW
  # this is not strictly necessary, though implementation of separate bounds is messy
  matches = []
  sample = sample[sample[sample.columns[2]] > 0]
  for i in targets.index:
    match = sample[(sample.mz <= targets.tmzupper[i]) & (sample.mz >= targets.tmzlower[i]) & (sample.rt <= targets.trtupper[i]) & (sample.rt >= targets.trtlower[i])]
    if len(match) > 0:
      match["id"] = targets.id[i]
      match["subid"] = targets.subid[i]
      match["trt"] = targets.trt[i]
      match["tmz"] = targets.tmz[i]
      matches.append(match)

  # handle case of no-match
  if len(matches) == 0:
    print(sample.columns[2], " no-match")
    return pd.DataFrame(columns=["id","subid","mz","rt",sample.columns[2],"tmz","trt"])

  # Concatenate matches and sort by rt, irt, id, subid
  # Get unique dataframe of matched mzRts and unique list of idSubids
  matches = pd.concat(matches, axis=0).sort_values(["rt","trt","id","subid"])
  # print(sample.columns[2], " initial matches: ", str(len(matches)))
  mzRts = matches[["mz","rt"]].drop_duplicates()
  idSubids = matches[["id","subid","trt"]].drop_duplicates()

  # Add mzRt pk and idSubid to mzRts
  # Add idSubid pk and mzRt to subIds
  # Add mzRt to matches
  mzRts["mzRt"] = range(0, len(mzRts))
  idSubids["idSubid"] = range(0, len(idSubids))

  # Add mzRt and idSubid to matches for quick indexing
  # Cost is approximately O(n^2), could we do better?
  # Sort may cost O(nlogn)
  matches = pd.merge(matches, mzRts, on=["mz","rt"])
  matches = pd.merge(matches, idSubids[["id","subid","idSubid"]], on=["id","subid"])
  mzRts["idSubid"] = [-1]*len(mzRts)
  idSubids["mzRt"] = [-1]*len(idSubids)
  matches["dtw"] = [0]*len(matches)
  matches["flag"] = [0]*len(matches) # not relevant yet

  # Get unique dataframe of ids
  # ids = idSubids[["id","trt"]].drop_duplicates()
  ids = idSubids[["id"]].drop_duplicates()

  # DTW Score based matching for initialization
  if dtw:
    matches, idSubids, mzRts = matchDTW(matches, idSubids, mzRts)

  # Cycle through first shot matching
  matchState = idSubids.mzRt[idSubids.mzRt != -1] # previously set to None
  matchCount = matchState.count() # previously set to 0
  drtBound = drtHypothesis = drtBoundLimit = 0
  boundTest = 0
  while True:
    # print(sample.columns[2] + " matchCount: " + str(matchCount) + " drtBound: " + str(drtBound) + " ids: " + str(len(idSubids.index)))
    if boundTest == boundTestLimit:
      matches, ids, idSubids, mzRts, drtBound, drtHypothesis = firstShot(matches, ids, idSubids, mzRts, drtBoundLimit)
    else:
      matches, ids, idSubids, mzRts, drtBound, drtHypothesis = firstShot(matches, ids, idSubids, mzRts)
      boundTest = boundTest + 1
      drtBoundLimit = drtBound
    newState = idSubids.mzRt[idSubids.mzRt != -1]
    newCount = newState.count()
    if newCount == matchCount:
      if matchCount == 0 or matchState.equals(newState):
        break
    matchState = newState
    matchCount = newCount

  # print("drtBound: " + str(drtBound))

  # Second shot matching
  matches, ids, idSubids, mzRts, drtHypothesis = secondShot(matches, ids, idSubids, mzRts, drtBound)
  # print("drtHypothesis Two: " + str(drtHypothesis))

  # Merge onto the original matches to filter final selection
  matches = pd.merge(matches, idSubids[["idSubid","mzRt","irt","drt"]], on=["idSubid","mzRt"], how="inner")
  # print(sample.columns[2], " done\n")
  return matches[["id","subid","mz","rt",sample.columns[2],"tmz","trt","irt","drt"]]

# Search directory for ADAP files
# NOTE: Currently depth functionality is defunct
def findADAP(dir, ADAP, depth):
  for element in os.listdir(dir):
    elementpath = os.path.join(dir, element)
    if os.path.isfile(elementpath):
      if len(element) > 4 and element[-4:] == ".csv":
        ADAP.put(element)
    elif depth < 0:
      ADAP = findADAP(elementpath, ADAP, depth+1)
  return ADAP

# Estimate overall shift of list
def shiftRt(targets, matches, trtSmallBound):
  ids = matches["id"].unique()
  shifts = []
  rtShift= 0
  for id in ids:
    rts = matches.loc[matches["id"] == id, ["rt","trt"]]
    if len(rts) >= 3:
      shifts = shifts + rts["rt"].sub(rts["trt"]).tolist()
  if len(shifts) > 0:
    rtShift = pd.Series(shifts).median()
    targets["trtupper"] = targets.trt.apply(lambda x: x+trtSmallBound+rtShift)
    targets["trtlower"] = targets.trt.apply(lambda x: x-trtSmallBound+rtShift)
  return [targets, rtShift]

# Read metabolic features
# Extract mz time information with label
# Search for targets and get best match, sample by sample
def runSample(targets, mzindex, rtindex, iindex, boundTestLimit, adapdir, ADAP, samples, lock, redundancy, dtw, shift, trtSmallBound):
  lock.acquire()
  if shift:
    ortBounds = targets[["trtupper", "trtlower"]].copy()
  while not ADAP.empty():
    adap = ADAP.get()
    lock.release()
    sample = pd.read_csv(adapdir + "/" + adap).iloc[:,[mzindex,rtindex,iindex]]
    sample.columns = ["mz","rt"]+[sample.columns[2].replace(redundancy, "")]
    # matches = targetSearch(targets, sample, 1)
    if shift:
      matches = targetSearch(targets, sample, boundTestLimit, dtw)
      targets, rtShift = shiftRt(targets, matches, trtSmallBound)
      print(adap + " " + str(rtShift))
      matches = targetSearch(targets, sample, boundTestLimit, dtw)
      targets[["trtupper","trtlower"]] = ortBounds[["trtupper", "trtlower"]]
    else:
      matches = targetSearch(targets, sample, boundTestLimit, dtw)
    samples.put(matches)
    lock.acquire()
  lock.release()
  return

def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "a:f:t:p:")
  except getopt.GetoptError as err:
    print(err)
    sys.exit(2)
  rawdir=None
  adapdir=None
  targetlist=None
  processors=1
  print(opts)
  for o, a in opts:
    if o == "-a":
      print(a)
      adapdir = os.path.abspath(a)
    elif o == "-f":
      print(a)
      featuredir = os.path.abspath(a)
    elif o == "-t":
      print(a)
      targetlist = os.path.abspath(a)
    elif o == "-p":
      print(a)
      processors = int(a)
    else:
      print(o)
      print(a)
      print("unhandled option")
      #sys.exit(2)
  if not (adapdir and featuredir and targetlist):
    print("missing options")
    sys.exit(2)

  # Get list of ADAP feature tables
  ADAP = Queue()
  ADAP=findADAP(adapdir, ADAP, 0)
  sampleNum = ADAP.qsize()

  # Set variables for reading file
  # Indexes are for the sample files, not the target list
  # NOTE: we need a good way to handle multiple target list formats
  mzindex = 0
  rtindex = 1
  iindex = 3
  redundancy = ".mzXML Peak area"

  # Set variables for processing
  boundTestLimit = 1 # this variable handles the number of times drtBound is optimized
  trtBound = 0.3
  trtSmallBound = 0.3
  # trt = "min" # apparently mzmine output is in minutes
  shift = True # parameter for whether we apply an automatic list shift
  dtw = True # parameter for whether we start with dtw

  # Read target list
  # calculate deltappm range
  targets = pd.read_csv(targetlist, sep=',')
  #targets.columns=["id","name","tmz","trt","monoisotopic","cas","subid"]
  targets.columns=["id","name","tmz","trt","monoisotopic","cas","subid","formula","concentration"]
  #targets.columns=["id","name","tmz","trt","monoisotopic","cas","subid","formula","concentration","note","realRt"]
  targets["tmzupper"] = targets.tmz.apply(lambda x: x+x*0.000006)
  targets["tmzlower"] = targets.tmz.apply(lambda x: x-x*0.000006)
  # if trt == "min":
  #   targets.trt = targets.trt*60.0
  targets["trtupper"] = targets.trt.apply(lambda x: x+trtBound)
  targets["trtlower"] = targets.trt.apply(lambda x: x-trtBound)

  # Set up and end processes for running samples
  samples = Queue()
  lock = Lock()
  for p in range(0,processors):
    worker = Process(target = runSample, args = (targets, mzindex, rtindex, iindex, boundTestLimit, adapdir, ADAP, samples, lock, redundancy, dtw, shift, trtSmallBound), daemon = True)
    worker.start()

  intensities = targets.copy()
  rtimes = targets.copy()
  masscharges = targets.copy()
  # I recognize that the triple merge is inefficient here
  # but I'm just gonna fix it later
  while sampleNum > 0:
    sampleNum = sampleNum - 1
    sample = samples.get()
    sampleId = sample.columns[4]
    intensities = pd.merge(intensities, sample.iloc[:,[0,1,4]],on=["id","subid"],how="left")
    sample.drop(columns=sampleId, inplace=True)
    sample.rename(columns={"rt":sampleId}, inplace=True)
    rtimes = pd.merge(rtimes, sample.iloc[:,[0,1,3]],on=["id","subid"],how="left")
    sample.drop(columns=sampleId, inplace=True)
    sample.rename(columns={"mz":sampleId}, inplace=True)
    masscharges = pd.merge(masscharges, sample.iloc[:,[0,1,2]],on=["id","subid"],how="left")

  intensities.fillna(0).to_csv(featuredir + "/" + "feature.sample.i.csv")
  rtimes.fillna(0).to_csv(featuredir + "/" + "feature.sample.rt.csv")
  masscharges.fillna(0).to_csv(featuredir + "/" + "feature.sample.mz.csv")

if __name__ == "__main__":
  main()

