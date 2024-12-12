import sys, os
import pandas as pd
import re

dir = sys.argv[1]
batchPattern = r'.*[Bb]atch([0-9]+).*\.csv$'

batches = []
for file in os.listdir(dir):
  path = os.path.join(dir, file)
  if os.path.isfile(path):
    isbatch = re.search(pattern, file)
    if isbatch:
      batch = pd.read_csv(path)
      batch["batch"] = int(isbatch.group(0))
      batches.append(batch)

batches = pd.concat(batches)
batch = batch[["Filename","Sample.ID","Batch"]]
batch.columns = ["sample","id","batch"]
batch["id"] = batch["id"].str.replace(r'_[0-9]+$','',regex = True)
batch["type"] = "ignore"
batch.loc[batch["sample"].str.fullmatch("^BL_.*"), "type"] = "subject"
batch.loc[batch["sample"].str.fullmatch("^wash.*"), "type"] = "wash"
batch.loc[batch["sample"].str.fullmatch("^NIST1958.*"), "type"] = "NIST1958"
batch.loc[batch["sample"].str.fullmatch("^NIST1975.*"), "type"] = "NIST1975"
batch.loc[batch["sample"].str.fullmatch("^BB.*"), "type"] = "BB"
batch.loc[batch["sample"].str.fullmatch("^BM.*"), "type"] = "BM"
batch.loc[batch["sample"].str.fullmatch("^RT.*"), "type"] = "RT"
batch.loc[batch["sample"].str.fullmatch("^Qstd.*"), "type"] = "Qstd"
batch.loc[batch["sample"].str.fullmatch("^BP1.*"), "type"] = "BP1"
batch.loc[batch["sample"].str.fullmatch("^BP2.*"), "type"] = "BP2"
batch.loc[batch["sample"].str.fullmatch("^BP3.*"), "type"] = "BP3"
batch.loc[batch["sample"].str.fullmatch("^BP4.*"), "type"] = "BP4"
batch.loc[batch["sample"].str.fullmatch("^BP5.*"), "type"] = "BP5"
batch.loc[batch["sample"].str.fullmatch("^CBP1.*"), "type"] = "CBP1"
batch.loc[batch["sample"].str.fullmatch("^CBP2.*"), "type"] = "CBP2"
batch.loc[batch["sample"].str.fullmatch("^CBP3.*"), "type"] = "CBP3"
batch.loc[batch["sample"].str.fullmatch("^CBP4.*"), "type"] = "CBP4"
batch.loc[batch["sample"].str.fullmatch("^CBP5.*"), "type"] = "CBP5"
batch.loc[batch["sample"].str.fullmatch("^QBP1.*"), "type"] = "QBP1"
batch.loc[batch["sample"].str.fullmatch("^QBP2.*"), "type"] = "QBP2"
batch.loc[batch["sample"].str.fullmatch("^QBP3.*"), "type"] = "QBP3"
batch.loc[batch["sample"].str.fullmatch("^QBP4.*"), "type"] = "QBP4"
batch.loc[batch["sample"].str.fullmatch("^QBP5.*"), "type"] = "QBP5"
batch["dilu"] = "NA"
batch.loc[batch["sample"].str.fullmatch("^BP.*"), "dilu"] = 1
batch.loc[batch["sample"].str.fullmatch("^CBP.*"), "dilu"] = 1
batch = batch[["sample","type","dilu","batch","id"]]
batch.to_csv("batch.csv",index=False)

