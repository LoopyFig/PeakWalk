#!/usr/bin/python3

import pandas as pd
import sys, os
import re

featDir = sys.argv[1]
fileType = sys.argv[2]

fileName = "feature." + fileType + ".csv"
libs = []
for libDir in os.listdir(featDir):
  if os.path.isdir(featDir + "/" + libDir):
    lib = pd.read_csv(featDir + "/" + libDir + "/" + fileName)
    lib.insert(0, "lib", re.search(r'.*?([^\/]+)\/?$', libDir).group(1))
    libs.append(lib)

libs = pd.concat(libs)
libs.sort_values(by = ["lib", "id", "subid"], inplace = True)
libs.to_csv(featDir + "/lib." + fileType + ".csv", index = False)

