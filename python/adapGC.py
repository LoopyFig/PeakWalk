#!/usr/bin/env python3

import getopt, sys, os

def main():
  batchFile = (
"""<?xml version="1.0" encoding="UTF-8"?><batch>
    <batchstep method="io.github.mzmine.modules.io.rawdataimport.RawDataImportModule">
        <parameter name="Raw data file names">
            <file>{mzXML}</file>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetectionModule">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scans">
            <ms_level>1</ms_level>
        </parameter>
        <parameter name="Mass detector" selected="Centroid">
            <module name="Centroid">
                <parameter name="Noise level">1000.0</parameter>
            </module>
            <module name="Exact mass">
                <parameter name="Noise level"/>
            </module>
            <module name="Local maxima">
                <parameter name="Noise level"/>
            </module>
            <module name="Recursive threshold">
                <parameter name="Noise level"/>
                <parameter name="Min m/z peak width"/>
                <parameter name="Max m/z peak width"/>
            </module>
            <module name="Wavelet transform">
                <parameter name="Noise level"/>
                <parameter name="Scale level"/>
                <parameter name="Wavelet window size (%)"/>
            </module>
        </parameter>
        <parameter name="Mass list name">masses</parameter>
        <parameter name="Output netCDF filename (optional)" selected="false"/>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_ADAPchromatogrambuilder.ADAPChromatogramBuilderModule">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scans">
            <ms_level>1</ms_level>
        </parameter>
        <parameter name="Mass list">masses</parameter>
        <parameter name="Min group size in # of scans">3</parameter>
        <parameter name="Group intensity threshold">1000.0</parameter>
        <parameter name="Min highest intensity">8000.0</parameter>
        <parameter name="m/z tolerance">
            <absolutetolerance>0.0</absolutetolerance>
            <ppmtolerance>5.0</ppmtolerance>
        </parameter>
        <parameter name="Suffix">chromatograms</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_smoothing.SmoothingModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
        <parameter name="Filename suffix">smoothed</parameter>
        <parameter name="Filter width">5</parameter>
        <parameter name="Remove original feature list">false</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.DeconvolutionModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
        <parameter name="Suffix">deconvoluted</parameter>
        <parameter name="Algorithm" selected="Baseline cut-off">
            <module name="Baseline cut-off">
                <parameter name="Min peak height">8000.0</parameter>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>1.0</max>
                </parameter>
                <parameter name="Baseline level">1000.0</parameter>
            </module>
            <module name="Noise amplitude">
                <parameter name="Min peak height"/>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="Amplitude of noise"/>
            </module>
            <module name="Savitzky-Golay">
                <parameter name="Min peak height"/>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="Derivative threshold level"/>
            </module>
            <module name="Local minimum search">
                <parameter name="Chromatographic threshold"/>
                <parameter name="Search minimum in RT range (min)"/>
                <parameter name="Minimum relative height"/>
                <parameter name="Minimum absolute height"/>
                <parameter name="Min ratio of peak top/edge"/>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
            </module>
            <module name="Wavelets (XCMS)">
                <parameter name="S/N threshold">10.0</parameter>
                <parameter name="Wavelet scales">
                    <min>0.25</min>
                    <max>5.0</max>
                </parameter>
                <parameter name="Peak duration range">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="Peak integration method">Use smoothed data</parameter>
                <parameter name="R engine">RCaller</parameter>
            </module>
            <module name="Wavelets (ADAP)">
                <parameter name="S/N threshold">10.0</parameter>
                <parameter name="S/N estimator">
                    <module name="Intensity window SN"/>
                    <module name="Wavelet Coeff. SN">
                        <parameter name="Peak width mult.">3.0</parameter>
                        <parameter name="abs(wavelet coeffs.)">true</parameter>
                    </module>
                </parameter>
                <parameter name="min feature height">10.0</parameter>
                <parameter name="coefficient/area threshold">110.0</parameter>
                <parameter name="Peak duration range">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="RT wavelet range">
                    <min>0.001</min>
                    <max>0.1</max>
                </parameter>
            </module>
        </parameter>
        <parameter measure="MEDIAN" name="m/z center calculation" weighting="NONE">CenterFunction</parameter>
        <parameter name="m/z range for MS2 scan pairing (Da)" selected="false"/>
        <parameter name="RT range for MS2 scan pairing (min)" selected="false"/>
        <parameter name="Remove original feature list">false</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.DeconvolutionModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
        <parameter name="Suffix">deconvoluted</parameter>
        <parameter name="Algorithm" selected="Wavelets (ADAP)">
            <module name="Baseline cut-off">
                <parameter name="Min peak height"/>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="Baseline level"/>
            </module>
            <module name="Noise amplitude">
                <parameter name="Min peak height"/>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="Amplitude of noise"/>
            </module>
            <module name="Savitzky-Golay">
                <parameter name="Min peak height"/>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="Derivative threshold level"/>
            </module>
            <module name="Local minimum search">
                <parameter name="Chromatographic threshold"/>
                <parameter name="Search minimum in RT range (min)"/>
                <parameter name="Minimum relative height"/>
                <parameter name="Minimum absolute height"/>
                <parameter name="Min ratio of peak top/edge"/>
                <parameter name="Peak duration range (min)">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
            </module>
            <module name="Wavelets (XCMS)">
                <parameter name="S/N threshold">10.0</parameter>
                <parameter name="Wavelet scales">
                    <min>0.25</min>
                    <max>5.0</max>
                </parameter>
                <parameter name="Peak duration range">
                    <min>0.0</min>
                    <max>10.0</max>
                </parameter>
                <parameter name="Peak integration method">Use smoothed data</parameter>
                <parameter name="R engine">RCaller</parameter>
            </module>
            <module name="Wavelets (ADAP)">
                <parameter name="S/N threshold">8.0</parameter>
                <parameter name="S/N estimator" selected="Intensity window SN">
                    <module name="Intensity window SN"/>
                    <module name="Wavelet Coeff. SN">
                        <parameter name="Peak width mult.">3.0</parameter>
                        <parameter name="abs(wavelet coeffs.)">true</parameter>
                    </module>
                </parameter>
                <parameter name="min feature height">8000.0</parameter>
                <parameter name="coefficient/area threshold">110.0</parameter>
                <parameter name="Peak duration range">
                    <min>0.0</min>
                    <max>1.0</max>
                </parameter>
                <parameter name="RT wavelet range">
                    <min>0.0</min>
                    <max>0.08</max>
                </parameter>
            </module>
        </parameter>
        <parameter measure="MEDIAN" name="m/z center calculation" weighting="NONE">CenterFunction</parameter>
        <parameter name="m/z range for MS2 scan pairing (Da)" selected="false"/>
        <parameter name="RT range for MS2 scan pairing (min)" selected="false"/>
        <parameter name="Remove original feature list">false</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.filter_duplicatefilter.DuplicateFilterModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
        <parameter name="Name suffix">filtered</parameter>
        <parameter name="Filter mode">OLD AVERAGE</parameter>
        <parameter name="m/z tolerance">
            <absolutetolerance>0.0</absolutetolerance>
            <ppmtolerance>5.0</ppmtolerance>
        </parameter>
        <parameter name="RT tolerance" type="absolute">0.05</parameter>
        <parameter name="Require same identification">false</parameter>
        <parameter name="Remove original peaklist">false</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.tools.sortpeaklists.SortPeakListsModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.align_join.JoinAlignerModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
        <parameter name="Feature list name">Aligned feature list</parameter>
        <parameter name="m/z tolerance">
            <absolutetolerance>0.0</absolutetolerance>
            <ppmtolerance>5.0</ppmtolerance>
        </parameter>
        <parameter name="Weight for m/z">3.0</parameter>
        <parameter name="Retention time tolerance" type="absolute">0.1</parameter>
        <parameter name="Weight for RT">1.0</parameter>
        <parameter name="Require same charge state">false</parameter>
        <parameter name="Require same ID">false</parameter>
        <parameter name="Compare isotope pattern" selected="false">
            <parameter name="Isotope m/z tolerance"/>
            <parameter name="Minimum absolute intensity"/>
            <parameter name="Minimum score"/>
        </parameter>
        <parameter name="Compare spectra similarity" selected="false">
            <parameter name="Mass list"/>
            <parameter name="Spectral m/z tolerance">
                <absolutetolerance>0.001</absolutetolerance>
                <ppmtolerance>10.0</ppmtolerance>
            </parameter>
            <parameter name="MS level">2</parameter>
            <parameter name="Compare spectra similarity">
                <module name="Weighted dot-product cosine">
                    <parameter name="Weights">MassBank (mz^2 * I^0.5)</parameter>
                    <parameter name="Minimum  cos similarity">0.7</parameter>
                    <parameter name="Remove unmatched signals">false</parameter>
                </module>
                <module name="Composite dot -product identity (similar to NIST search)">
                    <parameter name="Weights">MassBank (mz^2 * I^0.5)</parameter>
                    <parameter name="Minimum  cos similarity">0.7</parameter>
                    <parameter name="Remove unmatched signals">false</parameter>
                </module>
            </parameter>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.io.csvexport.CSVExportModule">
        <parameter name="Feature lists" type="BATCH_LAST_PEAKLISTS"/>
        <parameter name="Filename">
            <current_file>{output}</current_file>
            <last_file>{output}</last_file>
        </parameter>
        <parameter name="Field separator">,</parameter>
        <parameter name="Export common elements">
            <item>Export row m/z</item>
            <item>Export row retention time</item>
            <item>Export row identity (main ID)</item>
        </parameter>
        <parameter name="Export data file elements">
            <item>Peak area</item>
        </parameter>
        <parameter name="Export quantitation results and other information">false</parameter>
        <parameter name="Identification separator">;</parameter>
        <parameter name="Filter rows">ALL</parameter>
    </batchstep>
</batch>"""
    )

  try:
      opts, args = getopt.getopt(sys.argv[1:], "i:o:")
  except getopt.GetoptError as err:
    print(err)
    sys.exit(2)
  input=None
  output=None
  target=None
  print(opts)
  for o, a in opts:
    if o == "-i":
      print(a)
      input = os.path.abspath(a)
    elif o == "-o":
      print(a)
      output = os.path.abspath(a)
    else:
      print(o)
      print(a)
      print("unhandled option")
      #sys.exit(2)
  if not (input and output):
    print("missing options")
    sys.exit(2)

  mzXML=[]
  mzXML=findRaw(input, mzXML, 0)
  for sample in mzXML:
    samplepath = os.path.join(input, sample)
    xmlpath = os.path.join(output, sample[:-5] + "adap.xml")
    outputpath = os.path.join(output, sample[:-5] + "csv")
    f=open(output + "/" + sample[:-5] + "adap.xml", "w")
    f.write(batchFile.format(mzXML=samplepath, output=outputpath))
    f.close()

def findRaw(dir, mzXML, depth):
  for element in os.listdir(dir):
    elementpath = os.path.join(dir, element)
    if os.path.isfile(elementpath):
      if len(element) > 6 and element[-6:] == ".mzXML":
        mzXML.append(element)
    elif depth < 3:
      mzXML = findRaw(elementpath, mzXML, depth+1)
  return mzXML

if __name__ == "__main__":
  main()

