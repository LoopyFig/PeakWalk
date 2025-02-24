#!/usr/bin/env python3

import getopt, sys, os

batchOptions = [
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
</batch>""",
"""<?xml version="1.0" encoding="UTF-8"?><batch mzmine_version="4.4.3">
    <batchstep method="io.github.mzmine.modules.io.import_rawdata_mzxml.MzXMLImportModule" parameter_version="1">
        <parameter name="File names">
            <file>{mzXML}</file>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetectionModule" parameter_version="1">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scan filters" selected="true">
            <parameter name="Scan number"/>
            <parameter name="Base Filtering Integer"/>
            <parameter name="Retention time"/>
            <parameter name="Mobility"/>
            <parameter name="MS level filter" selected="All MS levels">1</parameter>
            <parameter name="Scan definition"/>
            <parameter name="Polarity">Any</parameter>
            <parameter name="Spectrum type">ANY</parameter>
        </parameter>
        <parameter name="Scan types (IMS)">All scan types</parameter>
        <parameter name="Denormalize fragment scans (traps)">false</parameter>
        <parameter name="Mass detector" selected_item="Auto">
            <module name="Factor of lowest signal">
                <parameter name="Noise factor">2.5</parameter>
            </module>
            <module name="Auto">
                <parameter name="Noise level">1000.0</parameter>
            </module>
            <module name="Centroid">
                <parameter name="Noise level">10000.0</parameter>
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
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_adapchromatogrambuilder.ModularADAPChromatogramBuilderModule" parameter_version="1">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scan filters" selected="true">
            <parameter name="Scan number"/>
            <parameter name="Base Filtering Integer"/>
            <parameter name="Retention time"/>
            <parameter name="Mobility"/>
            <parameter name="MS level filter" selected="MS1, level = 1">1</parameter>
            <parameter name="Scan definition"/>
            <parameter name="Polarity">Any</parameter>
            <parameter name="Spectrum type">ANY</parameter>
        </parameter>
        <parameter name="Minimum consecutive scans">3</parameter>
        <parameter name="Minimum intensity for consecutive scans">1000.0</parameter>
        <parameter name="Minimum absolute height">10000.0</parameter>
        <parameter name="m/z tolerance (scan-to-scan)">
            <absolutetolerance>0.0</absolutetolerance>
            <ppmtolerance>10.0</ppmtolerance>
        </parameter>
        <parameter name="Suffix">chromatograms</parameter>
        <parameter name="Allow single scan chromatograms"/>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_chromatogramdeconvolution.minimumsearch.MinimumSearchFeatureResolverModule" parameter_version="2">
        <parameter name="Feature lists" type="BATCH_LAST_FEATURELISTS"/>
        <parameter name="Suffix">resolved</parameter>
        <parameter name="Original feature list">REMOVE</parameter>
        <parameter name="MS/MS scan pairing" selected="false">
            <parameter name="MS1 to MS2 precursor tolerance (m/z)">
                <absolutetolerance>0.01</absolutetolerance>
                <ppmtolerance>10.0</ppmtolerance>
            </parameter>
            <parameter name="Retention time filter" selected="Use feature edges" unit="MINUTES">0.2</parameter>
            <parameter name="Minimum relative feature height" selected="false">0.25</parameter>
            <parameter name="Minimum required signals" selected="false">1</parameter>
            <parameter name="Limit by ion mobility edges">false</parameter>
            <parameter name="Merge MS/MS spectra (TIMS)">false</parameter>
            <parameter name="Minimum signal intensity (absolute, TIMS)" selected="false">250.0</parameter>
            <parameter name="Minimum signal intensity (relative, TIMS)" selected="false">0.01</parameter>
        </parameter>
        <parameter name="Dimension">Retention time</parameter>
        <parameter name="Chromatographic threshold">0.85</parameter>
        <parameter name="Minimum search range RT/Mobility (absolute)">0.05</parameter>
        <parameter name="Minimum relative height">0.1</parameter>
        <parameter name="Minimum absolute height">10000.0</parameter>
        <parameter name="Min ratio of peak top/edge">1.0</parameter>
        <parameter name="Peak duration range (min/mobility)">
            <min>0.0</min>
            <max>1.0</max>
        </parameter>
        <parameter name="Minimum scans (data points)">3</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.io.export_features_csv_legacy.LegacyCSVExportModule" parameter_version="1">
        <parameter name="Feature lists" type="BATCH_LAST_FEATURELISTS"/>
        <parameter name="Filename">
            <current_file>{output}</current_file>
            <last_file>{output}</last_file>
        </parameter>
        <parameter name="Field separator">,</parameter>
        <parameter name="Export common elements">
            <item>Export row m/z</item>
            <item>Export row retention time</item>
        </parameter>
        <parameter name="Export data file elements">
            <item>Feature duration time</item>
            <item>Peak area</item>
        </parameter>
        <parameter name="Export quantitation results and other information">false</parameter>
        <parameter name="Identification separator">;</parameter>
        <parameter name="Filter rows">ALL</parameter>
    </batchstep>
</batch>"""]


def main():
  try:
      opts, args = getopt.getopt(sys.argv[1:], "r:a:s:v:")
  except getopt.GetoptError as err:
    print(err)
    sys.exit(2)
  rawdir=None
  adapdir=None
  suffix=".mzXML"
  batchFile = batchOptions[1]
  print(opts)
  for o, a in opts:
    if o == "-r":
      print(a)
      rawdir = os.path.abspath(a)
    elif o == "-a":
      print(a)
      adapdir = os.path.abspath(a)
    elif o == "-s":
      print(a)
      suffix = "." + a
    elif o == "-v":
      if a == "2":
        print("mzmine v2")
        batchFile = batchOptions[0]
      elif a == "4":
        print("mzmine v4")
        batchFile = batchOptions[1]
    else:
      print(o)
      print(a)
      print("unhandled option")
      #sys.exit(2)
  if not (rawdir and adapdir):
    print("missing options")
    sys.exit(2)

  raws=[]
  raws=findRaw(rawdir, raws, suffix, 0)
  for sample in raws:
    samplepath = os.path.join(rawdir, sample)
    xmlpath = os.path.join(adapdir, sample[:-len(suffix)] + ".adap.xml")
    outputpath = os.path.join(adapdir, sample[:-len(suffix)] + ".csv")
    f=open(adapdir + "/" + sample[:-len(suffix)] + ".adap.xml", "w")
    f.write(batchFile.format(mzXML=samplepath, output=outputpath))
    f.close()

def findRaw(dir, raws, suffix, depth):
  for element in os.listdir(dir):
    elementpath = os.path.join(dir, element)
    if os.path.isfile(elementpath):
      if len(element) > len(suffix) and element[-len(suffix):] == suffix:
        raws.append(element)
    elif depth < 3:
      raws = findRaw(elementpath, raws, suffix, depth+1)
  return raws

if __name__ == "__main__":
  main()

