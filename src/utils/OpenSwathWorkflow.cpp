// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

// Files
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathTSVWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>

#include <cassert>
#include <limits>

// #define OPENSWATH_WORKFLOW_DEBUG

using namespace OpenMS;

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenSwathWorkflow OpenSwathWorkflow

  @brief Complete workflow to run OpenSWATH

  This implements the OpenSwath workflow as described in Rost and Rosenberger
  et al. (Nature Biotechnology, 2014) and provides a complete, integrated
  analysis tool without the need to run multiple tools consecutively.

  It executes the following steps in order:

  <ul>
    <li>Reading of input files, which can be provided as one single mzML or multiple "split" mzMLs (one per SWATH)</li>
    <li>Computing the retention time transformation using RT-normalization peptides</li>
    <li>Reading of the transition list</li>
    <li>Extracting the specified transitions</li>
    <li>Scoring the peak groups in the extracted ion chromatograms (XIC)</li>
    <li>Reporting the peak groups and the chromatograms</li>
  </ul>

  The overall execution flow for this tool is described in the TOPPOpenSwathWorkflow documentation.

  See below or have a look at the INI file (via "OpenSwathWorkflow -write_ini myini.ini") for available parameters and more functionality.

  <h3>Input: SWATH maps and transition list </h3>
  SWATH maps can be provided as mzML files, either as single file directly from
  the machine (this assumes that the SWATH method has 1 MS1 and then n MS2
  spectra which are ordered the same way for each cycle). E.g. a valid method
  would be MS1, MS2 [400-425], MS2 [425-450], MS1, MS2 [400-425], MS2 [425-450]
  while an invalid method would be MS1, MS2 [400-425], MS2 [425-450], MS1, MS2
  [425-450], MS2 [400-425] where MS2 [xx-yy] indicates an MS2 scan with an
  isolation window starting at xx and ending at yy. OpenSwathWorkflow will try
  to read the SWATH windows from the data, if this is not possible please
  provide a tab-separated list with the correct windows using the
  -swath_windows_file parameter (this is recommended). Note that the software
  expects extraction windows (e.g. which peptides to extract from
  which window) which cannot have overlaps, otherwise peptides will be
  extracted from two different windows.

  Alternatively, a set of split files (n+1 mzML files) can be provided, each
  containing one SWATH map (or MS1 map).

  Since the file size can become rather large, it is recommended to not load the
  whole file into memory but rather cache it somewhere on the disk using a
  fast-access data format. This can be specified using the -readOptions cache
  parameter (this is recommended!).

  <h3>Parameters</h3>
  The current parameters are optimized for 2 hour gradients on SCIEX 5600 /
  6600 TripleTOF instruments with a peak width of around 30 seconds using iRT
  peptides.  If your chromatography differs, please consider adjusting
  -Scoring:TransitionGroupPicker:min_peak_width  to allow for smaller or larger
  peaks and adjust the -rt_extraction_window to use a different extraction
  window for the retention time. In m/z domain, consider adjusting
  -mz_extraction_window to your instrument resolution, which can be in Th or
  ppm (using -ppm).

  Furthermore, if you wish to use MS1 information, use the -use_ms1_traces flag
  and provide an MS1 map in addition to the SWATH data.

  If you encounter issues with peak picking, try to disable peak filtering by
  setting -Scoring:TransitionGroupPicker:compute_peak_quality false which will
  disable the filtering of peaks by chromatographic quality. Furthermore, you
  can adjust the smoothing parameters for the peak picking, by adjusting
  -Scoring:TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length or using a
  Gaussian smoothing based on your estimated peak width. Adjusting the signal
  to noise threshold will make the peaks wider or smaller.

  <h3>Output: Feature list and chromatograms </h3>
  The output of the OpenSwathWorkflow is a feature list, either as FeatureXML
  or as tsv (use -out_features or -out_tsv) while the latter is more memory
  friendly. If you analyze large datasets, it is recommended to only use
  -out_tsv and not -out_features. For downstream analysis (e.g. using mProphet or pyProphet)
  also the -out_tsv format is recommended.

  The feature list generated by -out_tsv is a tab-separated file. It can be
  used directly as input to the mProphet or pyProphet (a Python
  re-implementation of mProphet) software tool, see Reiter et al (2011, Nature
  Methods).

  In addition, the extracted chromatograms can be written out using the
  -out_chrom parameter.

  <h4> Feature list output format </h4>

  The tab-separated feature output contains the following information:

<CENTER>
  <table>
    <tr>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Header row </td>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Format </td>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Description </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> transition_group_id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique id for the transition group (all chromatographic traces that are analyzed together)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> peptide_group_label </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique id for the peptide group (will be the same for each charge state and heavy/light status) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> run_id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> An identifier for the run (currently always 0)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> filename </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The input filename </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> RT </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peak group retention time </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique identifier for the peak group</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Sequence </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peptide sequence (no modifications) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> FullPeptideName </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Full peptide sequence including modifications in Unimod format</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Charge </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Int </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Assumed charge state</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> m/z </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Precursor m/z</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peak group intensity (sum of all transitions)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> ProteinName </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Name of the associated protein</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> decoy </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Whether the transition is decoy or not (0 = false, 1 = true) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> assay_rt </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The expected RT in seconds (based on normalized iRT value) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> delta_rt </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The difference between the expected RT and the peak group RT in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> leftWidth </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The start of the peak group (left side) in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> main_var_xx_swath_prelim_score </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Initial score </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> norm_RT </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The peak group retention time in normalized (iRT) space </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> nr_peaks </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Int </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The number of transitions used </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> peak_apices_sum </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The sum of all peak apices (may be used as alternative intensity) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> potentialOutlier </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Potential outlier transitions (or "none" if none was detected)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> rightWidth </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The end of the peak group (left side) in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> rt_score </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The raw RT score (unnormalized) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> sn_ratio </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The raw S/N ratio </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> total_xic </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The total XIC of the chromatogram </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> var_... </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> One of multiple sub-scores used by OpenSWATH to describe the peak group </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Peak_Area </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak area) of MS1 traces separated by semicolon </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Peak_Apex </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak apex) of MS1 traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Fragment_Annotation </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Annotation of MS1 traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Peak_Area </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak area) of fragment ion traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Peak_Apex </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak apex) of fragment ion traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Fragment_Annotation </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=2> Annotation of fragment ion traces separated by semicolon </td>
    </tr>


  </table>
</CENTER>

  <h3>Execution flow:</h3>

  The overall execution flow for this tool is described in the TOPPOpenSwathWorkflow documentation.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathWorkflow.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathWorkflow.html

*/

/** @brief Extended documentation on OpenSwath

  The overall execution flow for this tool is as follows:

    - Parameter validation
    - Transition loading: loads input transitions into OpenSwath::LightTargetedExperiment
    - SWATH file loading:
      - Load SWATH files (see loadSwathFiles())
      - Annotate SWATH-files with user-defined windows (see OpenMS::SwathWindowLoader::annotateSwathMapsFromFile() )
      - Sanity check: there should be no overlap between the windows:
    - Perform RT and m/z calibration (see performCalibration() and OpenMS::OpenSwathRetentionTimeNormalization)
    - Set up chromatogram file output
    - Set up peakgroup file output
    - Extract and score (see OpenMS::OpenSwathWorkflow or OpenMS::OpenSwathWorkflowSonar)

*/
class TOPPOpenSwathWorkflow
  : public TOPPOpenSwathBase 
{
public:

  TOPPOpenSwathWorkflow()
    : TOPPOpenSwathBase("OpenSwathWorkflow", "Complete workflow to run OpenSWATH", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML,sqMass"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML','tsv','pqp')");
    setValidFormats_("tr", ListUtils::create<String>("traML,tsv,pqp"));
    registerStringOption_("tr_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    setValidStrings_("tr_type", ListUtils::create<String>("traML,tsv,pqp"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt", "<file>", "", "transition file ('TraML')", false);
    setValidFormats_("tr_irt", ListUtils::create<String>("traML,tsv,pqp"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt_nonlinear", "<file>", "", "transition file ('TraML')", false);
    setValidFormats_("tr_irt_nonlinear", ListUtils::create<String>("traML,tsv,pqp"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.", false, true);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerInputFile_("swath_windows_file", "<file>", "", "Optional, tab separated file containing the SWATH windows for extraction: lower_offset upper_offset \\newline 400 425 \\newline ... Note that the first line is a header and will be skipped.", false, true);
    registerFlag_("sort_swath_maps", "Sort input SWATH files when matching to SWATH windows from swath_windows_file", true);

    registerFlag_("use_ms1_traces", "Extract the precursor ion trace(s) and use for scoring", true);
    registerFlag_("enable_uis_scoring", "Enable additional scoring of identification assays", true);

    // one of the following two needs to be set
    registerOutputFile_("out_features", "<file>", "", "output file", false);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_tsv", "<file>", "", "TSV output file (mProphet compatible TSV file)", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_osw", "<file>", "", "OSW output file (PyProphet compatible SQLite file)", false);
    setValidFormats_("out_osw", ListUtils::create<String>("osw"));

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms output in mzML (chrom.mzML) or sqMass (SQLite format)", false, true);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML,sqMass"));

    // misc options
    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the edge to still consider a precursor, in Thomson", false, true);
    registerFlag_("sonar", "data is scanning SWATH data");

    // RT, mz and IM windows
    registerDoubleOption_("rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).", false);
    registerDoubleOption_("extra_rt_extraction_window", "<double>", 0.0, "Output an XIC with a RT-window that by this much larger (e.g. to visually inspect a larger area of the chromatogram)", false, true);
    setMinFloat_("extra_rt_extraction_window", 0.0);
    registerDoubleOption_("ion_mobility_window", "<double>", -1, "Extraction window in ion mobility dimension (in milliseconds). This is the full window size, e.g. a value of 10 milliseconds would extract 5 milliseconds on either side.", false);
    registerDoubleOption_("mz_extraction_window", "<double>", 0.05, "Extraction window used (in Thomson, to use ppm see -ppm flag)", false);
    setMinFloat_("mz_extraction_window", 0.0);
    registerStringOption_("mz_extraction_window_unit", "<name>", "Th", "Unit for mz extraction", false, true);
    setValidStrings_("mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));

    // MS1 mz windows and ion mobility
    registerDoubleOption_("mz_extraction_window_ms1", "<double>", 0.05, "Extraction window used in MS1 (in ppm)", false);
    setMinFloat_("mz_extraction_window_ms1", 0.0);
    registerStringOption_("mz_extraction_window_ms1_unit", "<name>", "Th", "Unit of the MS1 m/z extraction window", false, true);
    setValidStrings_("mz_extraction_window_ms1_unit", ListUtils::create<String>("ppm,Th"));
    registerDoubleOption_("im_extraction_window_ms1", "<double>", -1, "Extraction window in ion mobility dimension for MS1 (in milliseconds).", false);

    registerStringOption_("use_ms1_ion_mobility", "<name>", "true", "Also perform precursor extraction using the same ion mobility window as for fragment ion extraction", false, true);
    setValidStrings_("use_ms1_ion_mobility", ListUtils::create<String>("true,false"));

    // iRT mz and IM windows
    registerDoubleOption_("irt_mz_extraction_window", "<double>", 0.05, "Extraction window used for iRT and m/z correction (in Thomson, use ppm use -ppm flag)", false, true);
    setMinFloat_("irt_mz_extraction_window", 0.0);
    registerDoubleOption_("irt_im_extraction_window", "<double>", -1, "Ion mobility extraction window used for iRT (in 1/K0 or milliseconds)", false, true);
    registerStringOption_("irt_mz_extraction_window_unit", "<name>", "Th", "Unit for mz extraction", false, true);
    setValidStrings_("irt_mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));


    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false, true);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false, true);

    registerFlag_("split_file_input", "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)", true);
    registerFlag_("use_elution_model_score", "Turn on elution model score (EMG fit to peak)", true);

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first. If you choose cache, make sure to also set tempDirectory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("normal,cache,cacheWorkingInMemory,workingInMemory"));

    registerStringOption_("mz_correction_function", "<name>", "none", "Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.", false, true);
    setValidStrings_("mz_correction_function", ListUtils::create<String>("none,regression_delta_ppm,unweighted_regression,weighted_regression,quadratic_regression,weighted_quadratic_regression,weighted_quadratic_regression_delta_ppm,quadratic_regression_delta_ppm"));

    registerStringOption_("tempDirectory", "<tmp>", "/tmp/", "Temporary directory to store cached files for example", false, true);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", ListUtils::create<String>("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 250, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 250-1000)", false, true);
    setMinInt_("batchSize", 0);
    registerIntOption_("outer_loop_threads", "<number>", -1, "How many threads should be used for the outer loop (-1 use all threads, use 4 to analyze 4 SWATH windows in memory at once).", false, true);

    registerIntOption_("ms1_isotopes", "<number>", 0, "The number of MS1 isotopes used for extraction", false, true);
    setMinInt_("ms1_isotopes", 0);

    registerIntOption_("min_ms1_chromatograms", "<number>", 1, "The minimal number of MS1 isotopes (including monoisotopic peak) required for scoring", false, true);
    setMinInt_("min_ms1_chromatograms", 1);

    registerIntOption_("min_transitions", "<number>", 6, "Minimal number of transitions used for scoring", false, true);
    setMinInt_("min_transitions", 3);

        registerIntOption_("max_transitions", "<number>", 6, "Maximum number of transitions used for scoring", false, true);
    setMinInt_("max_transitions", 3);

    registerSubsection_("Scoring", "Scoring parameters section");
    registerSubsection_("Library", "Library parameters section");

    registerSubsection_("RTNormalization", "Parameters for the RTNormalization for iRT petides. This specifies how the RT alignment is performed and how outlier detection is applied. Outlier detection can be done iteratively (by default) which removes one outlier per iteration or using the RANSAC algorithm.");
    registerSubsection_("Debugging", "Debugging");
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.remove("rt_extraction_window");
      feature_finder_param.setValue("stop_report_after_feature", 5);
      feature_finder_param.setValue("rt_normalization_factor", 100.0); // for iRT peptides between 0 and 100 (more or less)

      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:remove_overlapping_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:write_sn_log_messages", "false"); // no log messages
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerMRM:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width", 30.0);
      feature_finder_param.setValue("uis_threshold_sn",0);
      feature_finder_param.setValue("uis_threshold_peak_area",0);
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_bin_count");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:stop_after_feature");

      // EMG Scoring - turn off by default since it is very CPU-intensive
      feature_finder_param.remove("Scores:use_elution_model_score");
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");

      // remove these parameters
      feature_finder_param.remove("add_up_spectra");
      feature_finder_param.remove("spacing_for_spectra_resampling");
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
    else if (name == "RTNormalization")
    {
      Param p;

      p.setValue("alignmentMethod", "linear", "How to perform the alignment to the normalized RT space using anchor points. 'linear': perform linear regression (for few anchor points). 'interpolated': Interpolate between anchor points (for few, noise-free anchor points). 'lowess' Use local regression (for many, noisy anchor points). 'b_spline' use b splines for smoothing.");
      p.setValidStrings("alignmentMethod", ListUtils::create<String>("linear,interpolated,lowess,b_spline"));
      p.setValue("lowess:span", 2.0/3, "Span parameter for lowess");
      p.setMinFloat("lowess:span", 0.0);
      p.setMaxFloat("lowess:span", 1.0);
      p.setValue("b_spline:num_nodes", 5, "Number of nodes for b spline");
      p.setMinInt("b_spline:num_nodes", 0);

      p.setValue("outlierMethod", "iter_residual", "Which outlier detection method to use (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
      p.setValidStrings("outlierMethod", ListUtils::create<String>("iter_residual,iter_jackknife,ransac,none"));

      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.");
      p.setValidStrings("useIterativeChauvenet", ListUtils::create<String>("true,false"));

      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");

      p.setValue("estimateBestPeptides", "false", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).");
      p.setValidStrings("estimateBestPeptides", ListUtils::create<String>("true,false"));

      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else if (name == "Debugging")
    {
      Param p;
      p.setValue("irt_mzml", "", "Chromatogram mzML containing the iRT peptides");
      // p.setValidFormats_("irt_mzml", ListUtils::create<String>("mzML"));
      p.setValue("irt_trafo", "", "Transformation file for RT transform");
      // p.setValidFormats_("irt_trafo", ListUtils::create<String>("trafoXML"));
      return p;
    }
    else if (name == "Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
    }
  }

  ExitCodes main_(int, const char **) override
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");

    Param irt_detection_param = getParam_().copy("RTNormalization:", true);

    //tr_file input file type
    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));
    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = FileHandler::getType(tr_file);
      writeDebug_(String("Input file type (-tr): ") + FileTypes::typeToName(tr_type), 2);
    }

    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type for '-tr' !");
      return PARSE_ERROR;
    }

    String out = getStringOption_("out_features");
    String out_tsv = getStringOption_("out_tsv");
    String out_osw = getStringOption_("out_osw");

    String irt_tr_file = getStringOption_("tr_irt");
    String nonlinear_irt_tr_file = getStringOption_("tr_irt_nonlinear");
    String trafo_in = getStringOption_("rt_norm");
    String swath_windows_file = getStringOption_("swath_windows_file");

    String out_chrom = getStringOption_("out_chrom");
    bool split_file = getFlag_("split_file_input");
    bool use_emg_score = getFlag_("use_elution_model_score");
    bool force = getFlag_("force");
    bool sonar = getFlag_("sonar");
    bool sort_swath_maps = getFlag_("sort_swath_maps");
    bool use_ms1_traces = getFlag_("use_ms1_traces");
    bool enable_uis_scoring = getFlag_("enable_uis_scoring");
    int batchSize = (int)getIntOption_("batchSize");
    int outer_loop_threads = (int)getIntOption_("outer_loop_threads");
    int ms1_isotopes = (int)getIntOption_("ms1_isotopes");
    int min_ms1_chromatograms = (int)getIntOption_("min_ms1_chromatograms");
    int min_transitions = (int)getIntOption_("min_transitions");
    int max_transitions = (int)getIntOption_("max_transitions");
    Size debug_level = (Size)getIntOption_("debug");

    double min_rsq = getDoubleOption_("min_rsq");
    double min_coverage = getDoubleOption_("min_coverage");

    Param debug_params = getParam_().copy("Debugging:", true);

    String readoptions = getStringOption_("readOptions");
    String mz_correction_function = getStringOption_("mz_correction_function");
    String tmp = getStringOption_("tempDirectory");

    ///////////////////////////////////
    // Parameter validation
    ///////////////////////////////////

    bool load_into_memory = false;
    if (readoptions == "cacheWorkingInMemory")
    {
      readoptions = "cache";
      load_into_memory = true;
    }
    else if (readoptions == "workingInMemory")
    {
      readoptions = "normal";
      load_into_memory = true;
    }

    bool is_sqmass_input  = (FileHandler::getTypeByFileName(file_list[0]) == FileTypes::SQMASS);
    if (is_sqmass_input && !load_into_memory)
    {
      std::cout << "When using sqMass input files, it is highly recommended to use the workingInMemory option as otherwise data access will be very slow." << std::endl;
    }

    if (trafo_in.empty() && irt_tr_file.empty())
    {
      std::cout << "Since neither rt_norm nor tr_irt is set, OpenSWATH will " <<
        "not use RT-transformation (rather a null transformation will be applied)" << std::endl;
    }
    if ( int(!out.empty()) + int(!out_tsv.empty()) + int(!out_osw.empty()) != 1 )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Either out_features, out_tsv or out_osw needs to be set (but not two or three at the same time)");
    }
    if (!out_osw.empty() && tr_type != FileTypes::PQP)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "OSW output files can only be generated in combination with PQP input files (-tr).");
    }

    // Check swath window input
    if (!swath_windows_file.empty())
    {
      LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      LOG_INFO << "Read Swath maps file with " << swath_prec_lower.size() << " windows." << std::endl;
      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    bool use_ms1_im = getStringOption_("use_ms1_ion_mobility") == "true";

    ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.mz_extraction_window  = getDoubleOption_("mz_extraction_window");
    cp.ppm                   = getStringOption_("mz_extraction_window_unit") == "ppm";
    cp.rt_extraction_window  = getDoubleOption_("rt_extraction_window");
    cp.im_extraction_window  = getDoubleOption_("ion_mobility_window");
    cp.extraction_function   = getStringOption_("extraction_function");
    cp.extra_rt_extract      = getDoubleOption_("extra_rt_extraction_window");

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range for iRT measurements
    cp_irt.mz_extraction_window = getDoubleOption_("irt_mz_extraction_window");
    cp_irt.im_extraction_window = getDoubleOption_("irt_im_extraction_window");
    cp_irt.ppm                  = getStringOption_("irt_mz_extraction_window_unit") == "ppm";

    ChromExtractParams cp_ms1 = cp;
    cp_ms1.mz_extraction_window  = getDoubleOption_("mz_extraction_window_ms1");
    cp_ms1.ppm                   = getStringOption_("mz_extraction_window_ms1_unit") == "ppm";
    cp_ms1.im_extraction_window  = getDoubleOption_("im_extraction_window_ms1");

    Param feature_finder_param = getParam_().copy("Scoring:", true);
    Param tsv_reader_param = getParam_().copy("Library:", true);
    if (use_emg_score)
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "true");
    }
    else
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    }
    if (use_ms1_traces)
    {
      feature_finder_param.setValue("Scores:use_ms1_correlation", "true");
      feature_finder_param.setValue("Scores:use_ms1_fullscan", "true");
    }
    if (enable_uis_scoring)
    {
      feature_finder_param.setValue("Scores:use_uis_scores", "true");
    }

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp = loadTransitionList(tr_type, tr_file, tsv_reader_param);
    LOG_INFO << "Loaded " << transition_exp.getProteins().size() << " proteins, " <<
      transition_exp.getCompounds().size() << " compounds with " << transition_exp.getTransitions().size() << " transitions." << std::endl;

    if (tr_type == FileTypes::PQP)
    {
      remove(out_osw.c_str());
      if (!out_osw.empty())
      {
        std::ifstream  src(tr_file.c_str(), std::ios::binary);
        std::ofstream  dst(out_osw.c_str(), std::ios::binary);

        dst << src.rdbuf();
      }
    }

    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////
    boost::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;
    if (!loadSwathFiles(file_list, exp_meta, swath_maps, split_file, tmp, readoptions, 
                        swath_windows_file, min_upper_edge_dist, force,
                        sort_swath_maps, sonar))
    {
      return PARSE_ERROR;
    }

    ///////////////////////////////////
    // Get the transformation information (using iRT peptides)
    ///////////////////////////////////
    String irt_trafo_out = debug_params.getValue("irt_trafo");
    String irt_mzml_out = debug_params.getValue("irt_mzml");
    TransformationDescription trafo_rtnorm;
    if (nonlinear_irt_tr_file.empty())
    {
      trafo_rtnorm = performCalibration(trafo_in, irt_tr_file, swath_maps,
                                        min_rsq, min_coverage, feature_finder_param,
                                        cp_irt, irt_detection_param, mz_correction_function,
                                        debug_level, sonar, load_into_memory,
                                        irt_trafo_out, irt_mzml_out);
    }
    else
    {
      ///////////////////////////////////
      // First perform a simple linear transform, then do a second, nonlinear one
      ///////////////////////////////////

      Param linear_irt = irt_detection_param;
      linear_irt.setValue("alignmentMethod", "linear");
      trafo_rtnorm = performCalibration(trafo_in, irt_tr_file, swath_maps,
                                        min_rsq, min_coverage, feature_finder_param,
                                        cp_irt, linear_irt, "none",
                                        debug_level, sonar, load_into_memory,
                                        irt_trafo_out, irt_mzml_out);

      cp_irt.rt_extraction_window = 900; // extract some substantial part of the RT range (should be covered by linear correction)
      cp_irt.rt_extraction_window = 600; // extract some substantial part of the RT range (should be covered by linear correction)

      ///////////////////////////////////
      // Get the secondary transformation (nonlinear)
      ///////////////////////////////////
      OpenSwath::LightTargetedExperiment transition_exp_nl;
      transition_exp_nl = loadTransitionList(FileHandler::getType(nonlinear_irt_tr_file), nonlinear_irt_tr_file, tsv_reader_param);

      std::vector< OpenMS::MSChromatogram > chromatograms;
      OpenSwathCalibrationWorkflow wf;
      wf.setLogType(log_type_);
      wf.simpleExtractChromatograms_(swath_maps, transition_exp_nl, chromatograms,
                                    trafo_rtnorm, cp_irt, sonar, load_into_memory);

      trafo_rtnorm = wf.doDataNormalization_(transition_exp_nl, chromatograms, min_rsq,
                                        min_coverage, feature_finder_param, irt_detection_param,
                                        swath_maps, mz_correction_function,
                                        cp_irt.mz_extraction_window, cp_irt.ppm);

    }

    ///////////////////////////////////
    // Set up chromatogram output
    // Either use chrom.mzML or sqlite DB
    ///////////////////////////////////
    Interfaces::IMSDataConsumer * chromatogramConsumer;
    prepareChromOutput(&chromatogramConsumer, exp_meta, transition_exp, out_chrom);

    ///////////////////////////////////
    // Set up peakgroup file output
    ///////////////////////////////////
    FeatureMap out_featureFile;
    OpenSwathTSVWriter tsvwriter(out_tsv, file_list[0], use_ms1_traces, sonar, enable_uis_scoring); // only active if filename not empty
    OpenSwathOSWWriter oswwriter(out_osw, file_list[0], use_ms1_traces, sonar, enable_uis_scoring); // only active if filename not empty

    ///////////////////////////////////
    // Extract and score
    ///////////////////////////////////
    if (sonar)
    {
      OpenSwathWorkflowSonar wf(use_ms1_traces);
      wf.setLogType(log_type_);
      wf.performExtractionSonar(swath_maps, trafo_rtnorm, cp, cp_ms1, feature_finder_param, transition_exp,
          out_featureFile, !out.empty(), tsvwriter, oswwriter, chromatogramConsumer, batchSize, load_into_memory);
    }
    else
    {
      OpenSwathWorkflow wf(use_ms1_traces, use_ms1_im, outer_loop_threads);
      wf.setLogType(log_type_);
      wf.performExtraction(swath_maps, trafo_rtnorm, cp, cp_ms1, feature_finder_param, transition_exp,
          out_featureFile, !out.empty(), tsvwriter, oswwriter, chromatogramConsumer, batchSize, ms1_isotopes, min_ms1_chromatograms,
          max_transitions, min_transitions, load_into_memory);
    }

    if (!out.empty())
    {
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FeatureXMLFile().store(out, out_featureFile);
    }

    delete chromatogramConsumer;

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathWorkflow tool;
  return tool.main(argc, argv);
}

/// @endcond
