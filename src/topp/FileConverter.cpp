// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/EDTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FileConverter FileConverter

  @brief Converts between different MS file formats.

  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FileConverter \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_GenericWrapper (e.g. for calling external converters) </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> any tool operating on the output format</td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any vendor software exporting supported formats (e.g. mzXML) </td>
  </tr>
  </table>
  </CENTER>

  The main use of this tool is to convert data from external sources to the formats used by OpenMS/TOPP. 
  Maybe most importantly, data from MS experiments in a number of different formats can be converted to mzML, 
  the canonical file format used by OpenMS/TOPP for experimental data. (mzML is the PSI approved format and 
  supports traceability of analysis steps.)

  Many different format conversions are supported, and some may be more useful than others. Depending on the
  file formats involved, information can be lost during conversion, e.g. when converting	featureXML to mzData. 
  In such cases a warning is shown.

  The input and output file types are determined from	the file extensions or from the first few lines of the 
  files. If file type determination is not possible, the input or output file type has to be given explicitly.

  Conversion with the same output as input format is supported. In some cases, this can be helpful to remove
  errors from files, to update file formats to new versions, or to check whether information is lost upon 
  reading or writing.

  Some information about the supported input types:
  @ref OpenMS::MzMLFile "mzML"
  @ref OpenMS::MzXMLFile "mzXML"
  @ref OpenMS::MzDataFile "mzData"
  @ref OpenMS::MascotGenericFile "mgf"
  @ref OpenMS::DTA2DFile "dta2d"
  @ref OpenMS::DTAFile "dta"
  @ref OpenMS::FeatureXMLFile "featureXML"
  @ref OpenMS::ConsensusXMLFile "consensusXML"
  @ref OpenMS::MS2File "ms2"
  @ref OpenMS::XMassFile "fid/XMASS"
  @ref OpenMS::MsInspectFile "tsv"
  @ref OpenMS::SpecArrayFile "peplist"
  @ref OpenMS::KroenikFile "kroenik"
  @ref OpenMS::EDTAFile "edta"

  See @ref TOPP_IDFileConverter for similar functionality for protein/peptide identification file formats.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FileConverter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_FileConverter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileConverter :
  public TOPPBase
{
public:
  TOPPFileConverter() :
    TOPPBase("FileConverter", "Converts between different MS file formats.")
  {

  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file ");
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    String formats("mzData,mzXML,mzML,dta,dta2d,mgf,featureXML,consensusXML,ms2,fid,tsv,peplist,kroenik,edta");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));

    formats = "mzData,mzXML,mzML,dta2d,mgf,featureXML,consensusXML,edta";
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>(formats));
    registerStringOption_("out_type", "<type>", "", "output file type -- default: determined from file extension or content\n", false);
    setValidStrings_("out_type", ListUtils::create<String>(formats));
    registerFlag_("TIC_DTA2D", "Export the TIC instead of the entire experiment in mzML/mzData/mzXML -> DTA2D conversions.", true);
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input file names
    String in = getStringOption_("in");

    //input file type
    FileHandler fh;
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }


    //output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    bool TIC_DTA2D = getFlag_("TIC_DTA2D");

    writeDebug_(String("Output file type: ") + FileTypes::typeToName(out_type), 1);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    typedef MSExperiment<Peak1D> MSExperimentType;
    MSExperimentType exp;

    typedef MSExperimentType::SpectrumType SpectrumType;

    typedef FeatureMap<> FeatureMapType;

    FeatureMapType fm;
    ConsensusMap cm;

    writeDebug_(String("Loading input file"), 1);

    if (in_type == FileTypes::CONSENSUSXML)
    {
      ConsensusXMLFile().load(in, cm);
      cm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML))
      {
        // You you will lose information and waste memory. Enough reasons to issue a warning!
        writeLog_("Warning: Converting consensus features to peaks. You will lose information!");
        exp.set2DData(cm);
      }
    }
    else if (in_type == FileTypes::EDTA)
    {
      EDTAFile().load(in, cm);
      cm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML))
      {
        // You you will lose information and waste memory. Enough reasons to issue a warning!
        writeLog_("Warning: Converting consensus features to peaks. You will lose information!");
        exp.set2DData(cm);
      }
    }
    else if (in_type == FileTypes::FEATUREXML ||
             in_type == FileTypes::TSV ||
             in_type == FileTypes::PEPLIST ||
             in_type == FileTypes::KROENIK)
    {
      fh.loadFeatures(in, fm, in_type);
      fm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML))
      {
        // You will lose information and waste memory. Enough reasons to issue a warning!
        writeLog_("Warning: Converting features to peaks. You will lose information! Mass traces are added, if present as 'num_of_masstraces' and 'masstrace_intensity_<X>' (X>=0) meta values.");
        exp.set2DData<true>(fm);
      }
    }
    else
    {
      fh.loadExperiment(in, exp, in_type, log_type_);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    writeDebug_(String("Writing output file"), 1);

    if (out_type == FileTypes::MZML)
    {
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZML));
      MzMLFile f;
      f.setLogType(log_type_);
      ChromatogramTools().convertSpectraToChromatograms(exp, true);
      f.store(out, exp);
    }
    else if (out_type == FileTypes::MZDATA)
    {
      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZDATA));
      MzDataFile f;
      f.setLogType(log_type_);
      ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
      f.store(out, exp);
    }
    else if (out_type == FileTypes::MZXML)
    {
      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZXML));
      MzXMLFile f;
      f.setLogType(log_type_);
      ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
      f.store(out, exp);
    }
    else if (out_type == FileTypes::DTA2D)
    {
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 FORMAT_CONVERSION));
      DTA2DFile f;
      f.setLogType(log_type_);
      ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
      if (TIC_DTA2D)
      {
        // store the total ion chromatogram (TIC)
        f.storeTIC(out, exp);
      }
      else
      {
        // store entire experiment
        f.store(out, exp);
      }


    }
    else if (out_type == FileTypes::MGF)
    {
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 FORMAT_CONVERSION));
      MascotGenericFile f;
      f.setLogType(log_type_);
      f.store(out, exp);
    }
    else if (out_type == FileTypes::FEATUREXML)
    {
      if ((in_type == FileTypes::FEATUREXML) || (in_type == FileTypes::TSV) ||
          (in_type == FileTypes::PEPLIST) || (in_type == FileTypes::KROENIK))
      {
        fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      }
      else if (in_type == FileTypes::CONSENSUSXML || in_type == FileTypes::EDTA)
      {
        ConsensusMap::convert(cm, true, fm);
      }
      else // not loaded as feature map or consensus map
      {
        // The feature specific information is only defaulted. Enough reasons to issue a warning!
        writeLog_("Warning: Converting peaks to features will lead to incomplete features!");
        fm.clear();
        fm.reserve(exp.getSize());
        typedef FeatureMapType::FeatureType FeatureType;
        FeatureType feature;
        feature.setQuality(0, 1); // override default
        feature.setQuality(1, 1); // override default
        feature.setOverallQuality(1); // override default
        for (MSExperimentType::ConstIterator spec_iter = exp.begin();
             spec_iter != exp.end();
             ++spec_iter
             )
        {
          feature.setRT(spec_iter->getRT());
          for (SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
               peak1_iter != spec_iter->end();
               ++peak1_iter
               )
          {
            feature.setMZ(peak1_iter->getMZ());
            feature.setIntensity(peak1_iter->getIntensity());
            feature.setUniqueId();
            fm.push_back(feature);
          }
        }
        fm.updateRanges();
      }

      addDataProcessing_(fm, getProcessingInfo_(DataProcessing::
                                                FORMAT_CONVERSION));
      FeatureXMLFile().store(out, fm);
    }
    else if (out_type == FileTypes::CONSENSUSXML)
    {
      if ((in_type == FileTypes::FEATUREXML) || (in_type == FileTypes::TSV) ||
          (in_type == FileTypes::PEPLIST) || (in_type == FileTypes::KROENIK))
      {
        fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
        ConsensusMap::convert(0, fm, cm);
      }
      // nothing to do for consensus input
      else if (in_type == FileTypes::CONSENSUSXML || in_type == FileTypes::EDTA)
      {
      }
      else // experimental data
      {
        ConsensusMap::convert(0, exp, cm, exp.size());
      }

      addDataProcessing_(cm, getProcessingInfo_(DataProcessing::
                                                FORMAT_CONVERSION));
      ConsensusXMLFile().store(out, cm);
    }
    else if (out_type == FileTypes::EDTA)
    {
      if (fm.size() > 0 && cm.size() > 0)
      {
        LOG_ERROR << "Internal error: cannot decide on container (Consensus or Feature)! This is a bug. Please report it!";
        return INTERNAL_ERROR;
      }
      if (fm.size() > 0) EDTAFile().store(out, fm);
      else if (cm.size() > 0) EDTAFile().store(out, cm);
    }
    else
    {
      writeLog_("Unknown output file type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFileConverter tool;
  return tool.main(argc, argv);
}

/// @endcond