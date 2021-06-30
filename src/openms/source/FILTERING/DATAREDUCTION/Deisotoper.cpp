// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------


#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>

namespace OpenMS
{
// static
void Deisotoper::deisotopeWithAveragineModel(MSSpectrum& spec,
  double fragment_tolerance,
  bool fragment_unit_ppm,
  bool rem_low_intensity,
  int min_charge,
  int max_charge,
  bool keep_only_deisotoped,
  unsigned int min_isopeaks,
  unsigned int max_isopeaks,
  bool make_single_charged,
  bool use_averagine_model,
  bool annotate_charge,
  bool annotate_iso_peak_count,
  bool add_up_intensity,
  bool used_for_open_search)
{
  OPENMS_PRECONDITION(spec.isSorted(), "Spectrum must be sorted.");

  if (min_isopeaks < 2 || max_isopeaks < 2 || min_isopeaks > max_isopeaks)
  {
    throw Exception::IllegalArgument(__FILE__,
      __LINE__,
      OPENMS_PRETTY_FUNCTION,
      "Minimum/maximum number of isotopic peaks must be at least 2 (and min_isopeaks <= max_isopeaks).");
  }

  if (spec.empty()) { return; }

  // discard low-intensity peaks
  if (rem_low_intensity)
  { 
    UInt max_num_peaks = used_for_open_search ? 1000 : 5000;

    // remove 0 intensity peaks
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakSpectrum(spec);

    // only keep max_num_peaks highest peaks
    NLargest nlargest_filter = NLargest(max_num_peaks);
    nlargest_filter.filterPeakSpectrum(spec);
    
    spec.sortByPosition();
  }

  Size charge_index(0);
  Size iso_peak_count_index(0);

  // reserve integer data array to store charge of peaks
  if (annotate_charge)
  {
    // expand to hold one additional integer data array to hold the charge
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("charge");
    charge_index = spec.getIntegerDataArrays().size() - 1;
  }
  // reserve integer data array to store number of isotopic peaks for each isotopic pattern
  if (annotate_iso_peak_count)
  {
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("iso_peak_count");
    iso_peak_count_index = spec.getIntegerDataArrays().size() - 1;
  }

  // during discovery phase, work on a constant reference (just to make sure we do not modify spec)
  const MSSpectrum& old_spectrum = spec;

  std::vector<size_t> extensions;
  std::vector<std::vector<size_t>> clusters;
  std::vector<int> charges_of_extensions;
  std::vector<int> KL_of_extensions;
  // threshold for the averagine check with i peaks is in averagine_check_threshold[i] (0 and 1 peak are not checked)
  const float averagine_check_threshold[7] = {0.0f, 0.0f, 0.05f, 0.1f, 0.2f, 0.4f, 0.6f};

  std::vector<std::vector<size_t>> best_clusters;
  std::vector<double> best_clusters_KL;
  std::vector<int> charges_of_best_clusters;


  bool has_precursor_data(false);
  double precursor_mass(0);
  if (old_spectrum.getPrecursors().size() == 1)
  {
    has_precursor_data = true;
    int precursor_charge = old_spectrum.getPrecursors()[0].getCharge();
    precursor_mass = (old_spectrum.getPrecursors()[0].getMZ() * precursor_charge) - (Constants::PROTON_MASS * precursor_charge);
  }

  for (size_t current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
  {
    // Monoisotopic peaks with intensity 0.0 interfere with averagine check when normalizing the spectrum peaks height
    if (use_averagine_model && !rem_low_intensity && old_spectrum[current_peak].getIntensity() == 0.0) continue;

    const double current_mz = old_spectrum[current_peak].getMZ();
    clusters.clear();
    charges_of_extensions.clear();
    for (int q = max_charge; q >= min_charge; --q) // important: test charge hypothesis from high to low
    {      
      // try to extend isotopes from mono-isotopic peak
      // if appropriate extension larger than min_isopeaks possible:
      //   - save charge q in mono_isotopic_peak[]
      //   - annotate_charge all isotopic peaks with feature number

      bool has_min_isopeaks = true;
      const double tolerance_dalton = fragment_unit_ppm ? Math::ppmToMass(fragment_tolerance, current_mz) : fragment_tolerance;

      // do not bother testing charges q (and masses m) with: m/q > precursor_mass/q (or m > precursor_mass)
      if (has_precursor_data)
      {
        double current_theo_mass = (current_mz * q) - (Constants::PROTON_MASS * q);
        if (current_theo_mass > (precursor_mass + tolerance_dalton))
        {
          continue;
        }
      }

      extensions.clear();
      extensions.push_back(current_peak);

      float total_weight = old_spectrum[current_peak].getIntensity() * (old_spectrum[current_peak].getMZ() * q - (q - 1) * Constants::PROTON_MASS_U);
      float total_intensity = old_spectrum[current_peak].getIntensity();
      float current_KL = 0.0;
      for (unsigned int i = 1; i < max_isopeaks; ++i)
      {
        const double expected_mz = current_mz + static_cast<double>(i) * Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
        const int p = old_spectrum.findNearest(expected_mz, tolerance_dalton);
        if (p == -1)// test for missing peak
        {
          has_min_isopeaks = (i >= min_isopeaks);
          break;
        }
        else
        {
          // compare to averagine distribution
          if (use_averagine_model)
          {
            // compute average weight (weighted by intensity)
            total_weight += old_spectrum[p].getIntensity() * (old_spectrum[p].getMZ() * q - (q - 1) * Constants::PROTON_MASS_U);
            total_intensity += old_spectrum[p].getIntensity();

            // generate averagine distribution
            CoarseIsotopePatternGenerator gen(extensions.size() + 1);
            IsotopeDistribution distr = gen.estimateFromPeptideWeight(total_weight / total_intensity);

            // compute KL divergence (Sum over all x: P(x) * log(P(x) / Q(x));
            // normalize spectrum intensities as this is a density measure and the averagine distribution is also normalized to 1
            float KL = 0;
            for (unsigned int peak = 0; peak != extensions.size(); ++peak)
            {
              double Px = old_spectrum[extensions[peak]].getIntensity() / total_intensity;
              if (Px != 0.0)// Term converges to 0 for P(x) -> 0
              {
                KL += Px * log(Px / distr[peak].getIntensity());
              }
            }
            float curr_threshold = -1;
            if (distr.size() > extensions.size())
            {
              // also consider current peak
              double Px = old_spectrum[p].getIntensity() / total_intensity;
              if (Px != 0.0)
              {
                KL += Px * log(Px / distr[extensions.size()].getIntensity());
              }

              // choose threshold corresponding to cluster size (extensions do not include new peak yet)
              curr_threshold = (extensions.size() + 1 >= 6) ? averagine_check_threshold[6] : averagine_check_threshold[extensions.size() + 1];
            }

            // compare to threshold and stop extension if distribution does not fit well enough
            if (KL > curr_threshold)
            {
              has_min_isopeaks = (i >= min_isopeaks);
              break;
            }
            
            current_KL = KL;
          }
          // after model checks passed:
          extensions.push_back(p);
        }
      }

      if (has_min_isopeaks)
      {
        clusters.push_back(extensions);
        charges_of_extensions.push_back(q);
        KL_of_extensions.push_back(current_KL);
      }
    } // all charges tested, clusters complete

    // if current_peak is possible monoisotopic peak for a cluster, pick the best of its clusters, annotate peaks with a feature number
    if (!clusters.empty())
    {
      // pick cluster with largest size and highest charge (since all have the same monoisotopic peak)
      unsigned int best_idx = 0;
      Size largest_size = 0;
      int highest_charge = min_charge;
      for (unsigned int i = 0; i != clusters.size(); ++i)
      {
        if ((clusters[i].size() > largest_size) || ((clusters[i].size() == largest_size) && (charges_of_extensions[i] > highest_charge)))
        {
          largest_size = clusters[i].size();
          highest_charge = charges_of_extensions[i];
          best_idx = i;
        }
      }

      // save results
      best_clusters.push_back(clusters[best_idx]);
      best_clusters_KL.push_back(KL_of_extensions[best_idx]);
      charges_of_best_clusters.push_back(charges_of_extensions[best_idx]);
    }
  }
  
  // of all clusters found: of those conflicting, put the ones with not the best KL on a blacklist
  Size curr_cluster = 0;
  std::vector<Size> clusters_not_to_keep;
  for (auto it = best_clusters.cbegin(); it != best_clusters.cend(); ++it, ++curr_cluster)
  {
    // only treat current cluster if it is not already on the blacklist
    bool on_blacklist = false;
    for (auto rit = clusters_not_to_keep.rbegin(); rit != clusters_not_to_keep.rend() && *rit > curr_cluster; ++rit)
    {
      if (*rit == curr_cluster)
      {
        on_blacklist = true;
        break;
      }
    }
    if (on_blacklist) continue;

    // determine all clusters that have an overlap with the current cluster
    Size last_peak = it->back();
    std::vector<Size> curr_overlaps;
    for (Size temp_idx = curr_cluster + 1; temp_idx < best_clusters.size() && last_peak > best_clusters[temp_idx].front(); ++temp_idx)
    {
      curr_overlaps.push_back(temp_idx);
    }

    // if there is no overlap, nothing to do
    if (curr_overlaps.empty()) continue;

    // determine if current clusters KL is the one with highest difference to its threshold
    double curr_KL_dist = ((it->size() >= 6) ? averagine_check_threshold[6] : averagine_check_threshold[it->size()]) - best_clusters_KL[curr_cluster];
    bool curr_is_best = true;
    for (Size overlap_clusters_idx : curr_overlaps)
    {
      double temp_KL_dist = ((best_clusters[overlap_clusters_idx].size() >= 6) ? averagine_check_threshold[6] : averagine_check_threshold[best_clusters[overlap_clusters_idx].size()]) - best_clusters_KL[overlap_clusters_idx];
      if (temp_KL_dist > curr_KL_dist)
      {
        curr_is_best = false;
        break;
      }
    }

    // put either the current cluster or all overlaps on the blacklist
    if (curr_is_best)
    {
      for (Size overlap_clusters_idx : curr_overlaps)
      {
        clusters_not_to_keep.push_back(overlap_clusters_idx);
      }
    }
    else
    {
      clusters_not_to_keep.push_back(curr_cluster);
    }
  }

  int feature_number = 0;
  std::vector<Size> mono_isotopic_peak_charge(old_spectrum.size(), 0);
  std::vector<int> features(old_spectrum.size(), -1);
  std::vector<double> mono_iso_peak_intensity(old_spectrum.size(), 0);
  std::vector<Size> iso_peak_count(old_spectrum.size(), 1);

  for (Size curr_cluster = 0; curr_cluster < best_clusters.size(); ++curr_cluster)
  {
    // is curr_cluster on the blacklist?
    if (std::find(clusters_not_to_keep.begin(), clusters_not_to_keep.end(), curr_cluster) != clusters_not_to_keep.end()) continue;

    Size current_monoisotopic_peak = best_clusters[curr_cluster].front();
    mono_isotopic_peak_charge[current_monoisotopic_peak] = charges_of_best_clusters[curr_cluster];

    for (auto peak_it = best_clusters[curr_cluster].begin(); peak_it != best_clusters[curr_cluster].end(); ++peak_it)
    {
      features[*peak_it] = feature_number;
      if (add_up_intensity)
      {
        mono_iso_peak_intensity[*peak_it] += old_spectrum[*peak_it].getIntensity();
      }
    }

    if (annotate_iso_peak_count)
    {
      iso_peak_count[current_monoisotopic_peak] = best_clusters[curr_cluster].size();
    }

    ++feature_number;
  }
  

  // apply changes, i.e. select the indices which should survive
  std::vector<Size> select_idx;

  for (Size i = 0; i != spec.size(); ++i)
  {

    Size z = mono_isotopic_peak_charge[i];
    if (annotate_charge)
    {
      spec.getIntegerDataArrays()[charge_index].push_back((int)z);
    }
    if (annotate_iso_peak_count)
    {
      spec.getIntegerDataArrays()[iso_peak_count_index].push_back((int)iso_peak_count[i]);
    }
    if (add_up_intensity)
    {
      spec[i].setIntensity(mono_iso_peak_intensity[i]);
    }

    if (!keep_only_deisotoped)
    { // keep all unassigned peaks
      if (features[i] < 0)
      {
        select_idx.push_back(i);
        continue;
      }
    }

    if (z == 0) continue;

    // convert mono-isotopic peak with charge assigned by deisotoping
    if (make_single_charged)
    {
      spec[i].setMZ(spec[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
    }
    select_idx.push_back(i);
  }

  // properly subsets all datapoints (incl. dataArrays)
  spec.select(select_idx);
  spec.sortByPosition();
  return;
}

// static
void Deisotoper::deisotopeAndSingleCharge(MSSpectrum& spec,
                      double fragment_tolerance,
                      bool fragment_unit_ppm,
                      int min_charge,
                      int max_charge,
                      bool keep_only_deisotoped,
                      unsigned int min_isopeaks,
                      unsigned int max_isopeaks,
                      bool make_single_charged,
                      bool annotate_charge,
                      bool annotate_iso_peak_count,
                      bool use_decreasing_model,
                      unsigned int start_intensity_check,
                      bool add_up_intensity)
{
  OPENMS_PRECONDITION(spec.isSorted(), "Spectrum must be sorted.");

  if (min_isopeaks < 2 || max_isopeaks < 2 || min_isopeaks > max_isopeaks)
  {
    throw Exception::IllegalArgument(__FILE__,
		    __LINE__,
		    OPENMS_PRETTY_FUNCTION,
		    "Minimum/maximum number of isotopic peaks must be at least 2 (and min_isopeaks <= max_isopeaks).");
  }

  if (spec.empty()) { return; }

  Size charge_index(0);
  Size iso_peak_count_index(0);

  // reserve integer data array to store charge of peaks
  if (annotate_charge)
  {
    // expand to hold one additional integer data array to hold the charge
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("charge");
    charge_index = spec.getIntegerDataArrays().size()-1;
  }
  // reserve integer data array to store number of isotopic peaks for each isotopic pattern
  if (annotate_iso_peak_count)
  {
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("iso_peak_count");
    iso_peak_count_index = spec.getIntegerDataArrays().size()-1;
  }

  // during discovery phase, work on a constant reference (just to make sure we do not modify spec)
  const MSSpectrum& old_spectrum = spec;

  // determine charge seeds and extend them
  std::vector<size_t> mono_isotopic_peak(old_spectrum.size(), 0);
  std::vector<int> features(old_spectrum.size(), -1);
  std::vector<double> mono_iso_peak_intensity(old_spectrum.size(), 0);
  std::vector<Size> iso_peak_count(old_spectrum.size(), 1);
  int feature_number = 0;

  std::vector<size_t> extensions;

  bool has_precursor_data(false);
  double precursor_mass(0);
  if (old_spectrum.getPrecursors().size() == 1)
  {
    has_precursor_data = true;
    int precursor_charge = old_spectrum.getPrecursors()[0].getCharge();
    precursor_mass = (old_spectrum.getPrecursors()[0].getMZ() * precursor_charge) - (Constants::PROTON_MASS * precursor_charge);
  }

  for (size_t current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
  {
    const double current_mz = old_spectrum[current_peak].getMZ();
    if (add_up_intensity)
    {
      mono_iso_peak_intensity[current_peak] = old_spectrum[current_peak].getIntensity();
    }

    for (int q = max_charge; q >= min_charge; --q) // important: test charge hypothesis from high to low
    {
      // try to extend isotopes from mono-isotopic peak
      // if extension larger then min_isopeaks possible:
      //   - save charge q in mono_isotopic_peak[]
      //   - annotate_charge all isotopic peaks with feature number
      if (features[current_peak] == -1) // only process peaks which have no assigned feature number
      {
        bool has_min_isopeaks = true;
        const double tolerance_dalton = fragment_unit_ppm ? Math::ppmToMass(fragment_tolerance, current_mz) : fragment_tolerance;

        // do not bother testing charges q (and masses m) with: m/q > precursor_mass/q (or m > precursor_mass)
        if (has_precursor_data)
        {
          double current_theo_mass = (current_mz * q) - (Constants::PROTON_MASS * q);
          if (current_theo_mass > (precursor_mass + tolerance_dalton))
          {
            continue;
          }
        }

        extensions.clear();
        extensions.push_back(current_peak);
        for (unsigned int i = 1; i < max_isopeaks; ++i)
        {
          const double expected_mz = current_mz + static_cast<double>(i) * Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
          const int p = old_spectrum.findNearest(expected_mz, tolerance_dalton);
          if (p == -1) // test for missing peak
          {
            has_min_isopeaks = (i >= min_isopeaks);
            break;
          }
          else
          {
            // Possible improvement: include proper averagine model filtering
            // for now start at the peak with i = start_intensity_check to test hypothesis
            // if start_intensity_check = 0 or 1, start checking by comparing monoisotopic and second isotopic peak
            // if start_intensity_check = 2, start checking by comparing second isotopic peak with the third, etc.
            // Note: this is a common approach used in several other search engines
            if (use_decreasing_model && (i >= start_intensity_check) && (old_spectrum[p].getIntensity() > old_spectrum[extensions.back()].getIntensity()))
            {
              has_min_isopeaks = (i >= min_isopeaks);
              break;
            }
            // averagine check passed or skipped
            extensions.push_back(p);
            if (annotate_iso_peak_count)
            {
              iso_peak_count[current_peak] = i + 1; // with "+ 1" the monoisotopic peak is counted as well
            }
          }
        }

        if (has_min_isopeaks)
        {
          mono_isotopic_peak[current_peak] = q;
          for (unsigned int i = 0; i != extensions.size(); ++i)
          {
            features[extensions[i]] = feature_number;
            // monoisotopic peak intensity is already set above, add up the other intensities here
            if (add_up_intensity && (i != 0))
            {
              mono_iso_peak_intensity[current_peak] += old_spectrum[extensions[i]].getIntensity();
            }
          }
          ++feature_number;
        }
      }
    }
  }

  // apply changes, i.e. select the indices which should survive
  std::vector<Size> select_idx;

  for (size_t i = 0; i != spec.size(); ++i)
  {
    Size z = mono_isotopic_peak[i];
    if (annotate_charge)
    {
      spec.getIntegerDataArrays()[charge_index].push_back((int)z);
    }
    if (annotate_iso_peak_count)
    {
      spec.getIntegerDataArrays()[iso_peak_count_index].push_back((int)iso_peak_count[i]);
    }
    if (add_up_intensity)
    {
      spec[i].setIntensity(mono_iso_peak_intensity[i]);
    }

    if (!keep_only_deisotoped)
    { // keep all unassigned peaks
      if (features[i] < 0)
      {
        select_idx.push_back(i);
        continue;
      }
    }

    if (z == 0) continue;

    // convert mono-isotopic peak with charge assigned by deisotoping
    if (make_single_charged)
    {
      spec[i].setMZ(spec[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
    }
    select_idx.push_back(i);
  }

  // properly subsets all datapoints (incl. dataArrays)
  spec.select(select_idx);
  spec.sortByPosition();
  return;
}
} // namespace
