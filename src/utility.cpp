/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <math.h>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "utility.h"
#include "globals.h"
#include "Data.h"

namespace ranger {

void equalSplit(std::vector<uint>& result, uint start, uint end, uint num_parts) {

  result.reserve(num_parts + 1);

  // Return range if only 1 part
  if (num_parts == 1) {
    result.push_back(start);
    result.push_back(end + 1);
    return;
  }

  // Return vector from start to end+1 if more parts than elements
  if (num_parts > end - start + 1) {
    for (uint i = start; i <= end + 1; ++i) {
      result.push_back(i);
    }
    return;
  }

  uint length = (end - start + 1);
  uint part_length_short = length / num_parts;
  uint part_length_long = (uint) ceil(length / ((double) num_parts));
  uint cut_pos = length % num_parts;

  // Add long ranges
  for (uint i = start; i < start + cut_pos * part_length_long; i = i + part_length_long) {
    result.push_back(i);
  }

  // Add short ranges
  for (uint i = start + cut_pos * part_length_long; i <= end + 1; i = i + part_length_short) {
    result.push_back(i);
  }
}

void loadDoubleVectorFromFile(std::vector<double>& result, std::string filename) { // #nocov start

  // Open input file
  std::ifstream input_file;
  input_file.open(filename);
  if (!input_file.good()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  // Read the first line, ignore the rest
  std::string line;
  getline(input_file, line);
  std::stringstream line_stream(line);
  double token;
  while (line_stream >> token) {
    result.push_back(token);
  }
} // #nocov end

void drawWithoutReplacement(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    size_t num_samples) {
  if (num_samples < max / 10) {
    drawWithoutReplacementSimple(result, random_number_generator, max, num_samples);
  } else {
    //drawWithoutReplacementKnuth(result, random_number_generator, max, skip, num_samples);
    drawWithoutReplacementFisherYates(result, random_number_generator, max, num_samples);
  }
}

void drawWithoutReplacementSkip(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    const std::vector<size_t>& skip, size_t num_samples) {
  if (num_samples < max / 10) {
    drawWithoutReplacementSimple(result, random_number_generator, max, skip, num_samples);
  } else {
    //drawWithoutReplacementKnuth(result, random_number_generator, max, skip, num_samples);
    drawWithoutReplacementFisherYates(result, random_number_generator, max, skip, num_samples);
  }
}

void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    size_t num_samples) {

  result.reserve(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(max, false);

  std::uniform_int_distribution<size_t> unif_dist(0, max - 1);
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = unif_dist(random_number_generator);
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    const std::vector<size_t>& skip, size_t num_samples) {

  result.reserve(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(max, false);

  std::uniform_int_distribution<size_t> unif_dist(0, max - 1 - skip.size());
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = unif_dist(random_number_generator);
      for (auto& skip_value : skip) {
        if (draw >= skip_value) {
          ++draw;
        }
      }
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

void drawWithoutReplacementFisherYates(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max, size_t num_samples) {

  // Create indices
  result.resize(max);
  std::iota(result.begin(), result.end(), 0);

  // Draw without replacement using Fisher Yates algorithm
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t i = 0; i < num_samples; ++i) {
    size_t j = i + distribution(random_number_generator) * (max - i);
    std::swap(result[i], result[j]);
  }

  result.resize(num_samples);
}

void drawWithoutReplacementFisherYates(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max, const std::vector<size_t>& skip, size_t num_samples) {

  // Create indices
  result.resize(max);
  std::iota(result.begin(), result.end(), 0);

  // Skip indices
  for (size_t i = 0; i < skip.size(); ++i) {
    result.erase(result.begin() + skip[skip.size() - 1 - i]);
  }

  // Draw without replacement using Fisher Yates algorithm
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t i = 0; i < num_samples; ++i) {
    size_t j = i + distribution(random_number_generator) * (max - skip.size() - i);
    std::swap(result[i], result[j]);
  }

  result.resize(num_samples);
}

void drawWithoutReplacementWeighted(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    const std::vector<size_t>& indices, size_t num_samples, const std::vector<double>& weights) {

  result.reserve(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(indices.size(), false);

  std::discrete_distribution<> weighted_dist(weights.begin(), weights.end());
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = weighted_dist(random_number_generator);
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(indices[draw]);
  }
}

void drawWithoutReplacementWeighted(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max_index, size_t num_samples, const std::vector<double>& weights) {

  result.reserve(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(max_index + 1, false);

  std::discrete_distribution<> weighted_dist(weights.begin(), weights.end());
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = weighted_dist(random_number_generator);
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

double mostFrequentValue(const std::unordered_map<double, size_t>& class_count,
    std::mt19937_64 random_number_generator) {
  std::vector<double> major_classes;

  // Find maximum count
  size_t max_count = 0;
  for (auto& class_value : class_count) {
    if (class_value.second > max_count) {
      max_count = class_value.second;
      major_classes.clear();
      major_classes.push_back(class_value.first);
    } else if (class_value.second == max_count) {
      major_classes.push_back(class_value.first);
    }
  }

  if (major_classes.size() == 1) {
    return major_classes[0];
  } else {
    // Choose randomly
    std::uniform_int_distribution<size_t> unif_dist(0, major_classes.size() - 1);
    return major_classes[unif_dist(random_number_generator)];
  }
}

double computeConcordanceIndex(const Data& data, const std::vector<double>& sum_chf,
    const std::vector<size_t>& sample_IDs) {

  // Compute concordance index
  double concordance = 0;
  double permissible = 0;
  for (size_t i = 0; i < sum_chf.size(); ++i) {
    size_t sample_i = i;
    if (!sample_IDs.empty()) {
      sample_i = sample_IDs[i];
    }
    double time_i = data.get_y(sample_i, 0);
    double status_i = data.get_y(sample_i, 1);

    for (size_t j = i + 1; j < sum_chf.size(); ++j) {
      size_t sample_j = j;
      if (!sample_IDs.empty()) {
        sample_j = sample_IDs[j];
      }
      double time_j = data.get_y(sample_j, 0);
      double status_j = data.get_y(sample_j, 1);

      if (time_i < time_j && status_i == 0) {
        continue;
      }
      if (time_j < time_i && status_j == 0) {
        continue;
      }
      if (time_i == time_j && status_i == status_j) {
        continue;
      }

      permissible += 1;

      if (time_i < time_j && sum_chf[i] > sum_chf[j]) {
        concordance += 1;
      } else if (time_j < time_i && sum_chf[j] > sum_chf[i]) {
        concordance += 1;
      } else if (sum_chf[i] == sum_chf[j]) {
        concordance += 0.5;
      }

    }
  }

  return (concordance / permissible);

}

std::string uintToString(uint number) {
#if WIN_R_BUILD == 1
  std::stringstream temp;
  temp << number;
  return temp.str();
#else
  return std::to_string(number);
#endif
}

std::string beautifyTime(uint seconds) { // #nocov start
  std::string result;

  // Add seconds, minutes, hours, days if larger than zero
  uint out_seconds = (uint) seconds % 60;
  result = uintToString(out_seconds) + " seconds";
  uint out_minutes = (seconds / 60) % 60;
  if (seconds / 60 == 0) {
    return result;
  } else if (out_minutes == 1) {
    result = "1 minute, " + result;
  } else {
    result = uintToString(out_minutes) + " minutes, " + result;
  }
  uint out_hours = (seconds / 3600) % 24;
  if (seconds / 3600 == 0) {
    return result;
  } else if (out_hours == 1) {
    result = "1 hour, " + result;
  } else {
    result = uintToString(out_hours) + " hours, " + result;
  }
  uint out_days = (seconds / 86400);
  if (out_days == 0) {
    return result;
  } else if (out_days == 1) {
    result = "1 day, " + result;
  } else {
    result = uintToString(out_days) + " days, " + result;
  }
  return result;
} // #nocov end

// #nocov start
size_t roundToNextMultiple(size_t value, uint multiple) {

  if (multiple == 0) {
    return value;
  }

  size_t remainder = value % multiple;
  if (remainder == 0) {
    return value;
  }

  return value + multiple - remainder;
}
// #nocov end

void splitString(std::vector<std::string>& result, const std::string& input, char split_char) { // #nocov start

  std::istringstream ss(input);
  std::string token;

  while (std::getline(ss, token, split_char)) {
    result.push_back(token);
  }
} // #nocov end

void shuffleAndSplit(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all, size_t n_first,
    std::mt19937_64 random_number_generator) {

  // Reserve space
  first_part.resize(n_all);

  // Fill with 0..n_all-1 and shuffle
  std::iota(first_part.begin(), first_part.end(), 0);
  std::shuffle(first_part.begin(), first_part.end(), random_number_generator);

  // Copy to second part
  second_part.resize(n_all - n_first);
  std::copy(first_part.begin() + n_first, first_part.end(), second_part.begin());

  // Resize first part
  first_part.resize(n_first);
}

void shuffleAndSplitAppend(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all,
    size_t n_first, const std::vector<size_t>& mapping, std::mt19937_64 random_number_generator) {
  // Old end is start position for new data
  size_t first_old_size = first_part.size();
  size_t second_old_size = second_part.size();

  // Reserve space
  first_part.resize(first_old_size + n_all);
  std::vector<size_t>::iterator first_start_pos = first_part.begin() + first_old_size;

  // Fill with 0..n_all-1 and shuffle
  std::iota(first_start_pos, first_part.end(), 0);
  std::shuffle(first_start_pos, first_part.end(), random_number_generator);

  // Mapping
  for (std::vector<size_t>::iterator j = first_start_pos; j != first_part.end(); ++j) {
    *j = mapping[*j];
  }

  // Copy to second part
  second_part.resize(second_part.size() + n_all - n_first);
  std::vector<size_t>::iterator second_start_pos = second_part.begin() + second_old_size;
  std::copy(first_start_pos + n_first, first_part.end(), second_start_pos);

  // Resize first part
  first_part.resize(first_old_size + n_first);
}

std::string checkUnorderedVariables(const Data& data, const std::vector<std::string>& unordered_variable_names) { // #nocov start
  size_t num_rows = data.getNumRows();
  std::vector<size_t> sampleIDs(num_rows);
  std::iota(sampleIDs.begin(), sampleIDs.end(), 0);

  // Check for all unordered variables
  for (auto& variable_name : unordered_variable_names) {
    size_t varID = data.getVariableID(variable_name);
    std::vector<double> all_values;
    data.getAllValues(all_values, sampleIDs, varID, 0, sampleIDs.size());

    // Check level count
    size_t max_level_count = 8 * sizeof(size_t) - 1;
    if (all_values.size() > max_level_count) {
      return "Too many levels in unordered categorical variable " + variable_name + ". Only "
          + uintToString(max_level_count) + " levels allowed on this system.";
    }

    // Check positive integers
    if (!checkPositiveIntegers(all_values)) {
      return "Not all values in unordered categorical variable " + variable_name + " are positive integers.";
    }
  }
  return "";
} // #nocov end

bool checkPositiveIntegers(const std::vector<double>& all_values) { // #nocov start
  for (auto& value : all_values) {
    if (value < 1 || !(floor(value) == value)) {
      return false;
    }
  }
  return true;
} // #nocov end

double maxstatPValueLau92(double b, double minprop, double maxprop) {

  if (b < 1) {
    return 1.0;
  }

  // Compute only once (minprop/maxprop don't change during runtime)
  static double logprop = log((maxprop * (1 - minprop)) / ((1 - maxprop) * minprop));

  double db = dstdnorm(b);
  double p = 4 * db / b + db * (b - 1 / b) * logprop;

  if (p > 0) {
    return p;
  } else {
    return 0;
  }
}

double maxstatPValueLau94(double b, double minprop, double maxprop, size_t N, const std::vector<size_t>& m) {

  double D = 0;
  for (size_t i = 0; i < m.size() - 1; ++i) {
    double m1 = m[i];
    double m2 = m[i + 1];

    double t = sqrt(1.0 - m1 * (N - m2) / ((N - m1) * m2));
    D += 1 / M_PI * exp(-b * b / 2) * (t - (b * b / 4 - 1) * (t * t * t) / 6);
  }

  return 2 * (1 - pstdnorm(b)) + D;
}

double maxstatPValueUnadjusted(double b) {
  return 2 * pstdnorm(-b);
}

double dstdnorm(double x) {
  return exp(-0.5 * x * x) / sqrt(2 * M_PI);
}

double pstdnorm(double x) {
  return 0.5 * (1 + erf(x / sqrt(2.0)));
}

std::vector<double> adjustPvalues(std::vector<double>& unadjusted_pvalues) {
  size_t num_pvalues = unadjusted_pvalues.size();
  std::vector<double> adjusted_pvalues(num_pvalues, 0);

  // Get order of p-values
  std::vector<size_t> indices = order(unadjusted_pvalues, true);

  // Compute adjusted p-values
  adjusted_pvalues[indices[0]] = unadjusted_pvalues[indices[0]];
  for (size_t i = 1; i < indices.size(); ++i) {
    size_t idx = indices[i];
    size_t idx_last = indices[i - 1];

    adjusted_pvalues[idx] = std::min(adjusted_pvalues[idx_last],
        (double) num_pvalues / (double) (num_pvalues - i) * unadjusted_pvalues[idx]);
  }
  return adjusted_pvalues;
}

std::vector<double> logrankScores(const std::vector<double>& time, const std::vector<double>& status) {
  size_t n = time.size();
  std::vector<double> scores(n);

  // Get order of timepoints
  std::vector<size_t> indices = order(time, false);

  // Compute scores
  double cumsum = 0;
  size_t last_unique = -1;
  for (size_t i = 0; i < n; ++i) {

    // Continue if next value is the same
    if (i < n - 1 && time[indices[i]] == time[indices[i + 1]]) {
      continue;
    }

    // Compute sum and scores for all non-unique values in a row
    for (size_t j = last_unique + 1; j <= i; ++j) {
      cumsum += status[indices[j]] / (n - i);
    }
    for (size_t j = last_unique + 1; j <= i; ++j) {
      scores[indices[j]] = status[indices[j]] - cumsum;
    }

    // Save last computed value
    last_unique = i;
  }

  return scores;
}

void maxstat(const std::vector<double>& scores, const std::vector<double>& x, const std::vector<size_t>& indices,
    double& best_maxstat, double& best_split_value, double minprop, double maxprop) {
  size_t n = x.size();

  double sum_all_scores = 0;
  for (size_t i = 0; i < n; ++i) {
    sum_all_scores += scores[indices[i]];
  }

  // Compute sum of differences from mean for variance
  double mean_scores = sum_all_scores / n;
  double sum_mean_diff = 0;
  for (size_t i = 0; i < n; ++i) {
    sum_mean_diff += (scores[i] - mean_scores) * (scores[i] - mean_scores);
  }

  // Get smallest and largest split to consider, -1 for compatibility with R maxstat
  size_t minsplit = 0;
  if (n * minprop > 1) {
    minsplit = n * minprop - 1;
  }
  size_t maxsplit = n * maxprop - 1;

  // For all unique x-values
  best_maxstat = -1;
  best_split_value = -1;
  double sum_scores = 0;
  size_t n_left = 0;
  for (size_t i = 0; i <= maxsplit; ++i) {

    sum_scores += scores[indices[i]];
    n_left++;

    // Dont consider splits smaller than minsplit for splitting (but count)
    if (i < minsplit) {
      continue;
    }

    // Consider only unique values
    if (i < n - 1 && x[indices[i]] == x[indices[i + 1]]) {
      continue;
    }

    // If value is largest possible value, stop
    if (x[indices[i]] == x[indices[n - 1]]) {
      break;
    }

    double S = sum_scores;
    double E = (double) n_left / (double) n * sum_all_scores;
    double V = (double) n_left * (double) (n - n_left) / (double) (n * (n - 1)) * sum_mean_diff;
    double T = fabs((S - E) / sqrt(V));

    if (T > best_maxstat) {
      best_maxstat = T;

      // Use mid-point split if possible
      if (i < n - 1) {
        best_split_value = (x[indices[i]] + x[indices[i + 1]]) / 2;
      } else {
        best_split_value = x[indices[i]];
      }
    }
  }
}

std::vector<size_t> numSamplesLeftOfCutpoint(std::vector<double>& x, const std::vector<size_t>& indices) {
  std::vector<size_t> num_samples_left;
  num_samples_left.reserve(x.size());

  for (size_t i = 0; i < x.size(); ++i) {
    if (i == 0) {
      num_samples_left.push_back(1);
    } else if (x[indices[i]] == x[indices[i - 1]]) {
      ++num_samples_left[num_samples_left.size() - 1];
    } else {
      num_samples_left.push_back(num_samples_left[num_samples_left.size() - 1] + 1);
    }
  }

  return num_samples_left;
}

// #nocov start
std::stringstream& readFromStream(std::stringstream& in, double& token) {
  if (!(in >> token) && (std::fpclassify(token) == FP_SUBNORMAL)) {
    in.clear();
  }
  return in;
}
// #nocov end

double betaLogLik(double y, double mean, double phi) {

  // Avoid 0 and 1
  if (y < std::numeric_limits<double>::epsilon()) {
    y = std::numeric_limits<double>::epsilon();
  } else if (y >= 1) {
    y = 1 - std::numeric_limits<double>::epsilon();
  }
  if (mean < std::numeric_limits<double>::epsilon()) {
    mean = std::numeric_limits<double>::epsilon();
  } else if (mean >= 1) {
    mean = 1 - std::numeric_limits<double>::epsilon();
  }
  if (phi < std::numeric_limits<double>::epsilon()) {
    phi = std::numeric_limits<double>::epsilon();
  } else if (mean >= 1) {
    phi = 1 - std::numeric_limits<double>::epsilon();
  }

  return (lgamma(phi) - lgamma(mean * phi) - lgamma((1 - mean) * phi) + (mean * phi - 1) * log(y)
      + ((1 - mean) * phi - 1) * log(1 - y));
}

} // namespace ranger
