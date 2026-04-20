#pragma once

#include <vector>
#include <unordered_map>
//#include <limits>
#include "CompressedColorSets.hh"


using namespace std;

inline void process_word(uint64_t word, const uint64_t base, vector<int16_t> &results, const uint64_t freq)
{
    while (word)
    {
        // uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        results[base + bit] += freq;
        word &= word - 1; // clear lowest bit
    }
}

inline void read_bv(const uint64_t *data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<int16_t> &results)
{

    const uint64_t *ptr = data + (start * n_colors) / 64;
    uint64_t bit_offset = (start * n_colors) % 64;

    uint64_t color_id = 0;
    uint64_t bits_left = n_colors;

    if (bit_offset != 0)
    {
        // 1. Read the first word
        uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
        uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
        uint64_t word = (*ptr >> bit_offset) & mask;
        process_word(word, color_id, results, freq);

        ptr++;
        color_id += bits_to_read;
        bits_left -= bits_to_read;
    }

    // 2. Read aligned words in btw
    while (bits_left >= 64)
    {
        uint64_t word = *ptr++;
        process_word(word, color_id, results, freq);
        color_id += 64;
        bits_left -= 64;
    }

    // 3. Read the last word (if any)
    if (bits_left > 0)
    {
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;
        process_word(word, color_id, results, freq);
    }
}


void read_verydense(const int64_t pos, const uint64_t freq, const DeltaSet &EF, const vector<uint16_t> &L, const uint64_t n_colors, vector<int16_t> &results, uint16_t &dense)
{
    dense += freq;
    const size_t end = EF.get_start(pos); // exclusive end
    size_t start = (pos > 0) ? EF.get_start(pos - 1) : 0; // inclusive start
    while (start < end)
    {
        results[L[start++]] -= freq;
    }
}

void read_colors(const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results, const vector<pair<int64_t, uint64_t>> &fmin_v, uint16_t &dense)
{
    const sdsl::bit_vector &BV = CCS.getBV();
    const uint64_t *data = BV.data();
    const vector<uint16_t> &L = CCS.getL();
    const DeltaSet &EF = CCS.getEF();

    // pos < sparse_count; [sparse]
    // sparse_count <= pos < dense_count; [very dense]
    // pos >= dense; [bitmap]

    const uint64_t sparse_count = CCS.sparse_count;
    const uint64_t dense_count = CCS.dense_count;

    // Exploit the fact that the pos are sorted
    uint64_t i;
    for (i = 0; i < fmin_v.size(); i++)
    {
        const auto &[pos, freq] = fmin_v[i];
        // sparse
        if (pos < sparse_count && pos > 0)
        {
            // read from L - in sparse mode, start is always > 0
            const size_t end = EF.get_start(pos); // exclusive end
            size_t start = EF.get_start(pos - 1); // inclusive start
            while (start < end)
            {
                results[L[start++]] += freq;
            }
        }
        else
        {
            break;
        }
    }
    // very dense
    uint64_t j;
    for (j = i; j < fmin_v.size(); j++)
    {
        const auto &[pos, freq] = fmin_v[j];
        if (pos < dense_count)
        {
            // read complementary values from L
            read_verydense(pos, freq, EF, L, n_colors, results, dense);
        }
        else
        {
            break;
        }
    }
    // read from BV
    for (auto i = j; i < fmin_v.size(); i++)
    {
        const auto &[pos, freq] = fmin_v[i];
        uint64_t start = pos - dense_count;
        read_bv(data, start, freq, n_colors, results);
    }
}

void counting_sort(const vector<int16_t> &results, vector<pair<uint16_t, uint16_t>> &ans, const size_t found_fmin, const uint16_t n_colors)
{
    vector<uint16_t> counts(found_fmin + 1);

    for (size_t idx = 0; idx < n_colors; idx++)
    {
        counts[results[idx]]++;
    }

    // Cumulative Sums
    for (size_t c = 1; c < counts.size(); c++)
    {
        counts[c] += counts[c - 1];
    }

    ans.resize(n_colors);
    for (size_t idx = 0; idx < n_colors; idx++)
    {
        ans[counts[results[idx]] - 1] = {static_cast<uint16_t>(idx), results[idx]};
        counts[results[idx]]--;
    }
}

inline void pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results)
{
    // vector<int16_t> results(n_colors, 0);
    /* if (results.size() != n_colors) {
        results.assign(n_colors, 0);
    } else {
        std::fill(results.begin(), results.end(), 0);
    }
 */
    if (Fmin.empty())
    {
        return;
    } // not 0 as everything wuold be >=

    std::fill(results.begin(), results.end(), 0);

    std::sort(Fmin.begin(), Fmin.end());
    vector<pair<int64_t, uint64_t>> fmin_v;
    fmin_v.reserve(Fmin.size());

    for (size_t i = 0; i < Fmin.size();)
    {
        size_t j = i + 1;
        while (j < Fmin.size() && Fmin[j] == Fmin[i])
            ++j;
        fmin_v.emplace_back(Fmin[i], j - i);
        i = j;
    }

    /* // Count freq of each fmin
    std::unordered_map<int64_t, uint64_t> fmin_counts;
    for (auto v : Fmin) {
        fmin_counts[v]++;
    }

    // vector for sorted output so that it is possible to scan color_set_concat
    vector<pair<int64_t, uint64_t>> fmin_v(fmin_counts.begin(), fmin_counts.end());
    std::sort(fmin_v.begin(), fmin_v.end()); */
    uint16_t dense = 0;
    read_colors(CCS, n_colors, results, fmin_v, dense);
    // TODO keep track of which counters were incremented and set to zero only those
    // Add number of dense sets
    for (auto &r : results)
    {
        r += dense;
    }
    return;
}

inline int16_t pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results, const float t)
{ // vector<uint64_t>& results,
    if (Fmin.empty())
    {
        return 1;
    } // not 0 as everything wuold be >=

    std::fill(results.begin(), results.end(), 0);

    std::sort(Fmin.begin(), Fmin.end());
    vector<pair<int64_t, uint64_t>> fmin_v;
    fmin_v.reserve(Fmin.size());

    for (size_t i = 0; i < Fmin.size();)
    {
        size_t j = i + 1;
        while (j < Fmin.size() && Fmin[j] == Fmin[i])
            ++j;
        fmin_v.emplace_back(Fmin[i], j - i);
        i = j;
    }

    // Check the values above the minimum in search
    const size_t found_fmin = Fmin.size(); // # total finimizers
    int16_t T = found_fmin * t;
    uint16_t dense = 0;

    read_colors(CCS, n_colors, results, fmin_v, dense);

    // If we care about the number of matches, add number of dense sets
    for (auto& r: results){r+=dense;}

    // adjust T
    // T -= dense;
    //for (auto& r :results){ r+= dense;}
    return T;
}

inline void pseudoalignment_stats_sorted(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results,     vector<pair<uint16_t, uint16_t>> &ans)
{
    // vector<int16_t> results(n_colors, 0);
    /* if (results.size() != n_colors) {
        results.assign(n_colors, 0);
    } else {
        std::fill(results.begin(), results.end(), 0);
    }
 */
    if (Fmin.empty())
    {
        return;
    } // not 0 as everything wuold be >=

    std::fill(results.begin(), results.end(), 0);

    std::sort(Fmin.begin(), Fmin.end());
    vector<pair<int64_t, uint64_t>> fmin_v;
    fmin_v.reserve(Fmin.size());

    for (size_t i = 0; i < Fmin.size();)
    {
        size_t j = i + 1;
        while (j < Fmin.size() && Fmin[j] == Fmin[i])
            ++j;
        fmin_v.emplace_back(Fmin[i], j - i);
        i = j;
    }

    /* // Count freq of each fmin
    std::unordered_map<int64_t, uint64_t> fmin_counts;
    for (auto v : Fmin) {
        fmin_counts[v]++;
    }

    // vector for sorted output so that it is possible to scan color_set_concat
    vector<pair<int64_t, uint64_t>> fmin_v(fmin_counts.begin(), fmin_counts.end());
    std::sort(fmin_v.begin(), fmin_v.end()); */
    uint16_t dense = 0;
    read_colors(CCS, n_colors, results, fmin_v, dense);
    // TODO keep track of which counters were incremented and set to zero only those
    // Add number of dense sets
    for (auto &r : results)
    {
        r += dense;
    }
    const size_t found_fmin = Fmin.size(); // # total finimizers
    // Sort results so that the output is sorted
    counting_sort(results, ans, found_fmin, n_colors);
    return;
}


inline int16_t pseudoalignment_stats_sorted(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results, const float t, vector<pair<uint16_t, uint16_t>> &ans)
{
    if (Fmin.empty())
    {
        return 1;
    } // not 0 as everything wuold be >=

    std::fill(results.begin(), results.end(), 0);

    std::sort(Fmin.begin(), Fmin.end());
    vector<pair<int64_t, uint64_t>> fmin_v;
    fmin_v.reserve(Fmin.size());

    for (size_t i = 0; i < Fmin.size();)
    {
        size_t j = i + 1;
        while (j < Fmin.size() && Fmin[j] == Fmin[i])
            ++j;
        fmin_v.emplace_back(Fmin[i], j - i);
        i = j;
    }

    // Check the values above the minimum in search
    const size_t found_fmin = Fmin.size(); // # total finimizers
    int16_t T = found_fmin * t;
    uint16_t dense = 0;

    read_colors(CCS, n_colors, results, fmin_v, dense);

    // If we care about the number of matches, add number of dense sets
    for (auto& r: results){r+=dense;}

    // adjust T
    // T -= dense;
    //for (auto& r :results){ r+= dense;}
    //  Sort results so that the output is sorted
    counting_sort(results, ans, found_fmin, n_colors);
    return T;
}

