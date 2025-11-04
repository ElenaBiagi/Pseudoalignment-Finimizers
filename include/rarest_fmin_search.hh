#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

#include "BoundedDeque.hh"

#include "common.hh"
#include "bitsearch.hh"
#include "Buckets.hh"

#include "sdsl/bit_vectors.hpp"

/* struct Candidate {
    uint64_t len, _int, color, start;
    inline bool operator<(const Candidate& o) const {
        return std::tie(len, _int, color, start) < std::tie(o.len, o._int, o.color, o.start);
    }
}; */

// TODO Deal with empty SB
inline void FindShortFinimizer(const uint64_t int_sp_len, uint64_t int_sp, const std::unordered_map<uint32_t, pair<uint8_t, int64_t>> &sB, const uint64_t start, const uint64_t end, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t> &k_fmin)
{
    // 2. Prefix NOT found
    // Start from the longest possible prefix
    // if you find a real match, stop
    if (sB.empty())
    {
        return;
    }
    uint64_t sp_len = int_sp_len; // plen is here plen-1: sp_len must be < plen as the whole prefix was not found
    while (sp_len > 0)
    {
        auto it = sB.find(int_sp);
        if (it != sB.end() && sp_len == it->second.first)
        { // real match
            // all_fmin.insert(make_tuple(sp_len, int_sp, it->second.second, start));
            if ((start + sp_len - 1) > end)
            {
                next_candidates.push_back(make_tuple(sp_len + start - 1, int_sp, it->second.second, start));
            } // Sorted based on END
            else
            {
                tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {sp_len, int_sp, it->second.second, start};
                if (new_fmin < k_fmin)
                {
                    curr_candidates.clear();
                    k_fmin = new_fmin;
                    // curr_candidates.push_back(new_fmin);
                }
                else
                {
                    // push_monotone_increasing(curr_candidates, new_fmin);
                    while (curr_candidates.back() > new_fmin)
                    {
                        curr_candidates.pop_back();
                    }
                }
                curr_candidates.push_back(new_fmin);
            }
            return;
        }
        int_sp >>= 2;
        sp_len--;
    }
}

inline void FindPrefix(const Bucket &bucket_p, const uint64_t plen, const uint8_t s_len, const uint64_t int_s, const uint64_t int_p, const uint64_t start, const uint64_t end, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t> &k_fmin)
{
    auto [pos, len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len);
    if (pos > -1)
    {
        uint64_t f_int = (int_p << (len * 2)) | (int_s >> ((s_len - len) * 2)); //  Shift p_int to the left by len*2, and int_s to the right to remove the unused chars
        if ((start + plen + len - 1) > end)
        {
            next_candidates.push_back(make_tuple(plen + len + start - 1, f_int, bucket_p.color_set_ids[pos], start));
        } // Sorted based on END
        else
        {
            tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {plen + len, f_int, bucket_p.color_set_ids[pos], start};
            if (new_fmin < k_fmin)
            {
                curr_candidates.clear();
                k_fmin = new_fmin;
            }
            else
            {
                while (!curr_candidates.empty() && curr_candidates.back() > new_fmin)
                {
                    curr_candidates.pop_back();
                }
            }
            curr_candidates.push_back(new_fmin);
        }
    }
}

// The query is shorter than (k - plen)
inline void FindPrefix_short(const vector<optional<Bucket>> &buckets, const uint64_t plen, const uint8_t s_len, const uint64_t int_s, const uint64_t int_p, const uint64_t start, const uint64_t end, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t> &k_fmin)
{
    const Bucket &bucket_p = *buckets[int_p];
    auto [pos, len] = bitMagicSearch_short(bucket_p.tail_data, int_s, s_len);
    if (pos > -1)
    {
        uint64_t f_int = (int_p << (len * 2)) | (int_s >> ((s_len - len) * 2)); //  Shift p_int to the left by len*2, and int_s to the right to remove the unused chars
        if ((start + plen + len - 1) > end)
        {
            next_candidates.push_back(make_tuple(plen + len + start - 1, f_int, bucket_p.color_set_ids[pos], start));
        } // Sorted based on END
        else
        {
            tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {plen + len, f_int, bucket_p.color_set_ids[pos], start};
            if (new_fmin < k_fmin)
            {
                curr_candidates.clear();
                k_fmin = new_fmin;
            }
            else
            {
                while (curr_candidates.back() > new_fmin)
                {
                    curr_candidates.pop_back();
                }
            }
            curr_candidates.push_back(new_fmin);
        }
    }
}

void PickFinimizer(vector<int64_t> &Fmin, const uint64_t kmer_start, const uint64_t k, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t> &k_fmin)
{ //, const string& input){
    if (!curr_candidates.empty())
    {
        k_fmin = curr_candidates.front();

        // cout << input.substr(get<3>(k_fmin),get<0>(k_fmin)) << endl;
        Fmin.push_back(get<2>(k_fmin));
    }

    // 1. Check if this finimizer is good for the next k-mer (still in the window)
    while (!curr_candidates.empty() && get<3>(curr_candidates.front()) <= kmer_start)
    {
        curr_candidates.pop_front();
    }
    k_fmin = (curr_candidates.empty()) ? static_cast<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>(make_tuple(k + 1, 0, 0, kmer_start + 1)) : curr_candidates.front();

    // 2. Check if the NEXT finimizer would be good for the next k-mer
    if (!next_candidates.empty())
    {
        const auto &next_fmin = next_candidates.front(); // tuple<uint64_t, uint64_t, uint64_t, uint64_t>
        tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {get<0>(next_fmin) - get<3>(next_fmin) + 1, get<1>(next_fmin), get<2>(next_fmin), get<3>(next_fmin)};
        if (get<0>(next_fmin) <= kmer_start + k)
        { // end of the fmin is before end of next kmer
            if (new_fmin < k_fmin)
            { // always true if curr_candidates is empty
                curr_candidates.clear();
                k_fmin = new_fmin;
            }
            else
            {
                while (curr_candidates.back() > new_fmin)
                {
                    curr_candidates.pop_back();
                }
            }
            curr_candidates.push_back(new_fmin);
            next_candidates.pop_front();
        }
    }
}

void rarest_fmin_streaming_search(const string &input, const vector<Bucket> &buckets, const sdsl::rank_support_v5<1> &buckets_rs, const std::unordered_map<uint32_t, pair<uint8_t, int64_t>> &sB, const uint64_t plen, const uint64_t k, vector<int64_t> &Fmin)
{
    const int64_t str_len = input.size();

    uint64_t start = 0;
    uint64_t kmer_start = 0;

    BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> curr_candidates(k); // sort based on len, int (color, start)
    BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> next_candidates(k); // sort by end (start+len-1)
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> k_fmin = {k + 1, 0, 0, kmer_start};

    uint64_t int_p = prefix2int(input, start, plen); // start = 0
    uint8_t s_len = k - plen;
    uint64_t int_s = prefix2int(input, start + plen, s_len); // tail
    uint64_t int_sp;

    // Check the first characters
    auto correct_bucket_idx = buckets_rs.rank(int_p + 1);
    if (buckets_rs.rank(int_p) < correct_bucket_idx)
    {
        // 1. prefix found
        const Bucket &correct_bucket = buckets[correct_bucket_idx - 1];                                                             // 0 indexed
        FindPrefix(correct_bucket, plen, s_len, int_s, int_p, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
    }
    else
    {
        // 2. look for a shorter finimizer
        int_sp = int_p >> 2;
        FindShortFinimizer(plen - 1, int_sp, sB, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
    }

    // TODO Fix for very short queries

    // The first k-1 characters do not contail all possible finimizers for the first k-mer
    for (start = 1; start < k - 1 && start <= str_len - plen; start++)
    {
        int_p = stream_kmer(int_p, input[start + plen - 1], plen); // shorten by 1 at every loop iteration
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        s_len = (str_len >= start + k) ? k - plen : str_len - start - plen;
        int_s = prefix2int(input, start + plen, s_len);
        ; // tail
        auto correct_bucket_idx = buckets_rs.rank(int_p + 1);
        if (buckets_rs.rank(int_p) < correct_bucket_idx)
        {
            // 1. prefix found
            const Bucket &correct_bucket = buckets[correct_bucket_idx - 1];
            FindPrefix(correct_bucket, plen, s_len, int_s, int_p, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
        }
        else
        {
            int_sp = int_p >> 2;
            FindShortFinimizer(plen - 1, int_sp, sB, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
        }
    }

    uint64_t ss = start;
    s_len = k - plen; // Constant s_len
    for (start = ss; start <= str_len - k; start++)
    {
        int_p = stream_kmer(int_p, input[start + plen - 1], plen);
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        int_s = stream_kmer(int_s, input[start + plen + s_len - 1], s_len);

        auto correct_bucket_idx = buckets_rs.rank(int_p + 1);
        if (buckets_rs.rank(int_p) < correct_bucket_idx)
        {
            // 1. prefix found
            const Bucket &correct_bucket = buckets[correct_bucket_idx - 1];
            FindPrefix(correct_bucket, plen, s_len, int_s, int_p, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
        }
        else
        {
            int_sp = int_p >> 2;
            FindShortFinimizer(plen - 1, int_sp, sB, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
        }
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); // , input);
        kmer_start++;
    }
    // Shorter s_len
    ss = start;
    for (start = ss; start < str_len - plen + 1; start++)
    {
        int_p = stream_kmer(int_p, input[start + plen - 1], plen);
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        s_len = str_len - start - plen;
        int_s &= ((1ULL << (2 * s_len)) - 1); // Shorten int_s by 2 at the beginning

        auto correct_bucket_idx = buckets_rs.rank(int_p + 1);
        if (buckets_rs.rank(int_p) < correct_bucket_idx)
        {
            // 1. prefix found
            const Bucket &correct_bucket = buckets[correct_bucket_idx - 1];
            FindPrefix(correct_bucket, plen, s_len, int_s, int_p, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
        }
        else
        {
            int_sp = int_p >> 2;
            FindShortFinimizer(plen - 1, int_sp, sB, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input);
        }
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); // , input);
        kmer_start++;
    }

    // The last plen-1 characters cannot contain a prefix
    ss = start;
    uint64_t s_plen = plen;
    for (start = ss; start < str_len; start++)
    {
        s_plen--;
        int_p &= ((1ULL << (2 * s_plen)) - 1);                                                                      // Shorten int_p by 2
        FindShortFinimizer(s_plen, int_p, sB, start, kmer_start + k - 1, curr_candidates, next_candidates, k_fmin); // , input); // this shortens s_plen by 1 internally

        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); // , input);
        kmer_start++;
    }
    return;
}

void only_concat_read_colors(const uint64_t *data, const uint64_t n_colors, vector<uint64_t> &results, const vector<pair<int64_t, uint64_t>> &fmin_v)
{
    for (const auto &[start, freq] : fmin_v)
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

            while (word != 0)
            {
                uint64_t bit = __builtin_ctzll(word);
                results[bit] += freq;
                word &= word - 1;
            }

            ptr++;
            color_id += bits_to_read;
            bits_left -= bits_to_read;
        }

        // 2. Read aligned words in btw

        while (bits_left >= 64)
        {
            uint64_t word = *ptr++;
            for (uint64_t w = word; w != 0;)
            {
                uint64_t bit = __builtin_ctzll(w);
                results[color_id + bit] += freq;
                w &= w - 1;
            }
            color_id += 64;
            bits_left -= 64;
        }

        // 3. Read the last word (if any)
        if (bits_left > 0)
        {
            uint64_t mask = ((1ULL << bits_left) - 1);
            uint64_t word = *ptr & mask;

            while (word != 0)
            {
                uint64_t bit = __builtin_ctzll(word);
                results[color_id + bit] += freq;
                word &= word - 1;
            }
        }
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

// Not used now
// works for negative values in results
void counting_sort_neg(const vector<int16_t> &results, vector<pair<uint16_t, int16_t>> &ans, const size_t found_fmin, const uint16_t n_colors)
{
    if (results.empty())
        return;

    // min & max
    int16_t min_val = *min_element(results.begin(), results.end());
    int16_t max_val = *max_element(results.begin(), results.end());

    int16_t range = static_cast<size_t>(max_val - min_val + 1);
    vector<int16_t> counts(range, 0);

    for (uint16_t idx = 0; idx < n_colors; idx++)
    {
        counts[static_cast<int16_t>(results[idx] - min_val)]++; // offset by min_val
    }

    // Cumulative sums
    for (size_t i = 1; i < range; i++)
    {
        counts[i] += counts[i - 1];
    }

    for (uint16_t idx = n_colors; idx-- > 0;)
    {
        int16_t val = results[idx];
        size_t pos = counts[static_cast<int16_t>(val - min_val)] - 1;
        ans[pos] = {static_cast<uint16_t>(idx), val};
        counts[static_cast<int16_t>(val - min_val)]--;
    }
}
