#pragma once

#include <vector>
#include <optional>
#include <unordered_map>
#include <chrono>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>

#include "sdsl/bit_vectors.hpp"
#include "rarest_fmin_search.hh"
#include "Buckets.hh"
#include "Color_Set.hh"
#include "Color_Set_Storage.hh"

using namespace std;
using CCS_t = Color_Set_Storage<SDSL_Variant_Color_Set>;


static std::chrono::nanoseconds time_rarest_fmin(0);
static std::chrono::nanoseconds time_rarest_fmin_rc(0);
static std::chrono::nanoseconds time_combine(0);
static std::chrono::nanoseconds time_index_loading(0);
static std::chrono::nanoseconds time_output(0);

void print_search_timing_stats() {
    using namespace std::chrono;
    std::cerr << "Time to load the index: "
                << duration_cast<milliseconds>(time_index_loading).count() << " ms\n";

    std::cerr << "Time in rarest_fmin_streaming_search (fwd): "
                << duration_cast<milliseconds>(time_rarest_fmin).count() << " ms\n";
    std::cerr << "Time in combine_f_rc: "
                << duration_cast<milliseconds>(time_combine).count() << " ms\n";
    std::cerr << "Time to output results: "
                << duration_cast<milliseconds>(time_output).count() << " ms\n";
}

// Colored finimizers without much compression
class ColoredFinimizers {
public:

    vector<char> concat; // Finimizers concatenated in lexicographic order (ASCII characters).
    vector<uint8_t> lengths;
    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors.
    vector<uint8_t> lengths_by_freq;

    // Loads from the format output by the Rust CLI command `finimizer_matrix` with option --reverse.
    // That format has colexicographically sorted reverse finimizers. We reverse them to
    // get lex-sorted finimizers.
    void load(std::istream& in) {
        cerr << "Loading uncompressed tails" << endl;
        uint64_t n_finimizers;
        in.read(reinterpret_cast<char*>(&n_finimizers), sizeof(n_finimizers));

        uint64_t finimizer_total_length;
        in.read(reinterpret_cast<char*>(&finimizer_total_length), sizeof(finimizer_total_length));

        uint64_t n_colors;
        in.read(reinterpret_cast<char*>(&n_colors), sizeof(n_colors));

        lengths.resize(n_finimizers);
        in.read(reinterpret_cast<char*>(lengths.data()), n_finimizers * sizeof(uint8_t));

        concat.resize(finimizer_total_length);
        in.read(reinterpret_cast<char*>(concat.data()), finimizer_total_length * sizeof(char));

        int64_t n_bits = n_finimizers * n_colors;
        // The bits are in u64 Lsb format
        vector<uint64_t> words((n_bits + 63) / 64); // Ceil div by 64
        in.read(reinterpret_cast<char*>(words.data()), words.size() * sizeof(uint64_t));
        color_sets_concat.resize(n_bits);
        for(int64_t i = 0; i < words.size(); i++) {
            color_sets_concat.set_int(i*64, words[i]);
        }

        uint64_t lengths_by_freq_size;
        in.read(reinterpret_cast<char*>(&lengths_by_freq_size), sizeof(lengths_by_freq_size));
        lengths_by_freq.resize(lengths_by_freq_size);
        in.read(reinterpret_cast<char*>(lengths_by_freq.data()), lengths_by_freq_size * sizeof(uint8_t));

        cerr << "Reversing finimizer strings" << endl;
        int64_t start_in_concat = 0;
        for(uint64_t f_idx = 0; f_idx < lengths.size(); f_idx++) {
            int64_t s = start_in_concat;
            int64_t e = start_in_concat + lengths[f_idx];
            std::reverse(concat.begin() + s, concat.begin() + e);
            start_in_concat = e;
        }
        cerr << "Uncompressed index loaded" << endl;
    }

};



void true_or_crash(bool b, const char* error_message){
    if(!b){
        cerr << error_message << endl;
        exit(1);
    }
}

// Forward declarations for the combine_f_rc variants (implemented below)
void pseudoalignment_stats( vector<int64_t>& Fmin,
                  const class Color_Set_Storage<SDSL_Variant_Color_Set>& CCS,
                  const uint64_t n_colors,
                  std::vector<std::pair<uint16_t, uint16_t>>& ans);

uint16_t pseudoalignment_stats( vector<int64_t>& Fmin,
                      const class Color_Set_Storage<SDSL_Variant_Color_Set>& CCS,
                      const uint64_t n_colors,
                      std::vector<std::pair<uint16_t, uint16_t>>& ans,
                      const float t);

class CompressedColoredFinimizers {

private:
    sdsl::bit_vector unique_color_sets;

    // color_set_ids maps original finimizer index -> set id in storage
    vector<uint64_t> color_set_ids;

public:
    using CCS_t = Color_Set_Storage<SDSL_Variant_Color_Set>;
    CCS_t CCS; // storage backend (replaces old CompressedColorSets)

    vector<optional<Bucket>> buckets;
    unordered_map<uint32_t, pair<uint8_t, int64_t>> sB; // finimizers shorter than prefix length
    uint64_t n_colors = 0;
    uint64_t n_finimizers = 0;
    uint64_t plen = 0;
    uint64_t k = 0;

    int get_k() const { return k; }

    CompressedColoredFinimizers() = default;

    // Constructor: build the storage from deduplicated_cs (same style as before)
    CompressedColoredFinimizers(ColoredFinimizers&& cf, int64_t prefix_len, uint64_t kmer_size) {
        cerr << "Let's compress it!"<< endl;
        plen = prefix_len;
        cerr << "prefix length: "<< plen<< endl;
        k = kmer_size;
        cerr << "k-mer size: "<< k << endl;

        uint64_t n_buckets = (1ULL << (prefix_len * 2));
        buckets.resize(n_buckets);
        cerr << "total buckets: "<< (int)n_buckets << endl;

        n_finimizers = cf.lengths.size();
        cerr << "total finimizers: "<< (int)n_finimizers << endl;
        true_or_crash(n_finimizers > 0, "ERROR: 0 finimizers");

        true_or_crash(cf.color_sets_concat.size() % n_finimizers == 0, "ERROR: color set bitmap length not divisible by finimizer count");
        n_colors = cf.color_sets_concat.size() / n_finimizers;
        cerr << "n_colors: "<< (int)n_colors << endl;
        

        // TODO duplicate colorsets
        // Deduplicate color sets (same approach as before)
        cerr << "Deduplicate color sets" << endl;
        unordered_map<string, vector<size_t>> deduplicated_cs; // {cs:[fmin indices]}
        for (size_t i = 0; i < n_finimizers; i++) {
            sdsl::bit_vector bv(n_colors);
            for (size_t j = 0; j < n_colors; j++) {
                bv[j] = cf.color_sets_concat[i * n_colors + j];
            }
            string key((char*)bv.data(), ((n_colors + 63) / 64) * 8);
            deduplicated_cs[key].push_back(i);
        }

        this->color_set_ids.resize(n_finimizers);

        // Build CCS storage from deduplicated_cs
        // For each deduplicated pattern, create vector<int64_t> of set members and add_set to storage
        size_t next_set_id = 0;
        for (const auto& [key, old_offsets] : deduplicated_cs) {
            // reconstruct bitvector and push indices of 1-bits into vector<int64_t>
            sdsl::bit_vector bv(n_colors);
            // safe memcpy: key was formed from bv.data() earlier
            memcpy((char*)bv.data(), key.data(), key.size());

            vector<int64_t> members;
            members.reserve(sdsl::util::cnt_one_bits(bv));
            for (size_t c = 0; c < n_colors; ++c) {
                if (bv[c]) members.push_back((int64_t)c);
            }

            // Add set to storage
            CCS.add_set(members);

            // Assign color_set_ids for each original finimizer that had this pattern
            for (auto idx : old_offsets) {
                this->color_set_ids[idx] = next_set_id;
            }

            ++next_set_id;
        }

        // Finalize storage for queries
        CCS.prepare_for_queries();
        cerr << "Number of stored (deduplicated) color sets: " << next_set_id << endl;

        // Assign to final structure (done)

        // Deal with tails and buckets (same logic as before)
        cerr << "Deal with tails" << endl;
        int64_t first_nonegative_tail_idx = -1;
        int64_t f_start = 0;
        for(int64_t i = 0; i < (int64_t)n_finimizers; i++){
            if(cf.lengths[i] >= plen) {
                first_nonegative_tail_idx = i;
                break;
            }
            f_start += cf.lengths[i];
        }
        true_or_crash(first_nonegative_tail_idx >= 0, "ERROR: all tails shorter than prefix length");

        std::string_view cur_prefix(cf.concat.data() + f_start, plen);
        vector<std::string_view> cur_tails;
        vector<uint32_t> cur_color_set_ids;

        uint64_t p_int = prefix2int(cur_prefix, 0, plen);

        f_start = 0; // Go back to zero

        for (int64_t i = 0; i < (int64_t)n_finimizers; i++) {
            if(cf.lengths[i] < plen){
                std::string_view sprefix(cf.concat.data() + f_start, cf.lengths[i]);
                uint64_t sp_int = prefix2int(sprefix,0, cf.lengths[i]);
                sB[sp_int]= {cf.lengths[i],this->color_set_ids[i]};

            } else {
                std::string_view prefix(cf.concat.data() + f_start, plen);
                true_or_crash(f_start + plen <= (int64_t)cf.concat.size(),
                        "ERROR: out-of-bounds prefix access");

                if(prefix != cur_prefix) {
                    // Bucket changes -> encode currently collected tails
                    buckets[p_int]=Bucket(cur_tails, cur_color_set_ids, cf.lengths_by_freq);
                    p_int = prefix2int(prefix, 0, plen);
                    cur_tails.clear();
                    cur_color_set_ids.clear();
                }
                cur_tails.push_back(std::string_view(cf.concat.data() + f_start + plen, cf.lengths[i] - plen));
                cur_color_set_ids.push_back((uint32_t)this->color_set_ids[i]);
                cur_prefix = prefix;
            }
            f_start += cf.lengths[i];
        }

        if (!cur_tails.empty()) { // Last bucket
            buckets[p_int] = Bucket(cur_tails, cur_color_set_ids, cf.lengths_by_freq);
        }

    } // constructor end

    // Search interface (unchanged externally)
    void search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans) const {
        //cerr << "search"<< endl;
        const int64_t query_len = query.length();
        if (query_len < this->k) return;

        // Forward finimizer search
        vector<int64_t> Finimizers;
        Finimizers.reserve(query_len - k + 1);
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k, Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }

        {
            auto start = std::chrono::high_resolution_clock::now();
            pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, ans);
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
        //cerr << "search end" << endl;
    }

    // Threshold-based search: returns minimum value and fills ans
    uint16_t search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans, const float& t) const {
        //cerr << "search"<< endl;

        const int64_t query_len = query.length();
        if (query_len < this->k) return 0;

        // Forward finimizer search
        vector<int64_t> Finimizers;
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k, Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }

        uint16_t min_value;
        {
            auto start = std::chrono::high_resolution_clock::now();
            min_value = pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, ans, t);
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
        //cerr << "search end" << endl;

        return min_value;
    }

    // Serialize / load (keeps same file layout for buckets, sB, meta; color sets are saved via CCS.serialize)
    void serialize(const string& index_prefix) const {
        cerr << "Save the index" << endl;

        // Save color-set storage
        std::ofstream cs_out(index_prefix + ".colorsets.BIN", std::ios::binary);
        if (!cs_out) {
            std::cerr << "Error: Could not open colorsets file!" << std::endl;
            return;
        }
        CCS.serialize(cs_out);
        cs_out.close();

        // buckets (std::optional<Bucket>)
        std::ofstream buckets_out(index_prefix + ".buckets.BIN", std::ios::binary);
        if (!buckets_out) {
            std::cerr << "Error: Could not open buckets file!" << std::endl;
            return;
        }

        size_t num_buckets = buckets.size();
        buckets_out.write(reinterpret_cast<const char*>(&num_buckets), sizeof(num_buckets));
        for (const auto& bucket_opt : buckets) {
            bool present = bucket_opt.has_value();
            buckets_out.write(reinterpret_cast<const char*>(&present), sizeof(present));
            if (present) {
                bucket_opt->serialize(buckets_out);
            }
        }
        buckets_out.close();

        // sB
        std::ofstream sB_out(index_prefix + ".sB.BIN", std::ios::binary);
        if (!sB_out) {
            std::cerr << "Error: Could not open sB file!" << std::endl;
            return;
        }

        size_t map_size = sB.size();
        sB_out.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
        for (const auto& [key, val] : sB) {
            sB_out.write(reinterpret_cast<const char*>(&key), sizeof(uint32_t));
            sB_out.write(reinterpret_cast<const char*>(&val.first), sizeof(uint8_t));
            sB_out.write(reinterpret_cast<const char*>(&val.second), sizeof(int64_t));
        }
        sB_out.close();

        // metadata
        std::ofstream meta_out(index_prefix + ".meta", std::ios::binary);
        if (!meta_out) {
            std::cerr << "Error: Could not write metadata!" << std::endl;
            return;
        }
        meta_out.write(reinterpret_cast<const char*>(&n_colors), sizeof(n_colors));
        meta_out.write(reinterpret_cast<const char*>(&n_finimizers), sizeof(n_finimizers));
        meta_out.write(reinterpret_cast<const char*>(&plen), sizeof(plen));
        meta_out.write(reinterpret_cast<const char*>(&k), sizeof(k));
        meta_out.close();
        cerr << "DONE"<< endl;
    }

    void load(const string& index_prefix) {
        // color sets
        std::ifstream cs_in(index_prefix + ".colorsets.BIN", std::ios::binary);
        if (!cs_in) {
            std::cerr << "Error: Could not open colorsets file!" << std::endl;
            return;
        }
        CCS.load(cs_in);
        cs_in.close();

        // buckets
        std::ifstream buckets_in(index_prefix + ".buckets.BIN", std::ios::binary);
        if (!buckets_in) {
            std::cerr << "Error: Could not open buckets file!" << std::endl;
            return;
        }

        size_t num_buckets;
        buckets_in.read(reinterpret_cast<char*>(&num_buckets), sizeof(num_buckets));
        buckets.resize(num_buckets);
        for (size_t i = 0; i < num_buckets; ++i) {
            bool present;
            buckets_in.read(reinterpret_cast<char*>(&present), sizeof(present));
            if (present) {
                Bucket bucket;
                bucket.load(buckets_in);  // Make sure Bucket has a `load(std::istream&)` method
                buckets[i] = bucket;
            } else {
                buckets[i] = std::nullopt;
            }
        }
        buckets_in.close();

        // sB
        std::ifstream sB_in(index_prefix + ".sB.BIN", std::ios::binary);
        if (!sB_in) {
            std::cerr << "Error: Could not open sB file!" << std::endl;
            return;
        }

        size_t map_size;
        sB_in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
        sB.clear();
        for (size_t i = 0; i < map_size; ++i) {
            uint32_t key;
            std::pair<char,int64_t> val;
            sB_in.read(reinterpret_cast<char*>(&key), sizeof(uint32_t));
            sB_in.read(reinterpret_cast<char*>(&val.first), sizeof(uint8_t));
            sB_in.read(reinterpret_cast<char*>(&val.second), sizeof(int64_t));
            sB[key] = val;
        }
        sB_in.close();

        // metadata
        std::ifstream meta_in(index_prefix + ".meta", std::ios::binary);
        if (!meta_in) {
            std::cerr << "Error: Could not read metadata!" << std::endl;
            return;
        }
        meta_in.read(reinterpret_cast<char*>(&n_colors), sizeof(n_colors));
        meta_in.read(reinterpret_cast<char*>(&n_finimizers), sizeof(n_finimizers));
        meta_in.read(reinterpret_cast<char*>(&plen), sizeof(plen));
        meta_in.read(reinterpret_cast<char*>(&k), sizeof(k));
        meta_in.close();
    }

}; // class CompressedColoredFinimizers end


// ------------------------------ helper routines that use the storage abstraction ------------------------------

// Read from many single-set accesses (i_fmin_v: vector of {set_id, freq})
void read_colors(const CCS_t& CCS_storage, const uint64_t n_colors, vector<uint64_t>& results, const vector<pair<int64_t, uint64_t>>& fmin_v){
    //cerr << "read_colors (view-based)" << endl;

    for (const auto& [pos, freq] : fmin_v) {
        if (pos < 0) continue;
        // pos corresponds to deduplicated set id
        auto view = CCS_storage.get_color_set_by_id((int64_t)pos);
        // get colors and add
        vector<int64_t> colors = view.get_colors_as_vector();
        for (auto c : colors) {
            if ((uint64_t)c < n_colors) results[(size_t)c] += freq;
        }
    }
    //cerr << "end" << endl;
}


inline void pseudoalignment_stats(vector<int64_t>& Fmin, const CCS_t& CCS_storage, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
    vector<uint64_t> results(n_colors,0);
    std::sort(Fmin.begin(), Fmin.end()); 
    vector<pair<int64_t, uint64_t>> fmin_v;
    fmin_v.reserve(Fmin.size());

    for (size_t i = 0; i < Fmin.size();) {
        size_t j = i + 1;
        while (j < Fmin.size() && Fmin[j] == Fmin[i]) ++j;
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

    read_colors(CCS_storage, n_colors, results, fmin_v);

    const size_t found_fmin = Fmin.size(); // # total finimizers
    counting_sort(results, ans, found_fmin, n_colors);
    return;
}

inline uint16_t pseudoalignment_stats(vector<int64_t>& Fmin, const CCS_t& CCS_storage, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans,  const float t ){ //vector<uint64_t>& results,
    //if (Fmin.empty()){return 0;}
    
    vector<uint64_t> results(n_colors,0);
    /* if (results.size() != n_colors) {
        results.assign(n_colors, 0);  
    } else {
        std::fill(results.begin(), results.end(), 0);
    } */

 
    std::sort(Fmin.begin(), Fmin.end()); 
    vector<pair<int64_t, uint64_t>> fmin_v;
    fmin_v.reserve(Fmin.size());

    for (size_t i = 0; i < Fmin.size();) {
        size_t j = i + 1;
        while (j < Fmin.size() && Fmin[j] == Fmin[i]) ++j;
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

    read_colors(CCS_storage, n_colors, results, fmin_v);

    // Check the values above the minimum
    const size_t found_fmin = Fmin.size(); // # total finimizers

    counting_sort(results, ans, found_fmin, n_colors);
    
    // TODO use min value here?
    const uint64_t min_value = found_fmin * t;

    return min_value;
}