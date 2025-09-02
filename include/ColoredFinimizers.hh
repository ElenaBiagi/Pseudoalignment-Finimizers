#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <optional>
#include <bit>
#include <bitset>
#include <unordered_map>
#include <algorithm>

#include "sdsl/bit_vectors.hpp"
#include "common.hh"
#include "rarest_fmin_search.hh"
#include "Buckets.hh"

#include <chrono>


using namespace std;

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
    std::cerr << "Time in rarest_fmin_streaming_search (rev): "
                << duration_cast<milliseconds>(time_rarest_fmin_rc).count() << " ms\n";
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
        for(int64_t f_idx = 0; f_idx < lengths.size(); f_idx++) {
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

struct ColorSet {
    enum class ColorSetType : uint8_t { SPARSE, BITMAP, COMPLEMENT };
    // 3 different representations
    
    ColorSetType type;
    sdsl::bit_vector bitmap;
    vector<uint32_t> gaps;
    size_t size = 0;

    ColorSet() = default;

    ColorSet(const sdsl::bit_vector& bv) { // bv = color_set_concat_id

        size = sdsl::util::cnt_one_bits(bv);
        size_t n_colors = bv.size();
        size_t sparse_thr = n_colors / 4; //*0.25
        size_t dense_thr  = (3 * n_colors) / 4; // *0.75

        // Sparse
        if (size < sparse_thr) {
            type = ColorSetType::SPARSE;
            // TODO use sdsl::EliasFano?
            uint32_t prev = 0;
            for (size_t c = 0; c < n_colors; c++) {
                if (bv[c]) { // store 1s
                    gaps.push_back((uint32_t)(c - prev)); // Elias-Fano
                    prev = c + 1;
                }
            }
        } 
        // Dense
        else if (size > dense_thr) {
            type = ColorSetType::COMPLEMENT;
            // TODO use sdsl::EliasFano?
            uint32_t prev = 0;
            for (size_t c = 0; c < n_colors; c++) {
                if (!bv[c]) { // store 0s
                    gaps.push_back((uint32_t)(c - prev));
                    prev = c + 1;
                }
            }
        } else {
            type = ColorSetType::BITMAP;
            bitmap = bv;
        }
    }

    sdsl::bit_vector to_bitvector(size_t n_colors) const {
        sdsl::bit_vector out(n_colors); // all 0s
        if (type == ColorSetType::BITMAP) {
            return bitmap;
        } else if (type == ColorSetType::SPARSE) {
            size_t pos = 0;
            for (auto g : gaps) {
                pos += g;
                out[pos] = 1;
                pos++;
            }
        } else {
            std::fill(out.begin(), out.end(), 1);
            size_t pos = 0;
            for (auto g : gaps) {
                pos += g;
                out[pos] = 0;
                pos++;
            }
        }
        return out;
    }

    void serialize(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&type), sizeof(type));
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        if (type == ColorSetType::BITMAP) {
            sdsl::serialize(bitmap, out);
        } else {
            size_t len = gaps.size();
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            out.write(reinterpret_cast<const char*>(gaps.data()), len * sizeof(uint32_t));
        }
    }

    void load(std::istream& in, size_t n_colors) {
        in.read(reinterpret_cast<char*>(&type), sizeof(type));
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        if (type == ColorSetType::BITMAP) {
            sdsl::load(bitmap, in);
        } else {
            size_t len;
            in.read(reinterpret_cast<char*>(&len), sizeof(len));
            gaps.resize(len);
            in.read(reinterpret_cast<char*>(gaps.data()), len * sizeof(uint32_t));
        }
    }
};

void combine_f_rc(const std::vector<int64_t>& Fmin,
                  const std::vector<int64_t>& r_Fmin,
                  //const std::vector<ColorSet>& unique_color_sets,
                  //const std::vector<uint32_t>& color_set_ids,
                  const sdsl::bit_vector& color_sets_concat,
                  const uint64_t n_colors,
                  std::vector<std::pair<uint16_t, uint16_t>>& ans);

uint64_t combine_f_rc(const std::vector<int64_t>& Fmin,
                      const std::vector<int64_t>& r_Fmin,
                      //const std::vector<ColorSet>& unique_color_sets,
                      //const std::vector<uint32_t>& color_set_ids,
                      const sdsl::bit_vector& color_sets_concat,
                      const uint64_t n_colors,
                      std::vector<std::pair<uint16_t, uint16_t>>& ans,
                      const float& t);


class CompressedColoredFinimizers {
public:
    enum class ColorSetType : uint8_t { SPARSE, BITMAP, COMPLEMENT };
    

private:
    //vector<ColorSet> unique_color_sets;
    sdsl::bit_vector unique_color_sets;

    vector<uint32_t> color_set_ids;

public:
    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors. TODO: deduplicate.

    vector<optional<Bucket>> buckets; 
    unordered_map<uint32_t, pair<char, int64_t>> sB; // Create a hash table to store the finimizers shorter than the prefix length
    uint64_t n_colors;
    uint64_t n_finimizers;
    uint64_t plen;
    uint64_t k;

    int get_k() const{return k;} 

    CompressedColoredFinimizers() = default;

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

        /* unique_color_sets.reserve(deduplicated_cs.size());
        color_set_ids.resize(n_finimizers); // ids sorted based on the frequency of fmins length
        for (auto& [key, old_offsets] : deduplicated_cs) {
            sdsl::bit_vector bv(n_colors);
            memcpy((char*)bv.data(), key.data(), key.size()); // rebuild the bitvector from the string(hash key)
            unique_color_sets.emplace_back(bv); // bv stored as ColorSet
            uint32_t new_id = (uint32_t)(unique_color_sets.size() - 1); // -1 (zero indexed)
            // Assign to every id in the vector of indices its new color set id
            for (auto old_id : old_offsets) color_set_ids[old_id] = new_id;
        } */

        sdsl::bit_vector unique_color_sets(deduplicated_cs.size() * n_colors);

        color_set_ids.resize(n_finimizers); // ids sorted based on the frequency of fmins length
        uint64_t new_offset = 0;

        for (auto& [key, old_offsets] : deduplicated_cs) {


            sdsl::bit_vector bv(n_colors);
            memcpy((char*)bv.data(), key.data(), key.size()); // rebuild the bitvector from the string(hash key)
            // For every old offset write the new offset
            // Then write the color_id bitvector once, at new_offset  
            uint64_t new_color_id_offset = new_offset / n_colors;

            for (auto& c_id : old_offsets){
                color_set_ids[c_id] = new_color_id_offset;
                
            }
            uint64_t old_offset = old_offsets[0] * n_colors;

            // Copy the bv in unique_color_sets
            // TODO do this more efficiently
            for (size_t j = 0; j < n_colors; ++j) {
                unique_color_sets[new_offset + j] = cf.color_sets_concat[old_offset + j];
            }
            new_offset += n_colors;
        }
        // Assign to final structure
        this->color_sets_concat = std::move(unique_color_sets);

        cerr << "Deal with tails" << endl;
        int64_t first_nonegative_tail_idx = -1;
        int64_t f_start = 0;
        for(int64_t i = 0; i < n_finimizers; i++){
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

        for (int64_t i = 0; i < n_finimizers; i++) {
            if(cf.lengths[i] < plen){
                std::string_view sprefix(cf.concat.data() + f_start, cf.lengths[i]);
                uint64_t sp_int = prefix2int(sprefix,0, cf.lengths[i]);
                sB[sp_int]= {cf.lengths[i],color_set_ids[i]};

            } else {
                std::string_view prefix(cf.concat.data() + f_start, plen);
                true_or_crash(f_start + plen <= cf.concat.size(),
                        "ERROR: out-of-bounds prefix access");
                
                if(prefix != cur_prefix) {
                    // Bucket changes -> encode currently collected tails
                    buckets[p_int]=Bucket(cur_tails, cur_color_set_ids, cf.lengths_by_freq);
                    p_int = prefix2int(prefix, 0, plen);
                    cur_tails.clear();
                    cur_color_set_ids.clear();
                }
                cur_tails.push_back(std::string_view(cf.concat.data() + f_start + plen, cf.lengths[i] - plen));
                cur_color_set_ids.push_back(color_set_ids[i]);
                cur_prefix = prefix;
            }
            f_start += cf.lengths[i];
        }

        if (!cur_tails.empty()) { // Last bucket
            buckets[p_int] = Bucket(cur_tails, cur_color_set_ids, cf.lengths_by_freq);
        }

    }

    /* const ColorSet& get_color_set(size_t finimizer_id) const {
        return unique_color_sets[color_set_ids[finimizer_id]];
    } */

    void search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans) const {

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

        // Reverse complement search
        vector<int64_t> r_Finimizers;
        r_Finimizers.reserve(query_len - k + 1);
        string r_query = sbwt::get_rc(query);
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(r_query, this->buckets, this->sB, this->plen, this->k, r_Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin_rc += (end - start);
        }

        // Combine forward and reverse results
        {
            auto start = std::chrono::high_resolution_clock::now();
            // combine_f_rc(Finimizers, r_Finimizers, unique_color_sets, color_set_ids, this->n_colors, ans);
            combine_f_rc(Finimizers, r_Finimizers, this->color_sets_concat, this->n_colors, ans);
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
    }

    // Threshold-based search: returns minimum value and fills ans
    uint16_t search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans, const float& t) const {

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

        // reverse complement
        vector<int64_t> r_Finimizers;
        string r_query = sbwt::get_rc(query);
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(r_query, this->buckets, this->sB, this->plen, this->k, r_Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin_rc += (end - start);
        }

        // Combine the results of finimizers color ids for forward and reverse
        uint16_t min_value;
        {
            auto start = std::chrono::high_resolution_clock::now();
            min_value = combine_f_rc(Finimizers, r_Finimizers, this->color_sets_concat, this->n_colors, ans, t);
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }

        return min_value;
    }


    void serialize(const string& index_prefix) const {
        cerr << "Save the index" << endl;
        // color_sets_concat
        std::ofstream colors_out(index_prefix + ".colors.sdsl", std::ios::binary);
        if (!colors_out) {
            std::cerr << "Error: Could not open colors file!" << std::endl;
            return;
        }
        sdsl::serialize(color_sets_concat, colors_out);
        colors_out.close();
        /* ofstream colors_out(index_prefix + ".colors.hybrid", ios::binary);
        size_t num_sets = unique_color_sets.size();
        colors_out.write(reinterpret_cast<const char*>(&num_sets), sizeof(num_sets));
        for (auto& cs : unique_color_sets) {
            cs.serialize(colors_out);
        }
        colors_out.close(); */

        /* std::ofstream ids_out(index_prefix + ".color_set_ids.BIN", std::ios::binary);
        size_t ids_size = color_set_ids.size();
        ids_out.write(reinterpret_cast<const char*>(&ids_size), sizeof(ids_size));
        ids_out.write(reinterpret_cast<const char*>(color_set_ids.data()), ids_size * sizeof(uint32_t));
        ids_out.close(); */

        

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
            sB_out.write(reinterpret_cast<const char*>(&val.first), sizeof(char));
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
        // color_sets_concat
        std::ifstream colors_in(index_prefix + ".colors.sdsl", std::ios::binary);
        if (!colors_in) {
            std::cerr << "Error: Could not open colors file!" << std::endl;
            return;
        }
        sdsl::load(color_sets_concat, colors_in);
        colors_in.close();

        /* std::ifstream colors_in(index_prefix + ".colors.hybrid", ios::binary);
        size_t num_sets;
        colors_in.read(reinterpret_cast<char*>(&num_sets), sizeof(num_sets));
        unique_color_sets.resize(num_sets);
        for (size_t i = 0; i < num_sets; i++) {
            unique_color_sets[i].load(colors_in, n_colors);
        }
        colors_in.close(); */

        // color set ids are in buckets and sB
        std::ifstream ids_in(index_prefix + ".color_set_ids.BIN", std::ios::binary);
        size_t ids_size;
        ids_in.read(reinterpret_cast<char*>(&ids_size), sizeof(ids_size));
        color_set_ids.resize(ids_size);
        ids_in.read(reinterpret_cast<char*>(color_set_ids.data()), ids_size * sizeof(uint32_t));
        ids_in.close();

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
            sB_in.read(reinterpret_cast<char*>(&val.first), sizeof(char));
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

};

// TODO: how to sort these or exploit identical ones still keeping the pairs?
// two overlpaiing k-mers are likely to have the same fmin so they are likely to share the same fmin on both strands
// This should anyways keep the number of false pos low

void combine_f_rc(const vector<int64_t>& Fmin, const vector<int64_t>& r_Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
    // NEW pseudoaligment_stats
    vector<uint64_t> results;
    results.resize(n_colors, 0);

    const uint64_t* data = color_sets_concat.data();
    // for now assume they have the same size

    //option1 
    // check the ones that are the same and deal with only the others later
    // store the diff ones in a separate vector??
    // you can still sort Fmin -> correct results
    // How to deal with the reverse??? sort r_Fmin based on Fmin sorting

    
    const int n_fmin = Fmin.size(); // now this is input_len -k +1
    // 1. Make a vector of pairs and a vector of int64_t
    vector<pair<int64_t, int64_t>> p_Fmin; // pair of distinct color_sets ids
    p_Fmin.reserve(n_fmin);

    vector<int64_t> i_Fmin; // if colorset id is the same (identical) for forward and reverse
    i_Fmin.reserve(n_fmin);

    for (auto i = 0; i<n_fmin; i++){
        int64_t f = Fmin[i];
        int64_t r = r_Fmin[n_fmin - i - 1];
        // r string is the reverse complement
        if (f != r && f != -1 && r != -1){
            p_Fmin.push_back({f, r});
        }
        else if (f != -1){
            i_Fmin.push_back(f);
        } 
        else if (r != -1){
            i_Fmin.push_back(r);
        }
    }
    // cases:
    // == -1 -1 : nothing
    // == x x : store x
    // != -1 x : store x
    // != x -1 : store x
    // != x y : store x,y

    const size_t found_fmin = p_Fmin.size() + i_Fmin.size();
    //CSS.init(found_fmin, n_colors); // only here we know found_fmin

    // 2. Deal with the vector of int64_t: i_Fmin
        // 2a. Keep frequency 
        // Count freq of each i_fmin
    std::unordered_map<int64_t, uint64_t> i_fmin_counts;
    i_fmin_counts.reserve(i_Fmin.size());
    for (const auto& v: i_Fmin) {
        i_fmin_counts[v]++;
    }

    // 2b. sort i_fmin
    vector<pair<int64_t, uint64_t>> i_fmin_v(i_fmin_counts.begin(), i_fmin_counts.end());
    std::sort(i_fmin_v.begin(), i_fmin_v.end()); 
    read_colors(data, n_colors, results, i_fmin_v);

    // 3. Deal with the vector of pairs: p_Fmin
    struct pair_hash {
        size_t operator()(const pair<int64_t, int64_t>& p) const {
            return std::hash<int64_t>()(p.first) ^ (std::hash<int64_t>()(p.second) << 1);
        }
    };
    std::unordered_map<pair<int64_t,int64_t>, uint64_t, pair_hash> p_fmin_counts;
    p_fmin_counts.reserve(p_Fmin.size());

    for (const auto& v: p_Fmin) {
        p_fmin_counts[v]++;
    }
    vector<pair<pair<uint64_t, uint64_t>, uint64_t>> p_fmin_v(p_fmin_counts.begin(), p_fmin_counts.end()); // [{{f,r},counts},...]

    // TODO does it make sense to sort pairs??

    read_f_rc_colors(data, n_colors, results, p_fmin_v);

    counting_sort(results, ans, found_fmin, n_colors);
    return;
}

uint64_t combine_f_rc(const vector<int64_t>& Fmin, const vector<int64_t>& r_Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans, const float& t){

    // NEW pseudoaligment_stats
    vector<uint64_t> results;
    results.resize(n_colors, 0);

    const uint64_t* data = color_sets_concat.data();


    const int n_fmin = Fmin.size(); // now this is input_len -k +1
    // 1. Make a vector of pairs and a vector of int64_t
    vector<pair<int64_t, int64_t>> p_Fmin;
    p_Fmin.reserve(n_fmin);

    vector<int64_t> i_Fmin; // if colorset id is the same for forward and reverse
    i_Fmin.reserve(n_fmin);

    for (auto i = 0; i<n_fmin; i++){
        int64_t f = Fmin[i];
        int64_t r = r_Fmin[n_fmin - i - 1];
        if (f != r && f != -1 && r != -1){
            p_Fmin.push_back({f, r});
        }
        else if (f != -1){
            i_Fmin.push_back(f);
        } 
        else if (r != -1){
            i_Fmin.push_back(r);
        }
    }
    const size_t found_fmin = p_Fmin.size() + i_Fmin.size();
    //CSS.init(found_fmin, n_colors); // only here we know found_fmin


    // 2. Deal with the vector of int64_t: i_Fmin
        // 2a. Keep frequency 
        // Count freq of each i_fmin
    std::unordered_map<int64_t, uint64_t> i_fmin_counts;
    i_fmin_counts.reserve(i_Fmin.size());
    for (const auto& v: i_Fmin) {
        i_fmin_counts[v]++;
    }
        // 2b. sort i_fmin
    vector<pair<int64_t, uint64_t>> i_fmin_v(i_fmin_counts.begin(), i_fmin_counts.end());
    std::sort(i_fmin_v.begin(), i_fmin_v.end()); 
    read_colors(data, n_colors, results, i_fmin_v);


    // 3. Deal with the vector of pairs: p_Fmin
    struct pair_hash {
        size_t operator()(const pair<int64_t, int64_t>& p) const {
            return std::hash<int64_t>()(p.first) ^ (std::hash<int64_t>()(p.second) << 1);
        }
    };
    std::unordered_map<pair<int64_t,int64_t>, uint64_t, pair_hash> p_fmin_counts;
    p_fmin_counts.reserve(p_Fmin.size());

    for (const auto& v: p_Fmin) {
        p_fmin_counts[v]++;
    }
    vector<pair<pair<uint64_t, uint64_t>, uint64_t>> p_fmin_v(p_fmin_counts.begin(), p_fmin_counts.end());

    // TODO does it make sense to sort pairs??
    read_f_rc_colors(data, n_colors, results, p_fmin_v);

    counting_sort(results, ans, found_fmin, n_colors);

    const uint64_t min_value = found_fmin * t;

    return min_value;
}

// void combine_f_rc(const vector<int64_t>& Fmin, const vector<int64_t>& r_Fmin, const vector<ColorSet>& unique_color_sets, const vector<uint32_t>& color_set_ids, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
//     // NEW pseudoaligment_stats
//     vector<uint64_t> results;
//     results.resize(n_colors, 0);
//     // for now assume they have the same size

//     //option1 
//     // check the ones that are the same and deal with only the others later
//     // store the diff ones in a separate vector
    
//     const int n_fmin = Fmin.size(); // now this is input_len -k +1
//     // 1. Make a vector of pairs and a vector of int64_t
//     vector<pair<int64_t, int64_t>> p_Fmin; // pair of distinct color_sets ids
//     p_Fmin.reserve(n_fmin);

//     vector<int64_t> i_Fmin; // if colorset id is the same (identical) for forward and reverse
//     i_Fmin.reserve(n_fmin);

//     for (auto i = 0; i<n_fmin; i++){
//         int64_t f = Fmin[i];
//         int64_t r = r_Fmin[n_fmin - i - 1];
//         // r string is the reverse complement
//         if (f != r && f != -1 && r != -1){
//             p_Fmin.push_back({f, r});
//         }
//         else if (f != -1){
//             i_Fmin.push_back(f);
//         } 
//         else if (r != -1){
//             i_Fmin.push_back(r);
//         }
//     }
//     // cases:
//     // == -1 -1 : nothing
//     // == x x : store x
//     // != -1 x : store x
//     // != x -1 : store x
//     // != x y : store x,y

//     const size_t found_fmin = p_Fmin.size() + i_Fmin.size();

//     // 2. Deal with the vector of int64_t: i_Fmin
//         // 2a. Keep frequency 
//         // Count freq of each i_fmin
//     std::unordered_map<int64_t, uint64_t> i_fmin_counts;
//     i_fmin_counts.reserve(i_Fmin.size());
//     for (const auto& v: i_Fmin) {
//         i_fmin_counts[v]++;
//     }

//     // 2b. sort i_fmin
//     vector<pair<int64_t, uint64_t>> i_fmin_v(i_fmin_counts.begin(), i_fmin_counts.end());
//     std::sort(i_fmin_v.begin(), i_fmin_v.end()); 

//     // Read colors for single finimizers
//     for (auto& [idx, count] : i_fmin_v) {
//         const auto& cs = unique_color_sets[idx].to_bitvector(n_colors);
//         for (size_t j = 0; j < n_colors; j++) {
//             if (cs[j]) results[j] += count;
//         }
//     }

//     // 3. Deal with the vector of pairs: p_Fmin
//     struct pair_hash {
//         size_t operator()(const pair<int64_t, int64_t>& p) const {
//             return std::hash<int64_t>()(p.first) ^ (std::hash<int64_t>()(p.second) << 1);
//         }
//     };

//     std::unordered_map<pair<int64_t,int64_t>, uint64_t, pair_hash> p_fmin_counts;
//     p_fmin_counts.reserve(p_Fmin.size());
//     for (const auto& v: p_Fmin) {
//         p_fmin_counts[v]++;
//     }
    
//     vector<pair<pair<uint64_t, uint64_t>, uint64_t>> p_fmin_v(p_fmin_counts.begin(), p_fmin_counts.end()); // [{{f,r},counts},...]

//     /* vector<pair<pair<uint64_t, uint64_t>, uint64_t>> p_fmin_v;
//     for (auto& [key, count] : p_fmin_counts) {
//         p_fmin_v.push_back({{key.first, key.second}, count});
//     } */


//     // Read colors for pairs of finimizers
//     for (auto& [ids, count] : p_fmin_v) {
//         const auto& cs1 = unique_color_sets[ids.first].to_bitvector(n_colors);
//         const auto& cs2 = unique_color_sets[ids.second].to_bitvector(n_colors);
//         for (size_t j = 0; j < n_colors; j++) {
//             if (cs1[j] || cs2[j]) results[j] += count;
//         }
//     }

//     counting_sort(results, ans, found_fmin, n_colors);

// }
// uint64_t combine_f_rc(const vector<int64_t>& Fmin, const vector<int64_t>& r_Fmin, const vector<ColorSet>& unique_color_sets, const vector<uint32_t>& color_set_ids, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans, const float& t){

//     // NEW pseudoaligment_stats
//     vector<uint64_t> results;
//     results.resize(n_colors, 0);


//     const int n_fmin = Fmin.size(); // now this is input_len -k +1
//     // 1. Make a vector of pairs and a vector of int64_t
//     vector<pair<int64_t, int64_t>> p_Fmin;
//     p_Fmin.reserve(n_fmin);

//     vector<int64_t> i_Fmin; // if colorset id is the same for forward and reverse
//     i_Fmin.reserve(n_fmin);

//     for (auto i = 0; i<n_fmin; i++){
//         int64_t f = Fmin[i];
//         int64_t r = r_Fmin[n_fmin - i - 1];
//         if (f != r && f != -1 && r != -1){
//             p_Fmin.push_back({f, r});
//         }
//         else if (f != -1){
//             i_Fmin.push_back(f);
//         } 
//         else if (r != -1){
//             i_Fmin.push_back(r);
//         }
//     }
//     const size_t found_fmin = p_Fmin.size() + i_Fmin.size();
//     //cerr << "found fmin = "<< found_fmin << endl;
//     //CSS.init(found_fmin, n_colors); // only here we know found_fmin


//     // 2. Deal with the vector of int64_t: i_Fmin
//         // 2a. Keep frequency 
//         // Count freq of each i_fmin
//     std::unordered_map<int64_t, uint64_t> i_fmin_counts;
//     i_fmin_counts.reserve(i_Fmin.size());
//     for (const auto& v: i_Fmin) {
//         i_fmin_counts[v]++;
//     }
//         // 2b. sort i_fmin
//     vector<pair<int64_t, uint64_t>> i_fmin_v(i_fmin_counts.begin(), i_fmin_counts.end());
//     std::sort(i_fmin_v.begin(), i_fmin_v.end()); 
    
//     //cerr << "Read colors for single finimizers"<< endl;
//     for (auto& [idx, count] : i_fmin_v) {
//         const auto& cs = unique_color_sets[idx].to_bitvector(n_colors);
//         for (size_t j = 0; j < n_colors; j++) {
//             if (cs[j]) results[j] += count;
//         }
//     }

//     // 3. Deal with the vector of pairs: p_Fmin
//     struct pair_hash {
//         size_t operator()(const pair<int64_t, int64_t>& p) const {
//             return std::hash<int64_t>()(p.first) ^ (std::hash<int64_t>()(p.second) << 1);
//         }
//     };
//     std::unordered_map<pair<int64_t,int64_t>, uint64_t, pair_hash> p_fmin_counts;
//     p_fmin_counts.reserve(p_Fmin.size());

//     for (const auto& v: p_Fmin) {
//         p_fmin_counts[v]++;
//     }
//     vector<pair<pair<uint64_t, uint64_t>, uint64_t>> p_fmin_v(p_fmin_counts.begin(), p_fmin_counts.end());

//     // TODO does it make sense to sort pairs??
    
//     //cerr << "Read colors for pairs of finimizers" << endl;
//     for (auto& [ids, count] : p_fmin_v) {
//         const auto& cs1 = unique_color_sets[ids.first].to_bitvector(n_colors);
//         const auto& cs2 = unique_color_sets[ids.second].to_bitvector(n_colors);
//         for (size_t j = 0; j < n_colors; j++) {
//             if (cs1[j] || cs2[j]) results[j] += count;
//         }
//     }
//     counting_sort(results, ans, found_fmin, n_colors);

//     const uint64_t min_value = found_fmin * t;

//     return min_value;

// }

