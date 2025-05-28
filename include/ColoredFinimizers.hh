#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <optional>

#include "sdsl/bit_vectors.hpp"
#include "common.hh"
#include "rarest_fmin_search.hh"
#include "Buckets.hh"


using namespace std;

// Colored finimizers without much compression
class ColoredFinimizers {
public:

    vector<char> concat; // Finimizers concatenated in lexicographic order (ASCII characters).
    vector<uint8_t> lengths;
    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors.

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



void true_or_crash(bool b, char* error_message){
    if(!b){
        cerr << error_message << endl;
        exit(1);
    }
}

class CompressedColoredFinimizers {

public:

    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors. TODO: deduplicate.
    vector<optional<Bucket>> buckets; 

    std::unordered_map<uint32_t, pair<char, int64_t>> sB; // Create a hash table to store the finimizers shorter than the prefix length
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

        int64_t first_nonegative_tail_idx = -1;
        int64_t f_start = 0;
        for(int64_t i = 0; i < n_finimizers; i++){
            if(cf.lengths[i] >= prefix_len) {
                first_nonegative_tail_idx = i;
                break;
            }
            f_start += cf.lengths[i];
        }
        true_or_crash(first_nonegative_tail_idx >= 0, "ERROR: all tails shorter than prefix length");

        std::string_view cur_prefix(cf.concat.data() + f_start, prefix_len);
        vector<std::string_view> cur_tails;
        vector<uint32_t> cur_color_set_ids;

        uint64_t p_int = prefix2int(cur_prefix, 0, prefix_len);

        f_start = 0; // Go back to zero

        //set<string> vv = get_substrings((string)"CGCCGTGTCGATGGAGGCGCATTATAGGGAG");
        //set<string> vv = get_substrings((string)"AGTGGGTGCTCACCTTGACGGCTTTGCAGCG");
        
        
        /* for (auto v: vv){
            cerr << v << endl;
        }
        cerr << (int)vv.size()<< endl; */
        for(int64_t i = 0; i < n_finimizers; i++) {
            std::string_view fmin(cf.concat.data() + f_start, cf.lengths[i]);
            /* if (std::find(vv.begin(), vv.end(), std::string(fmin)) != vv.end()){
                cerr << fmin << " EXISTS !!"<< endl;
            } */
            if(cf.lengths[i] < prefix_len){
                std::string_view sprefix(cf.concat.data() + f_start, cf.lengths[i]);
                uint64_t sp_int = prefix2int(sprefix,0, cf.lengths[i]);
                //cerr << "Accessing sB (pos" << (int)p_int <<")...";
                sB[sp_int]= {cf.lengths[i],i}; // i= color_set_id
                //cerr << " ok"<< endl;

            } else {
                std::string_view prefix(cf.concat.data() + f_start, prefix_len);
                true_or_crash(f_start + prefix_len <= cf.concat.size(),
                        "ERROR: out-of-bounds prefix access");
                
                if(prefix != cur_prefix) {
                    // Bucket changes -> encode currently collected tails
                    buckets[p_int]=Bucket(cur_tails, cur_color_set_ids);
                    p_int = prefix2int(prefix, 0, prefix_len);
                    cur_tails.clear();
                    cur_color_set_ids.clear();
                }
                
                cur_tails.push_back(std::string_view(cf.concat.data() + f_start + prefix_len, cf.lengths[i] - prefix_len));
                cur_color_set_ids.push_back(i);
                cur_prefix = prefix;
            }
            f_start += cf.lengths[i];
        }

        if(cur_tails.size() > 0){ // Last bucket
            buckets[p_int]=Bucket(cur_tails, cur_color_set_ids);

        }
        color_sets_concat = std::move(cf.color_sets_concat);
    }

/*     // Transfer ownership of the index out of the builder
    unique_ptr<CompressedColoredFinimizers> get_index(){
        return std::move(this->index);
    } */
    void search(const std::string& query, vector<uint64_t>& results) const {
  
        const int64_t query_len = query.length();
      
        if (query.size() < this->k) return; 

        vector<uint64_t> Finimizers = rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k);
      
        // Check the colors for every finimizer found
        pseudoalignemnt_stats(Finimizers, this->color_sets_concat, this->n_colors, results);
        return;
    }

    void search(const std::string& query, vector<float>& results, const float& t) const {
        
        const int64_t query_len = query.length();

        if (query.size() < this->k) return; 

        vector<uint64_t> Finimizers = rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k);

        // Check the colors for every finimizer found
        pseudoalignemnt_stats(Finimizers, this->color_sets_concat, this->n_colors, results, t);
        return;
    }

    void serialize(const string& index_prefix) const {
        cerr << "Save the index"<< endl;
        // color_sets_concat
        std::ofstream colors_out(index_prefix + ".colors.sdsl", std::ios::binary);
        if (!colors_out) {
            std::cerr << "Error: Could not open colors file!" << std::endl;
            return;
        }
        sdsl::serialize(color_sets_concat, colors_out);
        colors_out.close();

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
