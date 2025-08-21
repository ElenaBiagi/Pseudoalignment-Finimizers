#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <optional>
#include <bit>
#include <bitset>

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



void true_or_crash(bool b, char* error_message){
    if(!b){
        cerr << error_message << endl;
        exit(1);
    }
}

class CompressedColoredFinimizers {
private:

    void subtract_and_filter(vector<uint8_t>& vec, uint64_t plen) {
        auto new_end = std::remove_if(vec.begin(), vec.end(), [plen](uint8_t& val) {
            if (val < plen) {return true;}     
            val -= plen;                    
            return false;// keep
        });
    }

    // TODO remove?
    /* vector<bit_vector> split_bitvector(const bit_vector& C, size_t n_colors) {
        size_t total_blocks = C.size() / n_colors;
        std::vector<bit_vector> result;
        result.reserve(total_blocks);

        for (size_t i = 0; i < total_blocks; ++i) {
            bit_vector block(n_colors);
            for (size_t j = 0; j < n_colors; ++j) {
                block[j] = C[i * n_colors + j];
            }
            result.push_back(std::move(block));
        }
        return result;
    } */

    struct BitVectorHash {
        std::size_t operator()(const sdsl::bit_vector& v) const noexcept {
            std::size_t h = 0;
            std::hash<uint64_t> hasher;
            // Access the packed 64-bit words directly
            const uint64_t* data = v.data();
            size_t n64 = (v.size() + 63) / 64;  // number of 64-bit words
            for (size_t i = 0; i < n64; ++i) {
                h ^= hasher(data[i]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            }
            return h;
        }
    };

    struct BitVectorEqual {
        bool operator()(const sdsl::bit_vector& a, const sdsl::bit_vector& b) const noexcept {
            if (a.size() != b.size()) return false;
            return std::equal(a.begin(), a.end(), b.begin());
        }
    };

    unordered_map<sdsl::bit_vector, vector<size_t>, BitVectorHash, BitVectorEqual> split_bitvector(const bit_vector& C, size_t n_colors) {
        
        // TODO split in a more efficient way
        size_t total_blocks = C.size() / n_colors;
        unordered_map<sdsl::bit_vector, vector<size_t>, BitVectorHash, BitVectorEqual> unique_color_ids;
        for (size_t i = 0; i < total_blocks; ++i) {
            bit_vector block(n_colors);
            uint64_t offset = i * n_colors;
            for (size_t j = 0; j < n_colors; ++j) {
                block[j] = C[offset + j];
            }
            unique_color_ids[block].push_back(i);
        }
        return unique_color_ids;
    }

    vector<sdsl::bit_vector> sort_bitvectors(unordered_map<sdsl::bit_vector, vector<size_t>, BitVectorHash, BitVectorEqual>& unique_color_ids) {

        // Sort the keys(bit_vectors) in a vector based on the number of 1s
        vector<sdsl::bit_vector> sorted_color_ids;
        sorted_color_ids.reserve(unique_color_ids.size());
        for (auto& K : unique_color_ids) {
            sorted_color_ids.push_back(K.first);
        }

        // Sort by number of 1s (descending) w popcount
        std::sort(sorted_color_ids.begin(), sorted_color_ids.end(),
            [](const sdsl::bit_vector& a, const sdsl::bit_vector& b) {
                return sdsl::util::cnt_one_bits(a) > sdsl::util::cnt_one_bits(b);
        });

        return sorted_color_ids;
    }

/* vector<unordered_map<string, vector<size_t>>> split_bitvector_and_count(const uint64_t* data, size_t total_bits, size_t n_colors) {
    if (n_colors == 0) throw runtime_error("n_colors must be > 0");
    if (total_bits % n_colors != 0) throw runtime_error("total_bits must be divisible by n_colors");

    size_t n_blocks = total_bits / n_colors;
    size_t total_words = (total_bits / 64);
    size_t words_per_block = (n_colors / 64) + 1;

    vector<unordered_map<string, vector<size_t>>> results(n_colors + 1);

    for (size_t start = 0; start < n_blocks; start++) {
        size_t count = 0;
        size_t base_bit = start * n_colors;
        size_t base_word = base_bit / 64;
        unsigned bit_offset = static_cast<unsigned>(base_bit % 64); // 0..63

        string block_bytes;
        block_bytes.reserve(words_per_block * sizeof(uint64_t));

        // Build each 64-bit chunk of the logical block by combining adjacent machine words when needed
        for (size_t k = 0; k < words_per_block; ++k) {
            size_t idx_lo = base_word + k;
            uint64_t lo = 0;
            uint64_t hi = 0;
            if (idx_lo < total_words) lo = data[idx_lo];
            if (idx_lo + 1 < total_words) hi = data[idx_lo + 1];

            uint64_t word;
            if (bit_offset == 0) {
                word = lo;
            } else {
                // shift/OR to assemble 64 aligned bits of the logical block
                word = (lo >> bit_offset) | (hi << (64 - bit_offset));
            }

            // mask the upper bits of the last chunk so identical logical blocks hash equal
            uint64_t bits_in_this_word = std::min<uint64_t>(64, n_colors - k * 64);
            uint64_t mask = (bits_in_this_word == 64) ? ~0ULL : ((1ULL << bits_in_this_word) - 1ULL);
            word &= mask;

            count += sdsl::bits::cnt(word);

            // append the full 8 bytes for consistent hashing (masked upper bits are zero)
            block_bytes.append(reinterpret_cast<const char*>(&word), sizeof(word));
        }

        // store occurrence
        results[count][block_bytes].push_back(start);
    }

    return results;
}   
 */

 // TODO remove?
 size_t count_ones(const sdsl::bit_vector& bv) {
        size_t count = 0;
        size_t n_bits = bv.size();
        size_t n_words = (n_bits + 63) / 64;

        const uint64_t* data = bv.data();

        // all full 64-bit words (except the last)
        for (size_t i = 0; i + 1 < n_words; ++i) {
            count += sdsl::bits::cnt(data[i]);
        }

        // last word
        if (n_bits % 64 != 0) {
            size_t last_bits = n_bits % 64;
            uint64_t mask = (uint64_t(1) << last_bits) - 1;
            uint64_t last_word = data[n_words - 1] & mask;
            count += sdsl::bits::cnt(last_word);
        } else if (n_words > 0) {
            // the last word can be full
            count += sdsl::bits::cnt(data[n_words - 1]);
        }

        return count;
    }

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
        cerr << "n_colors: "<< (int)n_colors << endl;

        subtract_and_filter(cf.lengths_by_freq, plen);

        // tlen frequency
        //Salmonella
        //vector<uint8_t> tlen_freq = {2, 3, 1, 4, 5,6,0,7,21,19,10,16,13,8,20,9,18,12,15,17,11,14};
        //Ecoli
        //vector<uint8_t> tlen_freq = {21, 19, 20, 16, 18, 17, 13, 15, 5, 4, 14, 6, 10, 12, 7, 11, 9, 8, 3, 2, 1};
        //TE
        //vector<uint8_t> tlen_freq = {4,3,5,19,21,16,20,6,13,18,7,17,10,15,14,2,12,11,9,8,1,0};



        // 1. Deduplicate color_sets
        cerr << "Deduplicate color_sets" << endl;

        // Create a vector of n_colors vectors to store in each the index at which a colorset id with that many colors appears
        
        //const uint64_t* data = color_sets_concat.data();
        //vector<unordered_map<string, vector<size_t>>> n_bits_set_to_1 = split_bitvector_and_count(reinterpret_cast<const uint64_t*>(color_sets_concat.data()), n_finimizers * n_colors, n_colors);
        
        auto Map_bits_set_to_1 = split_bitvector(cf.color_sets_concat, n_colors);

        auto Sorted_bits_set_to_1 = sort_bitvectors(Map_bits_set_to_1);

        // Store only unique
        size_t num_unique_blocks = Sorted_bits_set_to_1.size();
        sdsl::bit_vector unique_color_sets(num_unique_blocks * n_colors);

        // Store unique color sets ids per finimizer
        vector<uint32_t> color_set_ids(n_finimizers, 0); // One per finimizer: index into unique_color_sets

        // real offsets
        uint64_t new_offset = 0;

        for (auto& bv: Sorted_bits_set_to_1){
            vector<size_t>& old_offsets = Map_bits_set_to_1[bv];
            
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

        /* // Free up memory from now-unused vector
        Map_bits_set_to_1.clear();
        Map_bits_set_to_1.shrink_to_fit(); 
        Sorted_bits_set_to_1.clear();
        Sorted_bits_set_to_1.shrink_to_fit();*/

        // Assign to final structure
        this->color_sets_concat = std::move(unique_color_sets);

        cerr << "Deal with tails" << endl;
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

        for(int64_t i = 0; i < n_finimizers; i++) {
            
            if(cf.lengths[i] < prefix_len){
                std::string_view sprefix(cf.concat.data() + f_start, cf.lengths[i]);
                uint64_t sp_int = prefix2int(sprefix,0, cf.lengths[i]);
                sB[sp_int]= {cf.lengths[i],color_set_ids[i]}; // i= color_set_id

            } else {
                std::string_view prefix(cf.concat.data() + f_start, prefix_len);
                true_or_crash(f_start + prefix_len <= cf.concat.size(),
                        "ERROR: out-of-bounds prefix access");
                
                if(prefix != cur_prefix) {
                    // Bucket changes -> encode currently collected tails
                    buckets[p_int]=Bucket(cur_tails, cur_color_set_ids, cf.lengths_by_freq);
                    p_int = prefix2int(prefix, 0, prefix_len);
                    cur_tails.clear();
                    cur_color_set_ids.clear();
                }
                
                cur_tails.push_back(std::string_view(cf.concat.data() + f_start + prefix_len, cf.lengths[i] - prefix_len));
                cur_color_set_ids.push_back(color_set_ids[i]);
                cur_prefix = prefix;
            }
            f_start += cf.lengths[i];
        }

        if(cur_tails.size() > 0){ // Last bucket
            buckets[p_int]=Bucket(cur_tails, cur_color_set_ids, cf.lengths_by_freq);

        }
    }

    void search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans ) const{
  
        const int64_t query_len = query.length();
      
        if (query.size() < this->k) return; 

        // TODO: How to mark in Finimizers and r_Finimizers if a finimizer is not found????
        // now -1 if not found
        // option 2: store x + 64, (x=correct value) or 0 if not found 

        vector<int64_t> Finimizers;
        Finimizers.reserve(query_len - k +1);
        rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k, Finimizers);
        //cerr << Finimizers.size() << endl;
        
        // reverse complement
        //cerr << "reverse complement" << endl;
        vector<int64_t> r_Finimizers;
        r_Finimizers.reserve(query_len - k +1);
        string r_query = sbwt::get_rc(query);
        rarest_fmin_streaming_search(r_query, this->buckets, this->sB, this->plen, this->k, r_Finimizers);
        //cerr << r_Finimizers.size() << endl;

        // Combine the results of finimizers color ids for forward and reverse

        // Check the colors for every finimizer found
        //pseudoalignment_stats(Finimizers, this->color_sets_concat, this->n_colors, ans);// old
        combine_f_rc(Finimizers, r_Finimizers, this->color_sets_concat, this->n_colors, ans);

        return;
    }

    // TODO use min_value inside pseudoalignment stats ???
    uint16_t search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans, const float& t) const {
        
        const int64_t query_len = query.length();

        if (query.size() < this->k) return 0; 

        vector<int64_t> Finimizers;
        rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k, Finimizers);

        // reverse complement
        vector<int64_t> r_Finimizers;
        string r_query = sbwt::get_rc(query);
        rarest_fmin_streaming_search(r_query, this->buckets, this->sB, this->plen, this->k, r_Finimizers);

        // Combine the results of finimizers color ids for forward and reverse

        // Check the colors for every finimizer found
        //uint16_t min_value = pseudoalignment_stats(Finimizers, this->color_sets_concat, this->n_colors, ans,t);
        uint16_t min_value = combine_f_rc(Finimizers, r_Finimizers, this->color_sets_concat, this->n_colors, ans, t);

        return min_value;
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
