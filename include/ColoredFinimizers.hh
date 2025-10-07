#pragma once

#include <vector>
#include <optional>
#include <unordered_map>

#include "sdsl/bit_vectors.hpp"
#include "rarest_fmin_search.hh"
#include "CompressedColorSets.hh"
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

inline uint16_t pseudoalignment_stats( vector<int64_t>& Fmin, const CompressedColorSets& CCS, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans, const float t );

inline void pseudoalignment_stats( vector<int64_t>& Fmin, const CompressedColorSets& CCS, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans);


class CompressedColoredFinimizers {

private:
    sdsl::bit_vector unique_color_sets;

    vector<uint64_t> color_set_ids;

public:
    CompressedColorSets CCS; // L, EF, BV

    vector<optional<Bucket>> buckets; 
    unordered_map<uint32_t, pair<uint8_t, int64_t>> sB; // Create a hash table to store the finimizers shorter than the prefix length
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
        cerr << sdsl::util::cnt_one_bits(cf.color_sets_concat) << endl;

        cerr << "Deduplicate color sets" << endl;
        const uint64_t* data = cf.color_sets_concat.data();
        sdsl::bit_vector bv(n_colors);

        unordered_map<sdsl::bit_vector, vector<size_t>, BVHash, BVEqual> deduplicated_cs;

        for (size_t i = 0; i < n_finimizers; ++i) {
            read_colors_to_bv(data, n_colors, i, bv);

            auto it = deduplicated_cs.find(bv);
            if (it == deduplicated_cs.end()) {
                // First time seeing this color set
                deduplicated_cs.emplace(bv, vector<size_t>{i});
            } else {
                // Already seen
                it->second.push_back(i);
            }
        }

        cerr << "Unique color sets: " << deduplicated_cs.size() << endl;

        // Flatten into a single sdsl::bit_vector (unique_color_sets)
        unique_color_sets = sdsl::bit_vector(deduplicated_cs.size() * n_colors);
        uint64_t new_offset = 0;

        for (auto& kv : deduplicated_cs) {
            const sdsl::bit_vector& ucs = kv.first;
            for (size_t j = 0; j < n_colors; ++j) {
                unique_color_sets[new_offset + j] = ucs[j];
            }
            new_offset += n_colors;
        }


        this->color_set_ids.resize(n_finimizers); // ids sorted based on the frequency of fmins length
        CompressedColorSets CCS(deduplicated_cs, n_colors, this->color_set_ids);

        
        // Assign to final structure
        this->CCS = std::move(CCS);
        

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
                sB[sp_int]= {cf.lengths[i],this->color_set_ids[i]};

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
                cur_color_set_ids.push_back(this->color_set_ids[i]);
                cur_prefix = prefix;
            }
            f_start += cf.lengths[i];
        }

        if (!cur_tails.empty()) { // Last bucket
            buckets[p_int] = Bucket(cur_tails, cur_color_set_ids, cf.lengths_by_freq);
        }

    }

    void read_colors_to_bv(const uint64_t* data, const uint64_t n_colors, const uint64_t start, sdsl::bit_vector& bv){
        sdsl::util::set_to_value(bv, 0);
        const uint64_t* ptr = data + (start * n_colors) / 64;
        uint64_t bit_offset = (start * n_colors) % 64;

        uint64_t color_id = 0;
        uint64_t bits_left = n_colors;

        if (bit_offset != 0){
            // 1. Read the first word
            uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
            uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
            uint64_t word = (*ptr >> bit_offset) & mask;

            while (word != 0) {
                uint64_t bit = __builtin_ctzll(word);
                bv[bit]=1;
                word &= word - 1;
            }

            ++ptr;
            color_id += bits_to_read;
            bits_left -= bits_to_read;
        }

        // 2. Read aligned words in btw
        while (bits_left >= 64) {
            uint64_t word = *ptr++;
            for (uint64_t w = word; w != 0;) {
                uint64_t bit = __builtin_ctzll(w);
                bv[color_id + bit]=1;
                w &= w - 1;
            }
            color_id += 64;
            bits_left -= 64;
        }

        // 3. Read the last word (if any)
        if (bits_left > 0) {
            uint64_t mask = ((1ULL << bits_left) - 1);
            uint64_t word = *ptr & mask;

            while (word != 0) {
                uint64_t bit = __builtin_ctzll(word);
                bv[color_id + bit]=1;
                word &= word - 1;
            }
        }
    }

    

    void search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans, vector<int64_t>& Finimizers) const {
        const int64_t query_len = query.length();
        if (query_len < this->k) return;

        // Forward finimizer search
        //vector<int64_t> Finimizers;
        Finimizers.clear();
        Finimizers.reserve(query_len - k + 1);
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k, Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }
        // Color sets
        {
            auto start = std::chrono::high_resolution_clock::now();
            //combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans, results);
            pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, ans);
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
    }

    // Threshold-based search: returns minimum value and fills ans
    uint16_t search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans, const float t, vector<int64_t>& Finimizers) const {

        const int64_t query_len = query.length();
        if (query_len < this->k) return 0;

        // Forward finimizer search
        //vector<int64_t> Finimizers;
        Finimizers.clear();
        Finimizers.reserve(query_len - k + 1);
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(query, this->buckets, this->sB, this->plen, this->k, Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }

        // Combine the results of finimizers color ids for forward and reverse
        uint16_t min_value;
        {
            auto start = std::chrono::high_resolution_clock::now();
            //min_value = combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans, t, results);
            min_value = pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, ans, t);

            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
        return min_value;
    }


    void serialize(const string& index_prefix) const {
        cerr << "Save the index" << endl;

        CCS.serialize(index_prefix);

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

        CCS.load(index_prefix);

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
            std::pair<uint8_t,int64_t> val;
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

};

inline void process_word(uint64_t word, const uint64_t base, vector<uint64_t>& results, const uint64_t freq) {
    //cerr << "process_word" << endl;
    while (word) {
        //uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        results[base + bit] += freq;
        word &= word - 1; // clear lowest bit
    }
    //cerr << "end" << endl;

}

inline void read_bv(const uint64_t* data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<uint64_t>& results){
    //cerr << "read_bv" << endl;

    const uint64_t* ptr = data + (start * n_colors) / 64;
    uint64_t bit_offset = (start * n_colors) % 64;

    uint64_t color_id = 0;
    uint64_t bits_left = n_colors;


    if (bit_offset != 0){
        // 1. Read the first word
        uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
        uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
        uint64_t word = (*ptr >> bit_offset) & mask;

        /* while (word != 0) {
            uint64_t bit = __builtin_ctzll(word);
            results[bit]+=freq;
            word &= word - 1;
        } */
        process_word(word,color_id, results, freq);

        ptr++;
        color_id += bits_to_read;
        bits_left -= bits_to_read;
    }

    // 2. Read aligned words in btw
    while (bits_left >= 64) {
        uint64_t word = *ptr++;
        process_word(word, color_id, results, freq);
        /* for (uint64_t w = word; w != 0;) {
            
            uint64_t bit = __builtin_ctzll(w);
            results[color_id + bit]+=freq;
            w &= w - 1;
        } */
        color_id += 64;
        bits_left -= 64;
    }

    // 3. Read the last word (if any)
    if (bits_left > 0) {
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;

        /* while (word != 0) {
            uint64_t bit = __builtin_ctzll(word);
            results[color_id + bit]+=freq;
            word &= word - 1;
        } */
        process_word(word,color_id, results, freq);

    }
    //cerr << "end" << endl;
}

void read_colors(const CompressedColorSets& CCS, const uint64_t n_colors, vector<uint64_t>& results, const vector<pair<int64_t, uint64_t>>& fmin_v){
    const sdsl::bit_vector& BV = CCS.getBV();
    const vector<uint16_t>& L = CCS.getL();
    const sdsl::enc_vector<>& EF = CCS.getEF();

    /* cerr << "BV: "<< (BV.size()-63)/n_colors << endl;
    cerr << "EF: " << EF.size() << endl;
    cerr << "L: " << L.size() << endl; */

    const uint64_t* data = BV.data();
    const uint64_t L_size = EF.size();
    // Exploit the fact that the pos are sorted
    uint64_t i;
    for ( i=0; i< fmin_v.size(); i++) {
        const auto& [pos,freq] = fmin_v[i];
        if (pos < L_size){  // sparse -> read from L
            //if (pos >= EF.size()){ cerr << "pos = "<< pos << "== EF.size() = "<< EF.size() << endl;}
            //cerr << pos << ", " << EF.size() << endl;
            const size_t end = EF[pos]; // exclusive end
            size_t start = EF[pos-1]; // inclusive start
            while(start < end){ results[L[start++]]+=freq;}
        } else { break;}
    }
    for (auto j=i; j< fmin_v.size(); j++) {
        const auto& [pos,freq] = fmin_v[j];
        uint64_t start = pos - L_size;
        read_bv(data, start, freq, n_colors, results);
    }
}


inline void pseudoalignment_stats(vector<int64_t>& Fmin, const CompressedColorSets& CCS, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
    vector<uint64_t> results(n_colors,0);
    /* if (results.size() != n_colors) {
        results.assign(n_colors, 0);  
    } else {
        std::fill(results.begin(), results.end(), 0);
    }
 */
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

    read_colors(CCS, n_colors, results, fmin_v);

    const size_t found_fmin = Fmin.size(); // # total finimizers
    counting_sort(results, ans, found_fmin, n_colors);
    return;
}

inline uint16_t pseudoalignment_stats(vector<int64_t>& Fmin, const CompressedColorSets& CCS, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans,  const float t ){ //vector<uint64_t>& results,
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

    read_colors(CCS, n_colors, results, fmin_v);

    // Check the values above the minimum
    const size_t found_fmin = Fmin.size(); // # total finimizers

    counting_sort(results, ans, found_fmin, n_colors);
    
    // TODO use min value here?
    const uint64_t min_value = found_fmin * t;

    return min_value;
}
