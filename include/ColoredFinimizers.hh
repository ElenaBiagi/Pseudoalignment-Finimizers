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

void combine_f_rc(const std::vector<int64_t>& Fmin,
                  const std::vector<int64_t>& r_Fmin,
                  const CompressedColorSets& CCS,
                  const uint64_t n_colors,
                  std::vector<std::pair<uint16_t, uint16_t>>& ans);

uint64_t combine_f_rc(const std::vector<int64_t>& Fmin,
                      const std::vector<int64_t>& r_Fmin,
                      const CompressedColorSets& CCS,
                      const uint64_t n_colors,
                      std::vector<std::pair<uint16_t, uint16_t>>& ans,
                      const float& t);

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
        unordered_map<string, vector<size_t>> deduplicated_cs; // {cs:[fmin indices]}
        const uint64_t* data = cf.color_sets_concat.data();
        for (size_t i = 0; i < n_finimizers; i++) {
            sdsl::bit_vector bv = read_colors_to_bv(data, n_colors, i);
            string key((char*)bv.data(), ((n_colors + 63) / 64) * 8);
            
            /* // slightly slower
            std::ostringstream oss;
            sdsl::serialize(bv, oss);
            std::string key = oss.str(); */

            /* std::string key;
            key.resize((n_colors + 7) / 8, '\0');
            for (uint64_t j = 0; j < n_colors; j++) {
                if (bv[j]) key[j >> 3] |= (1 << (j & 7));
            } */

            deduplicated_cs[key].push_back(i);
        }        
        cerr << deduplicated_cs.size() << endl;

        this->color_set_ids.resize(n_finimizers); // ids sorted based on the frequency of fmins length
        CompressedColorSets CCS(deduplicated_cs, n_colors, this->color_set_ids, cf.color_sets_concat);

        
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

    sdsl::bit_vector read_colors_to_bv(const uint64_t* data, const uint64_t n_colors, const uint64_t start){
        sdsl::bit_vector bv(n_colors);

    
        const uint64_t* ptr = data + (start * n_colors) / 64;
        uint64_t bit_offset = (start * n_colors) % 64;

        uint64_t color_id = 0;

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

        // 2. Read aligned words in btw
        while (color_id + 64 <= n_colors) {
            uint64_t word = *ptr++;
            for (uint64_t w = word; w != 0;) {
                uint64_t bit = __builtin_ctzll(w);
                bv[color_id + bit]=1;
                w &= w - 1;
            }
            color_id += 64;
        }

        // 3. Read the last word (if any)
        uint64_t bits_left = n_colors - color_id;
        if (bits_left > 0) {
            uint64_t mask = ((1ULL << bits_left) - 1);
            uint64_t word = *ptr & mask;

            while (word != 0) {
                uint64_t bit = __builtin_ctzll(word);
                bv[color_id + bit]=1;
                word &= word - 1;
            }
        }
        return bv;
    }

    

    void search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans, vector<int64_t>& Finimizers, vector<int64_t>& r_Finimizers) const {
        //cerr << "search"<< endl;
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

        // Reverse complement search
        //vector<int64_t> r_Finimizers;
        r_Finimizers.clear();
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
            combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans);
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
        //cerr << "search end" << endl;
    }

    // Threshold-based search: returns minimum value and fills ans
    uint16_t search(const std::string& query, vector<pair<uint16_t, uint16_t>>& ans, const float& t, vector<int64_t>& Finimizers, vector<int64_t>& r_Finimizers) const {
        //cerr << "search"<< endl;

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

        // reverse complement
        //vector<int64_t> r_Finimizers;
        r_Finimizers.clear();
        r_Finimizers.reserve(query_len - k + 1);
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
            min_value = combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans, t);
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
        //cerr << "search end" << endl;

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

// TODO MOVE THIS TO COMPRESSED_COLOR_SETS ?

inline void read_bv(const uint64_t* data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<uint64_t>& results){
    //cerr << "read_bv" << endl;

    const uint64_t* ptr = data + (start * n_colors) / 64;
    uint64_t bit_offset = (start * n_colors) % 64;

    uint64_t color_id = 0;

    // 1. Read the first word
    uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
    uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
    uint64_t word = (*ptr >> bit_offset) & mask;

    while (word != 0) {
        uint64_t bit = __builtin_ctzll(word);
        results[bit]+=freq;
        word &= word - 1;
    }

    ptr++;
    color_id += bits_to_read;

    // 2. Read aligned words in btw
    while (color_id + 64 <= n_colors) {
        uint64_t word = *ptr++;
        for (uint64_t w = word; w != 0;) {
            uint64_t bit = __builtin_ctzll(w);
            results[color_id + bit]+=freq;
            w &= w - 1;
        }
        color_id += 64;
    }

    // 3. Read the last word (if any)
    uint64_t bits_left = n_colors - color_id;
    if (bits_left > 0) {
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;

        while (word != 0) {
            uint64_t bit = __builtin_ctzll(word);
            results[color_id + bit]+=freq;
            word &= word - 1;
        }
    }
    //cerr << "end" << endl;

}

void read_colors(const CompressedColorSets& CCS, const uint64_t n_colors, vector<uint64_t>& results, const vector<pair<int64_t, uint64_t>>& fmin_v){
    //cerr << "read_colors" << endl;

    const sdsl::bit_vector& BV = CCS.getBV();
    const vector<uint32_t>& L = CCS.getL();
    const vector<size_t>& EF = CCS.getEF();

    //cerr << "BV: "<< (BV.size()-63)/n_colors << endl;
    //cerr << "L: " << EF.size() << endl;

    const uint64_t* data = BV.data();
    const uint64_t L_size = EF.size();
    /* for (const auto& [pos,freq] : fmin_v) {
        if (pos >= L_size){  // Check how big start is (dense or sparse)

            uint64_t start = pos - L_size;
            read_bv(data, start, freq, n_colors, results);
        }
        else{ // read from L
            const size_t end = EF[pos]; // exclusive end
            size_t start = EF[pos-1]; // inclusive start
            while(start < end){ results[L[start++]]+=freq;}
        }
    } */
    // TODO exploit the fact that the pos are sorted
    uint64_t i;
    for ( i=0; i< fmin_v.size(); i++) {
        const auto& [pos,freq] = fmin_v[i];
        if (pos < L_size){  // sparse -> read from L
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

    //cerr << "end" << endl;

}

void combine_bv(const uint64_t* data, const uint64_t start1, const uint64_t start2, const uint64_t n_colors, vector<uint64_t>& results, const uint64_t freq) {
            //cerr << "combine_bv" << endl;

    const uint64_t bit_offset1 = start1 * n_colors;
    const uint64_t bit_offset2 = start2 * n_colors;

    uint64_t idx = 0;

    while (idx < n_colors) {
        
        uint64_t word1_idx = (bit_offset1 + idx) / 64;
        uint64_t word2_idx = (bit_offset2 + idx) / 64;

        uint64_t shift1 = (bit_offset1 + idx) % 64;
        uint64_t shift2 = (bit_offset2 + idx) % 64;

        uint64_t bits_to_process = std::min({n_colors - idx, 64 - shift1, 64 - shift2});

        uint64_t mask = (bits_to_process == 64) ? ~0ULL : ((1ULL << bits_to_process) - 1);

        uint64_t word1 = (data[word1_idx] >> shift1) & mask;
        uint64_t word2 = (data[word2_idx] >> shift2) & mask;
        uint64_t combined = word1 | word2;

        while (combined) {
            uint64_t bit = __builtin_ctzll(combined);
            uint64_t result_idx = idx + bit;
            //if (result_idx < n_colors) {
                results[result_idx] += freq;
            //} 
            combined &= combined - 1;
        }

        idx += bits_to_process;
    }
            //cerr << "end" << endl;

}

inline void process_word(uint64_t word, uint64_t base, const vector<uint32_t>& L, size_t& j, size_t end, vector<uint64_t>& results, uint64_t freq) {
    //cerr << "process_word" << endl;
    while (word) {
        uint64_t bit = __builtin_ctzll(word);
        uint64_t color = base + bit;

        // add smaller values in L
        while (j < end && L[j] < color) {
            results[L[j++]] += freq;
        }

        // skip L, write BV 
        if (j < end && L[j] == color) {
            j++;
        }

        // add BV values 
        results[color] += freq;
        word &= word - 1; // clear lowest bit
    }
    //cerr << "end" << endl;

}

void combine_bv_list(const uint64_t* data, const vector<uint32_t>& L, const vector<size_t>& EF, const uint64_t start1, const uint64_t pos2, const uint64_t n_colors, vector<uint64_t>& results, const uint64_t freq){
    
            //cerr << "combine_bv_lists" << endl;
// start2 is the list pos
    const size_t end2 = EF[pos2]; // exclusive end )
    size_t start2 = EF[pos2-1]; // inclusive start [

    // start1 is the bv pos
    const uint64_t* ptr = data + (start1 * n_colors) / 64;
    uint64_t bit_offset = (start1 * n_colors) % 64;

    uint64_t color_id = 0;

    // 1. Read the first word
    uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
    uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
    uint64_t word = (*ptr >> bit_offset) & mask;
    process_word(word, color_id, L, start2, end2, results, freq);
    ptr++;
    color_id += bits_to_read;

    // 2. Read aligned words in btw
    while (color_id + 64 <= n_colors) {
        uint64_t word = *ptr++;
        process_word(word, color_id, L, start2, end2, results, freq);
        color_id += 64;
    }

    // 3. Read the last word (if any)  
    if (color_id < n_colors) {
        uint64_t bits_left = n_colors - color_id;
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;
                process_word(word, color_id, L, start2, end2, results, freq);                
    }

    // L values left (if any)
    while (start2 < end2) {results[L[start2++]] += freq;}
        //cerr << "end" << endl;

}

inline void combine_lists(const vector<uint32_t>& L, const vector<size_t>& EF, const uint64_t pos1, const uint64_t pos2, vector<uint64_t>& results, const uint64_t freq){
        //cerr << "combine_lists" << endl;

    // AND
    const size_t end1 = EF[pos1]; // exclusive end
    const size_t start1 = EF[pos1-1]; // inclusive start

    const size_t end2 = EF[pos2]; // exclusive end
    const size_t start2 = EF[pos2-1]; 
            
    size_t i = start1, j = start2;
    while (i < end1 && j < end2) {
        if (L[i] < L[j]) {
            results[L[i++]] += freq;
        } else if (L[j] < L[i]) {
            results[L[j++]] += freq;
        } else {
            results[L[i]] += freq;
            i++; 
            j++;
        }
    }
    while (i < end1) results[L[i++]] += freq;
    while (j < end2) results[L[j++]] += freq;
        //cerr << "end" << endl;

}

void read_f_rc_colors(const CompressedColorSets& CCS, const uint64_t n_colors, vector<uint64_t>& results, const vector<pair<pair<uint64_t, uint64_t>, uint64_t>>& p_fmin_v){
    
    //cerr << "read_f_rc_colors" << endl;
    // option 1: compare f and r as you go
    // option 2: store f and r in 2 bitvectors and compare at the end
    //for (const auto& [[start, r_start], freq] : p_fmin_v) {
    const sdsl::bit_vector& BV = CCS.getBV();
    const vector<uint32_t>& L = CCS.getL();
    const vector<size_t>& EF = CCS.getEF();

    //cerr << "BV: "<< (BV.size()-63)/n_colors << endl;
    //cerr << "L: " << EF.size() << endl;

    const uint64_t* data = BV.data();
    const uint64_t L_size = EF.size();
    for (const auto& [key,freq] : p_fmin_v) {
        uint64_t start = key.first;
        uint64_t r_start = key.second;
        if (start >= L_size){
            start -= L_size;
            if (r_start >= L_size){
                r_start -= L_size;
                combine_bv(data, start, r_start, n_colors, results, freq);
                //uint64_t max = *std::max_element(results.begin(), results.end());
                //cerr << "Max value after combine_bv: " << max << std::endl;

            }
            else {
                combine_bv_list(data, L, EF, start, r_start, n_colors, results, freq);
                //uint64_t max = *std::max_element(results.begin(), results.end());
                //cerr << "Max value after combine_bv_list: " << max << std::endl;

            }
        } 
        else {
            if (r_start >= L_size){
                r_start -= L_size;
                combine_bv_list(data, L, EF, r_start, start, n_colors, results, freq);
                //uint64_t max = *std::max_element(results.begin(), results.end());
                //cerr << "Max value after combine_bv_list: " << max << std::endl;

            }
            else{
                combine_lists(L,EF, start, r_start, results, freq);
                //uint64_t max = *std::max_element(results.begin(), results.end());
                //cerr << "Max value after combine_list: " << max << std::endl;

            }
        }
    }
        //cerr << "end read_f_rc_colors" << endl;

}

// two overlpaiing k-mers are likely to have the same fmin so they are likely to share the same fmin on both strands
// This should anyways keep the number of false pos low

void combine_f_rc(const vector<int64_t>& Fmin, const vector<int64_t>& r_Fmin, const CompressedColorSets& CCS, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
    //cerr << "combine_f_rc" << endl;

    // NEW pseudoaligment_stats
    vector<uint64_t> results;
    results.resize(n_colors, 0);

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
    //CCS.init(found_fmin, n_colors); // only here we know found_fmin

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
    // TODO count sort max= max(EF)
    std::sort(i_fmin_v.begin(), i_fmin_v.end()); 

    read_colors(CCS, n_colors, results, i_fmin_v);

    //uint64_t max = *std::max_element(results.begin(), results.end());
    //cerr << "Max value after read_colors: " << max << std::endl;

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

    read_f_rc_colors(CCS, n_colors, results, p_fmin_v);
    //cerr << "results: " << results.size() << endl;
    //max = *std::max_element(results.begin(), results.end());
    //cerr << "Max value: " << max << std::endl;
    //cerr << "found_fmin: " << found_fmin << endl;
    //cerr << "ans: " << ans.size() << endl;
    //cerr << "n_colors: " << n_colors << endl;



    counting_sort(results, ans, found_fmin, n_colors);
    //cerr << "end combine_f_rc" << endl;

    return;
}

uint64_t combine_f_rc(const vector<int64_t>& Fmin, const vector<int64_t>& r_Fmin, const CompressedColorSets& CCS, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans, const float& t){
    //cerr << "combine_f_rc" << endl;
    // NEW pseudoaligment_stats
    vector<uint64_t> results;

    results.resize(n_colors, 0);

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
    //CCS.init(found_fmin, n_colors); // only here we know found_fmin


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
    read_colors(CCS, n_colors, results, i_fmin_v);


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
    read_f_rc_colors(CCS, n_colors, results, p_fmin_v);

    counting_sort(results, ans, found_fmin, n_colors);

    const uint64_t min_value = found_fmin * t;
    //cerr << "end" << endl;

    return min_value;
}
