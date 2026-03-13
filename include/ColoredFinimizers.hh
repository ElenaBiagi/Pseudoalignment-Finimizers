#pragma once

#include <vector>
#include <optional>
#include <unordered_map>
#include <map>

#include "sdsl/bit_vectors.hpp"
#include "CompressedColorSets.hh"

#include "Fluke8.hh"
#include "PrefTab.hh"

#include "batch_querying.hh"

#include <chrono>

using namespace std;

static std::chrono::nanoseconds time_rarest_fmin(0);
static std::chrono::nanoseconds time_rarest_fmin_rc(0);
static std::chrono::nanoseconds time_combine(0);
static std::chrono::nanoseconds time_index_loading(0);
static std::chrono::nanoseconds time_output(0);

void print_search_timing_stats()
{
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
class ColoredFinimizers
{
public:
    vector<char> concat; // Finimizers concatenated in lexicographic order (ASCII characters).
    vector<uint8_t> lengths;
    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors.
    vector<uint8_t> lengths_by_freq;

    struct VecConcat {
        std::vector<size_t> concat;
        std::vector<size_t> ends;
    };
    
    VecConcat sparse_colors;  // Sparse color representation (when meta=true)

    VecConcat load_sparse_colors(std::istream& file) {
        uint8_t is_sparse;
        file.read(reinterpret_cast<char*>(&is_sparse), 1);
        
        if (is_sparse != 1) {
            throw std::runtime_error("Expected sparse format");
        }
        
        uint64_t n_elements, n_sets;
        file.read(reinterpret_cast<char*>(&n_elements), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&n_sets), sizeof(uint64_t));
        
        std::vector<size_t> concat(n_elements);
        file.read(reinterpret_cast<char*>(concat.data()), n_elements * sizeof(size_t));
        
        std::vector<size_t> ends(n_sets);
        file.read(reinterpret_cast<char*>(ends.data()), n_sets * sizeof(size_t));
        
        return {concat, ends};
    }
    // Loads from the format output by the Rust CLI command `finimizer_matrix` with option --reverse.
    // That format has colexicographically sorted reverse finimizers. We reverse them to
    // get lex-sorted finimizers.
    void load(std::istream &in, bool meta)
    {
        cerr << "Loading uncompressed tails" << endl;
        uint64_t n_finimizers;
        in.read(reinterpret_cast<char *>(&n_finimizers), sizeof(n_finimizers));

        uint64_t finimizer_total_length;
        in.read(reinterpret_cast<char *>(&finimizer_total_length), sizeof(finimizer_total_length));

        uint64_t n_colors;
        in.read(reinterpret_cast<char *>(&n_colors), sizeof(n_colors));

        lengths.resize(n_finimizers);
        in.read(reinterpret_cast<char *>(lengths.data()), n_finimizers * sizeof(uint8_t));

        concat.resize(finimizer_total_length);
        in.read(reinterpret_cast<char *>(concat.data()), finimizer_total_length * sizeof(char));


        uint8_t is_sparse;
        in.read(reinterpret_cast<char *>(&is_sparse), 1);
        if (is_sparse == 1) {
            cerr << "meta" << endl;
            sparse_colors = load_sparse_colors(in);
        }
        else{
            int64_t n_bits = n_finimizers * n_colors;
            // The bits are in u64 Lsb format
            vector<uint64_t> words((n_bits + 63) / 64); // Ceil div by 64
            in.read(reinterpret_cast<char *>(words.data()), words.size() * sizeof(uint64_t));
            color_sets_concat.resize(n_bits);
            for (int64_t i = 0; i < words.size(); i++)
            {
                color_sets_concat.set_int(i * 64, words[i]);
            }
        }
        

        uint64_t lengths_by_freq_size;
        in.read(reinterpret_cast<char *>(&lengths_by_freq_size), sizeof(lengths_by_freq_size));
        lengths_by_freq.resize(lengths_by_freq_size);
        in.read(reinterpret_cast<char *>(lengths_by_freq.data()), lengths_by_freq_size * sizeof(uint8_t));

        cerr << "Reversing finimizer strings" << endl;
        int64_t start_in_concat = 0;
        for (uint64_t f_idx = 0; f_idx < lengths.size(); f_idx++)
        {
            int64_t s = start_in_concat;
            int64_t e = start_in_concat + lengths[f_idx];
            std::reverse(concat.begin() + s, concat.begin() + e);
            start_in_concat = e;
        }
        cerr << "Uncompressed index loaded" << endl;
    }

    vector<pair<int64_t, uint64_t>> fmin_stats()
    {
        cerr << "Printing finimizers stats:" << endl;
        std::sort(lengths.begin(), lengths.end());
        vector<pair<int64_t, uint64_t>> lengths_and_freq;
        lengths_and_freq.reserve(lengths_by_freq.size());

        for (size_t i = 0; i < lengths.size();)
        {
            size_t j = i + 1;
            while (j < lengths.size() && lengths[j] == lengths[i])
                ++j;
            lengths_and_freq.emplace_back(lengths[i], j - i);
            i = j;
        }

        return lengths_and_freq;
    }
};

void true_or_crash(bool b, const char *error_message)
{
    if (!b)
    {
        cerr << error_message << endl;
        exit(1);
    }
}

//inline int16_t pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results, const float t);

//inline void pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int16_t> &results);

class CompressedColoredFinimizers
{

public:
    CompressedColorSets CCS; // L, EF, BV

    uint64_t n_colors;
    uint64_t n_finimizers;
    uint64_t plen;
    uint64_t short_long_t;
    uint64_t k;
    bool meta;

    Fluke8 f8;
    PrefTab pt;

    vector<uint64_t> color_set_ids;

    int get_k() const { return k; }

    CompressedColoredFinimizers() = default;

    CompressedColoredFinimizers(ColoredFinimizers &&cf, int64_t prefix_len, int64_t short_long_t, uint64_t kmer_size, bool meta)
    {
        cerr << "Let's compress it!" << endl;
        this->plen = prefix_len;
        cerr << "Prefix length: " << plen << endl;
        this->short_long_t = short_long_t;
        cerr << "Threshold : " << short_long_t << endl;
        this->k = kmer_size;
        cerr << "k-mer size: " << k << endl;
        this->meta = meta;

        n_finimizers = cf.lengths.size();
        cerr << "total finimizers: " << (int)n_finimizers << endl;
        true_or_crash(n_finimizers > 0, "ERROR: 0 finimizers");

        if (meta) {
            // Sparse deduplication
            cerr << "Deduplicate sparse color sets" << endl;
            n_colors = cf.sparse_colors.ends.size();
            
            map<vector<size_t>, vector<size_t>> color_set_to_finimizers;
            
            size_t start = 0;
            for (size_t i = 0; i < n_finimizers; ++i)
            {
                size_t end = cf.sparse_colors.ends[i];
                
                //colors for this finimizer
                vector<size_t> color_set(cf.sparse_colors.concat.begin() + start, cf.sparse_colors.concat.begin() + end);
                
                color_set_to_finimizers[color_set].push_back(i);
                start = end;
            }
            
            cerr << "Unique color sets: " << color_set_to_finimizers.size() << endl;
            
            this->color_set_ids.resize(n_finimizers);
            vector<pair<vector<size_t>, vector<size_t>>> deduplicated_cs_sparse;
            uint64_t color_set_id = 0;
            for (auto &[color_set, finimizer_indices] : color_set_to_finimizers)
            {
                for (size_t fm_idx : finimizer_indices)
                {
                    this->color_set_ids[fm_idx] = color_set_id;
                }
                color_set_id++;
                
                deduplicated_cs_sparse.emplace_back(color_set, std::move(finimizer_indices));

            }
            CompressedColorSets CCS(deduplicated_cs_sparse, n_colors, this->color_set_ids);
            this->CCS = std::move(CCS);
        } else {
            true_or_crash(cf.color_sets_concat.size() % n_finimizers == 0, "ERROR: color set bitmap length not divisible by finimizer count");
            n_colors = cf.color_sets_concat.size() / n_finimizers;
            cerr << "n_colors: " << (int)n_colors << endl;
            //cerr << sdsl::util::cnt_one_bits(cf.color_sets_concat) << endl;

            cerr << "Deduplicate color sets" << endl;
            const uint64_t *data = cf.color_sets_concat.data();
            sdsl::bit_vector bv(n_colors);
            
            map<sdsl::bit_vector, vector<size_t>> color_set_to_finimizers;
            
            for (size_t i = 0; i < n_finimizers; ++i)
            {
                read_colors_to_bv(data, n_colors, i, bv);
                color_set_to_finimizers[bv].push_back(i);
            }
            
            vector<pair<sdsl::bit_vector, vector<size_t>>> deduplicated_cs;
            for (auto &[color_set, finimizer_indices] : color_set_to_finimizers)
            {
                deduplicated_cs.emplace_back(color_set, std::move(finimizer_indices));
            }

            cerr << "Unique color sets: " << deduplicated_cs.size() << endl;

            this->color_set_ids.resize(n_finimizers);
            CompressedColorSets CCS(deduplicated_cs, n_colors, this->color_set_ids);
            this->CCS = std::move(CCS);
        }

        // this->color_set_ids is in the right order for finimizers
        vector<uint64_t> short_color_set_ids;
        vector<uint64_t> long_color_set_ids;

        vector<std::string_view> short_fmins;
        vector<std::string_view> long_fmins;

        // TODO get rid of short_color_set_ids by overwriting color_set_ids
        //NB: Finimizers are assumed to be in ascending lex order

        int64_t f_start = 0;
        for (int64_t i = 0; i < n_finimizers; i++)
        {
            if (cf.lengths[i] <= short_long_t)
            {
                // fluke8
                short_color_set_ids.push_back(this->color_set_ids[i]);
                short_fmins.push_back(std::string_view(cf.concat.data() + f_start, cf.lengths[i]));
            }
            else
            {
                // prefix tab
                long_color_set_ids.push_back(this->color_set_ids[i]);
                long_fmins.push_back(std::string_view(cf.concat.data() + f_start, cf.lengths[i]));
            }
            f_start += cf.lengths[i];
        }
        uint64_t n_short_fmin = short_color_set_ids.size();
        
        short_color_set_ids.insert(short_color_set_ids.end(), long_color_set_ids.begin(), long_color_set_ids.end());
        this->color_set_ids = short_color_set_ids;
        

        // build fluke
        this->f8 = std::move(Fluke8(short_fmins, 0, short_long_t));

        // build prefix table using plen
        this->pt = std::move(PrefTab(long_fmins, plen));

    }

    void read_colors_to_bv(const uint64_t *data, const uint64_t n_colors, const uint64_t start, sdsl::bit_vector &bv)
    {
        sdsl::util::set_to_value(bv, 0);
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
                bv[bit] = 1;
                word &= word - 1;
            }

            ++ptr;
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
                bv[color_id + bit] = 1;
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
                bv[color_id + bit] = 1;
                word &= word - 1;
            }
        }
    }


    // void search(const std::string &query, vector<int16_t> &results, vector<int64_t> &Finimizers) const
    // {
    //     const int64_t query_len = query.length();
    //     if (query_len < this->k)
    //         return;

    //     // Forward finimizer search
    //     // vector<int64_t> Finimizers;
    //     Finimizers.clear();
    //     Finimizers.reserve(query_len - k + 1);
    //     {
    //         auto start = std::chrono::high_resolution_clock::now();
    //         //rarest_fmin_streaming_search(query, this->non_empty_buckets, this->non_empty_bv, this->non_empty_bv_rs, this->sB, this->plen, this->k, Finimizers);
    //         auto end = std::chrono::high_resolution_clock::now();
    //         time_rarest_fmin += (end - start);
    //     }
    //     // Color sets
    //     {
    //         auto start = std::chrono::high_resolution_clock::now();
    //         // combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans, results);
    //         pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, results);
    //         auto end = std::chrono::high_resolution_clock::now();
    //         time_combine += (end - start);
    //     }
    // }

    // Threshold-based search: returns minimum value and fills ans
    // uint16_t search(const std::string &query, vector<int16_t> &results, const float t, vector<int64_t> &Finimizers) const
    // {

    //     const int64_t query_len = query.length();
    //     if (query_len < this->k)
    //         return 0;

    //     // Forward finimizer search
    //     // vector<int64_t> Finimizers;
    //     Finimizers.clear();
    //     Finimizers.reserve(query_len - k + 1);
    //     {
    //         auto start = std::chrono::high_resolution_clock::now();
    //         //rarest_fmin_streaming_search(query, this->non_empty_buckets, this->non_empty_bv, this->non_empty_bv_rs, this->sB, this->plen, this->k, Finimizers);
    //         auto end = std::chrono::high_resolution_clock::now();
    //         time_rarest_fmin += (end - start);
    //     }

    //     // Combine the results of finimizers color ids for forward and reverse
    //     int16_t T;
    //     {
    //         auto start = std::chrono::high_resolution_clock::now();
    //         // min_value = combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans, t, results);
    //         T = pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, results, t);

    //         auto end = std::chrono::high_resolution_clock::now();
    //         time_combine += (end - start);
    //     }
    //     return T;
    // }

    void search_batch(const vector<std::string> &reads, vector<int16_t> &results, vector<int64_t> &Finimizers, const uint64_t batch_size, const uint64_t k, const float &t) const
    {   
        // It's not possible to reuse the same vector for every query as we now have a batch
        //Finimizers.clear();
        //Finimizers.reserve(query_len - k + 1);
        {
            auto start = std::chrono::high_resolution_clock::now();
            // what should this output?
            batch_querying(reads, this->f8, this->pt, batch_size, this->CCS, this->n_colors, this->color_set_ids, k, t);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }

        
        // {
        //     auto start = std::chrono::high_resolution_clock::now();
        //     // combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans, results);
        //     // pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, results);
        //     auto end = std::chrono::high_resolution_clock::now();
        //     time_combine += (end - start);
        // }
    }

    void search_batch(const vector<std::string> &reads, vector<int16_t> &results, vector<int64_t> &Finimizers, const uint64_t batch_size, const uint64_t k) const
    {   
        // It's not possible to reuse the same vector for every query as we now have a batch
        //Finimizers.clear();
        //Finimizers.reserve(query_len - k + 1);
        {
            auto start = std::chrono::high_resolution_clock::now();
            // what should this output?
            batch_querying(reads, this->f8, this->pt, batch_size, this->CCS, this->n_colors, this->color_set_ids, k);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }
    }

    void serialize(const std::string &index_prefix) const
    {
        string filename = index_prefix + ".fmin";
        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            cerr << "Error: Could not open file for writing: " << filename << endl;
            return;
        }

        cerr << "Save the index to " << filename << endl;

        CCS.serialize(out);

        uint64_t color_set_ids_size = color_set_ids.size();
        out.write(reinterpret_cast<const char *>(&color_set_ids_size), sizeof(uint64_t));
        out.write(reinterpret_cast<const char *>(color_set_ids.data()), color_set_ids_size * sizeof(uint64_t));

        f8.serialize(out);
        pt.serialize(out);

        // metadata
        out.write(reinterpret_cast<const char *>(&n_colors), sizeof(n_colors));
        out.write(reinterpret_cast<const char *>(&n_finimizers), sizeof(n_finimizers)); // TODO do we need this?
        out.write(reinterpret_cast<const char *>(&k), sizeof(k));

        out.close();
        cerr << "DONE" << endl;
    }
    void load(const std::string &index_prefix)
    {
        string filename = index_prefix + ".fmin";
        std::ifstream in(filename, std::ios::binary);
        if (!in)
        {
            cerr << "Error: Could not open file for reading: " << filename << endl;
            return;
        }

        cerr << "Loading index from " << filename << endl;

        CCS.load(in);

        uint64_t color_set_ids_size;
        in.read(reinterpret_cast<char *>(&color_set_ids_size), sizeof(uint64_t));
        color_set_ids.resize(color_set_ids_size);
        in.read(reinterpret_cast<char *>(color_set_ids.data()), color_set_ids_size * sizeof(uint64_t));

        f8.load(in);
        pt.load(in);

        // metadata
        in.read(reinterpret_cast<char *>(&n_colors), sizeof(n_colors));
        in.read(reinterpret_cast<char *>(&n_finimizers), sizeof(n_finimizers)); // TODO do we need this?
        //in.read(reinterpret_cast<char *>(&plen), sizeof(plen));
        in.read(reinterpret_cast<char *>(&k), sizeof(k));

        in.close();
        cerr << "DONE" << endl;
    }
};
