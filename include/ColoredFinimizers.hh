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

    // Loads from the format output by the Rust CLI command `finimizer_matrix` with option --reverse.
    // That format has colexicographically sorted reverse finimizers. We reverse them to
    // get lex-sorted finimizers.
    void load(std::istream &in)
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

        int64_t n_bits = n_finimizers * n_colors;
        // The bits are in u64 Lsb format
        vector<uint64_t> words((n_bits + 63) / 64); // Ceil div by 64
        in.read(reinterpret_cast<char *>(words.data()), words.size() * sizeof(uint64_t));
        color_sets_concat.resize(n_bits);
        for (int64_t i = 0; i < words.size(); i++)
        {
            color_sets_concat.set_int(i * 64, words[i]);
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

    vector<string> print_fmins()
    {   
        vector<string> fmins_vector;
        fmins_vector.reserve(lengths.size());
        size_t f_start = 0;
        for (size_t l : lengths){
            fmins_vector.emplace_back(concat.data() + f_start, l);
            f_start += l;
        }
        return fmins_vector; 
    }

    void save_fmins(const std::string &index_prefix)
    {
        std::string filename = index_prefix + ".fmins_only";
        std::ofstream out(filename, std::ios::binary);
        if (!out) {
            std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
            return;
        }

        std::cerr << "Saving fmins to " << filename << std::endl;

        // write size (optional but recommended)
        uint64_t size = concat.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        // write raw data
        out.write(concat.data(), concat.size());

        out.close();
    }

    void print_colors(const std::string &index_prefix)
    {
        std::string filename = index_prefix;
        std::ofstream out(filename);
        if (!out) {
            std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
            return;
        }

        std::cerr << "Saving colors to " << filename << std::endl;

        constexpr size_t NUM_COLORS = 3682;
        size_t total_bits = color_sets_concat.size();
        size_t num_fmins = total_bits / NUM_COLORS;

        for (size_t i = 0; i < num_fmins; ++i) {
            size_t offset = i * NUM_COLORS;
            for (size_t c = 0; c < NUM_COLORS; ++c) {
                out << color_sets_concat[offset + c];
            }
            out << '\n';
        }

        out.close();
    }

    void print_colors_as_list(const std::string &index_prefix)
    {
        std::string filename = index_prefix;
        std::ofstream out(filename);
        if (!out) {
            std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
            return;
        }

        std::cerr << "Saving colors to " << filename << std::endl;

        constexpr size_t NUM_COLORS = 3682;
        size_t total_bits = color_sets_concat.size();
        size_t num_fmins = total_bits / NUM_COLORS;

        for (size_t i = 0; i < num_fmins; ++i) {
            size_t offset = i * NUM_COLORS;
            for (size_t c = 0; c < NUM_COLORS-1; ++c) {
                if (color_sets_concat[offset + c] == 1){out << c <<", ";}
            }
            if (color_sets_concat[offset + NUM_COLORS-1] == 1){out << NUM_COLORS-1;}
            out << '\n';
        }

        out.close();
    }

    void save_colors(const std::string &index_prefix)
    {
        std::string filename = index_prefix + ".colors_only";
        std::ofstream out(filename, std::ios::binary);
        if (!out) {
            std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
            return;
        }
        sdsl::serialize(color_sets_concat, out);
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

inline int64_t pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, const float t);

inline void pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results);

inline void count_bases(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen);

inline void count_single_base(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen);

inline int64_t count_single_base(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen, const float t);


class CompressedColoredFinimizers
{

private:
    sdsl::bit_vector unique_color_sets;

    vector<uint64_t> color_set_ids;

    std::unordered_map<int, int> tlen_rank_map;

    void set_tlen_order(const std::vector<uint8_t> &ordered_tlens)
    {
        for (size_t i = 0; i < ordered_tlens.size(); ++i)
        {
            tlen_rank_map[ordered_tlens[i]] = static_cast<int>(i);
        }
    }

public:
    CompressedColorSets CCS; // L, EF, BV

    vector<Bucket> non_empty_buckets;
    sdsl::int_vector<1> non_empty_bv;
    sdsl::rank_support_v5<1> non_empty_bv_rs;           // try _v only ?
    unordered_map<uint32_t, pair<uint8_t, int64_t>> sB; // Create a hash table to store the finimizers shorter than the prefix length
    uint64_t n_colors;
    uint64_t n_finimizers;
    uint64_t plen;
    uint64_t k;

    int get_k() const { return k; }

    CompressedColoredFinimizers() = default;

    CompressedColoredFinimizers(ColoredFinimizers &&cf, int64_t prefix_len, uint64_t kmer_size)
    {
        cerr << "Let's compress it!" << endl;
        this->plen = prefix_len;
        cerr << "prefix length: " << plen << endl;
        this->k = kmer_size;
        cerr << "k-mer size: " << k << endl;

        uint64_t n_buckets = (1ULL << (prefix_len * 2));
        cerr << "total buckets: " << (int)n_buckets << endl;
        this->non_empty_bv = sdsl::int_vector<1>(n_buckets, 0);

        n_finimizers = cf.lengths.size();
        cerr << "total finimizers: " << (int)n_finimizers << endl;
        true_or_crash(n_finimizers > 0, "ERROR: 0 finimizers");

        true_or_crash(cf.color_sets_concat.size() % n_finimizers == 0, "ERROR: color set bitmap length not divisible by finimizer count");
        n_colors = cf.color_sets_concat.size() / n_finimizers;
        cerr << "n_colors: " << (int)n_colors << endl;
        cerr << sdsl::util::cnt_one_bits(cf.color_sets_concat) << endl;

        cerr << "Deduplicate color sets" << endl;
        const uint64_t *data = cf.color_sets_concat.data();
        sdsl::bit_vector bv(n_colors);

        unordered_map<sdsl::bit_vector, vector<size_t>, BVHash, BVEqual> deduplicated_cs;

        for (size_t i = 0; i < n_finimizers; ++i)
        {
            read_colors_to_bv(data, n_colors, i, bv);

            auto it = deduplicated_cs.find(bv);
            if (it == deduplicated_cs.end())
            {
                // First time seeing this color set
                deduplicated_cs.emplace(bv, vector<size_t>{i});
            }
            else
            {
                // Already seen
                it->second.push_back(i);
            }
        }

        cerr << "Unique color sets: " << deduplicated_cs.size() << endl;

        // Flatten into a single sdsl::bit_vector (unique_color_sets)
        unique_color_sets = sdsl::bit_vector(deduplicated_cs.size() * n_colors);
        uint64_t new_offset = 0;

        for (auto &kv : deduplicated_cs)
        {
            const sdsl::bit_vector &ucs = kv.first;
            for (size_t j = 0; j < n_colors; ++j)
            {
                unique_color_sets[new_offset + j] = ucs[j];
            }
            new_offset += n_colors;
        }

        this->color_set_ids.resize(n_finimizers); // ids sorted based on the frequency of fmins length
        CompressedColorSets CCS(deduplicated_cs, n_colors, this->color_set_ids);

        // Assign to final structure
        this->CCS = std::move(CCS);

        cerr << "Deal with tails" << endl;

        vector<uint8_t> real_tlen_freq;
        for (auto &f : cf.lengths_by_freq)
        {
            if (f > plen)
            {
                real_tlen_freq.push_back(f - plen);
            } // One could modify cf.lengths_by_freq directly if no value was <= plen
        }

        set_tlen_order(real_tlen_freq);

        int64_t first_nonegative_tail_idx = -1;
        int64_t f_start = 0;
        for (int64_t i = 0; i < n_finimizers; i++)
        {
            if (cf.lengths[i] >= plen)
            {
                first_nonegative_tail_idx = i;
                break;
            }
            f_start += cf.lengths[i];
        }
        true_or_crash(first_nonegative_tail_idx >= 0, "ERROR: all tails shorter than prefix length");

        std::string_view cur_prefix(cf.concat.data() + f_start, plen);
        vector<std::string_view> cur_tails;
        vector<uint64_t> cur_color_set_ids;

        uint64_t p_int = prefix2int(cur_prefix, 0, plen);

        f_start = 0; // Go back to zero

        non_empty_buckets.reserve(n_buckets);
        for (int64_t i = 0; i < n_finimizers; i++)
        {
            if (cf.lengths[i] < plen)
            {
                std::string_view sprefix(cf.concat.data() + f_start, cf.lengths[i]);
                uint64_t sp_int = prefix2int(sprefix, 0, cf.lengths[i]);
                sB[sp_int] = {cf.lengths[i], this->color_set_ids[i]};
            }
            else
            {
                std::string_view prefix(cf.concat.data() + f_start, plen);
                true_or_crash(f_start + plen <= cf.concat.size(),
                              "ERROR: out-of-bounds prefix access");

                if (prefix != cur_prefix)
                {
                    // Bucket changes -> encode currently collected tails
                    non_empty_bv[p_int] = 1; // mark non-empty buckets
                    non_empty_buckets.emplace_back(Bucket(cur_tails, cur_color_set_ids, tlen_rank_map));
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

        if (!cur_tails.empty())
        {                            // Last bucket
            non_empty_bv[p_int] = 1; // mark non-empty buckets
            non_empty_buckets.emplace_back(Bucket(cur_tails, cur_color_set_ids, tlen_rank_map));
        }
        non_empty_buckets.shrink_to_fit();

        // Check the density of non-empty buckets
        size_t marked = sdsl::util::cnt_one_bits(non_empty_bv);
        cerr << "Marked " << marked << " non-empty buckets out of " << n_buckets << endl;
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

    void buckets_stats()
    {
        cerr << plen << endl;
        cerr << "buckest sizes" << endl;
        for (const auto &bucket : non_empty_buckets)
        {
            cerr << bucket.color_set_ids.size() << " ";
        }
        cerr << endl;

        cerr << "ntails per tlen" << endl;
        int64_t pos = 0; // start from 0 now that we have a single vector

        int64_t tails_so_far = 0;

        uint64_t word_index = 0;
        uint8_t w_offset = 0;

        for (const auto &bucket : non_empty_buckets)
        {
            int64_t pos = 0; // start from 0 now that we have a single vector

            int64_t tails_so_far = 0;

            uint64_t word_index = 0;
            uint8_t w_offset = 0;

            const uint64_t *data = bucket.tail_data.data();
            while (pos < bucket.tail_data.size() - 128 - 4)
            { // break the loop once something is found
                // 1. check the length of the first tail, 5‐bit tlen
                word_index = pos / 64;
                w_offset = pos % 64;
                if (word_index >= bucket.tail_data.size())
                    continue;
                uint8_t tlen = (uint8_t)sdsl::bits::read_int(&data[word_index], w_offset, 5);
                if (tlen == 0)
                {
                    cerr << "[" << 0 << "," << 0 << "]";
                    continue;
                }
                pos += 5;

                // 2. check how many tails, vbyte #tails
                uint64_t ntails = 0;
                int shift = 0;
                uint8_t byte;
                do
                {
                    // if (shift >= 64) { throw std::runtime_error("Invalid VByte: too long");}
                    word_index = pos / 64; // >> 6
                    w_offset = pos % 64;   // & 63

                    byte = sdsl::bits::read_int(&data[word_index], w_offset, 8);
                    pos += 8;
                    ntails |= uint64_t(byte & 0x7F) << shift;
                    shift += 7;
                } while (byte & 0x80);

                pos += (tlen * ntails * 2);
                // 5. if fmin not found, add the tails seen so far
                // tails_so_far += ntails;
                cerr << "[" << (int)tlen << "," << ntails << "] ";
            }
            cerr << endl;
        }
    }
    // this->non_empty_buckets, this->non_empty_bv

    void search(const std::string &query, vector<int64_t> &results, vector<int64_t> &Finimizers, vector<int64_t> &last_seen, const bool count_bases) const
    {
        const int64_t query_len = query.length();
        if (query_len < this->k)
            return;

        // Forward finimizer search
        // vector<int64_t> Finimizers;
        Finimizers.clear();
        Finimizers.reserve(query_len - k + 1);
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(query, this->non_empty_buckets, this->non_empty_bv, this->non_empty_bv_rs, this->sB, this->plen, this->k, Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }
        // Color sets
        {
            auto start = std::chrono::high_resolution_clock::now();
            if (count_bases) {count_single_base(Finimizers, this->CCS, this->n_colors, results, last_seen);}
            else {pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, results);}
            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
    }

    // Threshold-based search: returns minimum value and fills ans
    uint64_t search(const std::string &query, vector<int64_t> &results, const float t, vector<int64_t> &Finimizers, vector<int64_t> &last_seen, const bool count_bases) const
    {

        const int64_t query_len = query.length();
        if (query_len < this->k)
            return 0;

        // Forward finimizer search
        // vector<int64_t> Finimizers;
        Finimizers.clear();
        Finimizers.reserve(query_len - k + 1);
        {
            auto start = std::chrono::high_resolution_clock::now();
            rarest_fmin_streaming_search(query, this->non_empty_buckets, this->non_empty_bv, this->non_empty_bv_rs, this->sB, this->plen, this->k, Finimizers);
            auto end = std::chrono::high_resolution_clock::now();
            time_rarest_fmin += (end - start);
        }

        // Color sets
        int64_t threshold;
        {
            auto start = std::chrono::high_resolution_clock::now();
            // min_value = combine_f_rc(Finimizers, r_Finimizers, this->CCS, this->n_colors, ans, t, results);
            // T = pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, results, t);
            if (count_bases) {threshold=count_single_base(Finimizers, this->CCS, this->n_colors, results, last_seen, t);}
            else {threshold=pseudoalignment_stats(Finimizers, this->CCS, this->n_colors, results, t);}

            auto end = std::chrono::high_resolution_clock::now();
            time_combine += (end - start);
        }
        return threshold;
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

        // non_empty_buckets
        size_t num_buckets = non_empty_buckets.size();
        out.write(reinterpret_cast<const char *>(&num_buckets), sizeof(num_buckets));
        for (const auto &bucket : non_empty_buckets)
        {
            bucket.serialize(out);
        }
        cerr << non_empty_buckets.size() << endl;

        // non_empty_bv
        non_empty_bv.serialize(out);

        // sB
        /* bool has_sB = !sB.empty();
        out.write(reinterpret_cast<const char *>(&has_sB), sizeof(has_sB));

        if (has_sB)
        { */
        size_t map_size = sB.size();
        out.write(reinterpret_cast<const char *>(&map_size), sizeof(map_size));
        for (const auto &[key, val] : sB)
        {
            out.write(reinterpret_cast<const char *>(&key), sizeof(uint32_t));
            out.write(reinterpret_cast<const char *>(&val.first), sizeof(uint8_t));
            out.write(reinterpret_cast<const char *>(&val.second), sizeof(int64_t));
        }
        //}

        // metadata
        out.write(reinterpret_cast<const char *>(&n_colors), sizeof(n_colors));
        out.write(reinterpret_cast<const char *>(&n_finimizers), sizeof(n_finimizers));
        out.write(reinterpret_cast<const char *>(&plen), sizeof(plen));
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

        // non_empty_buckets
        size_t num_buckets;
        in.read(reinterpret_cast<char *>(&num_buckets), sizeof(num_buckets));
        non_empty_buckets.resize(num_buckets);
        for (size_t i = 0; i < num_buckets; ++i)
        {
            non_empty_buckets[i].load(in);
        }

        // non_empty_bv
        non_empty_bv.load(in);
        sdsl::util::init_support(non_empty_bv_rs, &non_empty_bv); // build tables

        // sB
        /* bool has_sB = false;
        in.read(reinterpret_cast<char *>(&has_sB), sizeof(has_sB));

        sB.clear();
        if (has_sB)
        { */
        size_t map_size;
        in.read(reinterpret_cast<char *>(&map_size), sizeof(map_size));
        sB.clear();
        for (size_t i = 0; i < map_size; ++i)
        {
            uint32_t key;
            std::pair<uint8_t, int64_t> val;
            in.read(reinterpret_cast<char *>(&key), sizeof(uint32_t));
            in.read(reinterpret_cast<char *>(&val.first), sizeof(uint8_t));
            in.read(reinterpret_cast<char *>(&val.second), sizeof(int64_t));
            sB[key] = val;
        }
        //}

        // metadata
        in.read(reinterpret_cast<char *>(&n_colors), sizeof(n_colors));
        in.read(reinterpret_cast<char *>(&n_finimizers), sizeof(n_finimizers)); // TODO do we need this?
        in.read(reinterpret_cast<char *>(&plen), sizeof(plen));
        in.read(reinterpret_cast<char *>(&k), sizeof(k));

        in.close();
        cerr << "DONE" << endl;
    }
};

inline void process_word(uint64_t word, const uint64_t base, vector<uint64_t> &results, const uint64_t freq)
{
    while (word)
    {
        // uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        results[base + bit] += freq;
        word &= word - 1; // clear lowest bit
    }
}

inline void read_bv(const uint64_t *data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<uint64_t> &results)
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

inline void process_word(uint64_t word, const uint64_t base, vector<int64_t> &results, const uint64_t freq)
{
    while (word)
    {
        // uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        results[base + bit] += freq;
        word &= word - 1; // clear lowest bit
    }
}

inline void read_bv(const uint64_t *data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<int64_t> &results)
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

void read_verydense(const int64_t pos, const uint64_t freq, const DeltaSet &EF, const vector<uint64_t> &L, const uint64_t n_colors, vector<uint64_t> &results)
{
    const size_t end = EF.get_start(pos); // exclusive end
    size_t start = EF.get_start(pos - 1); // inclusive start
    for (size_t c = 0; c < n_colors; c++)
    {
        if (start < end && L[start] == c)
        {
            start++;
        }
        else
        {
            results[c] += freq;
        }
    }
}

void read_verydense(const int64_t pos, const uint64_t freq, const DeltaSet &EF, const vector<uint16_t> &L, const uint64_t n_colors, vector<int64_t> &results, uint64_t &dense)
{
    dense += freq;
    const size_t end = EF.get_start(pos); // exclusive end
    size_t start = EF.get_start(pos - 1); // inclusive start
    while (start < end)
    {
        results[L[start++]] -= freq;
    }
}

void read_colors(const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, const vector<pair<int64_t, uint64_t>> &fmin_v, uint64_t &dense)
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
        if (pos < sparse_count)
        {
            // read from L
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

inline void process_word_bases(uint64_t word, const uint64_t base, vector<int64_t> &results, const uint64_t freq, vector<int64_t> &last_seen, const uint64_t x,  const int64_t k)
{
    while (word)
    {
        // uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        
        int64_t diff = x + k - 2 - freq - last_seen[base + bit];
        results[base + bit] += min(diff,k) + freq -1;
        last_seen[base + bit] = x +k -2;

        word &= word - 1; // clear lowest bit
    }
    
                    
}
inline void read_bv_bases(const uint64_t *data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<int64_t> &results,  vector<int64_t> &last_seen, const uint64_t x, const int64_t k)
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
        process_word_bases(word, color_id, results, freq, last_seen, x, k);

        ptr++;
        color_id += bits_to_read;
        bits_left -= bits_to_read;
    }

    // 2. Read aligned words in btw
    while (bits_left >= 64)
    {
        uint64_t word = *ptr++;
        process_word_bases(word, color_id, results, freq, last_seen, x, k);
        color_id += 64;
        bits_left -= 64;
    }

    // 3. Read the last word (if any)
    if (bits_left > 0)
    {
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;
        process_word_bases(word, color_id, results, freq, last_seen, x, k);
    }
}

inline void read_colors_bases (const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen, const int64_t pos, uint64_t freq, const uint64_t x, const int64_t k ){
    freq=1;
    const sdsl::bit_vector &BV = CCS.getBV();
    const uint64_t *data = BV.data();
    const vector<uint16_t> &L = CCS.getL();
    const DeltaSet &EF = CCS.getEF();

    // pos < sparse_count; [sparse]
    // sparse_count <= pos < dense_count; [very dense]
    // pos >= dense; [bitmap]

    const uint64_t sparse_count = CCS.sparse_count;
    const uint64_t dense_count = CCS.dense_count;

    if ( pos < sparse_count)
        {
            // read from L
            const size_t end = EF.get_start(pos); // exclusive end
            size_t start = EF.get_start(pos - 1); // inclusive start
            while (start < end)
            {
                int64_t diff = (x + k - 2) - freq -last_seen[L[start]];
                results[L[start]] += min(diff,k) + freq -1;
                last_seen[L[start]] = x + k - 2;
                start++;
                
            }
        }
    else if (pos < dense_count)
        {
            // read complementary values from L
            //read_verydense(pos, freq, EF, L, n_colors, results);
            const size_t end = EF.get_start(pos); // exclusive end
            size_t start = EF.get_start(pos - 1); // inclusive start
            for (size_t c = 0; c < n_colors; c++)
            {
                if (start < end && L[start] == c)
                {
                    start++;
                }
                else
                {
                    int64_t diff = (x + k - 2) - freq - last_seen[c];
                    results[c] += min(diff,k) + freq -1;
                    last_seen[c] = x +k -2;
                }
            }
        }
         // BV
    else
        {
        uint64_t start = pos - dense_count;
        read_bv_bases(data, start, freq, n_colors, results, last_seen, x, k);
    }
        
}

inline void process_word_single_base(uint64_t word, const uint64_t base, vector<int64_t> &results, vector<int64_t> &last_seen, const uint64_t x,  const int64_t k)
{
    while (word)
    {
        // uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        
        int64_t diff = x + k - last_seen[base + bit];
        results[base + bit] += min(diff,k);
        last_seen[base + bit] = x + k;

        word &= word - 1; // clear lowest bit
    }
    
                    
}

inline void read_bv_single_base(const uint64_t *data, const uint64_t start, const uint64_t n_colors, vector<int64_t> &results,  vector<int64_t> &last_seen, const uint64_t x, const int64_t k)
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
        process_word_single_base(word, color_id, results, last_seen, x, k);

        ptr++;
        color_id += bits_to_read;
        bits_left -= bits_to_read;
    }

    // 2. Read aligned words in btw
    while (bits_left >= 64)
    {
        uint64_t word = *ptr++;
        process_word_single_base(word, color_id, results, last_seen, x, k);
        color_id += 64;
        bits_left -= 64;
    }

    // 3. Read the last word (if any)
    if (bits_left > 0)
    {
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;
        process_word_single_base(word, color_id, results, last_seen, x, k);
    }
}

inline void read_colors_single_base_covered (const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen, const int64_t pos, const uint64_t x){
    const int64_t k = 31;
    const sdsl::bit_vector &BV = CCS.getBV();
    const uint64_t *data = BV.data();
    const vector<uint16_t> &L = CCS.getL();
    const DeltaSet &EF = CCS.getEF();

    // pos < sparse_count; [sparse]
    // sparse_count <= pos < dense_count; [very dense]
    // pos >= dense; [bitmap]

    const uint64_t sparse_count = CCS.sparse_count;
    const uint64_t dense_count = CCS.dense_count;

    if (pos < sparse_count)
        {
            // read from L
            const size_t end = EF.get_start(pos); // exclusive end
            size_t start = EF.get_start(pos - 1); // inclusive start

            while (start < end)
            {
                int64_t diff = (x + k) -last_seen[L[start]];
                results[L[start]] += min(diff,k);
                last_seen[L[start]] = x + k;
                start++;
            }

        }
    else if (pos < dense_count)
        {
            // read complementary values from L
            //read_verydense(pos, freq, EF, L, n_colors, results);
            const size_t end = EF.get_start(pos); // exclusive end
            size_t start = EF.get_start(pos - 1); // inclusive start
            for (size_t c = 0; c < n_colors; c++)
            {
                if (start < end && L[start] == c)
                {
                    start++;
                }
                else
                {
                    int64_t diff = (x + k) - last_seen[c];
                    results[c] += min(diff,k);
                    last_seen[c] = x +k;
                }
            }
        }
         // BV
    else
        {
        uint64_t start = pos - dense_count;
        read_bv_single_base(data, start, n_colors, results, last_seen, x, k);
    }
    // Update the total tracker after processing all colors
    last_seen[last_seen.size()-1] = x + k;
    // Track total bases covered across all colors
    int64_t max_coverage = 0;
    for (size_t i = 0; i < results.size() - 1; i++) {
        max_coverage = std::max(max_coverage, results[i]);
    }
    results[results.size()-1] = max_coverage;
}

inline void count_single_base(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen)
{   
    if (Fmin.empty()) return;
    // read colors
    const int64_t k = 31;
    
    //vector<int64_t> last_seen(n_colors, -1);
    // Initialize last_seen to -k so that the first finimizer at x=0 covers exactly k bases
    std::fill(last_seen.begin(), last_seen.end(), -k);
    std::fill(results.begin(), results.end(), 0);

    // results is the bases counter
    uint64_t x = 0;
    
    while ( x < Fmin.size()){
        if (Fmin[x] != -1){read_colors_single_base_covered(CCS, n_colors, results, last_seen, Fmin[x], x);}
        x++;
    }
    return;
}


inline int64_t count_single_base(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen, const float t)
{   
    if (Fmin.empty()) return 0;
    // read colors
    const int64_t k = 31;
    
    //vector<int64_t> last_seen(n_colors, -1);
    // Initialize last_seen to -k so that the first finimizer at x=0 covers exactly k bases
    std::fill(last_seen.begin(), last_seen.end(), -k);
    std::fill(results.begin(), results.end(), 0);

    // results is the bases counter
    uint64_t x = 0;
    while ( x < Fmin.size()){
        if (Fmin[x] != -1){read_colors_single_base_covered(CCS, n_colors, results, last_seen, Fmin[x], x);}
        x++;
    }
    // return maximum number of bases covered 
    cerr << results[results.size()-1] << endl;
    return t*results[results.size()-1];
}

inline void count_bases(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, vector<int64_t> &last_seen)
{   
    if (Fmin.empty()) return;
    // read colors
    const int64_t k = 31;
    
    //vector<int64_t> last_seen(n_colors, -1);
    std::fill(last_seen.begin(), last_seen.end(), -1);
    std::fill(results.begin(), results.end(), 0);

    // results is the bases counter
    int64_t prev_f = Fmin[0];
    uint64_t freq = 1;
    int64_t pos;
    uint64_t x = 0;
    /* while ( x < Fmin.size() and Fmin[x]== -1){
        prev_f = Fmin[++x];
    }
    if (x==0){x++;} */
    //cerr << prev_f<< endl;
    while ( x < Fmin.size()){
        /* if (prev_f == Fmin[x]){
            freq++;
        }
        else{ */
            // read the colors and calculate this for each of them
            // NOW fmin results are NOT sorted
            //pos = prev_f; 
            pos = Fmin[x];
            if (pos != -1){read_colors_bases(CCS, n_colors, results, last_seen, pos, freq, x+1, k);}
            // reset
            freq = 1;
            //if (Fmin[x]!=-1) {
                prev_f = Fmin[x];
            //}
        //}
        x++;
    }
    // Check the last value(s)
    /* pos = Fmin[Fmin.size()-1];
    if (pos != -1){read_colors_bases(CCS, n_colors, results, last_seen, pos, freq, x, k);}
     */

    // Sort results so that the output is sorted
    //counting_sort(results, ans, found_fmin, n_colors);
    //cerr << "end count_bases"<< endl;
    return;
}

inline void pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results)
{
    // vector<int64_t> results(n_colors, 0);
    /* if (results.size() != n_colors) {
        results.assign(n_colors, 0);
    } else {
        std::fill(results.begin(), results.end(), 0);
    }
 */
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
    uint64_t dense = 0;
    read_colors(CCS, n_colors, results, fmin_v, dense);
    // TODO keep track of which counters were incremented and set to zero only those
    // Add number of dense sets
    for (auto &r : results)
    {
        r += dense;
    }
    const size_t found_fmin = Fmin.size(); // # total finimizers
    // Sort results so that the output is sorted
    // counting_sort(results, ans, found_fmin, n_colors);
    return;
}

inline int64_t pseudoalignment_stats(vector<int64_t> &Fmin, const CompressedColorSets &CCS, const uint64_t n_colors, vector<int64_t> &results, const float t)
{ // vector<uint64_t>& results,
    if (Fmin.empty())
    {
        return 1;
    } // not 0 as everything wuold be >=

    // vector<int64_t> results(n_colors, 0);
    /* if (results.size() != n_colors) {
        results.assign(n_colors, 0);
    } else {
        std::fill(results.begin(), results.end(), 0);
    } */

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

    // Check the values above the minimum in search
    const size_t found_fmin = Fmin.size(); // # total finimizers
    int64_t T = found_fmin * t;
    uint64_t dense = 0;

    read_colors(CCS, n_colors, results, fmin_v, dense);

    // If we care about the number of matches, add number of dense sets
    for (auto& r: results){r+=dense;}

    // adjust T
    // T -= dense;
    // for (auto& r :results){ r+= dense;}
    //  Sort results so that the output is sorted
    // counting_sort(results, ans, found_fmin, n_colors);
    return T;
}
