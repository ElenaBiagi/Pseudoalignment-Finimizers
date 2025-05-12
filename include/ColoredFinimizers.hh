#include <vector>
#include <iostream>
#include <fstream>
#include <optional>

#include "sdsl/bit_vectors.hpp"
#include "common.hh"

using namespace std;

namespace std {
    template<>
    struct hash<std::pair<uint32_t, char>> {
        std::size_t operator()(const std::pair<uint32_t, char>& p) const {
            return std::hash<uint32_t>{}(p.first) ^ (std::hash<char>{}(p.second) << 1);
        }
    };
}

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
    }
};

class Bucket {
public:

    struct Compact_tails{
        int tlen;
        uint64_t int_tail;
        uint32_t color_set_id;

        bool operator<(const Compact_tails& other) const {
            return tlen < other.tlen;  // sort by tlen // No need for them to be in lexicographic order
        }
    };
    
    // Compressed tail data
    // Bit layout: [tail length: u8][#tails: vbyte][concat of bitpacked tails]
    //             [tail length: u8][#tails: vbyte][concat of bitpacked tails] 
    //             [tail length: u8][#tails: vbyte][concat of bitpacked tails] 
    //             ...
    sdsl::int_vector<1> tail_data;

    // Color set ids for each tail
    vector<uint32_t> color_set_ids;

    Bucket(vector<std::string_view> tails, vector<uint32_t>& unsorted_color_set_ids) {
        
        // Convert tails to Compact_tails
        vector<Compact_tails> B_tails;
        B_tails.reserve(tails.size());
        for (size_t i = 0; i < tails.size(); i++) {
            Compact_tails ct;
            ct.tlen = tails[i].size();
            ct.int_tail = prefix2int(tails[i], 0, ct.tlen);
            ct.color_set_id = unsorted_color_set_ids[i];
            B_tails.push_back(ct);
        }
        std::sort(B_tails.begin(), B_tails.end()); 
        /* for (size_t i = 0; i < B_tails; i++)
            color_set_ids[i]=B_tails[i].color_set_id;
        } */

        // This permutes the color_set_ids
        color_set_ids.reserve(unsorted_color_set_ids.size());

        WriteTailsVector(B_tails);

    }

    void WriteTailsVector(const vector<Compact_tails>& B_tails){
        // tlen == 0 is a special case -> B_tails.size() == 1
        // Compressed tails
        //sdsl::int_vector<1> tail_data;
        
        // Color set ids for every tail
        //vector<uint32_t> color_set_ids;
        //color_set_ids.resize(B_tails.size()); 

        uint64_t total_bits = 0;
        int tlen;
        int cur_tlen = 0; // tlen cannot be 0 when compared to cur_tlen
        uint32_t ntails = 0;
        vector <int> tlens;
        unordered_map<int, uint8_t> m_ntails;
        
        // B_tails must be of length at least 0;
        tlen = B_tails[0].tlen;

        if (tlen == 0){
            total_bits += 5; // tlen
        }
        else{
            for (size_t i = 0; i < B_tails.size(); ++i) {
                tlen = B_tails[i].tlen;
                ntails++;
                if (tlen != cur_tlen){ 
                    tlens.push_back(tlen);
                    m_ntails[tlen]=ntails;
                    cur_tlen = tlen;
                    total_bits += 5; // tlen
                    total_bits += vbyte_encode(ntails).size() * 8; // vbyte encoded tail count
                    total_bits += ntails * tlen * 2; // each tail uses tlen*2 bits
                    ntails = 0;
                }
            }
        }
        
        tail_data.resize(total_bits); 
        uint64_t* data = tail_data.data();
        
        int64_t offset = 0;
        uint64_t word_index = 0;
        uint8_t w_offset = 0;
        
        // if tlen changed it was never 0        
        if (tlen == 0){
            color_set_ids[0]=B_tails[0].color_set_id; // nothing changed

            word_index = offset/64;
            w_offset = offset % 64;
            sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
            offset += 5;
        } 
        else {
            cur_tlen = 0; // tlen cannot be 0

            for (size_t i=0; i < B_tails.size(); i++){

                color_set_ids[i]=B_tails[i].color_set_id;

                tlen = B_tails[i].tlen;

                if (tlen != cur_tlen){
                    cur_tlen = tlen;
                    word_index = offset/64;
                    w_offset = offset % 64;
                    sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
                    offset += 5;
                        
                    const uint32_t tnumber = m_ntails[tlen];
                    auto vb = vbyte_encode(tnumber);
                    for (uint8_t b : vb) {
                        word_index = offset/64;
                        w_offset = offset % 64;
                        sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
                        offset += 8;
                    }
                } 
                // tails
                word_index = offset/64;
                w_offset = offset % 64;
                sdsl::bits::write_int(&data[word_index], B_tails[i].int_tail, w_offset, tlen*2);
                offset += (2 * tlen);
            }
        }

        return;
    }

    void serialize(std::ostream& out) const {
    sdsl::serialize(tail_data, out);
    
    size_t size = color_set_ids.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(color_set_ids.data()), size * sizeof(uint32_t));
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
    vector<optional<Bucket>> buckets; // TODO optional
    std::unordered_map<pair<uint32_t,char>, int64_t> sB; // Create a hash table to store the finimizers shorter than the prefix length
    uint64_t n_colors;
    uint64_t n_finimizers;
    uint64_t plen;
    uint64_t k;


    CompressedColoredFinimizers(ColoredFinimizers&& cf, int64_t prefix_len) {
        uint64_t n_buckets = (1ULL << (prefix_len * 2));
        buckets.resize(n_buckets);
        // TODO: What do we have in empty Buckets??

        n_finimizers = cf.lengths.size();
        true_or_crash(n_finimizers > 0, "ERROR: 0 finimizers");

        true_or_crash(cf.color_sets_concat.size() % n_finimizers == 0, "ERROR: color set bitmap length not divisible by finimizer count");
        n_colors = cf.color_sets_concat.size() / n_finimizers;

        vector<std::string_view> cur_bucket_nonnegative_tails;
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

        uint64_t p_int;
        for(int64_t i = 0; i < n_finimizers; i++) {
            if(cf.lengths[i] < prefix_len){
                // Skip negative tails. TODO: do something about them.
                //TODO we need to store the lengths (char) as well
                std::string_view sprefix(cf.concat.data() + f_start, cf.lengths[i]);

                sB[{prefix2int(sprefix,0, cf.lengths[i]), cf.lengths[i]}]=i; // i= color_set_id
            } else {
                std::string_view prefix(cf.concat.data() + f_start, prefix_len);
                p_int = prefix2int(prefix, 0, prefix_len);
                if(prefix != cur_prefix) {
                    // Bucket changes -> encode currently collected tails
                    buckets[p_int]=Bucket(cur_tails, cur_color_set_ids);
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

        //buckets.shrink_to_fit(); // we need exactly that many buckets
        color_sets_concat = std::move(cf.color_sets_concat);
    }

/*     // Transfer ownership of the index out of the builder
    unique_ptr<CompressedColoredFinimizers> get_index(){
        return std::move(this->index);
    } */

    void serialize(const string& index_prefix) const {
        
        // color_sets_concat
        std::ofstream colors_out(index_prefix + ".colors.sdsl", std::ios::binary);
        sdsl::serialize(color_sets_concat, colors_out);
        colors_out.close();

        // TODO buckets
        // WRONG
        // std::ofstream buckets_out(index_prefix + ".buckets.BIN", std::ios::binary);
        // size_t num_buckets = buckets.size();
        // buckets_out.write(reinterpret_cast<const char*>(&num_buckets), sizeof(num_buckets));

        // for (const auto& b : buckets) {
        //     b.serialize(buckets_out);
        // }

        // buckets_out.close();


        //sB

        std::ifstream sB_in(index_prefix + ".sB.BIN", std::ios::binary);
        size_t map_size;
        sB_in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
        for (size_t i = 0; i < map_size; ++i) {
            std::pair<uint32_t, char> key;
            int64_t val;
            sB_in.read(reinterpret_cast<char*>(&key.first), sizeof(uint32_t));
            sB_in.read(reinterpret_cast<char*>(&key.second), sizeof(char));
            sB_in.read(reinterpret_cast<char*>(&val), sizeof(int64_t));
            sB[key] = val;
        }
        sB_in.close(); 

        // n_colors, n_finimizers
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
        
        
    }

};

    

/*
int main(){
    UncompressedTails ut;
    ifstream in("out.bin");
    ut.load(in);
    for(auto c : ut.concat) cout << c; cout << endl;
    for(auto len : ut.lengths) cout << (int)len << " "; cout << endl;
    cout << ut.color_sets_concat << endl;
}
*/
