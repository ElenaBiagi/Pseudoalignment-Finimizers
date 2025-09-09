#pragma once

#include <vector>
#include <cstring>
#include <unordered_map>


#include "sdsl/bit_vectors.hpp"
#include <sdsl/enc_vector.hpp>


using namespace std;


class CompressedColorSets {
    private:
    // TODO Combine BV and L in a 
    // add a bitvector instead of branch 
    vector<uint32_t> L;
    sdsl::enc_vector<> EF;
    sdsl::bit_vector BV;

    public:
    const std::vector<uint32_t>& getL() const { return L; }
    const sdsl::enc_vector<>& getEF() const { return EF; }
    const sdsl::bit_vector& getBV() const { return BV; }

    // TODO compress colorset like Themisto
            // L + EF + BV
            // L = list of concatenated color sets
            // EF = Elias-Fano compressed ending pos in L, first pos is 0
            // BV = bitvector of concatenated color sets
            // color set ids x = if x < EF.size() {rank in L+1} (L[EF[x-1]...EF[x]])
            //                   if x >= EF.size() {x-EF.size()} (multiply by n_colors to get the index in BV[x*colors..(x+1)*colors])
            // BV indices will be added EF.size() at the end

    CompressedColorSets() = default;

    CompressedColorSets (const unordered_map<string, vector<size_t>>& deduplicated_cs, const uint64_t n_colors,  vector<uint64_t>& color_set_ids, const sdsl::bit_vector& color_sets_concat ){// cf.color_sets_concat
        // Fills in L, EF, BV
        //sdsl::bit_vector BV(deduplicated_cs.size() * n_colors); // this is too big
        vector<size_t> EF_v = {0}; // the first value has to be 0
        sdsl::bit_vector BV_v(deduplicated_cs.size() * n_colors,0); // ensure that it's all 0s
        BV.swap(BV_v);
        size_t BV_size = 0;
        sdsl::bit_vector BV_color_set_ids(color_set_ids.size(), 0);
        uint64_t new_offset = 0;
        const size_t sparse_thr = n_colors / 4; //*0.25
        //const size_t dense_thr  = (3 * n_colors) / 4; // *0.75

        const uint64_t* data = color_sets_concat.data();

        for (auto& [key, old_offsets] : deduplicated_cs) {
            const uint64_t start = old_offsets[0];

            // TODO access color_set_concat and save the value in a bv

            sdsl::bit_vector bv = read_colors_to_bv(data, n_colors, start);
            const size_t size = sdsl::util::cnt_one_bits(bv);
            // Sparse
            if (size < sparse_thr) {
                for (size_t c = 0; c < n_colors; c++) {
                    // store 1s explicitly
                    // TODO this could me more efficient
                    if (bv[c]) { L.push_back((uint32_t)(c)); } 
                }
                // color set ids = rank in L
                for (auto& c_id : old_offsets){ color_set_ids[c_id] = EF_v.size(); } // the minimum is 1
                EF_v.push_back(L.size()); // Keep track of ending pos // exclusive ends will be inclusive starts for the next interval
            } else {
                const uint64_t old_offset = old_offsets[0] * n_colors;
                //BV.resize(BV.size()+n_colors); // Resize BV every time.. not very efficient
                BV_size++; // augment every time a new color set is added
                for (size_t j = 0; j < n_colors; ++j) {
                    BV[new_offset + j] = color_sets_concat[old_offset + j];
                }
                
                // mark BV color set ids
                for (auto& c_id : old_offsets){ 
                    color_set_ids[c_id] = new_offset/n_colors;
                    BV_color_set_ids[c_id] = 1;
                }
                new_offset += n_colors;
            }
        }
        // Add EF.size() (after the loop) to the indices of BV
        for (auto b = 0; b < color_set_ids.size(); b++){
            if (BV_color_set_ids[b]){ color_set_ids[b]+= EF_v.size();}
        }
        BV.resize((BV_size * n_colors)+63);

        // Convert EF_v into real Elias-Fano econding 
        sdsl::enc_vector<> ef(EF_v);
        this->EF = std::move(ef);

        cerr << "BV: "<< (int)BV_size << endl;
        cerr << "L: " << EF_v.size() << endl;
        //uint64_t max = *std::max_element(L.begin(), L.end());
        //cerr << "Max value in L: " << max << std::endl;
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

    void serialize(const string& index_prefix) const {
        // BV
        std::ofstream BV_out(index_prefix + ".BV.sdsl", std::ios::binary);
        if (!BV_out) {
            std::cerr << "Error: Could not open BV file!" << std::endl;
            return;
        }
        sdsl::serialize(BV, BV_out);
        BV_out.close();

        // L
        std::ofstream L_out(index_prefix + ".L.BIN", std::ios::binary);
        size_t L_size = L.size();
        L_out.write(reinterpret_cast<const char*>(&L_size), sizeof(L_size));
        L_out.write(reinterpret_cast<const char*>(L.data()), L_size * sizeof(uint32_t));
        L_out.close();

        // EF
        std::ofstream EF_out(index_prefix + ".EF.sdsl", std::ios::binary);
        if (!EF_out) {
            std::cerr << "Error: Could not open EF file!" << std::endl;
            return;
        }
        sdsl::serialize(EF, EF_out);
        EF_out.close();
        cerr << "BV: "<< (BV.size()-63)/43 << endl;
        cerr << "EF: " << EF.size() << endl;

        cerr << "L: " << L.size() << endl;
    }

    void load(const string& index_prefix) {
        // BV
        std::ifstream colors_in(index_prefix + ".BV.sdsl", std::ios::binary);
        if (!colors_in) {
            std::cerr << "Error: Could not open colors file!" << std::endl;
            return;
        }
        sdsl::load(BV, colors_in);

        colors_in.close();

        // L
        std::ifstream L_in(index_prefix + ".L.BIN", std::ios::binary);
        size_t L_size;
        L_in.read(reinterpret_cast<char*>(&L_size), sizeof(L_size));
        L.resize(L_size);
        L_in.read(reinterpret_cast<char*>(L.data()), L_size * sizeof(uint32_t));
        L_in.close();

        // EF
        std::ifstream EF_in(index_prefix + ".EF.sdsl", std::ios::binary);
        if (!EF_in) {
            std::cerr << "Error: Could not open EF.sdsl !" << std::endl;
            return;
        }
        sdsl::load(EF, EF_in);
        EF_in.close();
    }

};
