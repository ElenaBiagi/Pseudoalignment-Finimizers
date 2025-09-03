#pragma once

#include <vector>
#include <cstring>
#include <unordered_map>


#include "sdsl/bit_vectors.hpp"


using namespace std;


class CompressedColorSets {
    private:
    vector<uint32_t> L;
    vector<size_t> EF = {0}; // the first value has to be 0
    sdsl::bit_vector BV;

    public:
    const std::vector<uint32_t>& getL() const { return L; }
    const std::vector<size_t>& getEF() const { return EF; }
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
        BV.resize(deduplicated_cs.size() * n_colors);
        size_t BV_size = 0;
        sdsl::bit_vector BV_color_set_ids(color_set_ids.size(), 0);
        uint64_t new_offset = 0;
        for (auto& [key, old_offsets] : deduplicated_cs) {

            sdsl::bit_vector bv(n_colors);
            memcpy((char*)bv.data(), key.data(), key.size());
            const size_t size = sdsl::util::cnt_one_bits(bv);
            const size_t sparse_thr = n_colors / 4; //*0.25
            //const size_t dense_thr  = (3 * n_colors) / 4; // *0.75

            // Sparse
            if (size < sparse_thr) {
                for (size_t c = 0; c < n_colors; c++) {
                    // store 1s explicitly
                    // TODO this could me more efficient
                    if (bv[c]) { L.push_back((uint32_t)(c)); } 
                }
                // color set ids = rank in L
                for (auto& c_id : old_offsets){ color_set_ids[c_id] = EF.size(); } // the minimum is 1
                EF.push_back(L.size()); // Keep track of ending pos // exclusive ends will be inclusive starts for the next interval
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
            if (BV_color_set_ids[b]){ color_set_ids[b]+= EF.size();}
        }
        BV.resize((BV_size * n_colors)+63);
        cerr << "BV: "<< (int)BV_size << endl;
        cerr << "L: " << EF.size() << endl;
    }

    void serialize(const string& index_prefix) const {
        std::ofstream BV_out(index_prefix + ".BV.sdsl", std::ios::binary);
        if (!BV_out) {
            std::cerr << "Error: Could not open BV file!" << std::endl;
            return;
        }
        sdsl::serialize(BV, BV_out);
        BV_out.close();

        std::ofstream L_out(index_prefix + ".L.BIN", std::ios::binary);
        size_t L_size = L.size();
        L_out.write(reinterpret_cast<const char*>(&L_size), sizeof(L_size));
        L_out.write(reinterpret_cast<const char*>(L.data()), L_size * sizeof(uint32_t));
        L_out.close();

         std::ofstream EF_out(index_prefix + ".EF.BIN", std::ios::binary);
        size_t EF_size = EF.size();
        EF_out.write(reinterpret_cast<const char*>(&EF_size), sizeof(EF_size));
        EF_out.write(reinterpret_cast<const char*>(EF.data()), EF_size * sizeof(size_t));
        EF_out.close();
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
        std::ifstream EF_in(index_prefix + ".EF.BIN", std::ios::binary);
        size_t EF_size;
        EF_in.read(reinterpret_cast<char*>(&EF_size), sizeof(EF_size));
        EF.resize(EF_size);
        EF_in.read(reinterpret_cast<char*>(EF.data()), EF_size * sizeof(size_t));
        EF_in.close();
    }

};
