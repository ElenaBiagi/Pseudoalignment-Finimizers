#pragma once

#include <vector>
#include <cstring>
#include <unordered_map>


#include "sdsl/bit_vectors.hpp"
#include <sdsl/enc_vector.hpp>

#include "deltaset.hpp"

using namespace std;

struct BVHash {
    size_t operator()(const sdsl::bit_vector& bv) const noexcept {
        const uint64_t* data = bv.data();
        size_t n64 = (bv.size() + 63) / 64;
        size_t h = 0;
        for (size_t i = 0; i < n64; ++i) {
            h ^= std::hash<uint64_t>{}(data[i]) 
                 + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

struct BVEqual {
    bool operator()(const sdsl::bit_vector& a, const sdsl::bit_vector& b) const noexcept {
        if (a.size() != b.size()) return false;
        const uint64_t* ad = a.data();
        const uint64_t* bd = b.data();
        size_t n64 = (a.size() + 63) / 64;
        for (size_t i = 0; i < n64; ++i) {
            if (ad[i] != bd[i]) return false;
        }
        return true;
    }
};

struct VectorHash {
    std::size_t operator()(const std::vector<uint64_t>& v) const noexcept {
        std::size_t h = 0;
        std::hash<uint64_t> hasher;
        for (auto x : v) {
            // Combine the hashes
            h ^= hasher(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};


class CompressedColorSets {
    private:
    // TODO Combine BV and L in a 
    // add a bitvector instead of branch 
    vector<uint16_t> L;
    DeltaSet EF;
    sdsl::bit_vector BV;

    public:
    const std::vector<uint16_t>& getL() const { return L; }
    const DeltaSet& getEF() const { return EF; }
    const sdsl::bit_vector& getBV() const { return BV; }
    uint64_t sparse_count = 0;
    uint64_t dense_count = 0;

    //TODO remove this
    float average(std::vector<float> const& v){
        if(v.empty()){
            return 0;
        }

        auto const count = static_cast<float>(v.size());
        return std::reduce(v.begin(), v.end()) / count;
    }

    CompressedColorSets() = default;
    
    CompressedColorSets (const unordered_map<vector<uint64_t>, vector<size_t>, VectorHash>& deduplicated_cs, const uint64_t n_colors,  vector<uint64_t>& color_set_ids){
        
        vector<float> BV_sizes;
        vector<float> L_sizes;
        vector<float> dense_sizes;
        
        if (n_colors == 0) throw runtime_error("n_colors must be > 0");
        // Fills in L, EF, BV
        //sdsl::bit_vector BV(deduplicated_cs.size() * n_colors); // this is too big
        
        vector<uint16_t> temp_cL;
        
        vector<size_t> temp_EF_v = {0}; // the first value has to be 0
        vector<size_t> temp_cEF_v = {}; // the first value is the last value of L
        sdsl::bit_vector temp_BV(deduplicated_cs.size() * ((n_colors * 15)/16),0); // ensure that it's all 0s
        //BV.swap(temp_BV);
        size_t BV_size = 0;
        sdsl::bit_vector BV_color_set_ids(color_set_ids.size(), 0);
        sdsl::bit_vector cL_color_set_ids(color_set_ids.size(), 0);
        uint64_t new_offset = 0;

        // TODO these could be computed once
        const size_t sparse_thr = n_colors / 16;
        const size_t dense_thr = (n_colors * 15) / 16;

        // TODO DO NOT add 0 to the beginning of ends
        // do while i<ends.size();
        // uint64t first = ends[0];
        // size = (first==0)? 1 : first;
        
        for (auto& [v, old_offsets] : deduplicated_cs) {
            const uint64_t start = old_offsets[0];
            const size_t size = v.size();

            // Sparse
            if (size < sparse_thr) {
                // TODO more efficient ?
                for (auto &value :v){
                    L.push_back(static_cast<uint16_t>(value)); 
                }
                // color set ids = rank in L
                uint64_t ef_index = temp_EF_v.size();

                for (auto& c_id : old_offsets){ color_set_ids[c_id] = ef_index; } // the minimum is 1

                
                temp_EF_v.push_back(L.size()); // Keep track of ending pos // exclusive ends will be inclusive starts for the next interval
                L_sizes.push_back((float)size);
            } else if (size < dense_thr){
            
                //const uint64_t old_offset = old_offsets[0] * n_colors;
                BV_size++; // augment every time a new color set is added

                for (auto &value :v){
                    temp_BV[new_offset + value] = 1;
                }
                
                // mark BV color set ids
                for (auto& c_id : old_offsets){ 
                    color_set_ids[c_id] = new_offset/n_colors;
                    BV_color_set_ids[c_id] = 1;
                }
                new_offset += n_colors;
                BV_sizes.push_back((float)size);

            } else {
                // very dense
                                
                uint64_t x = 0;
                for (auto &value :v){
                    while (x < value){
                        L.push_back(static_cast<uint16_t>(x));
                        x++;
                    }  
                }
                while (x < n_colors){
                    L.push_back(static_cast<uint16_t>(x));
                    x++;
                }  
                
                
                // color set ids = rank in temp_cL
                uint64_t cef_index = temp_cEF_v.size();
                for (auto& c_id : old_offsets){ 
                    color_set_ids[c_id] = cef_index;
                    cL_color_set_ids[c_id] = 1; 
                } // the minimum is 0 (+ sparse)
                temp_cEF_v.push_back(temp_cL.size()); // Keep track of ending pos // exclusive ends will be inclusive starts for the next interval
                dense_sizes.push_back((float)(n_colors - size));
            }
        }

        // counts
        this->sparse_count = temp_EF_v.size();
        
        // concatenate EF
        for (auto& v : temp_cEF_v){
            v+= temp_EF_v.back();
        }
        temp_EF_v.insert( temp_EF_v.end(), temp_cEF_v.begin(), temp_cEF_v.end() ); // temp_cEF_v[0] == temp_EF_v[-1]

        this->dense_count = temp_EF_v.size();

        // concatenate L
        L.insert( L.end(), temp_cL.begin(), temp_cL.end() );

        // shift color_set_ids 
        for (size_t b = 0; b < color_set_ids.size(); b++){
            // Add temp_EF_v.size() (after the loop) to the indices of cEF
            if (cL_color_set_ids[b]){ color_set_ids[b]+= sparse_count;}
            // Add temp_EF_v.size()+cEF.size() (after the loop) to the indices of BV
            if (BV_color_set_ids[b]){ color_set_ids[b]+= dense_count;}
        }
        temp_BV.resize((((BV_size * n_colors)+63)/64)*64); // only add the minimum number of bits to make it word aligned

        // Convert EF_v into a DeltaSet
        DeltaSet ef(temp_EF_v);

        this->EF = std::move(ef);
        this->BV = std::move(temp_BV);

        cerr << "BV: "<< (int)BV_size << endl;
        cerr << "sparse : " << sparse_count-1 << endl;
        cerr << "very dense:" << dense_count - (sparse_count-1) << endl;
        //uint64_t max = *std::max_element(L.begin(), L.end());
        //cerr << "Max value in L: " << max << std::endl;

        float BV_sizes_averge = average(BV_sizes);
        float L_sizes_averge = average(L_sizes);
        float dense_sizes_averge = average(dense_sizes);

        double max_BV = *std::max_element(BV_sizes.begin(), BV_sizes.end());
        double max_L = *std::max_element(L_sizes.begin(), L_sizes.end());
        double max_dense = *std::max_element(dense_sizes.begin(), dense_sizes.end());

        cerr << "BV sizes average = " << BV_sizes_averge << endl;
        cerr << "BV sizes max = " << max_BV << endl;

        cerr << "L sizes average = " << L_sizes_averge << endl;
        cerr << "L sizes max = " << max_L << endl;

        cerr << "dense sizes average = " << dense_sizes_averge << endl;
        cerr << "dense sizes max = " << max_dense << endl;
    }

    CompressedColorSets (const unordered_map<sdsl::bit_vector, vector<size_t>, BVHash, BVEqual>& deduplicated_cs, const uint64_t n_colors,  vector<uint64_t>& color_set_ids){
        
        vector<float> BV_sizes;
        vector<float> L_sizes;
        vector<float> dense_sizes;
        
        if (n_colors == 0) throw runtime_error("n_colors must be > 0");
        // Fills in L, EF, BV
        //sdsl::bit_vector BV(deduplicated_cs.size() * n_colors); // this is too big
        
        vector<uint16_t> temp_cL;
        
        vector<size_t> temp_EF_v = {0}; // the first value has to be 0
        vector<size_t> temp_cEF_v = {}; // the first value is the last value of L
        sdsl::bit_vector temp_BV(deduplicated_cs.size() * n_colors,0); // ensure that it's all 0s
        //BV.swap(temp_BV);
        size_t BV_size = 0;
        sdsl::bit_vector BV_color_set_ids(color_set_ids.size(), 0);
        sdsl::bit_vector cL_color_set_ids(color_set_ids.size(), 0);
        uint64_t new_offset = 0;
        const size_t sparse_thr = n_colors / 16;
        const size_t dense_thr = (n_colors * 15) / 16;

        for (auto& [bv, old_offsets] : deduplicated_cs) {
            const uint64_t start = old_offsets[0];

            // TODO access color_set_concat and save the value in a bv

            const size_t size = sdsl::util::cnt_one_bits(bv);
            // Sparse
            if (size < sparse_thr) {
                // store 1s explicitly
                // TODO more efficient ?
                for (size_t c = 0; c < n_colors; c++) {
                    if (bv[c]) {L.push_back(static_cast<uint16_t>(c));}    
                }
                // color set ids = rank in L
                uint64_t ef_index = temp_EF_v.size();
                for (auto& c_id : old_offsets){ color_set_ids[c_id] = ef_index; } // the minimum is 1
                temp_EF_v.push_back(L.size()); // Keep track of ending pos // exclusive ends will be inclusive starts for the next interval
                L_sizes.push_back((float)size);
            } else if (size < dense_thr){
            
                const uint64_t old_offset = old_offsets[0] * n_colors;
                //BV.resize(BV.size()+n_colors); // Resize BV every time.. not very efficient
                BV_size++; // augment every time a new color set is added
                for (size_t j = 0; j < n_colors; ++j) {
                    temp_BV[new_offset + j] = bv[j];
                }
                
                // mark BV color set ids
                for (auto& c_id : old_offsets){ 
                    color_set_ids[c_id] = new_offset/n_colors;
                    BV_color_set_ids[c_id] = 1;
                }
                new_offset += n_colors;
                BV_sizes.push_back((float)size);

            } else {
                // very dense
                for (size_t c = 0; c < n_colors; c++) {
                    if (!bv[c]) {temp_cL.push_back(static_cast<uint16_t>(c));}    
                }
                // color set ids = rank in temp_cL
                uint64_t cef_index = temp_cEF_v.size();
                for (auto& c_id : old_offsets){ 
                    color_set_ids[c_id] = cef_index;
                    cL_color_set_ids[c_id] = 1; 
                } // the minimum is 0 (+ sparse)
                temp_cEF_v.push_back(temp_cL.size()); // Keep track of ending pos // exclusive ends will be inclusive starts for the next interval
                dense_sizes.push_back((float)(n_colors - size));

            }
        }

        // counts
        this->sparse_count = temp_EF_v.size();
        
        // concatenate EF
        for (auto& v : temp_cEF_v){
            v+= temp_EF_v.back();
        }
        temp_EF_v.insert( temp_EF_v.end(), temp_cEF_v.begin(), temp_cEF_v.end() ); // temp_cEF_v[0] == temp_EF_v[-1]

        this->dense_count = temp_EF_v.size();

        // concatenate L
        L.insert( L.end(), temp_cL.begin(), temp_cL.end() );

        // shift color_set_ids 
        for (size_t b = 0; b < color_set_ids.size(); b++){
            // Add temp_EF_v.size() (after the loop) to the indices of cEF
            if (cL_color_set_ids[b]){ color_set_ids[b]+= sparse_count;}
            // Add temp_EF_v.size()+cEF.size() (after the loop) to the indices of BV
            if (BV_color_set_ids[b]){ color_set_ids[b]+= dense_count;}
        }
        temp_BV.resize((((BV_size * n_colors)+63)/64)*64); // only add the minimum number of bits to make it word aligned

        // Convert EF_v into a DeltaSet
        DeltaSet ef(temp_EF_v);

        this->EF = std::move(ef);
        this->BV = std::move(temp_BV);

        cerr << "BV: "<< (int)BV_size << endl;
        cerr << "sparse : " << sparse_count-1 << endl;
        cerr << "very dense:" << dense_count - (sparse_count-1) << endl;
        //uint64_t max = *std::max_element(L.begin(), L.end());
        //cerr << "Max value in L: " << max << std::endl;

        float BV_sizes_averge = average(BV_sizes);
        float L_sizes_averge = average(L_sizes);
        float dense_sizes_averge = average(dense_sizes);

        double max_BV = *std::max_element(BV_sizes.begin(), BV_sizes.end());
        double max_L = *std::max_element(L_sizes.begin(), L_sizes.end());
        double max_dense = *std::max_element(dense_sizes.begin(), dense_sizes.end());

        cerr << "BV sizes average = " << BV_sizes_averge << endl;
        cerr << "BV sizes max = " << max_BV << endl;

        cerr << "L sizes average = " << L_sizes_averge << endl;
        cerr << "L sizes max = " << max_L << endl;

        cerr << "dense sizes average = " << dense_sizes_averge << endl;
        cerr << "dense sizes max = " << max_dense << endl;
    }

    void serialize(std::ostream& out) const {

        // counts
        out.write(reinterpret_cast<const char*>(&sparse_count), sizeof(sparse_count));
        out.write(reinterpret_cast<const char*>(&dense_count), sizeof(dense_count));
        // BV
        sdsl::serialize(BV, out);
        // EF
        EF.serialize(out);
        // L
        size_t L_size = L.size();
        out.write(reinterpret_cast<const char*>(&L_size), sizeof(L_size));
        out.write(reinterpret_cast<const char*>(L.data()), L_size * sizeof(uint16_t));
    }

    void load(std::istream& in) {
        
        // counts
        in.read(reinterpret_cast<char*>(&sparse_count), sizeof(sparse_count));
        in.read(reinterpret_cast<char*>(&dense_count), sizeof(dense_count));
        // BV
        sdsl::load(BV, in);
        // EF
        EF.load(in);
        // L
        size_t L_size = 0;
        in.read(reinterpret_cast<char*>(&L_size), sizeof(L_size));
        L.resize(L_size);
        in.read(reinterpret_cast<char*>(L.data()), L_size * sizeof(uint16_t));
    }

};
