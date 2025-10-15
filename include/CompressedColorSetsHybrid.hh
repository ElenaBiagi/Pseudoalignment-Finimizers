#pragma once

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cassert>

#include "hybrid.hpp"  // from fulgor
#include "sdsl/bit_vectors.hpp"

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

namespace fulgor {

class CompressedColorSetsHybrid {
private:
    hybrid m_hybrid;                     // Stores all color sets in hybrid format
    vector<uint64_t> m_color_set_ids;    // Maps original color set IDs to hybrid indices
    uint64_t m_num_colors;               // Total number of colors

public:
    CompressedColorSetsHybrid() = default;

    
    // Build the hybrid color set structure from deduplicated color sets.
    // deduplicated_cs: unordered_map of bit_vector -> vector of original color set IDs
    // n_colors
    // color_set_ids: output vector mapping original color set IDs to hybrid indices    
    CompressedColorSetsHybrid(
        const unordered_map<sdsl::bit_vector, vector<size_t>, BVHash, BVEqual>& deduplicated_cs,
        uint64_t n_colors,
        vector<uint64_t>& color_set_ids
    ) {
        if (n_colors == 0) throw runtime_error("n_colors must be > 0");
        m_num_colors = n_colors;

        hybrid::builder hb(n_colors);
        color_set_ids.resize(0);

        size_t total_color_sets = 0;
        for (auto& [bv, ids] : deduplicated_cs) {
            total_color_sets += ids.size();
        }
        color_set_ids.resize(total_color_sets);

        uint64_t idx = 0; // global color set index
        for (auto& [bv, ids] : deduplicated_cs) {
            // Convert sdsl::bit_vector to sorted uint32_t array
            vector<uint32_t> colors;
            for (uint32_t c = 0; c < n_colors; ++c) {
                if (bv[c]) colors.push_back(c);
            }

            hb.encode_color_set(colors.data(), colors.size());

            // assign hybrid indices to original IDs
            for (auto id : ids) {
                color_set_ids[id] = idx;
                idx++;
            }
        }

        hb.build(m_hybrid);
        m_color_set_ids = color_set_ids;

        cerr << "Hybrid color sets built: " << total_color_sets << " sets, "
             << "total colors: " << n_colors << endl;
    }


    // Access a color set by its hybrid index.
    hybrid::forward_iterator get_color_set(uint64_t hybrid_idx) const {
        assert(hybrid_idx < m_hybrid.num_color_sets());
        return m_hybrid.color_set(hybrid_idx);
    }

    // Map original color set ID to hybrid index.
    uint64_t get_hybrid_index(uint64_t original_id) const {
        assert(original_id < m_color_set_ids.size());
        return m_color_set_ids[original_id];
    }

    uint64_t num_colors() const { return m_num_colors; }

    // total number of colors_sets
    uint64_t num_color_sets() const { return m_hybrid.num_color_sets(); }

/*     void serialize(const std::string& filename) const {
        // TODO
    }

    void load(const std::string& filename) {
        // TODO
    } */

    void print_stats() const {
        cerr << "Hybrid stats: " << endl;
        //cerr << "Number of colors: " << m_num_colors << endl;
        //cerr << "Number of color sets: " << num_color_sets() << endl;
        cerr << "Hybrid num bits: " << m_hybrid.num_bits() << endl;
    }
};

} // namespace fulgor
