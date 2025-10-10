#pragma once

#include <vector>
#include <cstring>
#include <unordered_map>


#include "sdsl/bit_vectors.hpp"
#include <sdsl/enc_vector.hpp>

#include "hybrid.hpp"

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


class CompressedColorSets {
    public:
 
    fulgor::hybrid m_hybrid;
    uint64_t n_colors;
    std::vector<uint64_t> color_set_ids;

    CompressedColorSets() = default;

    CompressedColorSets (const unordered_map<sdsl::bit_vector, vector<size_t>, BVHash, BVEqual>& deduplicated_cs, const uint64_t n_colors,  vector<uint64_t>& color_set_ids ){
        if (n_colors == 0) throw runtime_error("n_colors must be > 0");
        
        fulgor::hybrid::builder hb(n_colors);

        uint64_t next_id = 0;
        for (const auto& [bv, offsets] : deduplicated_cs) {
            vector<uint32_t> color_indices;
            for (size_t i = 0; i < bv.size(); ++i)
                if (bv[i]) {color_indices.push_back(static_cast<uint32_t>(i));}

            hb.encode_color_set(color_indices.data(), color_indices.size());

            for (auto oid : offsets) color_set_ids[oid] = next_id;
            next_id++;
        }

        hb.build(m_hybrid);
    }

    // Access a color set by its ID
    fulgor::hybrid::forward_iterator color_set(uint64_t id) const {
        return m_hybrid.color_set(id);
    }

    uint32_t num_colors() const { return m_hybrid.num_colors(); }
    uint64_t num_color_sets() const { return m_hybrid.num_color_sets(); }

    void print_stats() const {
        std::cout << "CompressedColorSets statistics:\n";
        std::cout << "  Number of colors: " << num_colors() << "\n";
        std::cout << "  Number of color sets: " << num_color_sets() << "\n";
        std::cout << "  Total bits used: " << m_hybrid.num_bits() << "\n";
    }

    void serialize(const std::string& prefix) const {
        // For simplicity, using SDsl to serialize offsets and color sets
        std::ofstream out(prefix + ".hybrid.sdsl", std::ios::binary);
        if (!out) throw std::runtime_error("Cannot open file for serialization.");
        sdsl::serialize(m_hybrid, out);
        out.close();
    }

    void load(const std::string& prefix) {
        std::ifstream in(prefix + ".hybrid.sdsl", std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open file for deserialization.");
        sdsl::load(m_hybrid, in);
        in.close();
    }

};
