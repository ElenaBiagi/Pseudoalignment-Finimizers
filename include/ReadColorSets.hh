#pragma once

#include <vector>
#include<bitset>
#include "common.hh"
#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>

inline void process_word(uint64_t word, const uint64_t base, vector<uint64_t>& results, const uint64_t freq, vector<uint64_t>& non_zero_count_indices) {
    while (word) {
        //uint64_t bit = std::countr_zero(word);
        uint64_t bit = __builtin_ctzll(word);
        results[base + bit] += freq;
        non_zero_count_indices.push_back(base+bit);
        word &= word - 1; // clear lowest bit
    }
}


inline void read_bv(const uint64_t* data, const uint64_t start, const uint64_t freq, const uint64_t n_colors, vector<uint64_t>& results, vector<uint64_t>& non_zero_count_indices){

    const uint64_t* ptr = data + (start * n_colors) / 64;
    uint64_t bit_offset = (start * n_colors) % 64;

    uint64_t color_id = 0;
    uint64_t bits_left = n_colors;


    if (bit_offset != 0){
        // 1. Read the first word
        uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
        uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
        uint64_t word = (*ptr >> bit_offset) & mask;
        process_word(word,color_id, results, freq, non_zero_count_indices);

        ptr++;
        color_id += bits_to_read;
        bits_left -= bits_to_read;
    }

    // 2. Read aligned words in btw
    while (bits_left >= 64) {
        uint64_t word = *ptr++;
        process_word(word, color_id, results, freq, non_zero_count_indices);
        color_id += 64;
        bits_left -= 64;
    }

    // 3. Read the last word (if any)
    if (bits_left > 0) {
        uint64_t mask = ((1ULL << bits_left) - 1);
        uint64_t word = *ptr & mask;
        process_word(word,color_id, results, freq, non_zero_count_indices);

    }
}

