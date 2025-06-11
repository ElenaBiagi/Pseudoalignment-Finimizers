
// GENERAL IDEA: create distinct int vectors of diff widths and combine them together into a 64bit int keeping track of number and width 
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <string.h>
#include <random>
#include <ctime> 
#include <chrono>
#include <unistd.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>


#include "common.hh"

using namespace std;

void printBinary(uint64_t v){ // prints in the reverse order
   for(int i=0;i<64;i++){
      cerr << ((v>>i)&1);
      //cerr << ((v>>(63-i))&1);
   }
   //cerr << '\n';
}



#define hasless(x,n) (((x)-~0UL/255*(n))&~(x)&~0UL/255*128)

#define haszero(v, W, mask2, mask3, Wmask) ((((v) - mask2) & ~(v) & mask3)&(~Wmask))


#define hasvaluesupply(x,W,mask, mask2, mask3, Wmask) \
(haszero((x ^ mask), W, mask2, mask3, Wmask))


// TODO: Adapt createMask to this as this is faster
constexpr uint64_t masks[10][3] = {
   {0,0,0}, //0th entry is not to be used!
   {0,0,0},
   {0,0,0},//{0b0101010101010101010101010101010101010101010110010101010101010101,0,0}, //2-bit patterns (32 of them) // at most 4 of them in our case
   {0,0,0},
   {0b0001000100010001000100010001000100010001000100010001000100010001,0,0}, //4-bit patterns (16 of them)
   {0b000010000100001000010000100001000010000100001000010000100001,0,0}, //5-bit (12 of them)
   {0b000001000001000001000001000001000001000001000001000001000001,0,0}, //6 (10)
   {0b000000100000010000001000000100000010000001000000100000010000001,0,0}, //7 (9)
   {0b0000000100000001000000010000000100000001000000010000000100000001,0,0}, //8 (8)
   {0b000000001000000001000000001000000001000000001000000001000000001,0,0}, //9 (7)
};

// TODO REMOVE ?
constexpr uint64_t rmasks[10][3] = {
   {0,0,0}, //0th entry is not to be used!
   {0,0,0},
   {0,0,0},//{0b1010101010101010101010101010101010101010101010101010101010101010,0,0}, //4-bit patterns (16 of them)
   {0,0,0},
   {0b1000100010001000100010001000100010001000100010001000100010001000,0,0}, //4-bit patterns (16 of them)
   {0b100001000010000100001000010000100001000010000100001000010000,0,0}, //5-bit (12 of them)
   {0b100000100000100000100000100000100000100000100000100000100000,0,0}, //6 (10)
   {0b100000010000001000000100000010000001000000100000010000001000000,0,0}, //7 (9)
   {0b1000000010000000100000001000000010000000100000001000000010000000,0,0}, //8 (8)
   {0b100000000100000000100000000100000000100000000100000000100000000,0,0}, //9 (7)
};

constexpr std::tuple<uint64_t, uint64_t, uint64_t> createMask(uint64_t key, int keyLen) {
    uint64_t kmask = 0;
    uint64_t mask2 = 0;
    uint64_t mask3 = 0;

    const uint64_t pattern2 = 1ULL;                    // 1 preceded by KeyLen -1 0s
    const uint64_t pattern3 = (1ULL << (keyLen - 1));  // 1 followed by KeyLen -1 0s

    int pos = 0;
    while (pos + keyLen <= 64) {
        kmask |= (key << pos);
        mask2 |= (pattern2 << pos);
        mask3 |= (pattern3 << pos);
        pos += keyLen;
    }

    return {kmask, mask2, mask3};
}

uint64_t read_unaligned_64bits(const uint64_t* data, size_t total_bits, size_t offset_bits) {
   //assert(offset_bits + 64 <= total_bits); // total_bits is useless

   size_t word_index = offset_bits / 64;
   size_t bit_offset = offset_bits % 64;

   if (bit_offset == 0) {
        return data[word_index];
   }

   uint64_t first_w = data[word_index] >> bit_offset;
   
   uint64_t second_w = data[word_index + 1] << (64-bit_offset);  // Padding is now 128 // ERROR heap-buffer-overflow

    return first_w | second_w;
}

inline int64_t slam(const sdsl::int_vector<1> &T, const uint64_t* data, const int64_t offset, const uint8_t W, const uint64_t key, const uint16_t ntails){
   // input T, offset at which the true tails start, W(tlen), key, #tails 

   // Look at 64 bits at a time starting from offset (skip tlen and ntail (5+8+?))
   // Look at tlen*ntail*2 bits
   // Mask all the bits after that

   auto masks = createMask(key, W);
   uint64_t mask  = std::get<0>(masks);
   uint64_t mask2 = std::get<1>(masks);
   uint64_t mask3 = std::get<2>(masks);

   uint64_t tails_per_word = std::min<uint64_t>(64 / W, ntails);

   uint64_t j = 0;
   for (uint64_t i = 0; i < ntails; i+=tails_per_word) {
      size_t bit_offset = offset + i * W; // size_t bit_offset = offset + i * 64;// size_t bit_offset = offset + (i - word_index) * 64;
      uint64_t w = read_unaligned_64bits(data, T.size(), bit_offset); 

      uint64_t tails_in_this_group = std::min(tails_per_word, ntails - i);

      uint64_t bits_used = tails_in_this_group * W;
      uint64_t Wmask = (bits_used < 64) ? (~0ULL << bits_used) : 0;

      uint64_t found = hasvaluesupply(w,W,mask,mask2, mask3, Wmask);
      
      if(found){
         uint64_t lz = __builtin_clzll(found);
         bool needsCorrection = (found>>(63-lz-W))&1;
         return (j*(tails_per_word))+((64-lz)/W)-1-needsCorrection;
      }
      j++; 
   // 1. I'm only looking at words that start at 0 -> no need to shift masks [OK]
   // 2. the word starts at 0 so no smaller tails -> no need to mask smaller characters [OK]
   // 3. it is easy to know where longer tails start -> We need to MASK LONGER TAILS [OK]
   }
   return -1;
}


// output: {pos in T (to get colors), tlen}
pair<int64_t, uint8_t> bitMagicSearch(const sdsl::int_vector<1> &T, const uint64_t s, const char slen){//, vector<uint64_t>& tailsSoFar){ // we know the width of the query
   // input: T, offset in T, len substring after prefix
   
   int64_t pos = 0;// start from 0 now that we have a single vector

   int64_t tails_so_far = 0;

   uint64_t word_index = 0;
   uint8_t w_offset = 0;

   const uint64_t* data = T.data();

   // Check first the shorter lengths. stop once a match is found
   while (pos < T.size()-128){ // break the loop once something is found
      // 1. check the length of the first tail, 5‐bit tlen
      word_index = pos/64;
      w_offset = pos %64;

      uint8_t tlen = (uint8_t)sdsl::bits::read_int(&data[word_index], w_offset, 5);
      if (tlen > slen){ return {-1,0};}
      if (tlen == 0){
         return {0,0};
      }
      pos += 5;

      // 2. check how many tailS, vbyte #tails
      uint64_t ntails = 0;
      int shift = 0;
      while (true) {
         word_index = pos/64;
         w_offset = pos %64;

         uint8_t byte = sdsl::bits::read_int(&data[word_index], w_offset, 8); 
         pos += 8;
         ntails |= uint64_t(byte & 0x7F) << shift;
         if ((byte & 0x80) == 0) break;
         shift += 7;
         if (shift >= 64) throw std::runtime_error("Invalid vbyte: too long");
      }      
      // 3. Extract substring
      // Starting at pos (64 - slen*2) extract the first 2*tlen bits of s
      uint64_t key = (s >> ((slen - tlen) * 2)) & ((1ULL << (tlen * 2)) - 1); // extract suffix of length tlen

      // 4. Look for substring where the tails of that length start 
      int64_t res = slam (T,data, pos, tlen*2, key, ntails); // bitwise operations 
      if (res!=-1) {
         // Add to the result the number of finimizers preceeding this. [DONE in rarest_fmin_streaming_search] 
         return {res+tails_so_far, tlen};}
      pos+= (tlen*ntails*2);
      // 5. if fmin not found, add the tails seen so far
      tails_so_far += ntails;
      //if (tlen == 21){break;} // there is no tail longer than 21 now
   }
   return {-1,0};
}