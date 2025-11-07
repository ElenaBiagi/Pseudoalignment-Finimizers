
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

// TODO remove
void printBinary(uint64_t v){ // prints in the reverse order
   for(int i=0;i<64;i++){
      cerr << ((v>>i)&1);
      //cerr << ((v>>(63-i))&1);
   }
   //cerr << '\n';
}

#define haszero(v, mask2, mask3, Wmask) ((((v) - mask2) & ~(v) & mask3)&(~Wmask))


#define hasvaluesupply(w, mask, mask2, mask3, Wmask) \
(haszero((w ^ mask), mask2, mask3, Wmask))


// TODO: this now works only for k=31 and plen=10
static constexpr uint64_t old_masks23[43][2] = {
   {0x0000000000000000, 0x0000000000000000}, // W = 0 (unused)
   {0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF}, // 1-bit
   {0x5555555555555555, 0xAAAAAAAAAAAAAAAA}, // 2-bit
   {0x1249249249249249, 0x4924924924924924}, // 3-bit
   {0x1111111111111111, 0x8888888888888888}, // 4-bit
   {0x0084210842108421, 0x0842108421084210}, // 5-bit
   {0x0041041041041041, 0x0820820820820820}, // 6-bit
   {0x0102040810204081, 0x4081020408102040}, // 7-bit
   {0x0101010101010101, 0x8080808080808080}, // 8-bit
   {0x0040201008040201, 0x4020100804020100}, // 9-bit
   {0x0004010040100401, 0x0802008020080200}, // 10-bit
   {0x0000100200400801, 0x0040080100200400}, // 11-bit
   {0x0001001001001001, 0x0800800800800800}, // 12-bit
   {0x0000008004002001, 0x0008004002001000}, // 13-bit
   {0x0000040010004001, 0x0080020008002000}, // 14-bit
   {0x0000200040008001, 0x0800100020004000}, // 15-bit
   {0x0001000100010001, 0x8000800080008000}, // 16-bit
   {0x0000000400020001, 0x0004000200010000}, // 17-bit
   {0x0000001000040001, 0x0020000800020000}, // 18-bit
   {0x0000004000080001, 0x0100002000040000}, // 19-bit
   {0x0000010000100001, 0x0800008000080000}, // 20-bit
   {0x0000040000200001, 0x4000020000100000}, // 21-bit
   {0x0000000000400001, 0x0000080000200000}, // 22-bit
   {0x0000000000800001, 0x0000200000400000}, // 23-bit
   {0x0000000001000001, 0x0000800000800000}, // 24-bit
   {0x0000000002000001, 0x0002000001000000}, // 25-bit
   {0x0000000004000001, 0x0008000002000000}, // 26-bit
   {0x0000000008000001, 0x0020000004000000}, // 27-bit
   {0x0000000010000001, 0x0080000008000000}, // 28-bit
   {0x0000000020000001, 0x0200000010000000}, // 29-bit
   {0x0000000040000001, 0x0800000020000000}, // 30-bit
   {0x0000000080000001, 0x2000000040000000}, // 31-bit
   {0x0000000100000001, 0x8000000080000000}, // 32-bit
   {0x0000000000000001, 0x0000000100000000}, // 33-bit
   {0x0000000000000001, 0x0000000200000000}, // 34-bit
   {0x0000000000000001, 0x0000000400000000}, // 35-bit
   {0x0000000000000001, 0x0000000800000000}, // 36-bit
   {0x0000000000000001, 0x0000001000000000}, // 37-bit
   {0x0000000000000001, 0x0000002000000000}, // 38-bit
   {0x0000000000000001, 0x0000004000000000}, // 39-bit
   {0x0000000000000001, 0x0000008000000000}, // 40-bit
   {0x0000000000000001, 0x0000010000000000}, // 41-bit
   {0x0000000000000001, 0x0000020000000000}, // 42-bit
};
static constexpr uint64_t masks23[21][2] = {
   {0x5555555555555555, 0xAAAAAAAAAAAAAAAA}, // 2-bit
   {0x1111111111111111, 0x8888888888888888}, // 4-bit
   {0x0041041041041041, 0x0820820820820820}, // 6-bit
   {0x0101010101010101, 0x8080808080808080}, // 8-bit
   {0x0004010040100401, 0x0802008020080200}, // 10-bit
   {0x0001001001001001, 0x0800800800800800}, // 12-bit
   {0x0000040010004001, 0x0080020008002000}, // 14-bit
   {0x0001000100010001, 0x8000800080008000}, // 16-bit
   {0x0000001000040001, 0x0020000800020000}, // 18-bit
   {0x0000010000100001, 0x0800008000080000}, // 20-bit
   {0x0000000000400001, 0x0000080000200000}, // 22-bit
   {0x0000000001000001, 0x0000800000800000}, // 24-bit
   {0x0000000004000001, 0x0008000002000000}, // 26-bit
   {0x0000000010000001, 0x0080000008000000}, // 28-bit
   {0x0000000040000001, 0x0800000020000000}, // 30-bit
   {0x0000000100000001, 0x8000000080000000}, // 32-bit
   {0x0000000000000001, 0x0000000200000000}, // 34-bit
   {0x0000000000000001, 0x0000000800000000}, // 36-bit
   {0x0000000000000001, 0x0000002000000000}, // 38-bit
   {0x0000000000000001, 0x0000008000000000}, // 40-bit
   {0x0000000000000001, 0x0000020000000000}, // 42-bit
};

static constexpr uint64_t tails[] = {
    32, // 2 
    16, // 4
    10, // 6
     8, // 8
     6, // 10
     5, // 12
     4, // 14
     4, // 16
     3, // 18
     3, // 20
     2, // 22
     2, // 24
     2, // 26
     2, // 28
     2, // 30
     2, // 32
     1, // 34
     1, // 36
     1, // 38
     1, // 40
     1 // 42
};


 uint64_t read_unaligned_64bits(const uint64_t* data, size_t offset_bits) {

   size_t word_index = offset_bits / 64;
   size_t bit_offset = offset_bits % 64;

   if (bit_offset == 0) {
        return data[word_index];
   }

   const uint64_t first_w = data[word_index] >> bit_offset;
   
   const uint64_t second_w = data[word_index + 1] << (64-bit_offset);  // Padding is now 128

   return first_w | second_w;
}

// TODO T is useless
inline int64_t SearchTail(const sdsl::int_vector<1> &T, const uint64_t* data, const int64_t offset, const uint8_t W, const uint64_t key, const uint64_t ntails){
   // input T, offset at which the true tails start, W(tlen), key, #tails 

   // Look at 64 bits at a time starting from offset (skip tlen and ntail (5+8+?))
   // Look at tlen*ntail*2 bits
   // Mask all the bits after that

   const uint64_t mask2 = masks23[(W/2)-1][0];
   const uint64_t mask3 = masks23[(W/2)-1][1];
   const uint64_t mask = mask2*key; //~0ULL/255 * key;

   const uint64_t tails_per_word = std::min<uint64_t>(tails[(W/2)-1], ntails);

   uint64_t j = 0;
   for (uint64_t i  = 0; i < ntails; i+=tails_per_word) {
      const size_t bit_offset = offset + i * W; // size_t bit_offset = offset + i * 64;// size_t bit_offset = offset + (i - word_index) * 64;
      const uint64_t w = read_unaligned_64bits(data, bit_offset);

      uint64_t tails_in_this_group = tails_per_word - ((i + tails_per_word - ntails) * (((ntails - i) / tails_per_word) == 0));

      uint64_t bits_used = tails_in_this_group * W;
      uint64_t Wmask = (~0ULL << bits_used) & -(bits_used < 64);

      uint64_t found = hasvaluesupply(w,mask,mask2, mask3, Wmask);
      
      if(found){
         uint64_t lz = __builtin_clzll(found);
         int needsCorrection = (found>>(63-lz-W))&1;
         return (j*(tails_per_word))+((64-lz)/W)-1-needsCorrection;
      }
      j++;
   // 1. I'm only looking at words that start at 0 -> no need to shift masks [OK]
   // 2. the word starts at 0 so no smaller tails -> no need to mask smaller characters [OK]
   // 3. it is easy to know where longer tails start -> We need to MASK LONGER TAILS [OK]
   }
   return -1;
}

inline uint64_t extract_tail(const uint64_t* data, const int64_t bit_offset, const uint8_t W){
   read_unaligned_64bits(data, bit_offset);
   const uint64_t w = read_unaligned_64bits(data, bit_offset);
   const uint64_t mask = (1ULL << W) - 1;
   return w & mask;
}

inline int64_t Tails_binary_search(const uint64_t* data, const int64_t offset, const uint8_t W, const uint64_t key, const uint64_t ntails){
   // BINARY SEARCH
   int64_t k = 0;
   // offset == where the tails start
   // offset + W = where the 2nd tail (at index 1) starts
   // offset + (b * W) = where the b-th (+1, index 0) tail starts
   int64_t new_offset;
   for (int64_t b = ntails/2; b >= 1; b /= 2) {
      new_offset = offset + ((k+b)*W); 
      while (k+b < ntails && extract_tail(data, new_offset, W) <= key) k += b;
   }
   if (extract_tail(data, new_offset, W) == key) {
      // x found at index k
      return k;
   }
}

inline int64_t SearchTail_BS(const sdsl::int_vector<1> &T, const uint64_t* data, const int64_t offset, const uint8_t W, const uint64_t key, const uint64_t ntails){
   const uint64_t mask2 = masks23[(W/2)-1][0];
   const uint64_t mask3 = masks23[(W/2)-1][1];
   const uint64_t mask = mask2*key; //~0ULL/255 * key;

   const uint64_t tails_per_word = std::min<uint64_t>(tails[(W/2)-1], ntails);

   uint64_t j = 0;
   for (uint64_t i  = 0; i < ntails; i+=tails_per_word) {
      const size_t bit_offset = offset + i * W; // size_t bit_offset = offset + i * 64;// size_t bit_offset = offset + (i - word_index) * 64;
      const uint64_t w = read_unaligned_64bits(data, bit_offset);

      uint64_t tails_in_this_group = tails_per_word - ((i + tails_per_word - ntails) * (((ntails - i) / tails_per_word) == 0));

      uint64_t bits_used = tails_in_this_group * W;
      uint64_t Wmask = (~0ULL << bits_used) & -(bits_used < 64);

      uint64_t found = hasvaluesupply(w,mask,mask2, mask3, Wmask);
      
      if(found){
         uint64_t lz = __builtin_clzll(found);
         int needsCorrection = (found>>(63-lz-W))&1;
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
inline pair<int64_t, uint8_t> bitMagicSearch(const sdsl::int_vector<1> &T, const uint64_t s, const uint8_t slen){//, vector<uint64_t>& tailsSoFar){ // we know the width of the query
   // input: T, offset in T, len substring after prefix
   // const uint8_t slen = k - plen; 
   int64_t pos = 0;// start from 0 now that we have a single vector

   int64_t tails_so_far = 0;

   uint64_t word_index = 0;
   uint8_t w_offset = 0;

   const uint64_t* data = T.data();
   // Check first the MOST FREQUENT lengths.
   while (pos < T.size()-128-4){ // break the loop once something is found
      // 1. check the length of the first tail, 5‐bit tlen
      word_index = pos/64;
      w_offset = pos %64;

      uint8_t tlen = (uint8_t)sdsl::bits::read_int(&data[word_index], w_offset, 5);
      //if (tlen > slen){ return {-1,0};}
      if (tlen == 0){return {0,0};}
      pos += 5;

      // 2. check how many tails, vbyte #tails
      uint64_t ntails = 0;
      int shift = 0; 
      uint8_t byte;
      do {
         //if (shift >= 64) { throw std::runtime_error("Invalid VByte: too long");}
         word_index = pos/64; // >> 6
         w_offset = pos %64; // & 63
         
         byte = sdsl::bits::read_int(&data[word_index], w_offset, 8);
         pos += 8;
         ntails |= uint64_t(byte & 0x7F) << shift;
         shift += 7;
      } while (byte & 0x80);       
      
      // 3. Extract substring
      // Starting at pos (64 - slen*2) extract the first 2*tlen bits of s
      uint64_t key = (s >> ((slen - tlen) * 2)) & ((1ULL << (tlen * 2)) - 1); // extract suffix of length tlen
         
      // 4. Look for substring where the tails of that length start 
      int64_t res = -1;
      if (ntails > 50 && tlen > 3){ // TODO select a proper tail lenght and tail number
         // BINARY SEARCH
         res = SearchTail_BS(T,data, pos, tlen*2, key, ntails); // bitwise operations 
      }
      else {            
         int64_t res = SearchTail(T,data, pos, tlen*2, key, ntails); // bitwise operations 
      }
      if (res!=-1) {
         // Add to the result the number of finimizers preceeding this. [DONE in rarest_fmin_streaming_search] 
         return {res+tails_so_far, tlen};
      }

      pos+= (tlen*ntails*2);
      // 5. if fmin not found, add the tails seen so far
      tails_so_far += ntails;
   }
   return {-1,0};
}

// output: {pos in T (to get colors), tlen}
// The query is shorter than (k - plen)
pair<int64_t, uint8_t> bitMagicSearch_short(const sdsl::int_vector<1> &T, const uint64_t s, const uint8_t slen){
   // input: T, offset in T, len substring after prefix
   
   int64_t pos = 0;// start from 0 now that we have a single vector

   int64_t tails_so_far = 0;

   uint64_t word_index = 0;
   uint8_t w_offset = 0;

   const uint64_t* data = T.data();
   // Check first the longer lengths.
   while (pos < T.size()-128-4){ // break the loop once something is found
      // 1. check the length of the first tail, 5‐bit tlen
      word_index = pos/64;
      w_offset = pos %64;

      uint8_t tlen = (uint8_t)sdsl::bits::read_int(&data[word_index], w_offset, 5);
      if (tlen == 0){return {0,0};}
      pos += 5;

      // 2. check how many tails, vbyte #tails
      uint64_t ntails = 0;
      int shift = 0; 
      uint8_t byte;
      do {
         //if (shift >= 64) { throw std::runtime_error("Invalid VByte: too long");}
         word_index = pos/64; // >> 6
         w_offset = pos %64; // & 63
         
         byte = sdsl::bits::read_int(&data[word_index], w_offset, 8);
         pos += 8;
         ntails |= uint64_t(byte & 0x7F) << shift;
         shift += 7;
      } while (byte & 0x80);       
      // 3. Extract substring
      // Starting at pos (64 - slen*2) extract the first 2*tlen bits of s
      if (tlen <= slen){
         uint64_t key = (s >> ((slen - tlen) * 2)) & ((1ULL << (tlen * 2)) - 1); // extract suffix of length tlen
         int64_t res = -1;
         if (ntails > 50 && tlen > 3){  // TODO
            // BINARY SEARCH
            res = SearchTail_BS(T,data, pos, tlen*2, key, ntails); // bitwise operations 
         }
         else {
         // PROCEED AS BEFORE
            // 4. Look for substring where the tails of that length start 
            res = SearchTail(T,data, pos, tlen*2, key, ntails); // bitwise operations 
         }
         if (res!=-1) {
            // Add to the result the number of finimizers preceeding this. [DONE in rarest_fmin_streaming_search] 
            return {res+tails_so_far, tlen};
         }
      }
      pos+= (tlen*ntails*2);
      // 5. if fmin not found, add the tails seen so far
      tails_so_far += ntails;
   }
   return {-1,0};
}

