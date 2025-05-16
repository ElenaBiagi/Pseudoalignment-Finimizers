
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

uint64_t read_unaligned_64bits(const uint64_t* data, size_t offset_bits) {
    size_t word_index = offset_bits / 64;
    size_t bit_offset = offset_bits % 64;

    uint64_t low = data[word_index];
    uint64_t high = data[word_index + 1]; // TODO: check as this is safe only if enough padding is added

    // 128 bits
    __uint128_t wide = ((__uint128_t)high << 64) | low;
    wide >>= bit_offset;

    return static_cast<uint64_t>(wide & 0xFFFFFFFFFFFFFFFFULL);
}

inline int64_t slam(const sdsl::int_vector<1> &T, const uint64_t* data, const int64_t offset, const uint8_t W, const uint64_t key, const uint16_t ntails){
   // input T, offset at which the true tails start, W(tlen), key, #tails 

   // Look at 64 bits at a time starting from offset (skip tlen and ntail (5+8+?))
   // Look at tlen*ntail*2 bits
   // Mask all the bits after that = what is not a tail of the correct size

   auto masks = createMask(key, W);
   uint64_t mask  = std::get<0>(masks);
   uint64_t mask2 = std::get<1>(masks);
   uint64_t mask3 = std::get<2>(masks);

   uint64_t tails_per_word = std::min<uint64_t>(64 / W, ntails);
   uint64_t total_groups = (ntails + tails_per_word - 1) / tails_per_word;
   
   uint64_t j = 0;
   for (uint64_t i = 0; i < total_groups; ++i) {
      size_t bit_offset = offset + i * tails_per_word * W; // size_t bit_offset = offset + i * 64;// size_t bit_offset = offset + (i - word_index) * 64; // size_t bit_offset = offset + i * W; //


      uint64_t w = read_unaligned_64bits(data, bit_offset);
      //printBinary(w); cerr << " w at bit_offset=" << bit_offset << endl;      

      uint64_t tails_in_this_group = std::min(tails_per_word, ntails - i * tails_per_word);

      uint64_t bits_used = tails_in_this_group * W;
      uint64_t Wmask = (bits_used < 64) ? (~0ULL << bits_used) : 0;
      //printBinary(Wmask); cerr << " = Wmask " << endl;

      uint64_t found = hasvaluesupply(w,W,mask,mask2, mask3, Wmask);
      
      if(found){
         uint64_t lz = __builtin_clzll(found);
         bool needsCorrection = (found>>(63-lz-W))&1;
         int64_t result = (j*(tails_per_word))+((64-lz)/W)-1-needsCorrection;

         /* cerr << "width = "<< (int)W << endl;
         cerr << "key = "<< key << endl;
         cerr << "ntails = " << ntails << endl;
         printBinary(w); cerr << " w" << endl;
         printBinary(mask); cerr << " mask for the key" << endl;

         printBinary((w) ^ mask); cerr << " (w) ^ mask" << endl;
         printBinary(mask2); cerr << " mask2" << endl;

         printBinary(((w ) ^ mask) - mask2); cerr << " (w ^ mask)- (mask2)" << endl;
         printBinary(~(w ^ mask)); cerr << " ~(w ^ mask)" << endl;

         printBinary(((w ^ mask) - mask2) & ~(w  ^ mask)); cerr << "((w  ^ mask) - (mask2) & ~(w ^ mask)" << endl;
         printBinary(mask3); cerr << " mask3" << endl;
         printBinary(((w  ^ mask) - mask2 & ~(w  ^ mask) & mask3)); cerr << " ((w ^ mask)- mask2) & ~(w  ^ mask) & mask3)" << endl;

         printBinary(Wmask); cerr << " mask wrong width" << endl;

         printBinary(((w  ^ mask) - mask2 & ~(w  ^ mask) & mask3) & (~Wmask)); cerr << " ((w ^ mask)- mask2) & ~(w  ^ mask) & mask3) & (Wmask)" << endl;

         printBinary(found); cerr << " found" << endl;
         
         cerr << "j = "<< j << endl;
         cerr << "lz = "<< lz<< endl;

         
         if (needsCorrection) cerr << "needsCorrection: " << needsCorrection << '\n';
                  
         cerr << "ok result = "<< result << endl;

         cerr << endl; */
         return result;
      }
      j++; 
   // 1. I'm only looking at words that start at 0 -> no need for shifting masks [OK]
   // 2. the word starts at 0 so no smaller tails -> no need to mask smaller characters [OK]
   // 3. it is easy to know where longer tails start -> We need to MASK LONGER TAILS [OK]

   }

   return -1;

}


// output: pos in T (to get colors), tlen
pair<int64_t, uint8_t> bitMagicSearch(const sdsl::int_vector<1> &T, uint64_t s, char slen){ // we know the width of the query
   // input: T, offset in T, string or substring after prefix
   // EVERY PREFIX HAS A DIFFERENT INT
   // T.size()= found prefixes THIS IS NOT TRUE!!
   
   int64_t pos = 0;// start from 0 now that we have a single vector // we know from where we need to look at the vector -> we start from here
   //BitReader br(T, pos);

   int64_t tails_so_far = 0;

   uint64_t word_index = 0;
   uint8_t w_offset = 0;
   bool firstTail = true;
   cerr << "T.size(): "<< T.size();
   const uint64_t* data = T.data();
   // Check first the shorter lengths. stop once a match is found
   while (pos < T.size()){ // break the loop once something is found
      // 1. check the length of the first tail, 5‐bit tlen
      //uint8_t char = read_bits(T, data, pos, 5);
      //pos += 5;

      word_index = pos/64;
      w_offset = pos %64;
      //print_bit_vector(T, pos);
      uint8_t tlen = (uint8_t)sdsl::bits::read_int(&data[word_index], w_offset, 5);
      //cerr << "tlen: " << (int)tlen << endl; 
      if (tlen == 0){
         if (firstTail){return {0,0};}
         return {-1,0};
      }
      firstTail = false; // be sure to report ) only if it's the first tail length read

      pos += 5;
      //cerr << "tlen = "<< (int)tlen << endl;

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
      //cerr << "ntails: "<< ntails << endl;
      
      // 3. Extract substring
      // TODO Satrting at pos (64 - slen*2) extract the first 2*tlen bits of s
      //uint64_t key = (s >> (64 - (64-(s_len*2)) - t_len*2)) & ((1ULL << (t_len*2)) - 1);
      uint64_t key = (s >> ((slen - tlen) * 2)) & ((1ULL << (tlen * 2)) - 1);

      //uint64_t key = prefix2int(s, 0, tlen); // TODO EXTRACT SUFFIX OF LENGTH TLEN; // The max length is (k-plen)*2= 42 if k=31 and plen=10, we need at least that many bits

      // 4. Look for substring where the tails of that length start (WHERE??)
      int64_t res = slam (T,data, pos, tlen*2, key, ntails); // TODO real bitwise operations 
      if (res!=-1) {return {res+tails_so_far, tlen};}
      pos+= (tlen*ntails*2);
      // 5. if fmin not found, add the tails seen so far
      tails_so_far += ntails;
      //cerr << "tails_so_far = " << tails_so_far << endl;
   }
   // TODO add to the result the number of finimizers preceeding this. [DONE in rarest_fmin_streaming_search] 
   return {-1,0};
}


