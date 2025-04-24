// To compile: (in Pseudoaligment-Finimizers) srun g++ nks/slam.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o slam -lz -Wno-deprecated-declarations

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

#include "common.hh"

using namespace std;

void printBinary(uint64_t v){ // prints in the reverse order
   for(int i=0;i<64;i++){
      cerr << ((v>>i)&1);
      //cerr << ((v>>(63-i))&1);
   }
   //cerr << '\n';
}

//__builtin_popcountll
//__builtin_clzll

#define hasless(x,n) (((x)-~0UL/255*(n))&~(x)&~0UL/255*128)

#define haszero(v) (((v) - 0x0101010101010101ULL) & ~(v) & 0x8080808080808080ULL)

#define hasvalue(x,n) \
(haszero((x) ^ (~0ULL/255 * (n))))

#define hasvaluesupply(x,mask) \
(haszero((x) ^ (mask)))
//(hasless((x) ^ (mask),1))

#define haszero2(v, W, mask2, mask3, Wmask) ((((v) - mask2) & ~(v) & mask3)&(~Wmask))
   //   (((v) - 0x1111111111111111ULL) & ~(v) & 0x8888888888888888ULL)


#define hasvaluesupply2(x,W,mask, mask2, mask3, Wmask) \
(haszero2((x ^ mask), W, mask2, mask3, Wmask))
//(hasless((x) ^ (mask),1))


constexpr uint64_t masks[10][3] = {
   {0,0,0}, //0th entry is not to be used!
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0b0001000100010001000100010001000100010001000100010001000100010001,0,0}, //4-bit patterns (16 of them)
   {0b000010000100001000010000100001000010000100001000010000100001,0,0}, //5-bit (12 of them)
   {0b000001000001000001000001000001000001000001000001000001000001,0,0}, //6 (10)
   {0b000000100000010000001000000100000010000001000000100000010000001,0,0}, //7 (9)
   {0b0000000100000001000000010000000100000001000000010000000100000001,0,0}, //8 (8)
   {0b000000001000000001000000001000000001000000001000000001000000001,0,0}, //9 (7)
};

constexpr uint64_t rmasks[10][3] = {
   {0,0,0}, //0th entry is not to be used!
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0b1000100010001000100010001000100010001000100010001000100010001000,0,0}, //4-bit patterns (16 of them)
   {0b100001000010000100001000010000100001000010000100001000010000,0,0}, //5-bit (12 of them)
   {0b100000100000100000100000100000100000100000100000100000100000,0,0}, //6 (10)
   {0b100000010000001000000100000010000001000000100000010000001000000,0,0}, //7 (9)
   {0b1000000010000000100000001000000010000000100000001000000010000000,0,0}, //8 (8)
   {0b100000000100000000100000000100000000100000000100000000100000000,0,0}, //9 (7)
};

void createMask(uint64_t key, int keyLen, int startPos, uint64_t &kmask, uint64_t &mask2, uint64_t &mask3) {
   kmask = 0;
   mask2 = 0;
   mask3 = 0;
   // all the patterns have the same length
   uint64_t pattern2 = 1ULL;                  // 1 preceded by KeyLen -1 0s
   uint64_t pattern3 = (1ULL << (keyLen - 1)); // 1 followed by KeyLen -1 0s

   int pos = startPos;

   while (pos + keyLen <= 64) {
       kmask |= (key << pos);
       mask2 |= (pattern2 << pos);
       mask3 |= (pattern3 << pos);
       pos += keyLen;
   }

   return;
}

tuple<vector<int>, vector<pair<int,int>>, vector<uint64_t>> count_itemsPerWord(int W, uint64_t key, uint64_t info, uint64_t &kmask, uint64_t &mask2, uint64_t &mask3){
   vector<uint64_t> Wmasks;

   // for every width != W, create a mask and combine & all the masks
   // not 1*width*bucketsize but it has to be at the right pos
   //uint64_t kmask, mask2, mask3;
   int smaller = 0;
   int larger = 0;
   int S = 0; // smaller items
   bool kmask_created = false;
   uint64_t *infop = &info;
   vector<int> itemsPerWord;
   vector<pair<int,int>> SLitems;
   vector<uint64_t> bucketSize; // could be uint32_t
   vector<uint64_t> lengths; 
   for(int i=0;i<8;i++){
      uint8_t x = *(((uint8_t*)(infop))+i);
      if ((uint64_t)x >0){
         if (i%2) {
            //cerr << "length = ";
            lengths.push_back((uint64_t)x);
            }
         else {
            //cerr << "number = ";
            bucketSize.push_back((uint64_t)x);
         }            
         cerr << (uint64_t)x << endl;
      }
   }
   // Check if there is more than one word
   //Not necessary now
   uint64_t i = 0; // How many buckets do we have? 
   int diffWitems=0;
   while (i < bucketSize.size()){
      uint64_t totl = 0;
      uint64_t new_totl = totl;
      uint64_t number = bucketSize[i];
      int items = 0; 
      int diffWItems=0; 
      int pos = 0; // pos at which the wrong width items start
      uint64_t mask = 0;
      while(i < bucketSize.size() and (new_totl + (lengths[i] * number)) <= 64){
         totl + (lengths[i] * number);
         if (lengths[i] != W){
            diffWItems += (int)number;
            if (lengths[i] < W){
               smaller+= diffWItems * lengths[i];
               S += diffWItems;

            }
            else {//if (lengths[i] > W){
               larger+= diffWItems * lengths[i];
            }
            //cerr << "pos = "<< pos << endl;
            //cerr << "diffWItems = "<< diffWItems * lengths[i] << endl;
            uint64_t mask_ = ((1ULL << (diffWItems* lengths[i])) - 1) << pos;
            //printBinary(mask_); cerr << " mask W = " << lengths[i] << endl;
            mask = mask | mask_;
            // reset for a new mask
            pos += diffWItems * lengths[i];
            
            diffWItems = 0;
         }else{
            // create the repated mask
            createMask(key, W, pos, kmask, mask2, mask3);
            //printBinary(kmask); cerr << "repeated mask" << endl;
            pos += lengths[i] * number; // the values before this pos contain ok items
         }         
      
         new_totl = totl;
         items += (int)number;
         i++;
         number = bucketSize[i];
      }
      new_totl = totl;
      while( i < bucketSize.size() and (new_totl + lengths[i]) <= 64){ 
         totl += lengths[i];
         new_totl = totl;
         items ++;
         if (lengths[i] != W){
            diffWItems ++;
            if (lengths[i] < W){
               smaller++;
               S++;
            } 
            else {//(lengths[i] > W){
               larger++;
            }
         } else{
            // create the repated mask if it has not been already created
            if (!kmask_created){ 
               createMask(key, W, pos, kmask, mask2, mask3);
               kmask_created = true;}
            pos += lengths[i]; // the values before this pos contain ok items
         }
         number --;
         if (number = 0){
            uint64_t mask_ = ((1ULL << (diffWItems* lengths[i])) - 1) << pos;
            //printBinary(mask_); cerr << " mask W = " << lengths[i] << endl;
            mask = mask | mask_;
            // reset for a new mask
            pos += diffWItems;
            diffWItems = 0;

            i++;
            number = bucketSize[i];
         }
      }
      itemsPerWord.push_back(items);
      SLitems.push_back({smaller,S});
      //cerr << "larger = " << larger << endl;
      //cerr << "smaller = " << smaller << endl;
      //cerr << "smaller items " << S<< endl;
      Wmasks.push_back(mask);
   }
   return {itemsPerWord,SLitems, Wmasks};
}

inline pair<int64_t,string> slam(const sdsl::int_vector<0> &T, int64_t start, char W, uint32_t key, uint16_t ntails){
   // input T, offset at which the true tails start, W(tlen), key, #tails 
   // TODO: bitwise operations
   // look at 64 bits at a time and mask what is not a tail of the correct size
   return {1,"1"};

}

uint64_t read_bits(const sdsl::int_vector<0>& vec, size_t offset, size_t bit_width){
    uint64_t value = 0;
    for (size_t i = 0; i < bit_width; ++i) {
        // Shift previous bits up by 1, then OR in the next bit
        value <<= 1;
        value |= vec[offset + i];
    }
    return value;
}


// Report the lenght of the found tail and the pos in t to get the colors
// output: pos in T (to get colors), tailsubstring
pair<int64_t, string> bitMagicSearch_new(const sdsl::int_vector<0> &T, int64_t start, string s){ // we know the width of the query
   // input: T, offset in T, string or substring after prefix
   // EVERY PREFIX HAS A DIFFERENT INT
   // T.size()= found prefixes THIS IS NOT TRUE!!
   
   
   // NEW METHOD FOR EVERY KEY: SLAM

   int64_t pos = start; // we know from where we need to look at the vector -> we start from here
   
   // Check first the shorter lengths. stop once a match is found
   while (true){ // break the loop once something is found
      // 1. check the length of the first tail, 5‐bit tlen
      char tlen = read_bits(T, pos, 5);
      pos += 5;

      // 2. check how many tailS, vbyte #tails
      uint64_t ntails = 0;
      int shift = 0;
      while (true) {
         uint8_t byte = read_bits(T, pos, 8);
         pos += 8;
         ntails |= uint64_t(byte & 0x7F) << shift;
         if ((byte & 0x80) == 0) break;
         shift += 7;
      }

      // 3. Extract substring
      uint32_t key = prefix2int(s, 0, tlen); // TODO EXTRACT SUFFIX OF LENGHT TLEN; // The max length is (k-plen)*2= 21 if k=31 and plen=10, we need at least that many bits

      // 4. Look for substring where the tails of that length start (WHERE??)
      pair<int64_t,string> res = slam (T, pos, tlen, key, ntails); // TODO real bitwise operations 
   
   }
   // TODO add to the result the number of finimizers preceeding this. [DONE in rarest_fmin_streaming_search]
   // option 1. store it in B = {prefix: {start, #prev tails}} [using this option now]
   // option 2. calculate it (might be very slow) 
   return {0,"0"};
}

inline int64_t bitMagicSearch2(uint64_t X, int W, uint64_t key, uint64_t info){ // we know the width of the query
   // We should knwo the number of tails with a given length
    
   cerr << W <<  " query width" << endl;
   cerr << key << " key" << endl;
   //printBinary(key); cerr << " key" << endl;

   //uint64_t mask = masks[W][0]*key; //~0ULL/255 * key; 
   //printBinary(mask); cerr << " mask" << endl;

   //uint64_t mask2 = masks[W][0]*1; // this works for different widths
   //uint64_t mask3 = rmasks[W][0]*1;

   uint64_t mask, mask2, mask3;

   // These cannot be set here as different words might have different number of items
   // Make sure that every work contains WHOLE items
   //uint64_t itemsPerWord = 64/W; 
   //uint64_t maxBitsPerWord = itemsPerWord*W; // this is not used

   tuple<vector<int>,vector<pair<int,int>>, vector<uint64_t>> results = count_itemsPerWord(W, key, info, mask, mask2, mask3); // one value for every word
   vector<int> itemsPerWord = get<0>(results);
   vector<pair<int,int>> SLitems = get<1>(results);
   
   vector<uint64_t> Wmasks = get<2>(results);
   //printBinary(mask); cerr << " mask" << endl;

   //uint64_t mask2 = get<3>(results);
   //uint64_t mask3 = get<4>(results);
   
   //uint64_t maskwW = maskWrongW(info);
   // TODO use another mask to mask the values of a different width !!!
   // based on info for each word
   // if lenght != width


   //uint64_t nbits = X.size()*W;
   //cerr << nbits << endl;
   //uint64_t fullWords = X.size()*W/64;
   //for(uint64_t i=0;i<fullWords;i++){
   //uint64_t dataSize = X.size()*W/64 + (X.size()*W%64 > 0); 
   uint64_t dataSize = 1;
   for(uint64_t i=0;i<dataSize;i++){

      int smaller = SLitems[i].first;
      int S = SLitems[i].second;

      //cerr << "smaller = " << smaller << endl;
      //cerr << "smaller items " << S<< endl;


      uint64_t w = X; //uint64_t w = X.data()[i];
      printBinary(w); cerr << " w " << endl;
      /* 
      printBinary(~Wmasks[i]); cerr << " mask wrong width" << endl;
      printBinary(mask); cerr << " mask" << endl;

      printBinary((w) ^ mask); cerr << " (w) ^ mask" << endl;
      printBinary(mask2); cerr << " mask2" << endl;
      //printBinary(0x1111111111111111ULL); cerr << " 0x1111111111111111ULL" << endl;

      printBinary(((w ) ^ mask) - mask2); cerr << " (w ^ mask)- (mask2)" << endl;
      printBinary(~(w ^ mask)); cerr << " ~(w ^ mask)" << endl;

      printBinary(((w ^ mask) - mask2) & ~(w  ^ mask)); cerr << "((w  ^ mask) - (mask2) & ~(w ^ mask)" << endl;
      printBinary(mask3); cerr << " mask3" << endl;

      printBinary(((w  ^ mask) - mask2 & ~(w  ^ mask) & mask3) & (~Wmasks[i])); cerr << " ((w ^ mask)- mask2) & ~(w  ^ mask) & mask3) & (~Wmasks[i])" << endl;
 */
      //uint64_t found = hasValueSupplyDebug(w,mask);//hasvaluesupply(w,mask); 
      uint64_t Wmask = Wmasks[i];
      uint64_t found = hasvaluesupply2(w,W,mask,mask2, mask3, Wmask);
      printBinary(found); cerr << " found" <<  endl;

      if(found){
         //printBinary(found); cerr << endl;

         uint64_t lz = __builtin_clzll(found);
         cerr << "lz = "<< lz<< endl;
         bool needsCorrection = (found>>(63-lz-W))&1;
         cerr << "needsCorrection: " << needsCorrection << '\n';
         uint64_t pos =  lz/W; //__builtin_clzll(found)/W;
         //cerr << key << ' ' << X[i*(itemsPerWord)+((itemsPerWord)-pos-1-needsCorrection)] << '\n';
         //for(int i=0;i<X.size();i++){
         //   cerr << X[i] << ' ';
         //}
         //cerr << '\n';
         //if(__builtin_popcountll(found) > 1){
         //   exit(1);
         //}
         //cerr << itemsPerWord<< endl;
         cerr << "pos " << pos << endl; // TODO FIX: pos is based on W as it depedens on mask 3 and this should be flipped 
     
         // items in previous words
         uint64_t result = 0;
         for (int j=0;j<i; j++){
            result += (itemsPerWord[j]);
         }
         // 64 - pos - C[W]=smaller values]
         /* cerr << "flip " << 64-lz<< endl;
         cerr << "remove smaller all " << (64-lz)- smaller << endl;
         cerr << "Divide by W " << ((64-lz)- smaller)/W << endl;
         cerr << "correct if needed " << (((64-lz)- smaller)/W) - needsCorrection << endl;
         cerr << "Add smaller items " << ((((64-lz)- smaller)/W) - needsCorrection) + S << endl;
         cerr << "Correct for 0 " << (((((64-lz)- smaller)/W) - needsCorrection) + S ) - 1 << endl; */
         cerr << "result = " << (((((64-lz)- smaller)/W) - needsCorrection) + S ) - 1 << endl;
         
         return result + (((((64-lz)- smaller)/W) - needsCorrection) + S ) - 1;
      }
   }
   return -1;
}

