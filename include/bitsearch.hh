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


inline uint64_t hasValueSupplyDebug(uint64_t x, uint64_t mask){
   printBinary(x);
   printBinary(mask);
   printBinary(x);
   printBinary(0x0101010101010101ULL);
   printBinary(x - 0x0101010101010101ULL);
   printBinary(~x);
   printBinary((x - 0x0101010101010101ULL)&~x);
   printBinary(((x - 0x0101010101010101ULL)&~x)&0x8080808080808080ULL);
   return (((x - 0x0101010101010101ULL)&~x)&0x8080808080808080ULL);
}

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

inline int64_t slam(const sdsl::int_vector<0> &T, int64_t start, char W, uint32_t key, uint16_t ntails){
   // input T, offset at which the true tails start, W(tlen), key, #tails 
   // TODO: bitwise operations
   // look at 64 bits at a time and mask what is not a tail of the correct size
   return 1;

}

uint64_t read_bits(const sdsl::int_vector<0>& vec,
                   size_t offset,
                   size_t bit_width){
    uint64_t value = 0;
    for (size_t i = 0; i < bit_width; ++i) {
        // Shift previous bits up by 1, then OR in the next bit
        value <<= 1;
        value |= vec[offset + i];
    }
    return value;
}

//TODO the following 2 are defined in common.hh
/* char get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}
uint32_t prefix2int(const std::string& s, uint64_t offset, char plen){
    uint32_t h = 0;
    for (uint64_t i = 0; i < (uint64_t)plen; i++) {
        uint32_t b = get_char_idx(s[offset + i]);
        h <<= 2;
        h |= b;
    }
    return h;
} */

// Report the lenght of the found tail and the pos in t to get the colors
pair<int64_t, char> bitMagicSearch_new(const sdsl::int_vector<0> &T, int64_t start, string s){ // we know the width of the query
   // input: T, offset in T, string or substring after prefix
   // EVERY PREFIX HAS A DIFFERENT INT
   // T.size()= found prefixes
   // we know from where we need to look at the vector -> we start from here
   // 5 bits for the length
   // vbyte for the number of tails 
   // Check first the shorter lengths. stop once a match is found

   // 1. check the length of the first tail
   // 2. check how many tails
   // 3. Extract substring
   // 4. Look for substring where the tails of that length start (WHERE??)
   
   // NEW METHOD FOR EVERY KEY: SLAM

   int64_t pos = start;
   
   while (true){ // break the loop once something is found
      // 1. 5‐bit tlen
      char tlen = read_bits(T, pos, 5);
      pos += 5;

      // 2. vbyte #tails
      uint64_t ntails = 0;
      int shift = 0;
      while (true) {
         uint8_t byte = read_bits(T, pos, 8);
         pos += 8;
         ntails |= uint64_t(byte & 0x7F) << shift;
         if ((byte & 0x80) == 0) break;
         shift += 7;
      }
      uint32_t key = prefix2int(s, 0, tlen); // TODO EXTRACT SUFFIX OF LENGHT TLEN; // The max length is (k-plen)*2= 21 if k=31 and plen=10, we need at least that many bits

      int64_t found = slam (T, pos, tlen, key, ntails); // TODO real bitwise operations 
   
   }
     // TODO add to the result the number of finimizers preceeding this.
   // option 1. store it in B = {prefix: {start, #prev tails}}
   // option 2. calculate it (might be very slow) 
   return {0,0};
   ////////////////// ///////////////////////// /////////////////////////// /////////////////////////////
   // 3. TODO actual tails
/*    for (uint64_t i = 0; i < ntails; ++i) {
      uint64_t tail = read_bits(T, pos, tlen * 2);
      pos += tlen * 2;
      // 
   } */
   
   /// TODO FIX THIS
   uint64_t info; // this is the first thing you encounter at offset
   int W;
   uint64_t key; // as to be etxracted from S depending on W

   cerr << W <<  " query width" << endl;
   cerr << key << " key" << endl;

   uint64_t mask, mask2, mask3;

   // These cannot be set here as different words might have different number of items
   // Make sure that every word contains WHOLE items
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


   //uint64_t nbits = T.size()*W;
   //cerr << nbits << endl;
   //uint64_t fullWords = T.size()*W/64;
   //for(uint64_t i=0;i<fullWords;i++){
   //uint64_t dataSize = T.size()*W/64 + (T.size()*W%64 > 0); 
   uint64_t dataSize = 1;
   for(uint64_t i=0;i<dataSize;i++){

      int smaller = SLitems[i].first;
      int S = SLitems[i].second;

      //cerr << "smaller = " << smaller << endl;
      //cerr << "smaller items " << S<< endl;


      uint64_t w;// = T; //uint64_t w = T.data()[i];
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
         //cerr << key << ' ' << T[i*(itemsPerWord)+((itemsPerWord)-pos-1-needsCorrection)] << '\n';
         //for(int i=0;i<T.size();i++){
         //   cerr << T[i] << ' ';
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
         
         return {-1,-1};//result + (((((64-lz)- smaller)/W) - needsCorrection) + S ) - 1;
      }
   }
   return {-1,-1};
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

inline int64_t bitMagicSearch2_old(uint64_t X, uint64_t key){
   uint64_t width = 8;
   uint64_t mask = masks[width][0]*key; //~0ULL/255 * key;
   uint64_t itemsPerWord = 64/width;
   uint64_t maxBitsPerWord = itemsPerWord*width;
   cerr << maxBitsPerWord << endl;
   //uint64_t fullWords = X.size()*X.width()/64;
   //for(uint64_t i=0;i<fullWords;i++){
   uint64_t dataSize = 1;
   for(uint64_t i=0;i<dataSize;i++){
      uint64_t w = X;
         printBinary(w); cerr << " w" << endl;
         printBinary(mask); cerr << " mask" << endl;
      //uint64_t found = hasValueSupplyDebug(w,mask);//hasvaluesupply(w,mask);
      uint64_t found = hasvaluesupply(w,mask);
      if(found){
         printBinary(found); cerr << " found" << endl;
         printBinary(w); cerr << " w" << endl;
         //cerr << w << '\n';
         printBinary(mask); cerr << " mask" << endl;
         //cerr << mask << '\n';
         printBinary(w ^ mask); cerr << " w ^ mask" << endl;
         printBinary(0x0101010101010101ULL); cerr << " 0x0101010101010101ULL" << endl;
         printBinary((w ^ mask)-0x0101010101010101ULL); cerr << " (w ^ mask)-0x0101010101010101ULL" << endl;
	      //printBinary((w ^ mask)-0X8080808080808080ULL); cerr << " (w ^ mask)-0x8080808080808080ULL" << endl;
         //printBinary( ~(w ^ mask)); cerr << " ~(w ^ mask)" << endl;
	      //printBinary(((w ^ mask)-0x8080808080808080ULL) & ~(w ^ mask)); cerr << "((w ^ mask)-8080808080808080ULL) & ~(w ^ mask)"<< endl;
         printBinary(((w ^ mask)-0x0101010101010101ULL) & ~(w ^ mask)); cerr << "((w ^ mask)-0x0101010101010101ULL) & ~(w ^ mask)" << endl;
	      //printBinary(((w ^ mask)-0x0101010101010101ULL) & ~(w ^ mask) & 0x0101010101010101ULL); cerr << endl;
         //printBinary(0x0101010101010101ULL); cerr << " 0x0101010101010101ULL" << endl;
         printBinary(0x8080808080808080ULL); cerr << " 0x8080808080808080ULL" << endl;
         //printBinary(((w ^ mask)-0x8080808080808080ULL) & ~(w ^ mask) & 0x0101010101010101ULL);
         printBinary(((w ^ mask)-0x0101010101010101ULL) & ~(w ^ mask) & 0x8080808080808080ULL); cerr << " ((w ^ mask)-0x0101010101010101ULL) & ~(w ^ mask) & 0x8080808080808080ULL)" << endl;

         uint64_t lz = __builtin_clzll(found);
         bool needsCorrection = (found>>(63-lz-width))&1;
         //cerr << "needsCorrection: " << needsCorrection << '\n';
         uint64_t pos =  lz/width; //__builtin_clzll(found)/X.width();
         //cerr << key << ' ' << X[i*(itemsPerWord)+((itemsPerWord)-pos-1-needsCorrection)] << '\n';
         //for(int i=0;i<X.size();i++){
         //   cerr << X[i] << ' ';
         //}
         //cerr << '\n';
         //if(__builtin_popcountll(found) > 1){
         //   exit(1);
         //}
         //cerr << itemsPerWord<< endl;
         cerr << pos << endl;
         cerr << needsCorrection << endl;
         cerr << "result = "<< i*(itemsPerWord)+((itemsPerWord)-pos-1-needsCorrection)<< endl;
         return i*(itemsPerWord)+((itemsPerWord)-pos-1-needsCorrection);
      }
   }
   return -1;
}


//TODO: Currently fails if the sought byte is an even number x and
//the number x+1 is also present. Is this a "feature" of the method
//and can it be corrected somehow?
//30.8.2023, 23:00, Update: Have corrected for the above case, now 
//most tests pass, but there still seems to be some weird behaviour
//--- probably another oddity I'm not understanding.
// correction needed also with x odd and x-1 (even)

// TODO this is now wrong but BitMagicSearch2 is correct
inline int64_t bitMagicSearch(sdsl::int_vector<> &X, int W, uint64_t key){ // we know the width of the query
   //int W = int(X.width());
   cerr << W <<  " width" << endl;
   cerr << int(X.width()) <<  " real width" << endl;

   uint64_t mask = masks[W][0]*key; //~0ULL/255 * key;
   cerr << key << " key" << endl;
   printBinary(key); cerr << " key" << endl;
   printBinary(mask); cerr << " mask" << endl;
   uint64_t mask2 = masks[W][0]*1; // this works for different widths
   uint64_t mask3 = rmasks[W][0]*1;
   uint64_t itemsPerWord = 64/W;
   uint64_t maxBitsPerWord = itemsPerWord*W;
   uint64_t nbits = X.size()*W;
   //cerr << nbits << endl;
   //uint64_t fullWords = X.size()*W/64;
   //for(uint64_t i=0;i<fullWords;i++){
   uint64_t dataSize = X.size()*W/64 + (X.size()*W%64 > 0); 
   for(uint64_t i=0;i<dataSize;i++){
      uint64_t w = X.data()[i];
      printBinary(w); cerr << " w " << endl;
      printBinary(w ^ mask); cerr << " w ^ mask" << endl;
      //printBinary(0x0101010101010101ULL); cerr << " 0x0101010101010101ULL" << endl;
      printBinary(mask2); cerr << " mask2" << endl;
      printBinary(0x1111111111111111ULL); cerr << " 0x1111111111111111ULL" << endl;

      printBinary((w ^ mask) - (0x1111111111111111ULL)); cerr << " (w ^ mask)- (mask2)" << endl;
      printBinary(~(w ^ mask)); cerr << " ~(w ^ mask)" << endl;

      printBinary(((w ^ mask) - (0x1111111111111111ULL)) & ~(w ^ mask)); cerr << "((w ^ mask) - (mask2) & ~(w ^ mask)" << endl;
      printBinary(0x8888888888888888ULL); cerr << " 0x8888888888888888ULL" << endl;

      printBinary(((w ^ mask) - (0x1111111111111111ULL) & ~(w ^ mask) & 0x8888888888888888ULL)); cerr << " ((w ^ mask)-0x1111111111111111ULL) & ~(w ^ mask) & 0x8888888888888888ULL)" << endl;


      //uint64_t found = hasValueSupplyDebug(w,mask);//hasvaluesupply(w,mask);
      uint64_t found = hasvaluesupply2(w,W,mask,mask2, mask3, mask); // TODO this is now wrong
      printBinary(found); cerr << " found" <<  endl;

      if(found){
         //printBinary(found); cerr << endl;
         //printBinary(w); cerr << endl;
         //cerr << w << '\n';
         //printBinary(mask);
         //cerr << mask << '\n';
         //printBinary(w ^ mask);
         //printBinary((w ^ mask)-0x0101010101010101ULL);
         //printBinary(((w ^ mask)-0x0101010101010101ULL) & ~(w ^ mask));
         //printBinary(0x8080808080808080ULL);
         //printBinary(((w ^ mask)-0x0101010101010101ULL) & ~(w ^ mask) & 0x8080808080808080ULL);

         uint64_t lz = __builtin_clzll(found);
         bool needsCorrection = (found>>(63-lz-W))&1;
         //cerr << "needsCorrection: " << needsCorrection << '\n';
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
         cerr << pos << endl;
         cerr << needsCorrection << endl;
         cerr << "result " << i*(itemsPerWord)+((itemsPerWord)-pos-1-needsCorrection)<< endl;
         return i*(itemsPerWord)+((itemsPerWord)-pos-1-needsCorrection);
      }
   }
   return -1;
}

inline int64_t binSearch(sdsl::int_vector<> &X, uint64_t key){
   int64_t mid = 0;
   int64_t low = 0;
   int64_t hi = X.size()-1;
   while (low <= hi) {
       mid = (low + hi) >> 1;
       uint64_t d = X[mid];
       //cerr << d << ' ' << key << '\n';
       if (d == key){
         //cerr << "Found\n";
         return mid;
       }else if (d > key){
         hi = mid - 1;
       }else{
         // This gets the insertion point right on the last loop.
         low = ++mid;
       }
   }
   //cerr << "Not found\n";
   return -1; //-mid - 1;
}

// /* int main(int argc, char **argv) {
//     if(argc != 4){
//        cerr << "Usage: " << argv[0] << "{int width} {# of ints} {# of queries}\n";
//        exit(1);
//     }

//     const uint64_t width = atoi(argv[1]);
//     if(!width || width > 16){
//        cerr << "Error: width should be less than 17 and greater than 0\n";
//        exit(1);
//     }

//     uint64_t maxint = (1<<width)-1;
//     uint64_t n = atoi(argv[2]);
//     uint64_t nqueries = atoi(argv[3]);

//     cerr << "width: " << width << '\n'; 
//     cerr << "n: " << n << '\n'; 
//     cerr << "nqueries: " << nqueries << '\n';

//     uint64_t X1 = 0;
//     uint64_t *X1p = &X1;
//     X1 |= (14ULL<<56); printBinary(X1);
//     X1 |= (15ULL<<48); printBinary(X1);
//     X1 |= (16ULL<<40); printBinary(X1);
//     X1 |= (17ULL<<32); printBinary(X1);
//     X1 |= (18ULL<<24); printBinary(X1);
//     X1 |= (19ULL<<16); printBinary(X1);
//     X1 |= (20ULL<<8); printBinary(X1);
//     X1 |= (21ULL<<0); printBinary(X1);
//     cerr << "X1: " << X1 << '\n';

//     for(int i=0;i<8;i++){
//        uint8_t x = *(((uint8_t*)(X1p))+i);
//        cerr << (uint64_t)x << ' '; 
//     }
//     cerr << '\n';
//     /* cerr << "OLD"<< endl;

//     bitMagicSearch2_old(X1,18);

//     cerr << "NEW"<< endl;
//    */
//    //bitMagicSearch2(X1,8,18);


//     uint64_t X2 = 0;
//     uint64_t *X2p = &X2;
//      X2 |= (21ULL<<56);
//      X2 |= (20ULL<<48); 
//      X2 |= (19ULL<<40); 
//      X2 |= (18ULL<<32); 
//      X2 |= (17ULL<<24); 
//      X2 |= (16ULL<<16); 
//      X2 |= (15ULL<<8); 
//      X2 |= (14ULL<<0);
//     cerr << "X2: " << X2 << '\n';

//     for(int i=0;i<8;i++){
//       uint8_t x = *(((uint8_t*)(X2p))+i);
//       cerr << (uint64_t)x << ' '; 
//    }
//    cerr << '\n';

//    cerr << '\n';
//    /* cerr << "OLD"<< endl;
//    bitMagicSearch2_old(X2,18);
//    cerr << "NEW"<< endl;
//    */
//    //bitMagicSearch2(X2,8,18);
    

//     cerr << endl;

//     uint64_t X3 = 0;
//     uint64_t *X3p = &X3;
//      X3 |= (21ULL<<56);
//      X3 |= (20ULL<<48); 
//      X3 |= (19ULL<<40); 
//      X3 |= (18ULL<<32); 
//      X3 |= (16ULL<<24); 
//      X3 |= (15ULL<<16); 
//      X3 |= (14ULL<<8); 
//      X3 |= (13ULL<<0);
//     cerr << "X3: " << X3 << '\n';

//     for(int i=0;i<8;i++){
//       uint8_t x = *(((uint8_t*)(X3p))+i);
//       cerr << (uint64_t)x << ' '; 
//    }
//    cerr << '\n';

//     /* cerr << '\n';
//     cerr << "OLD"<< endl;
//     bitMagicSearch2_old(X3,18);
    
//     cerr << "NEW"<< endl;
//  */   //bitMagicSearch2(X3,8,18);

//     cerr << endl;

//     uint64_t X4 = 0;
//     uint64_t *X4p = &X4;
//     // width = 4 for the first 3, then 
//     X4 |= (1ULL<<63); 
//      X4 |= (8ULL<<59); 
//      X4 |= (9ULL<<55);
//      X4 |= (10ULL<<51); 
//      X4 |= (11ULL<<47);
//      X4 |= (12ULL<<43);
//      X4 |= (13ULL<<39); // 1101
//      X4 |= (14ULL<<35);
//      X4 |= (15ULL<<31);

//      X4 |= (141ULL<<23); // 10001101
//      X4 |= (221ULL<<15);  // 11011101
//      //X4 |= (19ULL<<8); 
//      //X4 |= (18ULL<<0); 

//      X4 |= (3549ULL<<0); //110111011101
     
//     cerr << "X4: " << X4 << '\n';

//    uint64_t infoX4 = 0; // 8 bits will not be enough for some number of tails (max tails per bucket)
//    uint64_t *infoX4p = &infoX4;
// /*    infoX4 |= (4ULL<<56);
//    infoX4 |= (8ULL<<48); 
//    infoX4 |= (8ULL<<40); 
//    infoX4 |= (2ULL<<32);
//    infoX4 |= (15ULL<<24);
//    infoX4 |= (1ULL<<16);
//    infoX4 |= (1ULL<<8);
//    infoX4 |= (1ULL<<0); */

//    infoX4 |= (1ULL<<56);
//    infoX4 |= (1ULL<<48); 
//    infoX4 |= (4ULL<<40); //lenght
//    infoX4 |= (8ULL<<32);
//    infoX4 |= (8ULL<<24);
//    infoX4 |= (2ULL<<16);
//    infoX4 |= (15ULL<<8);
//    infoX4 |= (1ULL<<0);


//    /* for(int i=0;i<8;i++){
//       uint8_t x = *(((uint8_t*)(infoX4p))+i);
//       if (i%2) {cerr << "length = ";}
//       else {cerr << "number = ";}
//       cerr << (uint64_t)x << endl;
//    } */

// uint64_t val = *X4p;

//  // 4-bit values
// for (int i = 0; i < 8; ++i) {
//    uint8_t nibble = (val >> (59 - i * 4)) & 0xF;
//    cerr << (uint64_t)nibble << ' ';
// }

// // 8-bit values 
// for (int i = 0; i < 2; ++i) {
//    uint8_t byte = (val >> (15 + i * 8)) & 0xFF;
//    cerr << (uint64_t)byte << ' ';
// }
// // 15-bit values 
// for (int i = 0; i < 1; ++i) {
//    uint16_t chunk = (val >> (15 * (0 - i))) & 0x7FFF;
//    cerr << (uint64_t)chunk << ' ';
// }

// cerr << '\n';
//    cerr << "X4 "<< endl;
//    bitMagicSearch2(X4,4,13, infoX4);
//    cerr << "correct output = 5"<< endl;
//    cerr << endl;
//    cerr << "X4 "<< endl;
//    bitMagicSearch2(X4,4,15, infoX4);
//    cerr << "correct output = 3"<< endl;
//    cerr << endl;
//    bitMagicSearch2(X4,4,14, infoX4);
//    cerr << "correct output = 4"<< endl;
//    cerr << endl;
//    bitMagicSearch2(X4,8,18, infoX4);
//    cerr << "correct output = x"<< endl;
//    cerr << endl;
//    /* cerr << "X4 "<< endl;
//    bitMagicSearch2(X4,8,19, infoX4);
//    cerr << "correct output = 1" << endl; // ok
//    cerr << endl;  */
// /*    cerr << "X4 "<< endl;
//    bitMagicSearch2(X4,8,221, infoX4);
//    cerr << "correct output = 1" << endl;
//    cerr << endl;
//    cerr << "X4 "<< endl;
//    bitMagicSearch2(X4,15,3549, infoX4);
//    cerr << "correct output = 0" << endl;
//    cerr << endl;  */


//    uint64_t X6 = 0;
//    uint64_t *X6p = &X6;
//    // width = 4 for the first 3, then 
//    X6 |= (3549ULL<<49); //110111011101
   
//    X6 |= (141ULL<<41); // 10001101
//    X6 |= (221ULL<<33);  // 11011101

//     X6 |= (15ULL<<29); 
//     X6 |= (14ULL<<25); 
//     X6 |= (13ULL<<21); // 1101
//     X6 |= (12ULL<<17); 
//     X6 |= (11ULL<<13);
//     X6 |= (10ULL<<9);
//     X6 |= (9ULL<<5); 
//     X6 |= (8ULL<<1);
    
//     X6 |= (1ULL<<0);

    
//     //X6 |= (19ULL<<8); 
//     //X6 |= (18ULL<<0); 

    
//    cerr << "X6: " << X6 << '\n';

//   uint64_t infoX6 = 0; // 8 bits will not be enough for some number of tails (max tails per bucket)
//   uint64_t *infoX6p = &infoX6;
//   /* infoX6 |= (4ULL<<56);
//   infoX6 |= (8ULL<<48); 
//   infoX6 |= (8ULL<<40); 
//   infoX6 |= (2ULL<<32);
//   infoX6 |= (15ULL<<24);
//   infoX6 |= (1ULL<<16);
//   infoX6 |= (1ULL<<8);
//   infoX6 |= (1ULL<<0); */

//   infoX6 |= (15ULL<<56);
//   infoX6 |= (1ULL<<48); 
//   infoX6 |= (8ULL<<40); //lenght
//   infoX6 |= (2ULL<<32);
//   infoX6 |= (4ULL<<24);
//   infoX6 |= (8ULL<<16);
//   infoX6 |= (1ULL<<8);
//   infoX6 |= (1ULL<<0);


// cerr << 1 << " ";
//   // 4-bit values
//  for (int i = 0; i < 8; ++i) {
//     uint8_t nibble = (*X6p >> (29 - i * 4)) & 0xF;
//     cerr << (uint64_t)nibble << ' ';
//  }
 
//  // 8-bit values 
//  for (int i = 0; i < 2; ++i) {
//     uint8_t byte = (*X6p >> (41 + i * 8)) & 0xFF;
//     cerr << (uint64_t)byte << ' ';
//  }
//  // 15-bit values 
//  for (int i = 0; i < 1; ++i) {
//    uint16_t chunk = (*X6p >> 49) & 0x7FFF;
//    cerr << (uint64_t)chunk << ' ';
// }
//  cerr << endl;

//  cerr << "X6 "<< endl;
//  bitMagicSearch2(X6,4,13, infoX6);
//  cerr << "correct output = 6"<< endl;
//  cerr << endl;
//  cerr << "X4 "<< endl;
//  bitMagicSearch2(X6,4,15, infoX6);
//  cerr << "correct output = 8"<< endl;
//  cerr << endl;
//  bitMagicSearch2(X6,4,14, infoX6);
//  cerr << "correct output = 7"<< endl;
//  cerr << endl;
//  bitMagicSearch2(X6,8,18, infoX6);
//  cerr << "correct output = x"<< endl;
//  cerr << endl;
//  /* cerr << "X6 "<< endl;
//  bitMagicSearch2(X6,8,19, infoX6);
//  cerr << "correct output = 1" << endl; // ok
//  cerr << endl;  */
//  cerr << "X6 "<< endl;
//  bitMagicSearch2(X6,8,221, infoX6);
//  cerr << "correct output = 10" << endl;
//  cerr << endl;
//  cerr << "X6 "<< endl;
//  bitMagicSearch2(X6,15,3549, infoX6);
//  cerr << "correct output = 11" << endl;
//  cerr << endl; 


// /* 
//    cerr << "X3 "<< endl;

//    for(int i=0;i<8;i++){
//       uint8_t x = *(((uint8_t*)(X3p))+i);
//       cerr << (uint64_t)x << ' '; 
//    }
//    cerr << '\n';
//    uint64_t infoX3 = 0; // 8 bits will not be enough for some number of tails (max tails per bucket)
//    uint64_t *infoX3p = &infoX3;
//    infoX3 |= (8ULL<<56);
//    infoX3 |= (8ULL<<48);  
//    bitMagicSearch2(X3,8,16,infoX3);
//    cerr << "correct output = 3"<<endl;
//    cerr << endl;
   
//    cerr << "X3 "<< endl;
//    for(int i=0;i<8;i++){
//       uint8_t x = *(((uint8_t*)(X3p))+i);
//       cerr << (uint64_t)x << ' '; 
//    }
//    cerr << '\n';
//    bitMagicSearch2(X3,8,15,infoX3);
//    cerr << "correct output = 2"<<endl;
//    cerr << endl; */
   

//     uint64_t X5 = 0;
//     uint64_t *X5p = &X5;
//    // width = 4 for the first 3, then 
//     X5 |= (8ULL<<60); 
//     X5 |= (9ULL<<56);
//     X5 |= (10ULL<<52); 
//     X5 |= (11ULL<<48);
//     X5 |= (12ULL<<44);
//     X5 |= (13ULL<<40); 
//     X5 |= (14ULL<<36);
//     X5 |= (15ULL<<32);
//     X5 |= (7ULL<<28);

//     X5 |= (6ULL<<24); 
//     X5 |= (5ULL<<20);
//     X5 |= (4ULL<<16);

//     X5 |= (3ULL<<12);  
//     X5 |= (2ULL<<8); 
//     X5 |= (1ULL<<4);

//     X5 |= (0ULL<<0); 

//    /* cerr << "reverse" << endl;
//     for (int i = 0; i < 16; i++) {
//       uint64_t val4 = (X5 >> (60 - i * 4)) & 0xF;
//       cerr << val4 << ' ';
//   }
//   cerr << '\n';

//   uint64_t infoX5 = 0; // 8 bits will not be enough for some number of tails (max tails per bucket)
//    uint64_t *infoX5p = &infoX5;
//    infoX5 |= (4ULL<<56);
//    infoX5 |= (16ULL<<48);  
//    cerr << endl;
//   cerr << "new X5 "<< endl;

//     bitMagicSearch2(X5,4,13, infoX5);
//     cerr << "correct output = 10"<< endl;
//     cerr << endl;

//     bitMagicSearch2(X5,4,14, infoX5);
//     cerr << "correct output = 9"<< endl;
//     cerr << endl;
//  */

//     sdsl::int_vector<> Y(8,0,8);
//        Y[0] = 14;
//        Y[1] = 15;
//        Y[2] = 16;
//        Y[3] = 17;
//        Y[4] = 18;
//        Y[5] = 19;
//        Y[6] = 20;
//        Y[7] = 21;
//        /* cerr << Y.data()[0] << '\n';
//        cerr << Y.get_int(0,8) << '\n';
//        cerr << Y.get_int(8,8) << '\n';
//        cerr << Y.get_int(16,8) << '\n';
//        cerr << Y.get_int(24,8) << '\n';
//        cerr << Y.get_int(32,8) << '\n';
//        cerr << Y.get_int(40,8) << '\n';
//        cerr << Y.get_int(48,8) << '\n';
//        cerr << Y.get_int(56,8) << '\n';
//     cerr << '\n'; */
//     //bitMagicSearch(Y,int(Y.width()),18);

//     sdsl::int_vector<> YY(8,0,8);
//        Y[0] = 6;
//        Y[1] = 7;
//        Y[2] = 8;
//        Y[3] = 9;
//        Y[4] = 18;
//        Y[5] = 19;
//        Y[6] = 20;
//        Y[7] = 21;
//        /* cerr << Y.data()[0] << '\n';
//        cerr << Y.get_int(0,8) << '\n';
//        cerr << Y.get_int(8,8) << '\n';
//        cerr << Y.get_int(16,8) << '\n';
//        cerr << Y.get_int(24,8) << '\n';
//        cerr << Y.get_int(32,8) << '\n';
//        cerr << Y.get_int(40,8) << '\n';
//        cerr << Y.get_int(48,8) << '\n';
//        cerr << Y.get_int(56,8) << '\n';
//     cerr << '\n'; */
//     //bitMagicSearch(Y,int(Y.width()),18);


//     sdsl::int_vector<> Z(16,0,4); // (n,x,l), with n equals size, x default integer value, l width of integer 
//        Z[0] = 1; Z[1] = 7; Z[2] = 7; Z[3] = 3; Z[4] = 4; Z[5] = 5; Z[6] = 6; Z[7] = 7;
//        Z[8] = 8; Z[9] = 9; Z[10] = 10; Z[11] = 11; Z[12] = 12; Z[13] = 13; Z[14] = 14; Z[15] = 15;
//        //cerr << Z.data()[0] << '\n'; //empty
//        //cerr << Z.get_int(0) << '\n';
//        //cerr << int(Z.width()) << '\n';
      
//     cerr << '\n';

//     //uint64_t z = (uint64_t)-1;
//     //cerr << z << '\n'; //2^64 - 1
//  /*    int W = int(Z.width());
//     bitMagicSearch(Z,W,uint64_t(15));
//     bitMagicSearch(Z,W,uint64_t(8));
//     bitMagicSearch(Z,W,uint64_t(4)); */
    
//     cerr<< "done"<< endl;

//     exit(1);    

//     //generate a random permutation of [0..(2^width)-1], take the first n elements, sort them
//     vector<uint64_t> V;
//     for(int i=0;i<(1<<width);i++){ // this gives error
//        V.push_back(i);
//     }
//     std::srand(unsigned(std::time(0)));
//     std::random_shuffle(V.begin(),V.end());
//     for(uint64_t i=0;i<n;i++){
//        cerr << V[i] << ' ';
//     }
//     cerr << '\n';

//     std::sort(V.begin(),V.begin()+n);
//     for(uint64_t i=0;i<n;i++){
//        cerr << V[i] << ' ';
//     }
//     cerr << '\n';

//     //generate an SDSL int_vector of desired width containing the n elements above
//     sdsl::int_vector<> X(n,0,width);
//     cerr << "X.size() = " << X.size() << '\n';
//     for(uint64_t i=0;i<n;i++){
//        X[i] = V[i];
//     }
//     for(uint64_t i=0;i<n;i++){
//        cerr << X[i] << ' ';
//     }
//     cerr << '\n';

//     //generate queries, both mixed and positive

//     vector<uint16_t> randomQueries;
//     vector<uint16_t> positiveQueries;

//     std::random_device rd; // obtain a random number from hardware
//     std::mt19937 gen1(rd()); // seed the generator
//     std::uniform_int_distribution<> distr1(0, maxint); // define the range
//     std::mt19937 gen2(rd()); // seed the generator
//     std::uniform_int_distribution<> distr2(0, n-1); // define the range

//     for(uint64_t i=0; i<nqueries; i++){
//        uint64_t r1 = distr1(gen1);
//        randomQueries.push_back(r1);
//        //cerr << r1 << '\n';
//        uint64_t r2 = distr2(gen2);
//        positiveQueries.push_back(X[r2]);
//     }

//     cerr << randomQueries.size() << '\n';
//     cerr << positiveQueries.size() << '\n';

//     //bring it into cache
//     int64_t checksum = 0;
//     for(uint64_t i=0;i<n;i++){
//        checksum += X[i];
//     }
//     cerr << checksum << '\n';

//     cerr << "Searching by slow scan...\n";
//     auto start = chrono::steady_clock::now();
//     checksum = 0;
//     for(size_t i=0;i<randomQueries.size();i++){
//        int64_t res = -1;
//        for(uint64_t j=0;j<X.size();j++){
//           if(X[j] == randomQueries[i]){bitMagicSearch2
//              res = (int64_t)j;
//              break;
//           }
//        }
//        checksum += res;
//     }
//     auto end = chrono::steady_clock::now();
//     cout << "Time/query in nanoseconds: "
//         << chrono::duration_cast<chrono::nanoseconds>(end - start).count()/nqueries
//         << " ns" << endl;
//     cerr << checksum << '\n';

//     cerr << "Searching by binary search...\n";
//     start = chrono::steady_clock::now();
//     checksum = 0;
//     for(size_t i=0;i<randomQueries.size();i++){
//        int64_t res = binSearch(X,randomQueries[i]);
//        checksum += res;
//     }
//     end = chrono::steady_clock::now();
//     cout << "Time/query in nanoseconds: "
//         << chrono::duration_cast<chrono::nanoseconds>(end - start).count()/nqueries
//         << " ns" << endl;
//     cerr << checksum << '\n';

//     cerr << "Searching by bit magic search...\n";
//     start = chrono::steady_clock::now();
//     checksum = 0;
//     for(size_t i=0;i<randomQueries.size();i++){
//        int64_t res = bitMagicSearch(X,X.width(),randomQueries[i]);
//        checksum += res;
//     }
//     end = chrono::steady_clock::now();
//     cout << "Time/query in nanoseconds: "
//         << chrono::duration_cast<chrono::nanoseconds>(end - start).count()/nqueries
//         << " ns" << endl;
//     cerr << checksum << '\n';
// } */

