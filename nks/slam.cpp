// To compile: (in Pseudoaligment-Finimizers) srun g++ nks/slam.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o slam -lz -Wno-deprecated-declarations


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

void repeatKey(uint64_t key, int keyLen, int startPos, uint64_t &kmask, uint64_t &mask2, uint64_t &mask3) {
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

pair<vector<int>, vector<uint64_t>> count_itemsPerWord(int W, uint64_t key, uint64_t info, uint64_t &kmask, uint64_t &mask2, uint64_t &mask3){
   vector<uint64_t> Wmasks;

   // for every width != W, create a mask and combine & all the masks
   // not 1*width*bucketsize but it has to be at the right pos
   //uint64_t kmask, mask2, mask3;
   bool kmask_created = false;
   uint64_t *infop = &info;
   vector<int> itemsPerWord;
   vector<uint64_t> bucketSize; // could be uint32_t
   vector<uint64_t> lengths; 
   for(int i=0;i<8;i++){
      uint8_t x = *(((uint8_t*)(infop))+i);
      if ((uint64_t)x >0){
         if (i%2) {
            cerr << "length = ";
            lengths.push_back((uint64_t)x);
            }
         else {
            cerr << "number = ";
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
            repeatKey(key, W, pos, kmask, mask2, mask3);
            printBinary(kmask); cerr << "repeated mask" << endl;
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
         } else{
            // create the repated mask if it has not been already created
            if (!kmask_created){ 
               repeatKey(key, W, pos, kmask, mask2, mask3);
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
      Wmasks.push_back(mask);
      //Now we knwon how many of 
   }
   return {itemsPerWord,Wmasks};
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



inline int64_t bitMagicSearch2(uint64_t X, int W, uint64_t key, uint64_t info){ // we know the width of the query
   // We should knwo the number of tails with a given length
    
   cerr << W <<  " query width" << endl;
   cerr << key << " key" << endl;
   printBinary(key); cerr << " key" << endl;

   //uint64_t mask = masks[W][0]*key; //~0ULL/255 * key; 
   //printBinary(mask); cerr << " mask" << endl;

   //uint64_t mask2 = masks[W][0]*1; // this works for different widths
   //uint64_t mask3 = rmasks[W][0]*1;

   uint64_t mask, mask2, mask3;

   // These cannot be set here as different words might have different number of items
   // Make sure that every work contains WHOLE items
   //uint64_t itemsPerWord = 64/W; 
   //uint64_t maxBitsPerWord = itemsPerWord*W; // this is not used

   pair<vector<int>, vector<uint64_t>> results = count_itemsPerWord(W, key, info, mask, mask2, mask3); // one value for every word
   vector<int> itemsPerWord = results.first;
   vector<uint64_t> Wmasks = results.second;
   printBinary(mask); cerr << " mask" << endl;

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


      uint64_t w = X; //uint64_t w = X.data()[i];
      printBinary(w); cerr << " w " << endl;
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
         //if (pos>)
         cerr << needsCorrection << endl;
         cerr << "old items per word = "<< 64/W << endl;
         cerr << "new items per word = "<< itemsPerWord[0]<< endl;
         cerr << "old result " << i*(64/W)+((64/W)-pos-1-needsCorrection)<< endl;

         cerr << "result " << i*(itemsPerWord[i])+((itemsPerWord[i])-pos-1-needsCorrection)<< endl;
         return i*(itemsPerWord[i])+((itemsPerWord[i])-pos-1-needsCorrection);
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

int main(int argc, char **argv) {
    if(argc != 4){
       cerr << "Usage: " << argv[0] << "{int width} {# of ints} {# of queries}\n";
       exit(1);
    }

    const uint64_t width = atoi(argv[1]);
    if(!width || width > 16){
       cerr << "Error: width should be less than 17 and greater than 0\n";
       exit(1);
    }

    uint64_t maxint = (1<<width)-1;
    uint64_t n = atoi(argv[2]);
    uint64_t nqueries = atoi(argv[3]);

    cerr << "width: " << width << '\n'; 
    cerr << "n: " << n << '\n'; 
    cerr << "nqueries: " << nqueries << '\n';

    uint64_t X1 = 0;
    uint64_t *X1p = &X1;
    X1 |= (14ULL<<56); printBinary(X1);
    X1 |= (15ULL<<48); printBinary(X1);
    X1 |= (16ULL<<40); printBinary(X1);
    X1 |= (17ULL<<32); printBinary(X1);
    X1 |= (18ULL<<24); printBinary(X1);
    X1 |= (19ULL<<16); printBinary(X1);
    X1 |= (20ULL<<8); printBinary(X1);
    X1 |= (21ULL<<0); printBinary(X1);
    cerr << "X1: " << X1 << '\n';

    for(int i=0;i<8;i++){
       uint8_t x = *(((uint8_t*)(X1p))+i);
       cerr << (uint64_t)x << ' '; 
    }
    cerr << '\n';
    /* cerr << "OLD"<< endl;

    bitMagicSearch2_old(X1,18);

    cerr << "NEW"<< endl;
   */
   //bitMagicSearch2(X1,8,18);


    uint64_t X2 = 0;
    uint64_t *X2p = &X2;
     X2 |= (21ULL<<56);
     X2 |= (20ULL<<48); 
     X2 |= (19ULL<<40); 
     X2 |= (18ULL<<32); 
     X2 |= (17ULL<<24); 
     X2 |= (16ULL<<16); 
     X2 |= (15ULL<<8); 
     X2 |= (14ULL<<0);
    cerr << "X2: " << X2 << '\n';

    for(int i=0;i<8;i++){
      uint8_t x = *(((uint8_t*)(X2p))+i);
      cerr << (uint64_t)x << ' '; 
   }
   cerr << '\n';

   cerr << '\n';
   /* cerr << "OLD"<< endl;
   bitMagicSearch2_old(X2,18);
   cerr << "NEW"<< endl;
   */
   //bitMagicSearch2(X2,8,18);
    

    cerr << endl;

    uint64_t X3 = 0;
    uint64_t *X3p = &X3;
     X3 |= (21ULL<<56);
     X3 |= (20ULL<<48); 
     X3 |= (19ULL<<40); 
     X3 |= (18ULL<<32); 
     X3 |= (16ULL<<24); 
     X3 |= (15ULL<<16); 
     X3 |= (14ULL<<8); 
     X3 |= (13ULL<<0);
    cerr << "X3: " << X3 << '\n';

    for(int i=0;i<8;i++){
      uint8_t x = *(((uint8_t*)(X3p))+i);
      cerr << (uint64_t)x << ' '; 
   }
   cerr << '\n';

    /* cerr << '\n';
    cerr << "OLD"<< endl;
    bitMagicSearch2_old(X3,18);
    
    cerr << "NEW"<< endl;
 */   //bitMagicSearch2(X3,8,18);

    cerr << endl;

    uint64_t X4 = 0;
    uint64_t *X4p = &X4;
    // width = 4 for the first 3, then 
     X4 |= (8ULL<<60); 
     X4 |= (9ULL<<56);
     X4 |= (10ULL<<52); 
     X4 |= (11ULL<<48);
     X4 |= (12ULL<<44);
     X4 |= (13ULL<<40); // 1101
     X4 |= (14ULL<<36);
     X4 |= (15ULL<<32);

     X4 |= (141ULL<<24); // 10001101
     X4 |= (221ULL<<16);  // 11011101
     //X4 |= (19ULL<<8); 
     //X4 |= (18ULL<<0); 

     X4 |= (3549ULL<<1); //110111011101
     
    cerr << "X4: " << X4 << '\n';

   uint64_t infoX4 = 0; // 8 bits will not be enough for some number of tails (max tails per bucket)
   uint64_t *infoX4p = &infoX4;
   infoX4 |= (4ULL<<56);
   infoX4 |= (8ULL<<48); 
   infoX4 |= (8ULL<<40); 
   infoX4 |= (2ULL<<32);
   infoX4 |= (15ULL<<24);
   infoX4 |= (1ULL<<16);
   infoX4 |= (1ULL<<8);
   infoX4 |= (1ULL<<0);



   /* for(int i=0;i<8;i++){
      uint8_t x = *(((uint8_t*)(infoX4p))+i);
      if (i%2) {cerr << "length = ";}
      else {cerr << "number = ";}
      cerr << (uint64_t)x << endl;
   } */

   uint64_t val = *X4p;

// 4-bit values
for (int i = 0; i < 8; ++i) {
   uint8_t nibble = (val >> (60 - i * 4)) & 0xF;
   cerr << (uint64_t)nibble << ' ';
}

// 8-bit values 
for (int i = 0; i < 2; ++i) {
   uint8_t byte = (val >> (24 - i * 8)) & 0xFF;
   cerr << (uint64_t)byte << ' ';
}
// 16-bit values 
for (int i = 0; i < 1; ++i) {
   uint16_t chunk = (val >> (16 * (0 - i))) & 0xFFFF;
   cerr << (uint64_t)chunk << ' ';
}

cerr << '\n';
   cerr << "X4 "<< endl;
   bitMagicSearch2(X4,4,13, infoX4);
   cerr << "correct output = 5"<< endl;
   cerr << endl;
   cerr << "X4 "<< endl;
   bitMagicSearch2(X4,4,15, infoX4);
   cerr << "correct output = 3"<< endl;
   cerr << endl;
   bitMagicSearch2(X4,4,14, infoX4);
   cerr << "correct output = 4"<< endl;
   cerr << endl;
   bitMagicSearch2(X4,8,18, infoX4);
   cerr << "correct output = x"<< endl;
   cerr << endl;
   /* cerr << "X4 "<< endl;
   bitMagicSearch2(X4,8,19, infoX4);
   cerr << "correct output = 1" << endl; // ok
   cerr << endl;  */
  /*  cerr << "X4 "<< endl;
   bitMagicSearch2(X4,8,221, infoX4);
   cerr << "correct output = 1" << endl;
   cerr << endl;  */
   cerr << "X4 "<< endl;
   bitMagicSearch2(X4,15,3549, infoX4);
   cerr << "correct output = 0" << endl;
   cerr << endl; 

   /*cerr << "X3 "<< endl;

   for(int i=0;i<8;i++){
      uint8_t x = *(((uint8_t*)(X3p))+i);
      cerr << (uint64_t)x << ' '; 
   }
   cerr << '\n';
   uint64_t infoX3 = 0; // 8 bits will not be enough for some number of tails (max tails per bucket)
   uint64_t *infoX3p = &infoX3;
   infoX3 |= (8ULL<<56);
   infoX3 |= (8ULL<<48);  
   bitMagicSearch2(X3,8,16,infoX3);
   cerr << "correct output = 3"<<endl;
   cerr << endl;
   
   cerr << "X3 "<< endl;
   for(int i=0;i<8;i++){
      uint8_t x = *(((uint8_t*)(X3p))+i);
      cerr << (uint64_t)x << ' '; 
   }
   cerr << '\n';
   bitMagicSearch2(X3,8,15,infoX3);
   cerr << "correct output = 2"<<endl;
   cerr << endl; */
   

    uint64_t X5 = 0;
    uint64_t *X5p = &X5;
   // width = 4 for the first 3, then 
    X5 |= (8ULL<<60); 
    X5 |= (9ULL<<56);
    X5 |= (10ULL<<52); 
    X5 |= (11ULL<<48);
    X5 |= (12ULL<<44);
    X5 |= (13ULL<<40); 
    X5 |= (14ULL<<36);
    X5 |= (15ULL<<32);
    X5 |= (7ULL<<28);

    X5 |= (6ULL<<24); 
    X5 |= (5ULL<<20);
    X5 |= (4ULL<<16);

    X5 |= (3ULL<<12);  
    X5 |= (2ULL<<8); 
    X5 |= (1ULL<<4);

    X5 |= (0ULL<<0); 

   /* cerr << "reverse" << endl;
    for (int i = 0; i < 16; i++) {
      uint64_t val4 = (X5 >> (60 - i * 4)) & 0xF;
      cerr << val4 << ' ';
  }
  cerr << '\n';

  uint64_t infoX5 = 0; // 8 bits will not be enough for some number of tails (max tails per bucket)
   uint64_t *infoX5p = &infoX5;
   infoX5 |= (4ULL<<56);
   infoX5 |= (16ULL<<48);  
   cerr << endl;
  cerr << "new X5 "<< endl;

    bitMagicSearch2(X5,4,13, infoX5);
    cerr << "correct output = 10"<< endl;
    cerr << endl;

    bitMagicSearch2(X5,4,14, infoX5);
    cerr << "correct output = 9"<< endl;
    cerr << endl;
 */

    sdsl::int_vector<> Y(8,0,8);
       Y[0] = 14;
       Y[1] = 15;
       Y[2] = 16;
       Y[3] = 17;
       Y[4] = 18;
       Y[5] = 19;
       Y[6] = 20;
       Y[7] = 21;
       /* cerr << Y.data()[0] << '\n';
       cerr << Y.get_int(0,8) << '\n';
       cerr << Y.get_int(8,8) << '\n';
       cerr << Y.get_int(16,8) << '\n';
       cerr << Y.get_int(24,8) << '\n';
       cerr << Y.get_int(32,8) << '\n';
       cerr << Y.get_int(40,8) << '\n';
       cerr << Y.get_int(48,8) << '\n';
       cerr << Y.get_int(56,8) << '\n';
    cerr << '\n'; */
    //bitMagicSearch(Y,int(Y.width()),18);

    sdsl::int_vector<> YY(8,0,8);
       Y[0] = 6;
       Y[1] = 7;
       Y[2] = 8;
       Y[3] = 9;
       Y[4] = 18;
       Y[5] = 19;
       Y[6] = 20;
       Y[7] = 21;
       /* cerr << Y.data()[0] << '\n';
       cerr << Y.get_int(0,8) << '\n';
       cerr << Y.get_int(8,8) << '\n';
       cerr << Y.get_int(16,8) << '\n';
       cerr << Y.get_int(24,8) << '\n';
       cerr << Y.get_int(32,8) << '\n';
       cerr << Y.get_int(40,8) << '\n';
       cerr << Y.get_int(48,8) << '\n';
       cerr << Y.get_int(56,8) << '\n';
    cerr << '\n'; */
    //bitMagicSearch(Y,int(Y.width()),18);


    sdsl::int_vector<> Z(16,0,4); // (n,x,l), with n equals size, x default integer value, l width of integer 
       Z[0] = 1; Z[1] = 7; Z[2] = 7; Z[3] = 3; Z[4] = 4; Z[5] = 5; Z[6] = 6; Z[7] = 7;
       Z[8] = 8; Z[9] = 9; Z[10] = 10; Z[11] = 11; Z[12] = 12; Z[13] = 13; Z[14] = 14; Z[15] = 15;
       //cerr << Z.data()[0] << '\n'; //empty
       //cerr << Z.get_int(0) << '\n';
       //cerr << int(Z.width()) << '\n';
      
    cerr << '\n';

    //uint64_t z = (uint64_t)-1;
    //cerr << z << '\n'; //2^64 - 1
 /*    int W = int(Z.width());
    bitMagicSearch(Z,W,uint64_t(15));
    bitMagicSearch(Z,W,uint64_t(8));
    bitMagicSearch(Z,W,uint64_t(4)); */
    
    cerr<< "done"<< endl;

    exit(1);    

    //generate a random permutation of [0..(2^width)-1], take the first n elements, sort them
    vector<uint64_t> V;
    for(int i=0;i<(1<<width);i++){ // this gives error
       V.push_back(i);
    }
    std::srand(unsigned(std::time(0)));
    std::random_shuffle(V.begin(),V.end());
    for(uint64_t i=0;i<n;i++){
       cerr << V[i] << ' ';
    }
    cerr << '\n';

    std::sort(V.begin(),V.begin()+n);
    for(uint64_t i=0;i<n;i++){
       cerr << V[i] << ' ';
    }
    cerr << '\n';

    //generate an SDSL int_vector of desired width containing the n elements above
    sdsl::int_vector<> X(n,0,width);
    cerr << "X.size() = " << X.size() << '\n';
    for(uint64_t i=0;i<n;i++){
       X[i] = V[i];
    }
    for(uint64_t i=0;i<n;i++){
       cerr << X[i] << ' ';
    }
    cerr << '\n';

    //generate queries, both mixed and positive

    vector<uint16_t> randomQueries;
    vector<uint16_t> positiveQueries;

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen1(rd()); // seed the generator
    std::uniform_int_distribution<> distr1(0, maxint); // define the range
    std::mt19937 gen2(rd()); // seed the generator
    std::uniform_int_distribution<> distr2(0, n-1); // define the range

    for(uint64_t i=0; i<nqueries; i++){
       uint64_t r1 = distr1(gen1);
       randomQueries.push_back(r1);
       //cerr << r1 << '\n';
       uint64_t r2 = distr2(gen2);
       positiveQueries.push_back(X[r2]);
    }

    cerr << randomQueries.size() << '\n';
    cerr << positiveQueries.size() << '\n';

    //bring it into cache
    int64_t checksum = 0;
    for(uint64_t i=0;i<n;i++){
       checksum += X[i];
    }
    cerr << checksum << '\n';

    cerr << "Searching by slow scan...\n";
    auto start = chrono::steady_clock::now();
    checksum = 0;
    for(size_t i=0;i<randomQueries.size();i++){
       int64_t res = -1;
       for(uint64_t j=0;j<X.size();j++){
          if(X[j] == randomQueries[i]){
             res = (int64_t)j;
             break;
          }
       }
       checksum += res;
    }
    auto end = chrono::steady_clock::now();
    cout << "Time/query in nanoseconds: "
        << chrono::duration_cast<chrono::nanoseconds>(end - start).count()/nqueries
        << " ns" << endl;
    cerr << checksum << '\n';

    cerr << "Searching by binary search...\n";
    start = chrono::steady_clock::now();
    checksum = 0;
    for(size_t i=0;i<randomQueries.size();i++){
       int64_t res = binSearch(X,randomQueries[i]);
       checksum += res;
    }
    end = chrono::steady_clock::now();
    cout << "Time/query in nanoseconds: "
        << chrono::duration_cast<chrono::nanoseconds>(end - start).count()/nqueries
        << " ns" << endl;
    cerr << checksum << '\n';

    cerr << "Searching by bit magic search...\n";
    start = chrono::steady_clock::now();
    checksum = 0;
    for(size_t i=0;i<randomQueries.size();i++){
       int64_t res = bitMagicSearch(X,X.width(),randomQueries[i]);
       checksum += res;
    }
    end = chrono::steady_clock::now();
    cout << "Time/query in nanoseconds: "
        << chrono::duration_cast<chrono::nanoseconds>(end - start).count()/nqueries
        << " ns" << endl;
    cerr << checksum << '\n';
}

