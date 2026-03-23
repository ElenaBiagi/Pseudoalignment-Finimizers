#ifndef _PRED8_V4_H_
#define _PRED8_V4_H_

//A version of Pred8 that uses a bitvector to indicate empty/non-empty buckets and 
//uses rank on that bitvector to access bucket information (i.e. offset).

#include <stdio.h>
#include <stdlib.h>
#include <bit>
#include <bitset>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <ratio>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

#include "sdsl/bit_vectors.hpp"

class Pred8v4 {
   public:
    Pred8v4() {
    }
    ~Pred8v4() = default;
    Pred8v4(const vector<uint64_t> &data) {
       _n = data.size();
       _min = data[0];
       _u = data[data.size()-1];
       _nblocks = _u/256 + ((_u%256) > 1);

       cerr << "Pred8v4: _u _n _min _nblocks: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<'\n';
       _X_bv = sdsl::int_vector<1>(_nblocks+1);
       std::vector<uint32_t> X(_nblocks, 0);

       _nActiveBuckets = 0;
       for(uint64_t i=0;i<_n;i++){
          uint64_t v = data[i];
          if(!X[v>>8]) _nActiveBuckets++;
          X[v>>8]++;
          _X_bv[v>>8] = 1;
       }
       
       _X.resize(_nActiveBuckets);
       for(uint64_t i=0;i<_nActiveBuckets;i++) _X[i] = 0;

       _Y.resize(_n);

       cerr << "Pred8v4: _u _n _min _nblocks _nActiveBuckets: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<' '<<_nActiveBuckets<<'\n';
       cerr << "Pred8v4: sizeInBytes(): "<<sizeInBytes()<<'\n';
       uint64_t yi = 0;
       uint64_t xi = 0;
       for(uint64_t i=0;i<_n;){
          uint64_t v = data[i];
          uint64_t bcount = X[v >> 8];
          if(bcount){
             //in _X[xi] we will store for the x^th non-empty bucket: 1) starting pos in _Y; 2) bucket number; 3) # of elements
             _X[xi] = (yi << ((uint64_t)36)) | ((v>>(uint64_t)8) << (uint64_t)8) | (bcount-1);
             xi++;
             for(uint64_t j = 0; j < bcount; j++){
                v = data[i];
                _Y[yi] = v & 255;
                yi++;
                i++;
             }
          }//else{
          //   cerr << "Should never happen: "<<i<<'\n';
          //}
       }

       sdsl::util::init_support(this->_X_rs, &(this->_X_bv));

    }

    //  first is the index of the predecessor in the set
    //  second is the length of the prefix (in bits) that key has with its predecessor
    pair<uint64_t, uint64_t> inline getPredPrefix(int64_t key) const{
        _n_searches++;
        //cerr << "key: " << key << '\n';
        if (key < _min) {
            //cerr << "returning from 1\n";
            return {-1, 0};
        }
        cerr << "key: " << key << '\n';
        cerr << "key>>8: " << (key>>8) << '\n';
        uint32_t ex = _X_bv[key>>8];
        uint64_t px = _X_rs.rank(key>>8);
        cerr << "ex: " << (ex) << '\n';
        cerr << "px: " << (px) << '\n';
        //cerr << "key>>8: " << (key>>8) << "_nblocks: " << _nblocks << '\n';
        if(ex){
           //non-empty bucket
           _n_hard_searches++;
           uint64_t y = _X[px] >> 36; //this is where, in _Y, the bucket elements are located
           cerr << "y: "<<y<<'\n';
           //cerr << "Next bucket is empty? " << (1 - (_X[1+(key>>8)]&1)) << '\n';
           uint64_t bcount = (_X[px] & 255) + 1;
           cerr << "bcount: "<<bcount<<'\n';
           uint64_t k = key & 255;
           for(uint64_t j=0; j < bcount; j++){
              //cerr << "j: " << j << " _Y[y+j]: " << (uint64_t)_Y[y+j] << "k: " << k << '\n';
              if(_Y[y+j] >= k){
                 if(j > 0 || _Y[y+j] == k){
                    //cerr << "breaker: " << (uint64_t)_Y[y+j] << '\n';
                    uint64_t pred = _Y[y+j-(1-(_Y[y+j] == k))];
                    //cerr << "pred: " << pred << '\n'; 
                    //cerr << "k:    " << k << '\n';  
                    uint64_t lcp = countl_zero(~((pred) ^ (~k))); //lcp in bits, 64-bit oriented
                    //cerr << "lcp: "<<lcp<<'\n';
                    cerr << "returning from 2\n";
                    return {y+j-(1-(_Y[y+j] == k)),lcp};
                 }
                 //getting here means the key is smaller than everything in its bucket
                 //our predecessor will be the last element of the previous non-empty bucket
                 //Note: we wouldn't be here at all if key < _min, so there is a previous bucket that is non-empty
                 uint64_t bk = key>>8; //the bucket number of the key
                 uint64_t bp = (_X[px-1] << 28)>>36; //the bucket number of the previous non-empty bucket
                 cerr << "2+ bp bk:"<<bp<<' '<<bk<<'\n';
                 uint64_t lcp = countl_zero(~((bp<<8) ^ (~(bk<<8))));
                 cerr << "returning from 2+\n";
                 return {y-1, lcp};
              }
           }
           //getting here means that we did not find anything in the bucket that was >= k
           //this means the last element of the bucket is the predecessor of k, and it is not equal to k
           cerr << "returning from 3\n";
           uint64_t pred = _Y[y+bcount-1];
           uint64_t lcp = countl_zero(~((pred) ^ (~k))); //lcp in bits, 64-bit oriented
           //cerr << "lcp: "<<lcp<<'\n';
           return {y+bcount-1,lcp};
        }
        //the bucket that key belongs to is empty
        _n_easy_searches++;
        uint64_t bk = key>>8;
        uint64_t bp = (_X[px-1] << 28)>>36; //the bucket number of the previous non-empty bucket
        cerr << "4 bp bk:"<<bp<<' '<<bk<<'\n';
        uint64_t lcp = countl_zero(~((bp<<8) ^ (~(bk<<8))));
        cerr << "returning from 4\n";
        uint64_t y = _X[px] >> 36; //this is where, in _Y, the bucket elements are located
        uint64_t bcount = (_X[px] & 255) + 1;
        return {y+bcount-1, lcp};
    }
 

    size_t getu() const { return _u; }
    size_t getn() const { return _n; }

    uint64_t sizeInBytes() const{
        uint64_t sz = 5*sizeof(uint64_t) + (sizeof(uint32_t)*(_nblocks+1)) + _n;
        //TODO: add on the size of the bitvector and the rank support structure
        //sz += _X_bv.;
        //sz += _X_rs.;
        return sz;
    }

    int64_t serialize(std::ostream& os) const{
        cerr << "Pred8v4.serialize()...\n";
        int64_t written = 0;
        os.write((char *)&_u, sizeof(uint64_t));
        os.write((char *)&_n, sizeof(uint64_t));
        os.write((char *)&_min, sizeof(uint64_t));
        os.write((char *)&_nblocks, sizeof(uint64_t));
        os.write((char *)&_nActiveBuckets, sizeof(uint64_t));
        cerr << "_u _n _min _nblocks _nActiveBuckets: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<' '<<_nActiveBuckets<<'\n';
        os.write((char *)_X.data(),sizeof(uint64_t)*(_nActiveBuckets));
        os.write((char *)_Y.data(),sizeof(uint8_t)*_n);
        written += 4*sizeof(uint64_t) + (sizeof(uint32_t)*_nblocks) + _nblocks + _n;

        written += _X_bv.serialize(os);
        written += _X_rs.serialize(os);
        cerr << "Pred8.serialize() wrote " << written << " bytes.\n"; 
        return written;
    }

    void load(std::istream& is){
       cerr << "Pred8v4.load()...\n";
       is.read((char *)&_u, sizeof(uint64_t));
       is.read((char *)&_n, sizeof(uint64_t));
       is.read((char *)&_min, sizeof(uint64_t));
       is.read((char *)&_nblocks, sizeof(uint64_t));
       is.read((char *)&_nActiveBuckets, sizeof(uint64_t));
       cerr << "_u _n _min _nblocks _nActiveBuckets: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<' '<<_nActiveBuckets<<'\n';
       _X.resize(_nActiveBuckets);
       is.read((char *)_X.data(), sizeof(uint64_t)*(_nActiveBuckets));
       _Y.resize(_n);
       is.read((char *)_Y.data(), sizeof(uint8_t)*_n);

       _X_bv.load(is);
       _X_rs.load(is,&_X_bv);
    }
    
   //  Pred8v4(Pred8v4 &other){
   //     this->_u = other._u;
   //     this->_n = other._n;
   //     this->_min = other._min;
   //     this->_nblocks = other._nblocks;
   //     this->_nActiveBuckets = other._nActiveBuckets;
   //     this->_X = _X;
   //     this->_Y = _Y;
   //  }

    //stats
    mutable uint64_t _n_searches = 0;  
    mutable uint64_t _n_easy_searches = 0;  
    mutable uint64_t _n_hard_searches = 0;  

   private:

    uint64_t _u = 0;  // universe size
    uint64_t _n = 0;  // number of elements
    uint64_t _min = 0;  // value of the smallest element
    uint64_t _nblocks = 0;
    uint64_t _nActiveBuckets = 0;
    std::vector<uint64_t> _X;
    std::vector<uint8_t> _Y;
    sdsl::int_vector<1> _X_bv;
    sdsl::rank_support_v5<1> _X_rs;
};

#endif
