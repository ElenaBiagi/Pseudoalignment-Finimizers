#ifndef _PRED8_V2_H_
#define _PRED8_V2_H_

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

class Pred8v2 {
   public:
    Pred8v2() {
    }
    ~Pred8v2() = default;
    Pred8v2(const vector<uint64_t> &data) {
       _n = data.size();
       _min = data[0];
       _u = data[data.size()-1]; // - _min;
       _nblocks = _u/256 + ((_u%256) > 1);

       _X.resize(_nblocks+1); //+1 for a useful dummy at the end
       for(uint64_t i=0;i<_nblocks;i++) _X[i] = 0;
       std::vector<uint8_t> _C(_nblocks, 0);

       //NB: TODO currently assumes the 0th bucket is non-empty

       _nActiveBuckets = 0;
       for(uint64_t i=0;i<_n;i++){
          uint64_t v = data[i]; // - _min;
          //cerr << "data["<<i<<"]: "<< data[i] << '\n';
          //getc(stdin); 
          //if(v == 757315312){
          //   cerr << "\tHere it is, at i="<<i<<'\n';
          //   cerr << '\t' << (v>>8) << ' ' << _X[v>>8] << '\n';
          //}
          //if(v>>8 == 2958262){
          //   cerr << '\t' << (v>>8) << ' ' << _X[v>>8] << '\n';
          //}
          if(!_X[v>>8]) _nActiveBuckets++;
          _X[v>>8]++;
       }

       _Y.resize(_n);

       cerr << "Pred8v2: _u _n _min _nblocks _nActiveBuckets: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<' '<<_nActiveBuckets<<'\n';
       cerr << "Pred8v2: sizeInBytes(): "<<sizeInBytes()<<'\n';
       uint64_t yi = 0;
       for(uint64_t i=0;i<_n;){
          uint64_t v = data[i]; // - _min;
          uint64_t bcount = _X[v >> 8];
          _X[v >> 8] = (yi << 2) | (bcount > 0); //LSB==1 indicates bucket is non-empty, Note (bcount > 1) is always true here
          if(bcount){
             _C[v >> 8] = (uint8_t)(bcount-1);
             for(uint64_t j = 0; j < bcount; j++){
                v = data[i]; // - _min;
                _Y[yi] = v & 255;
                if((v >> 8) ==41931742 ){
                   cerr << "_Y["<<yi<<"]: "<< (uint64_t)_Y[yi] << '\n';
                   cerr << "data["<<i<<"]: " << data[i] << '\n';
                } 
                yi++;
                i++;
             }
          }//else{
          //   cerr << "Should never happen: "<<i<<'\n';
          //}
       }
       uint32_t lastNonEmptyX = _X[0] & 0xFFFFFFFE;
       uint32_t lastNonEmptyBucket = 0;
       uint32_t lastNonEmptyC = _C[0];
       for(uint64_t i=1;i<_nblocks;i++){
          if(!_X[i]){
             //Bucket i is an empty bucket
             if(lastNonEmptyBucket == i-1){
                //the previous bucket is non-empty
                _X[i] = (((lastNonEmptyX>>2) + lastNonEmptyC)<<2); //we want to point to it's last element
             }else{
                //store the bucket number of the previous non-empty bucket
                _X[i] = (lastNonEmptyBucket << 2) | 2; //two LSBs are 1 and 0 to indicate, resp. 1 last non-empty bucket is stored, and 0 this bucket is empty
             }
             //_X[i] = (((lastNonEmptyX>>1) + lastNonEmptyC)<<1); //we want to point to it's last element
          }else{
             lastNonEmptyX = _X[i] & 0xFFFFFFFE;
             lastNonEmptyC = _C[i];
             lastNonEmptyBucket = i;
          }
       }
       _X[_nblocks] = (((_X[_nblocks-1] >> 2) + _C[_nblocks-1])<<2);
       cerr << "_X[_nblocks]: " << _X[_nblocks] << '\n';
    }

    //  FIXME: as of 3 March 2026 this needs to be updated to not use _min, etc.
    //  p is the index of the predecessor in the set
    //  bool is 1 if the value at index p is equal to key, else 0
    pair<int64_t, bool> inline getPred(int64_t key) const{
        //cerr << "key: " << key << '\n';
        if (key < _min) {
            return {-1, false};
        }
        key = key - _min;
        uint32_t x = _X[key>>8];
        //cerr << "x: " << x << '\n';
        if(x&1){
           //non-empty bucket
           uint64_t y = x >> 1; //this is where, in _Y, the bucket elements are located
           uint64_t bcount = ((_X[1+(key>>8)]>>1) + (1 - (_X[1+(key>>8)]&1))) - y;//_C[key>>8] + 1;
           //uint64_t bcount = _C[key>>8] + 1;
           //cerr << "bcount: "<<bcount<<'\n';
           //cerr << "_C[key>>8]+1: "<<(_C[key>>8] + 1)<<'\n';
           uint64_t k = key & 255;
           for(uint64_t j=0; j < bcount; j++){
              if(_Y[y+j] >= k){
                 return {y+j-(1-(_Y[y+j] == k)),(_Y[y+j] == k)};
              }
           }
           //getting here means that we did not find anything in the bucket that was >= k
           //this means the last element of the bucket is the predecessor of k, and it is not equal to k
           return {y+bcount-1,false};
        }
        //the bucket that key belongs to is empty
        return {(x>>1), false};
    }

    int64_t rank(int64_t pos) const{
       pair<int64_t, bool> r = getPred(pos);
       return (r.first+1)-((uint64_t)(r.second));
    }

    //  first is the index of the predecessor in the set
    //  second is the length of the prefix (in bits) that key has with its predecessor
    pair<uint64_t, uint64_t> inline getPredPrefix(int64_t key) const{
        _n_searches++;
        //cerr << "key: " << key << '\n';
        if(key < _min) {
            //cerr << "returning from 1\n";
            return {-1, 0};
        }
        if(key > _u) {
           //TODO: double check the lcp is being computed correctly in this case
           uint64_t lux = countl_zero(~((_u) ^ (~key))); 
           return {_n-1,lux};
        }
        //cerr << "key: " << key << '\n';
        uint32_t x = _X[key>>8];
        //cerr << "x>>2: " << (x>>2) << '\n';
        //cerr << "key>>8: " << (key>>8) << "_nblocks: " << _nblocks << '\n';
        if(x&1){
           //non-empty bucket
           _n_hard_searches++;
           uint64_t y = x >> 2; //this is where, in _Y, the bucket elements are located
           //cerr << "Next bucket is empty? " << (1 - (_X[1+(key>>8)]&1)) << '\n';
           uint64_t bcount = ((_X[1+(key>>8)]>>2) + (1 - (_X[1+(key>>8)]&1))) - y;//_C[key>>8] + 1;
           //uint64_t bcount = _C[key>>8] + 1;
           //cerr << "bcount: "<<bcount<<'\n';
           //cerr << "_C[key>>8]+1: "<<(_C[key>>8] + 1)<<'\n';
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
                    //cerr << "returning from 2\n";
                    return {y+j-(1-(_Y[y+j] == k)),lcp};
                 }
                 //getting here means the key is smaller than everything in its bucket
                 //our predecessor will be the last element of the previous non-empty bucket
                 //TODO: triple check this case
                 uint64_t bk = key>>8;
                 uint64_t bp = bk-1; //assume this is an empty bucket preceded by a non-empty bucket
                 if(_X[bk] & 2){ //is that really the case?
                    bp = _X[bk]>>2; //no, so point to the bucket containing the predecessor
                 }
                 uint64_t lcp = countl_zero(~((bp<<8) ^ (~(bk<<8))));
                 //cerr << "returning from 2+\n";
                 return {(_X[bp+1]>>2), lcp};
              }
           }
           //getting here means that we did not find anything in the bucket that was >= k
           //this means the last element of the bucket is the predecessor of k, and it is not equal to k
           //cerr << "returning from 3\n";
           uint64_t pred = _Y[y+bcount-1];
           uint64_t lcp = countl_zero(~((pred) ^ (~k))); //lcp in bits, 64-bit oriented
           //cerr << "lcp: "<<lcp<<'\n';
           return {y+bcount-1,lcp};
        }
        //the bucket that key belongs to is empty
        _n_easy_searches++;
        uint64_t bk = key>>8;
        uint64_t bp = bk-1; //assume this is an empty bucket preceded by a non-empty bucket
        if(_X[bk] & 2){ //is that really the case?
           bp = _X[bk]>>2; //no, so point to the bucket containing the predecessor
        }
        uint64_t lcp = countl_zero(~((bp<<8) ^ (~(bk<<8))));
        //cerr << "returning from 4\n";
        return {(_X[bp+1]>>2), lcp};
    }
 

    size_t getu() const { return _u; }
    size_t getn() const { return _n; }

    uint64_t sizeInBytes() const{
        uint64_t sz = 5*sizeof(uint64_t) + (sizeof(uint32_t)*(_nblocks+1)) + _n;
        return sz;
    }

    int64_t serialize(std::ostream& os) const{
       cerr << "Pred8.serialize()...\n";
        int64_t written = 0;
        os.write((char *)&_u, sizeof(uint64_t));
        os.write((char *)&_n, sizeof(uint64_t));
        os.write((char *)&_min, sizeof(uint64_t));
        os.write((char *)&_nblocks, sizeof(uint64_t));
        os.write((char *)_X.data(),sizeof(uint32_t)*(_nblocks+1));
       cerr << "_u _n _min _nblocks _nActiveBuckets: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<' '<<_nActiveBuckets<<'\n';
        //os.write((char *)_C,sizeof(uint8_t)*_nblocks);
        os.write((char *)&_nActiveBuckets, sizeof(uint64_t));
        os.write((char *)_Y.data(),sizeof(uint8_t)*_n);
        written += 4*sizeof(uint64_t) + (sizeof(uint32_t)*_nblocks) + _nblocks + _n;
        return written;
    }

    void load(std::istream& is){
       cerr << "Pred8.load()...\n";
       is.read((char *)&_u, sizeof(uint64_t));
       is.read((char *)&_n, sizeof(uint64_t));
       is.read((char *)&_min, sizeof(uint64_t));
       is.read((char *)&_nblocks, sizeof(uint64_t));
       cerr << "_u _n _min _nblocks _nActiveBuckets: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<' '<<_nActiveBuckets<<'\n';
       _X.resize(_nblocks+1);
       is.read((char *)_X.data(), sizeof(uint32_t)*(_nblocks+1));
       //_C = new uint8_t[_nblocks];
       //is.read((char *)_C, sizeof(uint8_t)*_nblocks);
       is.read((char *)&_nActiveBuckets, sizeof(uint64_t));
       _Y.resize(_n);
       is.read((char *)_Y.data(), sizeof(uint8_t)*_n);
       //_X[_nblocks] = (((_X[_nblocks-1] >> 1) + _C[_nblocks-1])<<1);
    }
    
   //  Pred8v2(Pred8v2 &other){
   //     this->_u = other._u;
   //     this->_n = other._n;
   //     this->_min = other._min;
   //     this->_nblocks = other._nblocks;
   //     this->_nActiveBuckets = other._nActiveBuckets;
   //     this->_X = _X;
   //     //this->_C = _C;
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
    std::vector<uint32_t> _X;
    std::vector<uint8_t> _Y;
};

#endif
