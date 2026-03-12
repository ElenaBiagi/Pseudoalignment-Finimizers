#ifndef _PREF_TAB_H_
#define _PREF_TAB_H_

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <ratio>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

#include "Util.hh"


//#define BIN_SEARCH_THRESHOLD 65536
#define BIN_SEARCH_THRESHOLD 128

class PrefTab {
   public:
    PrefTab() {
    }

    ~PrefTab() {
       delete [] _X;
       delete [] _F;
    }

    PrefTab(const vector<std::string_view> &finstrings, uint64_t preflen) {

       vector<uint64_t> data;
       vector<uint8_t> lengths;
       for(const string_view& s : finstrings){ 
          uint64_t fini = 0;
          for(int i=0;i<s.length();i++){
             //pack a binary representation of the finimizer in the lower bits of fini
             fini |= ((uint64_t)char2bits(s[i])) << ((uint64_t)(2*(31-i)));
          }
          data.push_back(fini);
          lengths.push_back(s.length());
       }

       _nfins = data.size();
       _nprefs = ((uint64_t)1)<<(preflen<<1);
       uint64_t shift = 64-(preflen<<1);
       _preflen = preflen;

       _X = new uint32_t[_nprefs+1]; //+1 so we can eventually store _n in _X[_nprefs] and save a branch 
       for(uint64_t i=0;i<_nprefs+1;i++) _X[i] = 0;
       //Count number of items in each bucket
       for(uint64_t i=0;i<data.size();i++){
          _X[data[i]>>shift]++;
       }

       //Collect some statistics
       uint64_t maxcount = 0, emptycount = 0, nbigs = 0, nvbigs = 0;
       for(uint64_t i=0;i<_nprefs;i++){
          if(_X[i] > maxcount){ maxcount = _X[i]; }
          if(_X[i] > 128){ nbigs++; }
          if(_X[i] > 1024){ nvbigs++; }
          if(_X[i] == 0){ emptycount++; }
       }
       cerr << "_nprefs: "<<_nprefs<< '\n';
       cerr << "maxcount emptycount: "<< maxcount << ' ' << emptycount << '\n';
       cerr << "avg. items per bucket: "<<data.size()/(_nprefs-emptycount)<<'\n';
       cerr << "nbigs nvbigs: "<< nbigs << ' ' << nvbigs << '\n';

       //Prefix sum the counts in _X, creating pointers into _F.
       //Empty buckets are treated in a special way.
       uint32_t sum = 0;
       uint32_t lastnonempty = 0;
       uint32_t lastnonemptysz = 0;
       for(uint64_t i=0;i<_nprefs+1;i++){
          uint32_t v = _X[i];
          if(v){
             _X[i] = sum; //point to the first item in _X that's in this bucket
             lastnonempty = sum;
             lastnonemptysz = v;
             _X[i] = _X[i] << 1; //0 in LSB indicates non-empty
          }else{
             //empty bucket
             _X[i] = sum + v - 1; //point to the predecessor (in _F) of anything that lands in the bucket
             _X[i] = (_X[i] << 1) | 1; //1 in LSB indicates empty
          }
          sum += v;
       }
       _X[_nprefs] = (_nfins << 1) | 1;

       //Prepare the array of finimizers in _F, removing there prefixes and tagging them with their lengths in the least significant 5 bits
       _F = new uint64_t[_nfins];
       for(uint64_t i=0;i<_nfins;i++){
          //shift off prefix and tag with len
          _F[i] = (data[i] << (preflen<<1)) | lengths[i];
       }

       //run some tests
       //first for present elements
       for(uint64_t i=0;i<data.size();i++){
          pair<int64_t,bool> x = getPred(data[i]);
          if(i != x.first && x.second != true){
             cerr << i << ' ' << x.first << ' ' << x.second << '\n';
          }
       }
       cerr << "Test 1 passed.\n";
/*
       //now test for elements not in the set
       for(uint64_t i=1;i<data.size();i++){
          uint64_t missing = data[i-1]+1;
          while(missing < data[i]){
             pair<int64_t,bool> x = getPred(missing);
             if(i-1 != x.first && x.second != false){
                cerr << i << ' ' << x.first << ' ' << x.second << '\n';
             }
             missing++;
          }
       }
       cerr << "Test 2 passed.\n";
*/
    }

    //  p is the index of the predecessor in the set
    //  bool is 1 if the value at index p is equal to key, else 0
    pair<int64_t, bool> inline getPred(uint64_t key) const{
        uint64_t preflenbits = _preflen<<1;
        uint64_t shift = 64-preflenbits; 
        uint64_t p = key>>shift; //p has the prefix of the key of length preflenbits in its least significant bits
        key = key << preflenbits; //shift the prefix off the key
        uint32_t v = _X[p];
        if(v&1){
           //cerr << "Empty bucket\n";
           //empty bucket
           if(p!=0){
              return {v>>1, false};
           }
           return {-1, false};
        }
        //non-empty bucket
        v = v >> 1;
        //cerr << "Non-empty bucket\n";
        uint32_t blen = (_X[p+1]>>1) + (_X[p+1]&1) - v; //compute # items in bucket
        //cerr << "blen: " << blen << '\n';
        //cerr << "v: " << v << '\n';
        uint32_t i = v;
        for(;i<v+blen;i++){
           //cerr << "key: " << key << " curr: " << ((_F[i]>>preflenbits)<<preflenbits) <<'\n';
           if(((_F[i]>>preflenbits)<<preflenbits)>=key){
              break;
           }
        }
        return {i-(((_F[i]>>preflenbits)<<preflenbits)==key), (((_F[i]>>preflenbits)<<preflenbits)==key)};
    }

    pair<int64_t, uint64_t> finiLookup(uint64_t key, uint64_t keylen){ const{
       _n_searches++;
       uint64_t preflenbits = _preflen<<1;
       //cerr << "preflenbits: " << preflenbits << '\n';
       uint64_t shift = 64-preflenbits;
       uint64_t p = key>>shift; //p has the prefix of the key of length preflenbits in its least significant bits
       //cerr << "p: " << p << '\n';
       key = key << preflenbits; //shift the prefix off the key
       uint32_t v = _X[p];
       if(v&1){
          //cerr << "Empty bucket\n";
          //empty bucket
          _n_easy_searches++;
          _n_failed_searches++;
          return {-1,0};
       }
       //non-empty bucket
       _n_hard_searches++;
       v = v >> 1;
       //cerr << "Non-empty bucket\n";
       uint32_t blen = (_X[p+1]>>1) + (_X[p+1]&1) - v; //compute # items in bucket
       //cerr << "blen: " << blen << '\n';
       //cerr << "v: " << v << '\n';
       uint32_t i = v;
       //cerr << "v v+blen : " << v << ' ' << v+blen << '\n';
       
       if(blen < BIN_SEARCH_THRESHOLD){
          //cerr << "Scanning search: " << blen << " items \n";
          for(;i<v+blen;i++){
             //cerr << "Item: " << (i-v) << ", key: " << key << " curr: " << ((_F[i]>>preflenbits)<<preflenbits) <<'\n';
             if(((_F[i]>>preflenbits)<<preflenbits)>=key){
                break;
             }
          }
          //cerr << "Normal i: "<<i<<'\n';
          //cerr << (((_F[i]>>preflenbits)<<preflenbits)==key);
       }
       else{
          //cerr << "Binary search...\n";
          uint32_t low = v;
          uint32_t high = v+blen-1;
          while(low <= high){
             uint32_t mid = (low + high) >> 1;
             uint64_t midVal = ((_F[mid]>>preflenbits)<<preflenbits);
             if (midVal < key){
                low = mid + 1;
                i = low; //key not found, low is the index of the first element greater than key
             }else if (midVal > key){
                high = mid - 1;
                i = low; //key not found, low is the index of the first element greater than key
             }else{
                i = mid; break; // key found at mid
             }
          }
          //cerr << "Binary i: "<<i<<'\n';
       }
       //cerr << "i: " << i << '\n';
       i -= (1 - (((_F[i]>>preflenbits)<<preflenbits)==key)); //i now points to the predecessor
       //cerr << "Adjust i: " << i << '\n';
       //compute the LCP between the query and the predecessor
       uint64_t el = (_F[i] & (((uint64_t)1)<<preflenbits)-((uint64_t)1)); //the length of the finimizer
       uint64_t pred = ((_F[i]>>preflenbits)<<preflenbits);
       //uint64_t lcp = (__builtin_clzll((_F[i]>>preflenbits) ^ (~(key >> preflenbits))) - preflenbits)/2;
       uint64_t lcp = countl_zero(~((pred) ^ (~key)));
       lcp = _preflen + (lcp >> 1); //convert from bit matches to symbol matches
       //cerr << "pred el lcp: " << i << ' ' << el << ' ' << lcp << ' ' << ((lcp >= el) ? '*' : ' ') << '\n';
       if(lcp >= el && el <= keylen){
          return {i,el};
       }
       _n_failed_searches++;
       return {-1,0};
       //{i-(((_F[i]>>preflenbits)<<preflenbits)==key), (((_F[i]>>preflenbits)<<preflenbits)==key)};
    }

    uint64_t sizeInBytes() const{
       uint64_t sz = 0;
       sz += 3*sizeof(uint64_t);
       sz += _nprefs*sizeof(uint32_t);
       sz += _nfins*sizeof(uint64_t);
       return sz;
    }

    int64_t serialize(std::ostream& os) const{
       cerr << "PrefTab.serialize()...";
       int64_t written = 0;
       os.write((char *)&_preflen, sizeof(uint64_t));
       os.write((char *)&_nprefs, sizeof(uint64_t));
       os.write((char *)&_nfins, sizeof(uint64_t));
       os.write((char *)_X, (_nprefs+1)*sizeof(uint32_t));
       os.write((char *)_F, _nfins*sizeof(uint64_t));
       written += 3*sizeof(uint64_t);
       written += _nprefs*sizeof(uint32_t);
       written += _nfins*sizeof(uint64_t);
       for(int i=0;i<10;i++){
          cerr << _X[i] << ' ';
       }
       cerr << '\n';
       for(int i=0;i<10;i++){
          cerr << _F[i] << ' ';
       }
       cerr << '\n';
       return written;
    }

    void load(std::istream& is){
       cerr << "PrefTab.load()...";
       is.read((char *)&_preflen, sizeof(uint64_t));
       is.read((char *)&_nprefs, sizeof(uint64_t));
       is.read((char *)&_nfins, sizeof(uint64_t));
       cerr << "preflen nprefs nfins: "<<_preflen<<' '<<_nprefs<<' '<<_nfins<<'\n';
       _X = new uint32_t[_nprefs+1];
       is.read((char *)_X, (_nprefs+1)*sizeof(uint32_t));
       uint64_t nempties = 0;
       for(int i=0;i<_nprefs;i++){
          nempties += (_X[i] & 1);
       }
       cerr << "# empty buckets: "<<nempties<<'\n';
       _F = new uint64_t[_nfins];
       is.read((char *)_F, _nfins*sizeof(uint64_t));
       for(int i=0;i<10;i++){
          cerr << _X[i] << ' ';
       }
       cerr << '\n';
       for(int i=0;i<10;i++){
          cerr << _F[i] << ' ';
       }
       cerr << '\n';
    }
   
    //instrumenting statistics
    uint64_t _n_searches = 0;
    uint64_t _n_easy_searches = 0;
    uint64_t _n_hard_searches = 0;
    uint64_t _n_failed_searches = 0;
 
   private:
    uint64_t _preflen; //prefix length used for lookup table
    uint64_t _nprefs; //the number of prefixes == 1<<(_preflen<<1);
    uint64_t _nfins; //the number of finimizers == |_F|
    uint32_t *_X; //the lookup table _X[0.._nprefs-1]
    uint64_t *_F; //the finimizers, minus their prefixes, bit-packed and tagged with their lengths
};

#endif
