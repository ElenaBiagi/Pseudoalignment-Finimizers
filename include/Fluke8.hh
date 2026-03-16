#ifndef _FLUKE8_H_
#define _FLUKE8_H_

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
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

#include "Pred8v2.hh"
#include "Util.hh"

using namespace std;
using namespace std::chrono;

class Fluke8 {
   public:
    Fluke8() {
    }

    ~Fluke8() = default;
        
    Fluke8(const std::string filename, const uint64_t lowerlen, const uint64_t upperlen){
       cerr << "Fluke8, constructing from file...\n";
       std::ifstream is(filename);
       if (!is.is_open()) {
           std::cerr << "Error: Unable to open file." << std::endl;
           exit(1);
       }
       std::string line;
       uint64_t nfins = 0;
       cerr << "\tMaking first pass...\n";
       while (std::getline(is, line)) {
           if(line.length() >= lowerlen && line.length() <= upperlen){
              nfins++;
              if(nfins % 100000000 == 0) cerr << nfins << '\n';
           }
       }
       cerr << "\t" << nfins << " finimizers in range...\n";
       is.clear();
       is.seekg(0);

       cerr << "\tMaking second pass...\n";
       _n = nfins;
       vector<uint64_t> data;
       data.reserve(_n);
       while (std::getline(is, line)) {
          if(line.length() >= lowerlen && line.length() <= upperlen){
             int64_t fini = 0;
             for(int i=0;i<line.length();i++){
                //pack a binary representation of the finimizer in the lower bits of fini
                fini |= ((uint64_t)char2bits(line[i])) << ((uint64_t)(2*(upperlen-1-i)));
             }
             data.push_back(fini); 
          }
       }
       _min = data[0];
       _u = data[data.size()-1];
       _lowerlen = lowerlen;
       _upperlen = upperlen;
       _pred8 = Pred8v2(data);
       data.clear();

       is.clear();
       is.seekg(0);

       cerr << "\tMaking third pass...\n";
       _lengths.resize(_n);
       uint64_t i=0;
       while (std::getline(is, line)) {
          if(line.length() >= lowerlen && line.length() <= upperlen){
             _lengths[i++] = line.length();
          }
       }
       is.close();
    }


    Fluke8(const vector<std::string_view> &finstrings, const uint64_t lowerlen, const uint64_t upperlen) {
       vector<uint64_t> data;
       vector<uint8_t> lengths;
       for(const string_view& s : finstrings){
          uint64_t fini = 0;
          for(int i=0;i<s.length();i++){
             //pack a binary representation of the finimizer in the lower bits of fini
             fini |= ((uint64_t)char2bits(s[i])) << ((uint64_t)(2*(upperlen-1-i)));
          }
          data.push_back(fini);
          lengths.push_back(s.length());
       }

       _n = data.size();
       _min = data[0];
       _u = data[data.size()-1];
       _lowerlen = lowerlen;
       _upperlen = upperlen;
       _pred8 = Pred8v2(data);
       _lengths.resize(_n);
       for(uint64_t i=0;i<_n;i++){
          _lengths[i] = lengths[i];
       }
    }

    pair<int64_t,uint64_t> finiLookup(uint64_t key, uint64_t len) const{
       //cerr << "key: " << key << '\n';
       uint64_t cle = key>>(((uint64_t)64) - ((uint64_t)2)*_upperlen);
       //cerr << "cle: " << cle << '\n';
       pair<uint64_t, uint64_t> p = _pred8.getPredPrefix(cle);
       //cerr << p.first << ' ' << p.second << '\n';
       if(p.second){
          uint64_t lcp = p.second;
          lcp -= (((uint64_t)64)-(_upperlen<<1));
          lcp = lcp >> ((uint64_t)1); //lcp is now in symbols
          //cerr << "key: " << key << " lcp: "<<lcp<<" _lengths[p.first]: "<<((uint64_t)_lengths[p.first])<<" len: "<<len<<'\n';
          //cerr << ((lcp >= _lengths[p.first]) ? '*' : ' ') << '\n';
          if(lcp >= _lengths[p.first] && _lengths[p.first] <= len){
             return {p.first,_lengths[p.first]};
          }
       }
       _n_failed_searches++;
       return {-1,0};
    }

    pair<int64_t, bool> inline getPred(uint64_t key) const{
       uint64_t cle = key>>(((uint64_t)64) - ((uint64_t)2)*_upperlen);
       return _pred8.getPred(cle);
    }

    int64_t rank(int64_t pos) const{
       return _pred8.rank(pos);
    } 

    uint64_t getu() const { return _u; }
    uint64_t getn() const { return _n; }

    uint64_t sizeInBytes() const{ 
        uint64_t sz = 5*sizeof(uint64_t);
        sz += _n; //for lengths array
        sz += _pred8.sizeInBytes();
        return sz;
    }

    int64_t serialize(std::ostream& os) const{
        cerr << "Fluke8.serialize()...\n";
        int64_t written = 0;
        os.write((char *)&_u, sizeof(uint64_t));
        os.write((char *)&_n, sizeof(uint64_t));
        os.write((char *)&_min, sizeof(uint64_t));
        os.write((char *)&_upperlen, sizeof(uint64_t));
        os.write((char *)&_lowerlen, sizeof(uint64_t));
        cerr << "_u _n _min _upperlen _lowerlen: "<<_u<<' '<<_n<<' '<<_min<<' '<<_upperlen<<' '<<_lowerlen<<'\n';
        _pred8.serialize(os);
        os.write((char *)_lengths.data(), _n*sizeof(uint8_t));
        written += 5*sizeof(uint64_t) + _pred8.sizeInBytes() + _n;
        return written;
    }

    void load(std::istream& is){
       cerr << "Fluke8.load()...\n";
       is.read((char *)&_u, sizeof(uint64_t));
       is.read((char *)&_n, sizeof(uint64_t));
       is.read((char *)&_min, sizeof(uint64_t));
       is.read((char *)&_upperlen, sizeof(uint64_t));
       is.read((char *)&_lowerlen, sizeof(uint64_t));
       cerr << "_u _n _min _upperlen _lowerlen: "<<_u<<' '<<_n<<' '<<_min<<' '<<_upperlen<<' '<<_lowerlen<<'\n';
       _pred8.load(is);
       _lengths.resize(_n);
       is.read((char *)_lengths.data(), _n*sizeof(uint8_t));
    }
    
//    Fluke8(Fluke8 &other){
//       this->_u = other._u;
//       this->_n = other._n;
//       this->_min = other._min;
//       this->_pred8 = other._pred8;
//       this->_lengths = other._lengths;
//    }
//


    uint64_t getUpperlen(){
       return _upperlen;
    }

    uint64_t getLowerlen(){
       return _lowerlen;
    }

    Pred8v2 _pred8;
    mutable uint64_t _n_failed_searches = 0;
   private:
    uint64_t _u = 0;  // universe size
    uint64_t _n = 0;  // number of elements
    uint64_t _min = 0;  // value of the smallest element
    uint64_t _upperlen = 0;  // universe size
    uint64_t _lowerlen = 0;  // number of elements
    std::vector<uint8_t> _lengths;
};

#endif
