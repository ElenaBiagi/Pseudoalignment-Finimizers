
#include <cstdint>
#include <iostream>
#include <vector>

//find/replace 64 wih 32 to speed up // DONE
//find/replace uint8_t with uint16_t if needed

class DeltaSet{
private:
   uint64_t *_prefix_sums;
   uint16_t *_diffs;
   uint64_t _period;
   uint64_t _n;
public:
    DeltaSet() = default;

   DeltaSet(std::vector<uint64_t> &starts){
      _n = (uint64_t)starts.size();
      _prefix_sums = new uint64_t[_n/32];
      _diffs = new uint16_t[_n];
      _diffs[0] = starts[0];
      uint64_t pi = 0;
      for(uint64_t i=1;i<_n;i++){
         if((i%32)==0){
            _prefix_sums[pi++] = starts[i];
         }
         uint64_t diff = starts[i] - starts[i-1]; // TODO could check this value fits in 8 bits
         _diffs[i] = (uint16_t)diff;
      }
   }
   ~DeltaSet(){
      delete [] _diffs;
      delete [] _prefix_sums;
   }
   uint64_t get_start(uint64_t i) const {
      uint64_t p = _prefix_sums[i/32];
      for(uint64_t j = 32*(i/32); j < i; j++){
         p += _diffs[j];
      }
      return p;
   }

   size_t size() const {
        return _n;
   }

   void serialize(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&_n), sizeof(_n));
        out.write(reinterpret_cast<const char*>(&_period), sizeof(_period));

        size_t prefix_count = _n / _period + 1;
        out.write(reinterpret_cast<const char*>(_prefix_sums), prefix_count * sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(_diffs), _n * sizeof(uint16_t));
    }

    void load(std::istream& in) {
        delete[] _diffs;
        delete[] _prefix_sums;

        in.read(reinterpret_cast<char*>(&_n), sizeof(_n));
        in.read(reinterpret_cast<char*>(&_period), sizeof(_period));

        size_t prefix_count = _n / _period + 1;
        _prefix_sums = new uint64_t[prefix_count];
        _diffs = new uint16_t[_n];

        in.read(reinterpret_cast<char*>(_prefix_sums), prefix_count * sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(_diffs), _n * sizeof(uint16_t));
    }
};
