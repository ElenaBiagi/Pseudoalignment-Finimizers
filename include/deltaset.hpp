
#include <cstdint>
#include <iostream>
#include <vector>

//find/replace 64 wih 32 to speed up // DONE
//find/replace uint8_t with uint16_t if needed // DONE but also added a check to see if it's necessary
using namespace std;
class DeltaSet{
private:
    // arrays -> vectors, for easier memory management (avoid issues with serialize and load)
    std::vector<uint64_t> _prefix_sums;
    std::vector<uint16_t> _diffs;
    uint64_t _period = 32;
    uint64_t _n;
public:
    DeltaSet() = default;

    DeltaSet(std::vector<uint64_t> &starts){
        _n = (uint64_t)starts.size();
        if (_n == 0 ){return;}
        _period = 32;
        _prefix_sums.assign( (_n/_period)+1,0);
        _diffs.assign(_n,0);
        _prefix_sums[0] = starts[0];
        uint64_t pi = 1;
        for(uint64_t i=1;i<_n;i++){
            if((i%_period)==0){
                _prefix_sums[pi++] = starts[i];
            }
            uint64_t diff = starts[i] - starts[i-1]; 
            //if (diff > UINT8_MAX){ throw std::overflow_error("diff too large for uint8_t");} // check this value would fit in 8 bits
            _diffs[i-1] = (uint16_t)diff;
        }
    }

    uint64_t get_start(uint64_t i) const {
        uint64_t p = _prefix_sums[i/_period];
        for(uint64_t j = _period*(i/_period); j < i; j++){
            p += _diffs[j];
        }
        return p;
    }

/*     void read_all(uint64_t start, uint64_t end, uint64_t freq, vector<uint64_t>& results){
        uint64_t p = 0;
        for (auto i = start; i<end; i++){
            // no need to use a prefix sum right?
            //uint64_t p = _prefix_sums[i/32];
            //for(uint64_t j = 32*(i/32); j < i; j++){
            p += _diffs[i];
            results[p]+= freq;
        }
    } */

   size_t size() const {
        return _n;
   }

void serialize(std::ostream& out) const {
    if (_period == 0)
        throw std::runtime_error("EF::serialize(): invalid period (0)");

    out.write(reinterpret_cast<const char*>(&_n), sizeof(_n));
    out.write(reinterpret_cast<const char*>(&_period), sizeof(_period));

    const size_t prefix_count = (_n / _period) + 1;

    if (_prefix_sums.size() < prefix_count)
        throw std::runtime_error("EF::serialize(): prefix_sums smaller than expected");

    if (prefix_count > 0)
        out.write(reinterpret_cast<const char*>(_prefix_sums.data()),
                  prefix_count * sizeof(uint64_t));

    if (_diffs.size() < _n)
        throw std::runtime_error("EF::serialize(): diffs smaller than expected");

    if (_n > 0)
        out.write(reinterpret_cast<const char*>(_diffs.data()),
                  static_cast<size_t>(_n) * sizeof(uint16_t));
}

void load(std::istream& in) {
    in.read(reinterpret_cast<char*>(&_n), sizeof(_n));
    in.read(reinterpret_cast<char*>(&_period), sizeof(_period));

    if (!in)
        throw std::runtime_error("EF::load(): failed to read n/period");
    if (_period == 0)
        throw std::runtime_error("EF::load(): invalid period (0)");
    if (_n > 1e8) // sanity bound; adjust to your expected size
        throw std::runtime_error("EF::load(): suspiciously large n (possible corruption)");

    const size_t prefix_count = (_n / _period) + 1;

    _prefix_sums.assign(prefix_count, 0);
    _diffs.assign(_n, 0);

    if (prefix_count > 0) {
        in.read(reinterpret_cast<char*>(_prefix_sums.data()),
                prefix_count * sizeof(uint64_t));
        if (!in)
            throw std::runtime_error("EF::load(): failed to read prefix_sums");
    }

    if (_n > 0) {
        in.read(reinterpret_cast<char*>(_diffs.data()),
                static_cast<size_t>(_n) * sizeof(uint16_t));
        if (!in)
            throw std::runtime_error("EF::load(): failed to read diffs");
    }
}
};
