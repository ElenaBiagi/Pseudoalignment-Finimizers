#include <vector>
#include <tuple>
#include <cstdint>


using MyTuple = std::tuple<uint8_t, uint64_t, uint32_t, uint64_t>; // {f_len, f_int, f_color, start}

class CBuffer {
    std::vector<MyTuple> buffer;
    size_t head = 0;
    size_t size = 0;
    size_t capacity;

public:
    CBuffer(size_t k) : buffer(k), capacity(k) {}

    void insert(const MyTuple& t) {
        buffer[head] = t;
        head = (head + 1) % capacity;
        if (size < capacity) ++size;
    }

    template <typename Func>
    void for_each_recent(Func f) const {
        for (size_t i = 0; i < size; ++i) {
            size_t idx = (head + capacity - 1 - i) % capacity;
            if (!f(buffer[idx])) {
                break;
            }
        }
    }

    size_t current_size() const {
        return size;
    }
};
