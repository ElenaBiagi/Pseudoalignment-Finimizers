#pragma once
#include <vector>
#include <cstdint>
#include <iostream>
#include <stdexcept>

template <typename T>
class BoundedDeque {

private:
    std::vector<T> buf;       // Circular buffer
    int64_t front_idx;        // Index before the logical front
    int64_t back_idx;         // Index where the next element will be inserted at the back
    int64_t n_elements;

    int64_t mod_increment(int64_t i) const {
        return (i + 1) % buf.size();
    }

    int64_t mod_decrement(int64_t i) const {
        return (i - 1 + buf.size()) % buf.size();
    }

    int64_t physical_index(int64_t logical_index) const {
        return (mod_increment(front_idx) + logical_index) % buf.size();
    }

public:
    BoundedDeque(int64_t max_size)
        : buf(max_size), front_idx(max_size - 1), back_idx(0), n_elements(0) {}

    const T& back() const {
        return buf[mod_decrement(back_idx)];
    }

    const T& front() const {
        return buf[mod_increment(front_idx)];
    }

    int64_t size() const {
        return n_elements;
    }

    bool empty() const {
        return n_elements == 0;
    }

    void push_back(const T& x) {
        buf[back_idx] = x;
        back_idx = mod_increment(back_idx);
        n_elements++;
    }

    void push_front(const T& x) {
        buf[front_idx] = x;
        front_idx = mod_decrement(front_idx);
        n_elements++;
    }

    void pop_front() {
        front_idx = mod_increment(front_idx);
        n_elements--;
    }

    void pop_back() {
        back_idx = mod_decrement(back_idx);
        n_elements--;
    }

    void clear() {
        n_elements = 0;
        front_idx = buf.size() - 1;
        back_idx = 0;
    }

    // Indexing operators
    T& operator[](size_t index) {
        if (index >= static_cast<size_t>(n_elements)) {
            throw std::out_of_range("BoundedDeque index out of range");
        }
        return buf[physical_index(index)];
    }

    const T& operator[](size_t index) const {
        if (index >= static_cast<size_t>(n_elements)) {
            throw std::out_of_range("BoundedDeque index out of range");
        }
        return buf[physical_index(index)];
    }
};
