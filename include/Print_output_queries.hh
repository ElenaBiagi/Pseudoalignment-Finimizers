#pragma once

#include <string>
#include <cstring>

#include <filesystem>
#include <cstdio>
//#include <optional>
//#include <variant>
#include <sstream>

#include "SeqIO.hh"

using namespace std;

// template <typename out_stream_t>
// void write_out(const char *data, int64_t data_length, out_stream_t &output_writer, vector<char> &buffer, const size_t flush_t = 8 * 1024 * 1024)
// { // TODO find a better way to write threshold
//     for (int64_t i = 0; i < data_length; i++)
//     {
//         buffer.push_back(data[i]);
//         if (buffer.size() > flush_t && data[i] == '\n')
//         {
//             output_writer.write(buffer.data(), buffer.size());
//             buffer.clear(); // Let's hope this keeps the reserved capacity of the vector intact
//         }
//     }
// }

// uint16_t fast_int_to_string(uint16_t x, char *buffer)
// {
//     uint16_t i = 0;
//     // Write the digits in reverse order (reversed back at the end)
//     do
//     {
//         buffer[i++] = '0' + (x % 10);
//         x /= 10;
//     } while (x > 0);
//     std::reverse(buffer, buffer + i);
//     buffer[i] = '\0';
//     return i;
// }

// void print_out_queries(const vector<int16_t> results, const uint16_t q_id, const float t){
//     constexpr size_t flush_t = 8 * 1024 * 1024; // 8 MB //1 << 20; // 1MB

//     std::vector<char> output_buffer;
//     output_buffer.reserve(flush_t * 2);

//     // For output writing
//     char int_buf[32]; // Enough space for a 64-bit integer in ascii

//     if (t > 0){

//     }
//     else
//     {
//         auto start = std::chrono::high_resolution_clock::now();
//         for (auto idx = 0; idx < results.size(); idx++)
//         {
//             const auto &count = results[idx];
//             if (count > 0)
//             {
//                 uint16_t idx_len = fast_int_to_string(idx, int_buf);
//                 write_out(int_buf, idx_len, out, output_buffer, flush_t);
//                 // write_out(" ", 1, out, output_buffer, flush_t);

//                 // print the number of kmers/matches found
//                 write_out(":", 1, out, output_buffer, flush_t);

//                 uint16_t count_len = fast_int_to_string(count, int_buf);
//                 write_out(int_buf, count_len, out, output_buffer, flush_t);

//                 write_out(" ", 1, out, output_buffer, flush_t);
//             }
//         }
//         /* for (int a = static_cast<int>(ans.size()) - 1; a >= 0; a--)
//         {
//             const auto &[idx, count] = ans[a];
//             if (count >= 0)
//             {
//                 uint16_t idx_len = fast_int_to_string(idx, int_buf);
//                 write_out(int_buf, idx_len, out, output_buffer, flush_t);
//                 write_out(" ", 1, out, output_buffer, flush_t);
//             }
//         } */

//         write_out("\n", 1, out, output_buffer, flush_t);

//         auto end = std::chrono::high_resolution_clock::now();
//         //time_output += (end - start);
//     }
    
    

// }

void print_cout_queries(const vector<int16_t> results, const uint16_t q_id, const uint64_t T){
    
    cout << q_id << " ";

    //auto start = std::chrono::high_resolution_clock::now();
    for (auto idx = 0; idx < results.size(); idx++)
    {
        const auto &count = results[idx];
        if (count > T)
        {   
            cout << idx << ":" << count << " ";
        }
    }
    cout << "\n";

    //auto end = std::chrono::high_resolution_clock::now();

}

void print_cout_queries(const vector<int16_t> results, const uint16_t q_id){
    
    cout << q_id << " ";

    for (auto idx = 0; idx < results.size(); idx++)
    {
        const auto &count = results[idx];
        if (count > 0)
        {   
            cout << idx << ":" << count << " ";
        }
    }
    cout << "\n";

}

