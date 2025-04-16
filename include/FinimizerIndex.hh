#pragma once

#include <string>
#include <cstring>
#include <algorithm>
#include <bitset>
#include "sbwt/cxxopts.hpp"
#include "sbwt/globals.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/SubsetWT.hh"
#include "sbwt/stdlib_printing.hh"
#include "sbwt/SeqIO.hh"
#include "sbwt/SubsetMatrixRank.hh"
#include "sbwt/buffered_streams.hh"
#include "sbwt/variants.hh"
#include "sbwt/commands.hh"
#include <filesystem>
#include <cstdio>
#include <optional>

#include "sbwt/throwing_streams.hh"
#include "PackedStrings.hh"
#include "SeqIO.hh"
#include "common.hh"

class FinimizerIndex{

public:
// TODO remove
    struct QueryResult{
        vector<pair<int64_t, int64_t>> local_offsets; // Unitig id, distance from the start of the unitig
        int64_t n_found = 0;
    };

private:
    // Forbid copying because we have pointers to our internal data structures
    FinimizerIndex(const FinimizerIndex& other) = delete;
    FinimizerIndex& operator=(const FinimizerIndex& other) = delete;

public:

    // Note: if you add members, update size_in_bytes(), serialize(), and load()
    unique_ptr<plain_matrix_sbwt_t> sbwt; // These are smart pointers because they are passed in to the constructor
    unique_ptr<sdsl::int_vector<>> LCS; // These are smart pointers because they are passed in to the constructor
    std::unordered_map<std::string, std::set<int>> hashTable; // Create a hash table to store the list of colors for each Finimizer


    FinimizerIndex() {}

    //QueryResult 
    void search(const std::string& query, unordered_map<int, uint64_t>& results) const {
  
        //std::cerr << "Searching " << query << std::endl;

        const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();
        const vector<int64_t>& C = sbwt.get_C_array();
        const int64_t query_len = query.length();


        if (query.size() < k) return; 

        unordered_map<string, uint64_t> Finimizers = rarest_fmin_streaming_search(sbwt, *LCS, query);
        // Check the colors for every finimizer found
        pseudoalignemnt_stats(Finimizers, this->hashTable, results);

        return;
    }

    void search(const std::string& query, vector<pair<int, float>>& results, const float& t) const {
  
        const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();
        const vector<int64_t>& C = sbwt.get_C_array();
        const int64_t query_len = query.length();


        if (query.size() < k) return; 

        unordered_map<string, uint64_t> Finimizers = rarest_fmin_streaming_search(sbwt, *LCS, query);
        // Check the colors for every finimizer found
        pseudoalignemnt_stats(Finimizers, this->hashTable, results, t);

        return;
    }


    void serialize_HashTable(const std::unordered_map<std::string, std::set<int>>& hashTable, const std::string& hashTableName) const{
        std::ofstream hashTable_out(hashTableName, std::ios::binary);
        if (!hashTable_out) {
            std::cerr << "Error: Could not open file for writing!" << std::endl;
            return;
        }

        // Write the number of elements in the hash table
        size_t hashTableSize = hashTable.size();
        hashTable_out.write(reinterpret_cast<const char*>(&hashTableSize), sizeof(hashTableSize));

        for (const auto& [key, value] : hashTable) {
            // Write the size of the key and the key itself
            size_t keySize = key.size();
            hashTable_out.write(reinterpret_cast<const char*>(&keySize), sizeof(keySize));
            hashTable_out.write(key.data(), keySize);

            // Write the size of the value set and the set itself
            size_t valueSize = value.size();
            hashTable_out.write(reinterpret_cast<const char*>(&valueSize), sizeof(valueSize));
            for (const int elem : value) {
                hashTable_out.write(reinterpret_cast<const char*>(&elem), sizeof(elem));
            }
        }
        hashTable_out.close();
    }
 
    std::unordered_map<std::string, std::set<int>> load_HashTable(const std::string& hashTableName) {
    //std::unordered_map<std::string, std::set<int>> hashTable;
    hashTable = this->hashTable;

    std::ifstream inFile(hashTableName, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error: Could not open file for reading!" << std::endl;
        return hashTable;
    }

    // Read the number of elements in the hash table
    size_t hashTableSize;
    if (!inFile.read(reinterpret_cast<char*>(&hashTableSize), sizeof(hashTableSize))) {
        std::cerr << "Error: Failed to read hash table size!" << std::endl;
        return hashTable;
    }

    for (size_t i = 0; i < hashTableSize; ++i) {
        // Read the key
        size_t keySize;

        if (!inFile.read(reinterpret_cast<char*>(&keySize), sizeof(keySize))) {
            std::cerr << "Error: Failed to read key size!" << std::endl;
            return hashTable;
        }

        std::string key(keySize, '\0');
        if (!inFile.read(&key[0], keySize)) {
            std::cerr << "Error: Failed to read key!" << std::endl;
            return hashTable;
        }

        // Read the value set
        size_t valueSize;

        if (!inFile.read(reinterpret_cast<char*>(&valueSize), sizeof(valueSize))) {
            std::cerr << "Error: Failed to read value size!" << std::endl;
            return hashTable;
        }

        std::set<int> value;
        for (size_t j = 0; j < valueSize; ++j) {
            int elem;
            if (!inFile.read(reinterpret_cast<char*>(&elem), sizeof(elem))) {
                std::cerr << "Error: Failed to read value element!" << std::endl;
                return hashTable;
            }
            value.insert(elem);
        }

        // Insert the key-value pair into the hash table
        hashTable[key] = value;

    }

    inFile.close();
    return hashTable;
}

    void serialize(const string& index_prefix) const {

        std::ofstream LCS_out(index_prefix + ".LCS.sdsl");
        sdsl::serialize(*LCS, LCS_out);

        sbwt->serialize(index_prefix + ".sbwt");

        serialize_HashTable(hashTable, index_prefix + ".ht.BIN");
    }

    void load(const string& index_prefix) {

        LCS = make_unique<sdsl::int_vector<>>();
        ifstream LCS_in(index_prefix + ".LCS.sdsl");
        sdsl::load(*LCS, LCS_in);
        std::cerr<< "LCS_file loaded"<<std::endl;

        sbwt = make_unique<plain_matrix_sbwt_t>();
        sbwt->load(index_prefix + ".sbwt");
        std::cerr << "SBWT matrix loaded" << std::endl;

        load_HashTable(index_prefix + ".ht.BIN");
        std::cerr << "hashTable loaded" << std::endl;
        //printHashTable(hashTable);
    }

    // TODO: add hash table (later. not relevant now)
    // This also includes the rank structures which are not serialized
    int64_t size_in_bytes() const{
        int64_t total = 0;
        total += sdsl::size_in_bytes(*LCS);

        sbwt::SeqIO::NullStream ns;
        total += sbwt->serialize(ns);
        return total;
    }
};


class FinimizerIndexBuilder{
public:

    unique_ptr<plain_matrix_sbwt_t> sbwt;
    unique_ptr<sdsl::int_vector<>> LCS;

    unique_ptr<FinimizerIndex> index;



    // Takes ownership of sbwt and LCS
    FinimizerIndexBuilder(unique_ptr<plain_matrix_sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, const vector<string>& incolors) {
        index = make_unique<FinimizerIndex>();
        this->sbwt = move(sbwt); // Take ownership
        this->LCS = move(LCS); // Take ownership

        int64_t n_nodes = this->sbwt->number_of_subsets();

        //TODO REPLACE hashTable
        std::unordered_map<std::string, std::set<int>> hashTable; // Create a hash table to store the list of colors for each Finimizer
        
        //helper
        std::unordered_map<uint32_t, std::map<char, set< pair<uint32_t, set<int> >> >> B;

        // REAL DATA STR
        //std::unordered_map<uint32_t, uint32_t> B; // Create a hash table to store the prefixes of each bucket and a pointer to the start of the tails in the sdsl int vector
        std::unordered_map<uint32_t, std::set<int> > sB; // Create a hash table to store the finimizers shorter than the prefix lenght
        int_vector<0> T; // I can decide the width and modify every entry

        // Scan the genomes to get the list of colors for each finimizer using the hash table
        typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_colors_gzip;
        typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_colors_no_gzip;
        for (int i=0; i< incolors.size(); i++){ // colors are likely not a huge number
            //std::cerr << "scanning color " << i << std::endl;

            bool gzip_colors = SeqIO::figure_out_file_format(incolors[i]).gzipped;
    
            if (gzip_colors){
                run_colors_file<in_colors_gzip>(incolors[i], hashTable, i);
            }
            else{
                run_colors_file<in_colors_no_gzip>(incolors[i], hashTable, i);
            }
            std::cerr << "DONE"<< std::endl;
        }
        // TODO extract statistics
        get_stats(hashTable);
        
        index->sbwt = std::move(this->sbwt); // Transfer ownership
        index->LCS = std::move(this->LCS); // Transfer ownership 
        index->hashTable = std::move(hashTable); // Transfer ownership
        index->B = std::move(B); // Transfer ownership

    }

    // TODO fix return type
    template<typename reader_t>
    int64_t run_colors_file(const string& infile, unordered_map<std::string, std::set<int>>& hashTable, int& i){
        reader_t reader(infile);
        write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
        return from_reader_to_seq(reader, hashTable, i);
    }

    // TODO fix return type
    template<typename reader_t>
    int from_reader_to_seq(reader_t& reader, unordered_map<std::string, std::set<int>>& hashTable, int& i) {
        
        const int64_t k = sbwt->get_k();
        while(true){
            int64_t len = reader.get_next_read_to_buffer();
            if(len == 0) [[unlikely]] break;

            //const std::string& seq = remove_N_from_string(reader.read_buf);
            vector<string> seq_vector = split_by_N(reader.read_buf, k);
            for (const string &seq : seq_vector){
                scan_color(seq, hashTable, i);
            }
            std::cerr << "done seq"<< std::endl;
            
            for (const string &seq : seq_vector){
                const string reverse = sbwt::get_rc(seq);
                scan_color(reverse, hashTable, i);
            }
            std::cerr << "done rev"<< std::endl;

        }
        return 1;
    }
    // TODO fix return type
    int scan_color(const std::string& seq, unordered_map<std::string, std::set<int>>& hashTable, int& i) {
        //std::cerr << seq << std::endl;
        const int64_t n_nodes = sbwt->number_of_subsets();
        const int64_t k = sbwt->get_k();
        const vector<int64_t>& C = sbwt->get_C_array();

        int64_t freq;
        BoundedDeque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin(seq.size());
        const int64_t str_len = seq.size();
        tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, end}

        int64_t kmer = 0;
        int64_t start = 0;
        int64_t end;
        pair<int64_t, int64_t> I = {0, n_nodes - 1};
        int64_t I_start;
        tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
        char c;
        
        for (end = 0; end < str_len; end++) {
            c = static_cast<char>(seq[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
            //update the sbwt INTERVAL
            I = this->sbwt->update_sbwt_interval(&c, 1, I);
            // TODO REMOVE CHECK
             if (I.first ==-1){
                std::cerr << "This should be impossible!, pos " << end << ", char " << c << ", len " << end - start + 1 << " " << seq.substr(start, end - start + 1 ) << std::endl;
                return 0;
            }
            freq = (I.second - I.first + 1);
            I_start = I.first;
            if (freq == 1){ // 1. rarest 
                while (freq == 1) {  //2. shortest
                    curr_substr = {freq, end - start + 1, I_start, end};
                    // (2) drop the first char
                    // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                    start++;
                    I = drop_first_char(end - start + 1, I, *(this->LCS), n_nodes);
                    freq = (I.second - I.first + 1);
                    I_start = I.first;
                }
                if (w_fmin > curr_substr) {
                    all_fmin.clear();
                    w_fmin = curr_substr;
                } else{
                    while (all_fmin.back() > curr_substr) {all_fmin.pop_back();}
                }
                all_fmin.push_back(curr_substr);
            }
            if (end >= k -1 ){
                // Add the color to the finimizer
                string F = seq.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin));
                if (hashTable.find(F) != hashTable.end()) {
                    hashTable[F].insert(i);
                } else {
                    hashTable[F]={i};
                }
                
                kmer++;

                // Check if the current finimizer is still in this window
                while (get<3>(w_fmin)- get<1>(w_fmin)+1 < kmer) { // start
                    all_fmin.pop_front();
                    w_fmin = (all_fmin.size()==0) ? tuple<int64_t, int64_t, int64_t, int64_t> {n_nodes,k+1,kmer+1,kmer+k} : all_fmin.front();
                }
            }
        } 
        return 1;
    }


    uint64_t prefix2int(const string& s, uint64_t offset, char plen){
        uint64_t h = 0;
        for(uint64_t i=0; i<(uint64_t)plen; i++){
           uint64_t b = get_char_idx(s[i+offset]);
           h |= (b << (i<<1));
        }
        //cerr << h << '\n';
        return h;
    }

    uint64_t suffix2int(const std::string& s, uint64_t offset, char slen) {
        uint64_t h = 0;
        for (uint64_t i = 0; i < (uint64_t)slen; i++) {
            uint64_t b = get_char_idx(s[offset + slen - 1 - i]);
            h |= (b << (i << 1)); 
        }
        return h;
    }

    // TODO
    // DIVIDE INTO BUCKETS
    void Buckets (unordered_map<std::string, std::set<int>>& hashTable,unordered_map<uint32_t, std::map<char, set< pair<uint32_t, set<int> >> >>& B, unordered_map<uint32_t, std::set<int> >& sB, char plen){
        // create a hash table with all the possible strings of length plen

        // B={prefix:{tlen1:{{tail1,colors1},...}, tlen2:{{tail1,colors1},...},... }}
        uint32_t nbuckets = 1<<(plen<<1);
        for (uint32_t p=0; p<nbuckets; p++){
            B[p]={};
        }
        for (auto &f : hashTable){
            string fmin = f.first;
            char flen = fmin.length();
            if (flen < plen){ 
                uint32_t sfmin = prefix2int(fmin,0,flen);
                sB[sfmin]= f.second; // TODO this points you directly to the colors
            }
            else{
                uint32_t pfmin = prefix2int(fmin,0,plen);
                char tlen = flen - plen;
                // TODO extract tails
                uint32_t tail = 0; //TODO SUFFIX2int(fmin,0,plen);
                if (B.find(pfmin) != B.end()) {
                    //B[pfmin]=1; // TODO this should contain a pointer to the right place in the intvector of all the tails
                    if (B[pfmin].find(tlen) != B[pfmin].end()){
                        B[pfmin][tlen].insert({tail, f.second}); 
                    }
                    else{
                        B[pfmin][tlen]={{tail, f.second}}; 
                    }
                } else {
                    cerr << pfmin << " not found in B!";
                }
            }
        }
    }

    void write_bits(int_vector<0>& vec, uint64_t value, size_t offset, size_t bit_width) {
        for (size_t i = 0; i < bit_width; ++i) {
            vec[offset + i] = (value >> (bit_width - 1 - i)) & 1;
        }
    }

    vector<uint8_t> vbyte_encode(uint64_t x) {
        vector<uint8_t> bytes;
        do {
            uint8_t byte = x & 0x7F;
            x >>= 7;
            if (x != 0) byte |= 0x80;
            bytes.push_back(byte);
        } while (x != 0);
        return bytes;
    }

    void storeTails(unordered_map<uint32_t, std::map<char, set< pair<uint32_t, set<int> >> >>& B, int_vector<0>& T){
        size_t total_bits = 0;

        for (auto &prefix : B){
            for (auto &tails: prefix.second){
                // this should in tails order
                char tlen = tails.first; // 5 bits
                int32_t tnumber = tails.second.size(); // vbyte
                total_bits += 5; // tlen: 5 bits
                auto vb = vbyte_encode(tnumber);
                total_bits += vb.size() * 8; // vbyte: 8 bits per byte
            }

        }
        T = int_vector<0>(total_bits, 0);

    }

    void pack_all(const vector<pair<uint8_t, vector<uint64_t>>>& groups, int_vector<0>& out) {
        size_t total_bits = 0;
    
        // First, calculate total bits needed
        for (const auto& [tlen, numbers] : groups) {
            total_bits += 5; // tlen: fixed 5 bits
            auto vb = vbyte_encode(numbers.size());
            total_bits += vb.size() * 8; // vbyte: 8 bits per byte
            total_bits += numbers.size() * (2 * tlen); // data // !!!!!!!!!!!!!!!!!!!What is this???
        }
    
        out = int_vector<0>(total_bits, 0);
    
        // Now, write the data
        size_t offset = 0;
        for (const auto& [tlen, numbers] : groups) {
            // Write tlen (5 bits)
            write_bits(out, tlen, offset, 5);
            offset += 5;
    
            // Write tnumber (vbyte)
            auto vb = vbyte_encode(numbers.size());
            for (uint8_t b : vb) {
                write_bits(out, b, offset, 8);
                offset += 8;
            }
    
            // Write the tnumber integers using 2*tlen bits each
            for (uint64_t val : numbers) {
                write_bits(out, val, offset, 2 * tlen);
                offset += 2 * tlen;
            }
        }
    }
    

    // Transfer ownership of the index out of the builder
    unique_ptr<FinimizerIndex> get_index(){
        return std::move(this->index);
    }
};
