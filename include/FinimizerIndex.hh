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
#include "rarest_fmin_search.hh"

class FinimizerIndex{

private:
    // Forbid copying because we have pointers to our internal data structures
    FinimizerIndex(const FinimizerIndex& other) = delete;
    FinimizerIndex& operator=(const FinimizerIndex& other) = delete;

    int k = 31; // TODO FIX THIS
    uint8_t plen = 10;

public:

    // Note: if you add members, update size_in_bytes(), serialize(), and load()
    unique_ptr<plain_matrix_sbwt_t> sbwt; // These are smart pointers because they are passed in to the constructor
    unique_ptr<sdsl::int_vector<>> LCS; // These are smart pointers because they are passed in to the constructor

    std::unordered_map<uint32_t, pair<int64_t,int64_t>> B; // Create a hash table to store the prefixes of each bucket and a pointer to the start of the tails in the sdsl int vector
    std::unordered_map<uint32_t, int64_t> sB; // Create a hash table to store the finimizers shorter than the prefix length
    //unique_ptr<int_vector<1>> T; // tails
    std::vector<int_vector<1>> T;
    vector<set<int>> C;
    //TODO K and PLEN must be part of a structure
    //uint8_t plen; 
    //int k;

    FinimizerIndex() {}
    
    int get_k() const{return k;}


    void search(const std::string& query, unordered_map<int, uint64_t>& results) const {
  
        //std::cerr << "Searching " << query << std::endl;

        //const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        //const int64_t n_nodes = sbwt.number_of_subsets();
        //const int64_t k = sbwt.get_k();
        //const vector<int64_t>& C = sbwt.get_C_array();

        const int64_t query_len = query.length();
        // B is not empty!
      
        if (query.size() < this->k) return; 

        unordered_map<int64_t, uint64_t> Finimizers = rarest_fmin_streaming_search(query, this->B, this->sB, this->T, this->plen, this->k);
      
        // TODO Check the colors for every finimizer found
        //pseudoalignemnt_stats(Finimizers, this->hashTable, results);
        pseudoalignemnt_stats(Finimizers, this->C, results);
        //print_results(results);
        return;
    }

    // TODO check THIS ONCE THE ABOVE IS FIXED
    void search(const std::string& query, vector<pair<int, float>>& results, const float& t) const {
  
        /* const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();
        const vector<int64_t>& C = sbwt.get_C_array(); */
        
        const int64_t query_len = query.length();

        if (query.size() < this->k) return; 

        unordered_map<int64_t, uint64_t> Finimizers = rarest_fmin_streaming_search(query, this->B, this->sB, this->T, this->plen, this->k);

        // Check the colors for every finimizer found
        //pseudoalignemnt_stats(Finimizers, this->hashTable, results, t);
        pseudoalignemnt_stats(Finimizers, this->C, results, t);
        return;
    }

    // TODO remove
    /* void serialize_HashTable(const std::unordered_map<std::string, std::set<int>>& hashTable, const std::string& hashTableName) const{
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
    } */

    void serialize_sB(const std::unordered_map<uint32_t, int64_t >& sB, const std::string& sBName) const {
        std::ofstream sB_out(sBName, std::ios::binary);
        if (!sB_out) {
            std::cerr << "Error: Could not open file for writing!" << std::endl;
            return;
        }
    
        // Write the number of elements in the hash table
        size_t sBSize = sB.size();
        sB_out.write(reinterpret_cast<const char*>(&sBSize), sizeof(sBSize));
    
        // Write each key-value pair
        for (const auto& [key, value] : sB) {
            sB_out.write(reinterpret_cast<const char*>(&key), sizeof(key));
            sB_out.write(reinterpret_cast<const char*>(&value), sizeof(value));
        }
    
        sB_out.close();
    }

    void serialize_B(
        const std::unordered_map<uint32_t, std::pair<int64_t,int64_t>>& B, const std::string& BName) const{
        std::ofstream B_out(BName, std::ios::binary);
        if (!B_out) {
            std::cerr << "Error: Could not open file for writing!\n";
            return;
        }

        // Number of entries
        size_t BSize = B.size();
        B_out.write(reinterpret_cast<const char*>(&BSize), sizeof(BSize));

        // Write each key and its pair<first, second>
        for (const auto& [key, pr] : B) {
            int64_t v1 = pr.first;
            int64_t v2 = pr.second;
            B_out.write(reinterpret_cast<const char*>(&key), sizeof(key));
            B_out.write(reinterpret_cast<const char*>(&v1),  sizeof(v1));
            B_out.write(reinterpret_cast<const char*>(&v2),  sizeof(v2));
        }
        B_out.close();
    }

    void serialize_Colors(const std::vector<std::set<int>>& C, const std::string& filename) const {
        std::ofstream outFile(filename, std::ios::binary);
        if (!outFile) {
            std::cerr << "Error: Could not open colors file for writing!" << std::endl;
            return;
        }
    
        // Write the number of sets
        size_t numColors = C.size();
        outFile.write(reinterpret_cast<const char*>(&numColors), sizeof(numColors));
    
        // Write each set
        for (const auto& colorSet : C) {
            size_t setSize = colorSet.size();
            outFile.write(reinterpret_cast<const char*>(&setSize), sizeof(setSize));
    
            for (int val : colorSet) {
                outFile.write(reinterpret_cast<const char*>(&val), sizeof(val));
            }
        }
    
        outFile.close();
    }
 
    // TODO remove
    /* std::unordered_map<std::string, std::set<int>> load_HashTable(const std::string& hashTableName) {
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
 */

std::unordered_map<uint32_t, int64_t> load_sB(const std::string& sBName) {
    std::unordered_map<uint32_t, int64_t> sB;

    std::ifstream inFile(sBName, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error: Could not open file for reading!" << std::endl;
        return sB;
    }

    // Read the number of elements in the hash table
    size_t sBSize;
    if (!inFile.read(reinterpret_cast<char*>(&sBSize), sizeof(sBSize))) {
        std::cerr << "Error: Failed to read hash table size!" << std::endl;
        return sB;
    }

    for (size_t i = 0; i < sBSize; ++i) {
        uint32_t key;
        int64_t value;

        // Read key
        if (!inFile.read(reinterpret_cast<char*>(&key), sizeof(key))) {
            std::cerr << "Error: Failed to read key!" << std::endl;
            return sB;
        }

        // Read value
        if (!inFile.read(reinterpret_cast<char*>(&value), sizeof(value))) {
            std::cerr << "Error: Failed to read value!" << std::endl;
            return sB;
        }

        sB[key] = value;
    }

    inFile.close();
    return sB;
}


std::unordered_map<uint32_t, std::pair<int64_t,int64_t>> load_B(const std::string& BName) {
    std::unordered_map<uint32_t, std::pair<int64_t,int64_t>> B;

    std::ifstream inFile(BName, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error: Could not open file for reading!\n";
        return B;
    }

    // Read number of entries
    size_t BSize;
    if (!inFile.read(reinterpret_cast<char*>(&BSize), sizeof(BSize))) {
        std::cerr << "Error: Failed to read hash table size!\n";
        return B;
    }

    for (size_t i = 0; i < BSize; ++i) {
        uint32_t key;
        int64_t v1, v2;

        if (!inFile.read(reinterpret_cast<char*>(&key), sizeof(key))) {
            std::cerr << "Error: Failed to read key!\n";
            return B;
        }
        if (!inFile.read(reinterpret_cast<char*>(&v1), sizeof(v1))) {
            std::cerr << "Error: Failed to read first value!\n";
            return B;
        }
        if (!inFile.read(reinterpret_cast<char*>(&v2), sizeof(v2))) {
            std::cerr << "Error: Failed to read second value!\n";
            return B;
        }

        B[key] = {v1, v2};
    }
    inFile.close();
    return B;
}


    std::vector<std::set<int>> load_Colors(const std::string& filename) {
        std::vector<std::set<int>> C;
    
        std::ifstream inFile(filename, std::ios::binary);
        if (!inFile) {
            std::cerr << "Error: Could not open colors file for reading!" << std::endl;
            return C;
        }
    
        // Read number of sets
        size_t numColors;
        if (!inFile.read(reinterpret_cast<char*>(&numColors), sizeof(numColors))) {
            std::cerr << "Error: Failed to read number of color sets!" << std::endl;
            return C;
        }
        C.resize(numColors);
    
        for (size_t i = 0; i < numColors; ++i) {
            size_t setSize;
            if (!inFile.read(reinterpret_cast<char*>(&setSize), sizeof(setSize))) {
                std::cerr << "Error: Failed to read set size!" << std::endl;
                return C;
            }
    
            std::set<int> colorSet;
            for (size_t j = 0; j < setSize; ++j) {
                int val;
                if (!inFile.read(reinterpret_cast<char*>(&val), sizeof(val))) {
                    std::cerr << "Error: Failed to read color value!" << std::endl;
                    return C;
                }
                colorSet.insert(val);
            }
    
            C[i] = std::move(colorSet);
        }
    
        inFile.close();
        return C;
    }
    

    // TODO finimizerindex does not require sbwt nor LCS
    // TODO add plen and k
    void serialize(const string& index_prefix) const {

        // k and plen
        std::ofstream meta_out(index_prefix + ".meta", std::ios::binary);
        if (!meta_out) {
            std::cerr << "Error: Could not write metadata!" << std::endl;
            return;
        }
        meta_out.write(reinterpret_cast<const char*>(&k), sizeof(k));
        meta_out.write(reinterpret_cast<const char*>(&plen), sizeof(plen));
        meta_out.close();
        cerr << "k = " << k<< endl;
        /* std::ofstream LCS_out(index_prefix + ".LCS.sdsl");
        sdsl::serialize(*LCS, LCS_out);

        sbwt->serialize(index_prefix + ".sbwt"); */

        //serialize_HashTable(hashTable, index_prefix + ".ht.BIN");

        serialize_sB(sB, index_prefix + ".sB.BIN");

        serialize_B(B, index_prefix + ".B.BIN");

        //std::ofstream T_out(index_prefix + ".T.sdsl");
        //sdsl::serialize(*T.get(), T_out);
        std::ofstream T_out(index_prefix + ".T.sdsl", std::ios::binary);
        sdsl::serialize(T.size(), T_out); // Serialize the number of vectors
        for (const auto& vec : T) {
            sdsl::serialize(vec, T_out);
        }
        

        serialize_Colors(C, index_prefix + ".C.BIN");
    }

    // TODO finimizerindex does not require sbwt nor LCS
    // TODO add plen and k
    void load(const string& index_prefix) {
        // k and plen
        std::ifstream meta_in(index_prefix + ".meta", std::ios::binary);
        if (!meta_in) {
            std::cerr << "Error: Could not read metadata!" << std::endl;
            return;
        }
        meta_in.read(reinterpret_cast<char*>(&k), sizeof(k));
        meta_in.read(reinterpret_cast<char*>(&plen), sizeof(plen));
        meta_in.close();

        /* LCS = make_unique<sdsl::int_vector<>>();
        ifstream LCS_in(index_prefix + ".LCS.sdsl");
        sdsl::load(*LCS, LCS_in);
        std::cerr<< "LCS_file loaded"<<std::endl;

        sbwt = make_unique<plain_matrix_sbwt_t>();
        sbwt->load(index_prefix + ".sbwt");
        std::cerr << "SBWT matrix loaded" << std::endl;
        */
        /* hashTable=load_HashTable(index_prefix + ".ht.BIN");
        std::cerr << "hashTable loaded" << std::endl;
        //printHashTable(hashTable);
        */
        sB = load_sB(index_prefix + ".sB.BIN");
        std::cerr << "sB loaded" << std::endl;

        B = load_B(index_prefix + ".B.BIN");
        std::cerr << "B loaded" << std::endl;

        /* T = make_unique<sdsl::int_vector<1>>();
        ifstream T_in(index_prefix + ".T.sdsl", std::ios::binary);
        sdsl::load(*T, T_in); */
        std::ifstream T_in(index_prefix + ".T.sdsl", std::ios::binary);
        size_t num_vectors;
        sdsl::load(num_vectors, T_in);
        T.resize(num_vectors);
        for (auto& vec : T) {
            sdsl::load(vec, T_in);
        }
        std::cerr<< "Tails loaded"<<std::endl;

        C = load_Colors(index_prefix + ".C.BIN");
        std::cerr << "Colors loaded" << std::endl;
    }

    int64_t size_in_bytes() const {
        int64_t total = 0;

        total += sizeof(k);
        total += sizeof(plen);

        /* // LCS
        if (LCS) {total += sdsl::size_in_bytes(*LCS);}

        // SBWT 
        if (sbwt) {
            sbwt::SeqIO::NullStream ns;
            total += sbwt->serialize(ns);
        } */

        // T
        //if (T) {total += sdsl::size_in_bytes(*T);}
        for (const auto& vec : T) {
            total += sdsl::size_in_bytes(vec);
        }

        // B
        total += sizeof(std::pair<uint32_t, int64_t>) * B.size();
        total += sizeof(B); // Approximation for internal structure overhead

        // sB
        total += sizeof(std::pair<uint32_t, int64_t>) * sB.size();
        total += sizeof(sB);

        // C (vector<set<int>>)
        total += sizeof(C); // vector overhead
        for (const auto& s : C) {
            total += sizeof(std::set<int>);
            total += sizeof(int) * s.size(); // actual values
        }
        return total;
    }

};


class FinimizerIndexBuilder{
private:
    int k = 31;
    uint8_t plen = 10; // Prefix length

public:

    unique_ptr<plain_matrix_sbwt_t> sbwt;
    unique_ptr<sdsl::int_vector<>> LCS;


    //  TODO should B, sB, T, C and plen be part of this?
    /* // REAL DATA STR
    std::unordered_map<uint32_t, int64_t> B; // Create a hash table to store the prefixes of each bucket and a pointer to the start of the tails in the sdsl int vector
    std::unordered_map<uint32_t, int64_t> sB; // Create a hash table to store the finimizers shorter than the prefix length
    unique_ptr<int_vector<1>> T; // tails // width 0 so that I can decide the width and modify every entry
    vector<set<int>> C;
    uint8_t plen; // plen has to be given in input
 */

    unique_ptr<FinimizerIndex> index;

    int get_k() const{return k;}


    // Takes ownership of sbwt and LCS
    // TODO this should contain plen
    // k is now taken from the sbwt but must be linked to the new data str
    FinimizerIndexBuilder(unique_ptr<plain_matrix_sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, const vector<string>& incolors){ 
        index = make_unique<FinimizerIndex>();
        this->sbwt = move(sbwt); // Take ownership
        this->LCS = move(LCS); // Take ownership

        int64_t n_nodes = this->sbwt->number_of_subsets();

        // REAL DATA STR
        /* this->B = move(B);
        this->sB = move(sB);
        this->T = move(T);
        this->C = move(C);
        this->plen = move(plen); */

        std::unordered_map<uint32_t, pair<int64_t,int64_t>> B; // Create a hash table to store the prefixes of each bucket and a pointer to the start of the tails in the sdsl int vector
        std::unordered_map<uint32_t, int64_t> sB; // Create a hash table to store the finimizers shorter than the prefix length
        //int_vector<1> T; // width 0 so that I can decide the width and modify every entry
        std::vector<int_vector<1>> T;
        vector<set<int>> C;
        this->k = 31;
        //this->k = (int)index->sbwt->get_k(); //TODO FIX THIS

        //uint8_t plen; // plen has to be given in input

        //helpers
        std::unordered_map<uint32_t, std::map<uint8_t, set< pair<uint32_t, set<int> >> >> helperB;


        //TODO we are still using hashTable, but we could get rid of it
        std::unordered_map<std::string, std::set<int>> hashTable; // Create a hash table to store the list of colors for each Finimizer
        
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
        // TODO remove or modify to avoid using hash table
        // Extract statistics
        // get_stats(hashTable);

        // TODO 
        Buckets(hashTable, helperB, B, sB, C, this->plen);
        
        // TODO 
        storeTails(helperB, B, T, C, sB.size());
        

        // remove sbwt and LCS
        index->sbwt = std::move(this->sbwt); // Transfer ownership
        index->LCS = std::move(this->LCS); // Transfer ownership 
        //index->hashTable = std::move(hashTable); // Transfer ownership

        index->B = std::move(B); // Transfer ownership
        index->sB = std::move(sB); // Transfer ownership
        //index->T = std::move(T); // Transfer ownership
        index->T = std::move(T);
        index->C = std::move(C); // Transfer ownership

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
        cerr << "k= " << k << endl;
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
    
    

    void Buckets (unordered_map<std::string, std::set<int>>& hashTable,unordered_map<uint32_t, std::map<uint8_t, set< pair<uint32_t, set<int> >> >>& helperB, unordered_map<uint32_t, pair<int64_t,int64_t> >& B,  unordered_map<uint32_t, int64_t >& sB, vector<set<int>>& C, uint8_t plen){
        cerr << "Create buckets" << endl;
        // create a hash table with all the possible strings of length plen

        // helperB={prefix:{tlen1:{{tail1,colors1},...}, tlen2:{{tail1,colors1},...},... }}
        // helperB= prefix
        //             └─tlen
        //                └─ {tail,{colors}}
        
        // TODO THIS SHOULD BE IN INPUT OR DEPENDENT ON K
        //plen = 10;
        
        int64_t index_sp = 0;

        uint32_t nbuckets = 1<<(plen<<1);
        for (uint32_t p=0; p<nbuckets; p++){
            B[p]= {-1,-1};
            // INSERT EVERY POSSIBLE PREFIX IN B
            // helperB only contains the existing prefixes
        }
        for (auto &f : hashTable){
            string fmin = f.first;
            set<int> colors = f.second;
            uint8_t flen = fmin.length();
            if (flen < plen){ 
                uint32_t sfmin = prefix2int(fmin,0,flen);
                sB[sfmin]= index_sp++;  // ADD THE COLORS at the beginning of the colors vector<set<int>> and keep track of sB.size()
                C.push_back(colors);
            }
            else{
                uint32_t pfmin = prefix2int(fmin,0,plen);
                
                // Extract tails
                uint8_t tlen = flen - plen;
                
                if (tlen>0){
                    uint32_t tail = suffix2int(fmin,plen,tlen);// TODO deal with an EMPTY PREFIX
                    
                    if (helperB.find(pfmin) != helperB.end() and helperB[pfmin].find(tlen) != helperB[pfmin].end()) {
                        // The prefix is already there and also the correct length, add {tail, colors}
                            helperB[pfmin][tlen].insert({tail, colors}); 
                    } else {
                        // the prefix is not there yet, insert it!
                        // or the correct length is not there yet, insert it!
                        helperB[pfmin][tlen]={{tail, colors}};                     
                    }
                } else {
                    helperB[pfmin][tlen]={{0, colors}}; // here 0 will be interpreted as A*plen !!!!!!!!!!!!!!!
                }
            }
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

    // TODO deal with tlen=0
    void storeTails(unordered_map<uint32_t, std::map<uint8_t, set< pair<uint32_t, set<int> >> >>& helperB, std::unordered_map<uint32_t, pair<int64_t,int64_t>>& B, std::vector<int_vector<1>>& T, vector<set<int>>& C, int64_t sB_size){
        cerr << "Store Tails"<< endl;
        T.clear();
        
        // helperB={prefix:{tlen1:{{tail1,colors1},...}, tlen2:{{tail1,colors1},...},... }}
        // 1. Count the number of bits needed
/*         size_t total_bits = 0;
        for (auto &prefix : helperB){
            for (auto &tails: prefix.second){
                // this should be in tails order
                uint8_t tlen = tails.first; // 5 bits // tail length
                total_bits += 5; // tlen: 5 bits

                if (tlen > 0){
                    int32_t tnumber = tails.second.size(); // vbyte
                    auto vb = vbyte_encode(tnumber);
                    total_bits += vb.size() * 8; // vbyte: 8 bits per byte // number of tails
                    total_bits += tnumber * (2 * tlen); // actual tails 2bits/char
                }
            }
        }
        T = int_vector<1>(total_bits, 0); // is this what I want to do??
 */
        for (auto &prefix : helperB){
            size_t prefix_bits = 0;
            for (auto &tails: prefix.second){
                // this should be in tails order
                uint8_t tlen = tails.first; // 5 bits // tail length
                prefix_bits += 5; // tlen: 5 bits

                if (tlen > 0){
                    int32_t tnumber = tails.second.size(); // vbyte
                    auto vb = vbyte_encode(tnumber);
                    prefix_bits += vb.size() * 8; // vbyte: 8 bits per byte // number of tails
                    prefix_bits += tnumber * (2 * tlen); // actual tails 2bits/char
                }
            }
            T.push_back(int_vector<1>(prefix_bits, 0)); // bits for each prefix
        }


        // 2. Write data
        int64_t tails_so_far = sB_size;
        // offset >> tails_so_far

        int64_t i = 0;
        for (auto &prefix : helperB){        
            int64_t offset = 0;
            // add a pointer to the prefix in B
            B[prefix.first]={i, tails_so_far}; // store the index of the int_vector in T and not the offset anymore // we could store a pointer
            //cerr << "Current offset for prefix " << prefix.first << ": " << offset << endl;
            uint64_t* data = T[i].data();

            for (auto &tails: prefix.second){    // this should be in tails order
               // tails is std::pair<const char, std::set<std::pair<uint32_t, std::set<int>>>>
                // {tlen: [{tail1, {colors1}},..]}
                // tail length
                uint8_t tlen = tails.first;
                uint64_t word_index = offset/64;
                uint8_t w_offset = offset %64;
                sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
                
                // TODO remov this as it is useful only as a safety check
                uint8_t rxlen = sdsl::bits::read_int(&data[word_index], w_offset, 5); 
                if (rxlen != tlen){
                cerr << "Written tlen = " << (int)tlen << " at offset " << offset << endl;
                cerr << "rxlen at B[prefix].first: " << (int)rxlen << endl; // WRONG WHEN OFFSET IS A MULTIPLE OF 64
                }
                
                offset += 5;
                // deal with tlen=0
                if (tlen >0){
                    word_index = offset/64;
                    w_offset = offset %64;
                    // number of tails of length tlen
                    std::set<std::pair<uint32_t, std::set<int>>> tc = tails.second; //{tail1,{colors1}}
                    uint32_t tnumber = tc.size();
                    auto vb = vbyte_encode(tnumber);
                    for (uint8_t b : vb) {
                        sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
                        offset += 8;
                    }
                    tails_so_far += tnumber;

                    word_index = offset/64;
                    w_offset = offset %64;
                    for (const auto &p : tc){ //pair<uint32_t, set<int>
                        // tails
                        sdsl::bits::write_int(&data[word_index], p.first, w_offset, tlen*2);
                        offset += 2 * tlen;
                        
                        // 3. Write colors
                        // these are in the same order in which we are writing tails
                        // can be accessed with the offset obtained by searching tails (REMEMBER TO ADD sB.size()) // sB.size() has now been added to B directly
                        C.push_back(p.second);
                    }
                }
                else{
                    tails_so_far++;
                }
                 // if tlen=0, no need to write anything
                  //if the first 5 bits are zero you should know you are done // HOW DO YOU KNOW!??!?!?!?!   
            }
            i++;       
        }        
    }    

    // Transfer ownership of the index out of the builder
    unique_ptr<FinimizerIndex> get_index(){
        return std::move(this->index);
    }

};
