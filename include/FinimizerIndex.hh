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
    // PackedStrings unitigs;
    std::unordered_map<std::string, std::set<int>> hashTable; // Create a hash table to store the list of colors for each Finimizer


    FinimizerIndex() {}

    //QueryResult 
    void search(const std::string& query, vector<pair<int,uint64_t>>& results, set<int>& intersection, const float& t) const {
  
        //std::cerr << "Searching " << query << std::endl;

        const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();
        const vector<int64_t>& C = sbwt.get_C_array();
        const int64_t query_len = query.length();


        if (query.size() < k) return; 

        //TODO this is still excluding some kmers if a subtring is not found
        vector<string> Finimizers = rarest_fmin_streaming_search(sbwt, *LCS, query);
        // Check the colors for every finimizer found
        pseudoalignemnt_stats(Finimizers, this->hashTable, results, intersection, t);

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
        /* 
        std::ofstream packed_unitigs_out(index_prefix + ".packed_unitigs.sdsl");
        sdsl::serialize(unitigs.concat, packed_unitigs_out);
        
        std::ofstream unitig_endpoints_out(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::serialize(unitigs.ends, unitig_endpoints_out);
        */

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

        /*
        std::ifstream packed_unitigs_in(index_prefix + ".packed_unitigs.sdsl");
        sdsl::load(unitigs.concat, packed_unitigs_in);
        std::cerr << "unitigs loaded" << std::endl;

        std::ifstream unitig_endpoints_in(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::load(unitigs.ends, unitig_endpoints_in);
        std::cerr << "unitig endpoints loaded" << std::endl;
        */

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
        /* 
        total += sdsl::size_in_bytes(unitigs.concat);
        total += sdsl::size_in_bytes(unitigs.ends);
        */

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
    template<typename reader_t>
    FinimizerIndexBuilder(unique_ptr<plain_matrix_sbwt_t> sbwt, unique_ptr<sdsl::int_vector<>> LCS, reader_t& reader, const vector<string>& incolors) {
        index = make_unique<FinimizerIndex>();
        this->sbwt = move(sbwt); // Take ownership
        this->LCS = move(LCS); // Take ownership

        int64_t n_nodes = this->sbwt->number_of_subsets();

        std::unordered_map<std::string, std::set<int>> hashTable; // Create a hash table to store the list of colors for each Finimizer

        
        /* 
        //TODO if we want to keep this we must delete everything that is not necessary 
        pair<PackedStrings, sdsl::bit_vector> unitig_data = permute_unitigs(*(this->sbwt), reader);
        PackedStrings& unitigs = unitig_data.first;

        set<tuple<int64_t, int64_t, int64_t>>  finimizers;
        int64_t total_len = 0;
        vector<char> unitig_buf;
        for(int64_t i = 0; i < unitigs.number_of_strings(); i++){
            int64_t len = unitigs.get(i, unitig_buf);
            // we cannot iterate here and get finimizers colors as these are unitigs and not genomes
            set<tuple<int64_t, int64_t, int64_t>> new_search = add_sequence(unitig_buf.data(), fmin_bv, fmin_found, global_offsets, total_len, hashTable);
            total_len += len; 
            finimizers.insert(new_search.begin(), new_search.end());
        }
     
        print_finimizer_stats(finimizers, this->sbwt->number_of_kmers(), this->sbwt->number_of_subsets(), 1);
*/

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

        // TODO extact statistics
        get_stats(hashTable);
        
        index->sbwt = std::move(this->sbwt); // Transfer ownership
        index->LCS = std::move(this->LCS); // Transfer ownership 
        /* index->unitigs = std::move(unitigs); // Transfer ownership */
        index->hashTable = std::move(hashTable);
    }

    // TODO fix return type
    template<typename reader_t>
    int64_t run_colors_file(const string& infile, unordered_map<std::string, std::set<int>>& hashTable, int& i){
        reader_t reader(infile);
        write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
        return from_reader_to_seq(reader, hashTable, i);
    }

    /* set<tuple<int64_t, int64_t, int64_t>> add_sequence(const std::string& seq, sdsl::bit_vector& fmin_bv, sdsl::int_vector<>& fmin_found, vector<uint64_t>& global_offsets, const int64_t unitig_start, unordered_map<std::string, std::set<int>>& hashTable) {
        const int64_t n_nodes = sbwt->number_of_subsets();
        const int64_t k = sbwt->get_k();
        const vector<int64_t>& C = sbwt->get_C_array();
        int64_t freq;
        BoundedDeque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin(seq.size());
        const int64_t str_len = seq.size();
        tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, start}
        set<tuple<int64_t,int64_t, int64_t>> count_all_w_fmin;

        int64_t kmer = 0;
        int64_t start = 0;
        int64_t end;
        pair<int64_t, int64_t> I = {0, n_nodes - 1};
        int64_t I_start;
        tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
        char c;
        char char_idx;
        // the idea is to start from the first pos which is i and move until finding something of ok freq
        // then drop the first char keeping track of which char you are starting from
        // Start is always < k as start <= end and end <k
        // if start == end than the frequency higher than t
        for (end = 0; end < str_len; end++) {
            c = static_cast<char>(seq[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
            //update the sbwt INTERVAL
            I = this->sbwt->update_sbwt_interval(&c, 1, I);
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
                count_all_w_fmin.insert({get<1>(w_fmin),get<0>(w_fmin), get<2>(w_fmin) });// (length,freq,colex) freq = 1 thus == (freq, length,colex)
                // TODO need to insert input string back!!
                //Finimizer = input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin))

                // Save the finimizer in the hash table if it's the first time you encounter it
                if (fmin_found[get<2>(w_fmin)] == 0){
                    hashTable[seq.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin))] = {};
                }

                if (fmin_found[get<2>(w_fmin)] == 0 or fmin_found[get<2>(w_fmin)]< get<3>(w_fmin)){ // if the finimizer has been found before in a full kmer
                    fmin_bv[get<2>(w_fmin)]=1;
                    fmin_found[get<2>(w_fmin)] = get<3>(w_fmin);

                    if ((unitig_start + get<3>(w_fmin))> UINT64_MAX){
                        std::cerr<< "ISSUE: global offset exceedes the allowed bit range." << std::endl;
                    }
                    global_offsets[get<2>(w_fmin)]= unitig_start + get<3>(w_fmin);
                }
                // write_fasta({input.substr(kmer,k) + ' ' + to_string(get<0>(w_fmin)),input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin))},writer);
                kmer++;
                // Check if the current minimizer is still in this window
                while (get<3>(w_fmin)- get<1>(w_fmin)+1 < kmer) { // start
                    all_fmin.pop_front();
                    w_fmin = (all_fmin.size()==0) ? tuple<int64_t, int64_t, int64_t, int64_t> {n_nodes,k+1,kmer+1,kmer+k} : all_fmin.front();
                }
            }
        }
        return count_all_w_fmin;
    }
 */
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
        // this is the same as add_sequence but with genomes(colors) instead of unitigs
        //std::cerr << "i= " << i << endl;
        //std::cerr << seq<< std::endl;

        /* std::unordered_set<char> distinct_chars;
        // Insert each character into the set
        for (char c : seq) { distinct_chars.insert(c);}
        // Output the count of distinct characters
        if (distinct_chars.size() > 4){
            std::cerr << "Number of distinct characters: " << distinct_chars.size() << std::endl;
            return 0;
        } */

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
/*             int64_t char_idx = get_char_idx(c);
            if (char_idx == -1) [[unlikely]]{
                cerr << "Error: unknown character: " << c << endl;
                cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
                return {};
            } */
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

                    // Check if the current minimizer is still in this window
                    while (get<3>(w_fmin)- get<1>(w_fmin)+1 < kmer) { // start
                         all_fmin.pop_front();
                        w_fmin = (all_fmin.size()==0) ? tuple<int64_t, int64_t, int64_t, int64_t> {n_nodes,k+1,kmer+1,kmer+k} : all_fmin.front();
                    }
                }
        } 
        return 1;
    }

    // Transfer ownership of the index out of the builder
    unique_ptr<FinimizerIndex> get_index(){
        return std::move(this->index);
    }
};