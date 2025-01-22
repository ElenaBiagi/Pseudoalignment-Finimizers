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

    struct QueryResult{
        vector<pair<int64_t, int64_t>> local_offsets; // Unitig id, distance from the start of the unitig
        int64_t n_found = 0;
    };

private:
    // Forbid copying because we have pointers to our internal data structures
    FinimizerIndex(const FinimizerIndex& other) = delete;
    FinimizerIndex& operator=(const FinimizerIndex& other) = delete;

    void add_to_query_result(int64_t global_kmer_end, QueryResult& answer) const{
        int64_t global_kmer_start = global_kmer_end - sbwt->get_k() + 1;
        pair<int64_t, int64_t> local_start = unitigs.global_offset_to_local_offset(global_kmer_start);
        answer.local_offsets.push_back(local_start);
        answer.n_found++;
    }

    void walk_in_unitigs(const std::string& query, const PackedStrings& unitigs, int64_t global_kmer_end, QueryResult& answer, int64_t& kmer_end, const int64_t k) const{
        //cout << "take a nice unitig walk" << endl;
        int64_t unitig_id = answer.local_offsets.back().first;
        int64_t u_end = unitigs.ends[unitig_id]; // exlusive end
        int64_t max_match = std::min(u_end - global_kmer_end-1, (int64_t)(query.length()-kmer_end-1));

        if (global_kmer_end > u_end or max_match <= 0){return;}
        kmer_end++;
        //convert part of the query to an int_vector<2>
        sdsl::int_vector<2> query_v(max_match);
        
        for (int64_t i = 0; i < max_match;) {
            char c = query[kmer_end+i];
            switch(c){
                case 'A': query_v[i++] = 0; break;
                case 'C': query_v[i++] = 1; break;
                case 'G': query_v[i++] = 2; break;
                case 'T': query_v[i++] = 3; break;
                default: throw std::runtime_error("Invalid character: " + c);
            }
        }

        int64_t word_start = 0;
        int64_t global_kmer_end_copy = global_kmer_end;

        while(max_match > 0){
            int word_len = std::min({(int64_t)32, max_match});
            uint64_t query_word = query_v.get_int(word_start, (word_len)*2); // word start = least significant bit (word_len+1)*2-1
            uint64_t unitig_word = unitigs.concat.get_int(word_start + (global_kmer_end_copy+1)*2, (word_len)*2);
            //cerr << std::bitset<8 * sizeof(int64_t)>(query_word) << endl;
            //cerr << std::bitset<8 * sizeof(int64_t)>(unitig_word) << endl;


            int64_t result = query_word ^ unitig_word;
            //cerr << std::bitset<8 * sizeof(int64_t)>(result) << endl;

            if (result){
                int trailing_zeros = __builtin_ctzll(result);
                for (int i = 0; i < trailing_zeros/2; i++){
                    global_kmer_end++;
                    add_to_query_result(global_kmer_end, answer);
                    kmer_end++; // keep track of how many kmers were found    
                }
                break; // No need to check further. A mismatch has been found.
            }
            int trailing_zeros = (word_len)*2;
            for (int i = 0; i < trailing_zeros/2; i++){
                global_kmer_end++;
                add_to_query_result(global_kmer_end, answer);
                kmer_end++; // keep track of how many kmers were found    
            }
            max_match -= word_len;
            word_start += word_len*2;
        }
        kmer_end--; // will be updated later
    }


public:

    // Note: if you add members, update size_in_bytes(), serialize(), and load()
    unique_ptr<plain_matrix_sbwt_t> sbwt; // These are smart pointers because they are passed in to the constructor
    unique_ptr<sdsl::int_vector<>> LCS; // These are smart pointers because they are passed in to the constructor
    PackedStrings unitigs;
    sdsl::bit_vector fmin;
    sdsl::rank_support_v5<> fmin_rs;
    sdsl::int_vector<> global_offsets;
    sdsl::bit_vector Ustart;
    sdsl::rank_support_v5<> Ustart_rs;
    std::unordered_map<std::string, std::unordered_set<int>> hashTable; // Create a hash table to store the list of colors for each Finimizer


    FinimizerIndex() {}

    QueryResult search(const std::string& query) const {
        // For each k-mer S that is known to be in the SBWT
        //   - Find the finimizer x.
        //   - Walk forward in the SBWT from the colex rank of x (singleton interval), to the end of S,
        //     recording the rightmost branch point, and the distance from the end of x to this branch
        //     point.
        //   - If a branch point was found, the k-mer containing x is at the unitig after the rightmost
        //     branch. Otherwise, the k-mer containing x is at the unitig that x points to.
        //   - If there was a branch, the k- mer endpoint in that unitig is k + the number of steps taken after the last branch.
        //     

        const plain_matrix_sbwt_t& sbwt = *(this->sbwt.get());
        const int64_t n_nodes = sbwt.number_of_subsets();
        const int64_t k = sbwt.get_k();
        const vector<int64_t>& C = sbwt.get_C_array();
        const int64_t query_len = query.length();


        QueryResult answer{{}, 0};

        if(query.size() < k) answer; 


        //Find kmers and Finimizers together
        vector<string> Finimizers = rarest_fmin_streaming_search(sbwt, *LCS, query);
        //TODO Check the colors for every finimizer found
        pseudoalignemnt_stats(Finimizers, hashTable);

        return answer;
    }

    void serialize_HashTable(const std::unordered_map<std::string, std::unordered_set<int>>& hashTable, const std::string& hashTableName) const{
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

            // Write the size of the value vector and the vector itself
            size_t valueSize = value.size();
            hashTable_out.write(reinterpret_cast<const char*>(&valueSize), sizeof(valueSize));
            for (const int elem : value) {
                hashTable_out.write(reinterpret_cast<const char*>(&elem), sizeof(elem));
            }
        }

        hashTable_out.close();
    }

    std::unordered_map<std::string, std::unordered_set<int>> load_HashTable(const std::string& hashTableName) {
        std::unordered_map<std::string, std::unordered_set<int>> hashTable;

        std::ifstream inFile(hashTableName, std::ios::binary);
        if (!inFile) {
            std::cerr << "Error: Could not open file for reading!" << std::endl;
            return hashTable;
        }

        // Read the number of elements in the hash table
        size_t hashTableSize;
        inFile.read(reinterpret_cast<char*>(&hashTableSize), sizeof(hashTableSize));

        for (size_t i = 0; i < hashTableSize; ++i) {
            // Read the key
            size_t keySize;
            inFile.read(reinterpret_cast<char*>(&keySize), sizeof(keySize));
            std::string key(keySize, '\0');
            inFile.read(&key[0], keySize);

            // Read the value vector
            size_t valueSize;
            inFile.read(reinterpret_cast<char*>(&valueSize), sizeof(valueSize));
            std::unordered_set<int> value(value);
            for (size_t j = 0; j < valueSize; ++j) {
                int elem;
                inFile.read(reinterpret_cast<char*>(&elem), sizeof(elem));
                value.insert(elem);
            }

        hashTable[key] = value;

            // insert the key-value into the hash table
            hashTable[key] = value;
        }

        inFile.close();
        //std::cout << "Hash table loaded from " << hashTableName << std::endl;
        return hashTable;
    }

    // TODO: add hash table
    void serialize(const string& index_prefix) const {
        std::ofstream global_offsets_out(index_prefix + ".O.sdsl");
        sdsl::serialize(global_offsets, global_offsets_out);

        std::ofstream fmin_out(index_prefix + ".FBV.sdsl");
        sdsl::serialize(fmin, fmin_out);

        std::ofstream packed_unitigs_out(index_prefix + ".packed_unitigs.sdsl");
        sdsl::serialize(unitigs.concat, packed_unitigs_out);
        
        std::ofstream unitig_endpoints_out(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::serialize(unitigs.ends, unitig_endpoints_out);

        std::ofstream Ustart_out(index_prefix + ".Ustart.sdsl");
        sdsl::serialize(Ustart, Ustart_out);

        std::ofstream LCS_out(index_prefix + ".LCS.sdsl");
        sdsl::serialize(*LCS, LCS_out);

        sbwt->serialize(index_prefix + ".sbwt");

        serialize_HashTable(hashTable, index_prefix + "ht.BIN");
    }

    void load(const string& index_prefix) {

        LCS = make_unique<sdsl::int_vector<>>();
        ifstream LCS_in(index_prefix + ".LCS.sdsl");
        sdsl::load(*LCS, LCS_in);
        std::cerr<< "LCS_file loaded"<<std::endl;

        ifstream fmin_bv_in(index_prefix + ".FBV.sdsl");
        sdsl::load(fmin, fmin_bv_in);
        std::cerr<< "fmin_bv_file loaded"<<std::endl;
        sdsl::util::init_support(fmin_rs, &fmin);

        ifstream global_offsets_in(index_prefix + ".O.sdsl");
        sdsl::load(global_offsets, global_offsets_in);
        std::cerr<< "offsets loaded"<<std::endl;

        std::ifstream packed_unitigs_in(index_prefix + ".packed_unitigs.sdsl");
        sdsl::load(unitigs.concat, packed_unitigs_in);
        std::cerr << "unitigs loaded" << std::endl;

        std::ifstream unitig_endpoints_in(index_prefix + ".unitig_endpoints.sdsl");
        sdsl::load(unitigs.ends, unitig_endpoints_in);
        std::cerr << "unitig endpoints loaded" << std::endl;

        std::ifstream Ustart_in(index_prefix + ".Ustart.sdsl");
        sdsl::load(Ustart, Ustart_in);
        sdsl::util::init_support(Ustart_rs, &Ustart);
        std::cerr << "Ustart loaded" << std::endl;

        sbwt = make_unique<plain_matrix_sbwt_t>();
        sbwt->load(index_prefix + ".sbwt");
        std::cerr << "SBWT matrix loaded" << std::endl;

        load_HashTable(index_prefix + "ht.BIN");
        std::cerr << "hashTable loaded" << std::endl;
    }

    // TODO: add hash table (later. not relevant now)
    // This also includes the rank structures which are not serialized
    int64_t size_in_bytes() const{
        int64_t total = 0;
        total += sdsl::size_in_bytes(*LCS);
        total += sdsl::size_in_bytes(fmin);
        total += sdsl::size_in_bytes(fmin_rs);
        total += sdsl::size_in_bytes(global_offsets);
        total += sdsl::size_in_bytes(unitigs.concat);
        total += sdsl::size_in_bytes(unitigs.ends);
        total += sdsl::size_in_bytes(Ustart);
        total += sdsl::size_in_bytes(Ustart_rs);

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
        sdsl::bit_vector fmin_bv(n_nodes, 0); // Finimizer marks
        //sdsl::bit_vector fmin_found(n_nodes, 0);
        sdsl::int_vector fmin_found(n_nodes, 0);
        
        vector<uint64_t> global_offsets;
        global_offsets.reserve(n_nodes);
        global_offsets.resize(n_nodes, 0);

        std::unordered_map<std::string, std::unordered_set<int>> hashTable; // Create a hash table to store the list of colors for each Finimizer

        
        pair<PackedStrings, sdsl::bit_vector> unitig_data = permute_unitigs(*(this->sbwt), reader);
        PackedStrings& unitigs = unitig_data.first;
        sdsl::bit_vector& Ustart = unitig_data.second;

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

        sdsl::int_vector<> packed_global_offsets(finimizers.size(), 0, 64 - __builtin_clzll(*std::max_element(global_offsets.begin(), global_offsets.end())));
        
        int64_t global_offsets_idx = 0;
        for(int64_t i = 0; i < global_offsets.size(); i++){
            if(fmin_bv[i]) packed_global_offsets[global_offsets_idx++] = global_offsets[i];
        }
    
        print_finimizer_stats(finimizers, this->sbwt->number_of_kmers(), this->sbwt->number_of_subsets(), 1);

        // Scan the genomes to get the list of colors for each finimizer using the hash table
        typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_colors_gzip;
        typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_colors_no_gzip;
        for (int64_t i=0; i< incolors.size(); i++){
            std::cerr << "scanning color " << i << std::endl;

            bool gzip_colors = SeqIO::figure_out_file_format(incolors[i]).gzipped;

            if (gzip_colors){
                run_colors_file<in_colors_gzip>(incolors[i], hashTable, i);
                //scan_color<in_colors_gzip>(incolors[i], hashTable, i);
            }
            else{
                std::cerr << "running colors not gizipped" << std::endl;
                run_colors_file<in_colors_no_gzip>(incolors[i], hashTable, i);
                //scan_color<in_colors_no_gzip>(incolors[i], hashTable, i);
            }
            std::cerr << "DONE"<< std::endl;
            // TODO add reverse complement
            //const string reverse = sbwt::get_rc(incolors[i]);
            //std::cerr << "string reversed" << std::endl;
            //scan_color(reverse, hashTable, i);
        }

        index->sbwt = std::move(this->sbwt); // Transfer ownership
        index->LCS = std::move(this->LCS); // Transfer ownership 
        index->unitigs = std::move(unitigs); // Transfer ownership
        index->fmin = std::move(fmin_bv); // Transfer ownership
        index->fmin_rs = sdsl::rank_support_v5<>(&(index->fmin));
        index->global_offsets = std::move(packed_global_offsets); // Transfer ownership
        index->Ustart = std::move(Ustart); // Transfer ownership
        index->Ustart_rs = sdsl::rank_support_v5<>(&(index->Ustart));
        index->hashTable = std::move(hashTable);

    }

    // TODO fix return type
    template<typename reader_t>
    int64_t run_colors_file(const string& infile, unordered_map<std::string, std::unordered_set<int>> hashTable, int64_t i){
        reader_t reader(infile);
        write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
        return from_reader_to_seq(reader, hashTable, i);
    }

    set<tuple<int64_t, int64_t, int64_t>> add_sequence(const std::string& seq, sdsl::bit_vector& fmin_bv, sdsl::int_vector<>& fmin_found, vector<uint64_t>& global_offsets, const int64_t unitig_start, unordered_map<std::string, std::unordered_set<int>> hashTable) {
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

    // TODO fix return type
    // TODO remove hashTable
    template<typename reader_t>
    int from_reader_to_seq(reader_t& reader, unordered_map<std::string, std::unordered_set<int>> hashTable, int64_t i) {
        int64_t j =0;
        while(true){
            int64_t len = reader.get_next_read_to_buffer();
            if(len == 0) [[unlikely]] break;

            const std::string& seq =reader.read_buf;
            scan_color(seq, hashTable, j);
            std::cerr << seq<< std::endl;

            const string reverse = sbwt::get_rc(reader.read_buf);
            scan_color(reverse, hashTable, j);
            j++;
        }
        return 1;
    }
    // TODO fix return type
    // TODO remove hashTable
    //template<typename reader_t>
    int scan_color(const std::string& seq, unordered_map<std::string, std::unordered_set<int>> hashTable, int64_t i) {
        // this is the same as add_sequence but with genomes(colors) instead of unitigs

        std::cerr << seq<< std::endl;
        const int64_t n_nodes = sbwt->number_of_subsets();
        const int64_t k = sbwt->get_k();
        const vector<int64_t>& C = sbwt->get_C_array();
        int64_t freq;
        BoundedDeque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin(seq.size());
        const int64_t str_len = seq.size();
        tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len}; // {freq, len, I start, start}

        int64_t kmer = 0;
        int64_t start = 0;
        int64_t end;
        pair<int64_t, int64_t> I = {0, n_nodes - 1};
        int64_t I_start;
        tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;
        char c;
        char char_idx;
        
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
                // add the color to the finimizer
                hashTable[seq.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin))].insert(i);
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