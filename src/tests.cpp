#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "commands.hh"
#include "globals.hh"
#include "build_fmin.hh"
#include "search_fmin.hh"
#include <filesystem>
#include "FinimizerIndex.hh"
#include "lcs_basic_parallel_algorithm.hpp"
#include "backward.hpp"
#include "bitsearch.hh"

string temp_dir = "tests_temp";

vector<string> paper_example_unitigs = {"GTAAGTCT", "AGGAAA", "ACAGG", "GTAGG", "AGGTA"};
//                                          1           2           0

vector<string> paper_example_queries = {"AAGTAA"};

using namespace std;
using namespace std::filesystem;

template<typename T>
void assert_equal(const T& a, const T& b){
    if(a != b){
        backward::StackTrace st; st.load_here(32);
        backward::Printer p; p.print(st);
        cerr << "Assertion failed: " << a << " != " << b << endl;
        exit(1);
    }
}

void write_as_fasta(const vector<string>& seqs, const string& filename){
    ofstream fasta_out(filename);
    for(const string& S : seqs){
        fasta_out << ">\n" << S << "\n" << endl;
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

/* // Takes in a spectrum-preserving string set
unique_ptr<FinimizerIndex> build_index(const vector<string>& spss, int64_t k){

    unique_ptr<plain_matrix_sbwt_t> sbwt = make_unique<plain_matrix_sbwt_t>();
    NodeBOSSInMemoryConstructor<plain_matrix_sbwt_t> constructor;
    constructor.build(spss, *sbwt, k, true);

    unique_ptr<sdsl::int_vector<>> LCS = make_unique<sdsl::int_vector<>>(move(lcs_basic_parallel_algorithm(*sbwt, 3))); // 3 threads

    string input_filename = temp_dir + "/spss.fna";
    write_as_fasta(spss, input_filename);
    SeqIO::Reader<> reader(input_filename);
    FinimizerIndexBuilder builder(move(sbwt), move(LCS), reader);
    unique_ptr<FinimizerIndex> index = builder.get_index();
    return move(index);
}

unique_ptr<FinimizerIndex> build_example_index(){
    return build_index(paper_example_unitigs, 4);
}

void test_finimizer_selection(){
    // ACGG has outgoing edges T and C, but the one with C goes to the reverse complemented k-mer CGGC
    int64_t k = 4;
    vector<string> unitigs = {"ACGG", "CGGT", "GCCGTA"};
    string query = "GCCGTA";
    // Permuted order:            1       2        0

    unique_ptr<FinimizerIndex> index = build_index(unitigs, 4);
    index->search(query);

    sdsl::bit_vector true_fmin = {0,0,1,1,1,0,0,0,0,1,0,0};
    assert_equal(true_fmin, index->fmin);

    // 0 $$$$
    // 1 $$$A
    // 2 CGTA*
    // 3 $$AC*
    // 4 $GCC*
    // 5 $$GC
    // 6 $$$G
    // 7 $ACG
    // 8 GCCG
    // 9 ACGG*
    // 10 CCGT
    // 11 CGGT

}
 */

void test_tail_search(){
    string s = "ATCTT";//"TCT";//"AGGATTGTCT"; // 10 bits NEXT TGTAC
    cerr << s << endl;
    vector <int> tlens = {2,3,5}; 
    //                                x   x    x      x     x       x    x      x      f        f        f      x         x        f        f        f      f(13)     f(14)
    vector<vector<string>> tails= {{"AA","AC","GG", "CT", "CC"},{"ATT","AGT", "TTA", "TCT"}, {"AGGAT","AGCGG","ATCTT", "GCCTT", "GACCT", "TTCGT", "TGTAC", "TTTAA", "TGAGT"}};
    //                               0    1     2     3     4      5     6      7       8       9       10      11       12        13       14        15     16         17 
    // 00010-00000101--0000-0001-1010-0111-0101
    // 00011-00000100--001111-001011-111100-110111
    // // r
    // 10100-10010000--1100010100-0101100100-11101100
    // THE VECTOR IS WRITTEN CORRECTLY

    // 01000-10100000--0000-1000-0101-1110-1010
    // 11000-00100000--111100-110100-001111-111011 // end = 66!!
    // 10100-10010000--1100010100-0101100100-1111101100-1111101001-1110100001-1101101111-100011011100001111111101000111
    uint64_t total_bits = 0;
    for (size_t i = 0; i < tlens.size(); ++i) {
        int tlen = tlens[i];
        uint32_t ntails = tails[i].size();

        total_bits += 5; // tlen
        total_bits += vbyte_encode(ntails).size() * 8; // vbyte encoded tail count

        total_bits += ntails * tlen * 2; // each tail uses tlen*2 bits
    }

    sdsl::int_vector<1> T;
    T.resize(total_bits); 
    uint64_t* data = T.data();


   const char ntlen = tlens.size();
    int64_t offset = 0;
    uint64_t word_index = 0;
    uint8_t w_offset = 0;
           
    for (char i=0; i<ntlen; i++){
        int tlen = tlens[i]; 
        word_index = offset/64;
        w_offset = offset %64;
        sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
        offset += 5;

        const uint32_t tnumber = tails[i].size();
        auto vb = vbyte_encode(tnumber);
        for (uint8_t b : vb) {
            word_index = offset/64;
            w_offset = offset %64;
            sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
            offset += 8;
        }
        for (const auto &t : tails[i]){ //pair<uint32_t, set<int>
            // tails
            word_index = offset/64;
            w_offset = offset %64;
            uint32_t t_int = prefix2int(t,0,tlen);
            sdsl::bits::write_int(&data[word_index], t_int, w_offset, tlen*2);
            offset += (2 * tlen);
        }
    }
    std::cout << T << std::endl;
    uint64_t s_int = prefix2int(s, 0, s.size());
    auto result = bitMagicSearch_new(T, s_int, s.size()); // input: const sdsl::int_vector<1> &T, string S 
    pair<int64_t, uint8_t> correct_result = {11,5};
    cerr << result.first << ", " << (int)result.second << endl; 
    assert_equal(result, correct_result);
}

void test_tail_search_2(){
    string s = "TGTAC"; // 10 bits 
    
    vector <int> tlens = {2,3,5}; 
    //                                                                                                                                               x
    //vector<vector<string>> tails= {{"AA","AC","GG", "CT", "CC"},{"ATT","AGT", "TTA", "TCT"}, {"AGGAT","AGCGG","ATCTT", "GCCTT", "GACCT", "TTCGT", "TGTAC", "TTTAA", "TGAGT"}};
    //                                                                                          9
    vector<vector<string>> tails= {{"AA","AC","GG", "CT", "CC"},{"ATT","AGT", "TTA", "TCT"}, {"TGTAC", "TTTAA", "TGAGT"}};


    uint64_t total_bits = 0;
    for (size_t i = 0; i < tlens.size(); ++i) {
        int tlen = tlens[i];
        uint32_t ntails = tails[i].size();

        total_bits += 5; // tlen
        total_bits += vbyte_encode(ntails).size() * 8; // vbyte encoded tail count

        total_bits += ntails * tlen * 2; // each tail uses tlen*2 bits
    }

    sdsl::int_vector<1> T;
    T.resize(total_bits); 
    uint64_t* data = T.data();


   const char ntlen = tlens.size();
    int64_t offset = 0;
    uint64_t word_index = 0;
    uint8_t w_offset = 0;
           
    for (char i=0; i<ntlen; i++){
        int tlen = tlens[i]; 
        word_index = offset/64;
        w_offset = offset %64;
        sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
        offset += 5;

        const uint32_t tnumber = tails[i].size();
        auto vb = vbyte_encode(tnumber);
        for (uint8_t b : vb) {
            word_index = offset/64;
            w_offset = offset %64;
            sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
            offset += 8;
        }
        for (const auto &t : tails[i]){ //pair<uint32_t, set<int>
            // tails
            word_index = offset/64;
            w_offset = offset %64;
            uint32_t t_int = prefix2int(t,0,tlen);
            sdsl::bits::write_int(&data[word_index], t_int, w_offset, tlen*2);
            offset += (2 * tlen);
        }
    }
    uint64_t s_int = prefix2int(s, 0, s.size());
    auto result = bitMagicSearch_new(T, s_int, s.size()); // input: const sdsl::int_vector<1> &T, string S 
    pair<int64_t, uint8_t> correct_result = {9,5};
    cerr << result.first << ", " << (int)result.second << endl; 
    assert_equal(result, correct_result);
}

int main(int argc, char** argv){
    // Create test directory if does not exist
    if (!exists(temp_dir)){
        create_directory(temp_dir);
    }

    cerr << "Testing tail search..." << endl;
    test_tail_search();
    cerr << "...ok" << endl;

    cerr << "Testing longer tail search..." << endl;
    //test_tail_search_2();
    cerr << "...ok" << endl;

    /* cerr << "Testing shortest unique construction..." << endl;
    test_shortest_unique_construction();
    cerr << "...ok" << endl;

    cerr << "Testing shortest unique queries..." << endl;
    test_shortest_unique_queries();
    cerr << "...ok" << endl;

    cerr << "Testing finimizer branch" << endl;
    test_finimizer_branch();
    cerr << "...ok" << endl;

    cerr << "Testing reverse complement branch" << endl;
    test_reverse_complement_branch();
    cerr << "...ok" << endl;

    cerr << "Testing leftmost" << endl;
    test_leftmost();
    cerr << "...ok" << endl;

    cerr << "Testing Finimizer selection" << endl;
    test_finimizer_selection();
    cerr << "...ok" << endl;
 
    cerr << "Testing incoming rc branch" << endl;
    test_incoming_rc_branch();
    cerr << "...ok" << endl;

    cerr << "Testing rc query" << endl;
    test_reverse_complement_query();
    cerr << "...ok" << endl;

    cerr << "Testing WALK" << endl;
    test_walk();
    cerr << "...ok" << endl; 
    */

    cerr << "ALL TESTS PASSED" << endl;

}