#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "commands.hh"
#include "globals.hh"
#include "build_fmin.hh"
#include "search_fmin.hh"
#include <filesystem>
#include "backward.hpp"
#include "bitsearch.hh"
#include "ColoredFinimizers.hh"

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

//TODO TEST EVERYTHING


void test_tail_search(string& s, pair<int64_t, uint8_t>& correct_result, sdsl::int_vector<1> T){
    cerr << s << endl;
    // 00010-00000101--0000-0001-1010-0111-0101
    // 00011-00000100--001111-001011-111100-110111
    // // r
    // 10100-10010000--1100010100-0101100100-11101100
    // THE VECTOR IS WRITTEN CORRECTLY

    // 01000-10100000--0000-1000-0101-1110-1010
    // 11000-00100000--111100-110100-001111-111011 // end = 66!!
    // 10100-10010000--1100010100-0101100100-1111101100-1111101001-1110100001-1101101111-100011011100001111111101000111

    std::cout << T << std::endl;
    uint64_t s_int = prefix2int(s, 0, s.size());
    auto result = bitMagicSearch(T, s_int, s.size()); // input: const sdsl::int_vector<1> &T, string S 
    
    cerr << result.first << ", " << (int)result.second << endl; 
    assert_equal(result, correct_result);
}

int main(int argc, char** argv){
    // Create test directory if does not exist
    if (!exists(temp_dir)){
        create_directory(temp_dir);
    }

    cerr << "Testing tail search..." << endl;
    vector<vector<string>> tails= {{"AA","AC","GG", "CT", "CC"},{"ATT","AGT", "TTA", "TCT"}, {"AGGAT","AGCGG","ATCTT", "GCCTT", "GACCT", "TTCGT", "TGTAC", "TTTAA", "TGAGT"}};
    //                               0    1     2     3     4      5     6      7       8       9       10      11       12        13       14        15     16         17 
    vector <int> tlens = {2,3,5}; 
    vector<int> tails_so_far = {0,5,9};
    sdsl::int_vector<1> T = WriteTailsVector(tails,tlens);

    for (int j=0; j < tails.size(); j++){
        for (int i=0; i < tails[j].size(); i++){
            pair<int64_t, uint8_t> correct_result = {i+tails_so_far[j],tlens[j]};
            test_tail_search(tails[j][i], correct_result, T);
            cerr << "...ok" << endl << endl;
        }
    }
    cerr << "...ok" << endl;

    cerr << "Testing again tail search..." << endl;
    tails= {{"AA","AC","GG", "CT", "CC"},{"ATT","AGT", "TTA", "TCT"}, {"TGTAC", "TTTAA", "TGAGT"}};
    //        0    1     2     3     4      5     6      7       8       9        10        11      
    tlens = {2,3,5}; 
    tails_so_far = {0,5,9};
    T = WriteTailsVector(tails,tlens);

    for (int j=0; j < tails.size(); j++){
        for (int i=0; i < tails[j].size(); i++){
            pair<int64_t, uint8_t> correct_result = {i+tails_so_far[j],tlens[j]};
            test_tail_search(tails[j][i], correct_result, T);
            cerr << "...ok" << endl << endl;
        }
    }
    cerr << "ALL TESTS PASSED" << endl;

}