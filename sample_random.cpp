// To compile: g++ sample_random.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o sample_random -Wno-deprecated-declarations  -lz

#include "sbwt/SeqIO.hh"
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <memory>
#include <filesystem>

using namespace std;

// From Themisto codebase
vector<string> split_at_non_ACGT(const char* S, int64_t S_size) {
    vector<string> parts;
    int64_t start = 0;
    for (int64_t end = 0; end <= S_size; end++) { // Exclusive end point
        if (end == S_size || (S[end] != 'A' && S[end] != 'C' && S[end] != 'G' && S[end] != 'T')) {
            if (end > start) {
                parts.push_back(string(S + start, S + end));
            }
            start = end + 1;
        }
    }
    return parts;
}

int main(int argc, char** argv) {
    //srand(time(0));

    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <sequence_file> <number_of_samples> <sample_read_length> <engine_number>" << endl;
        return 1;
    }

    string seqfile = argv[1];
    int64_t howmany = stoll(argv[2]);
    int64_t sample_read_length = stoll(argv[3]);
    int i = std::stoi(argv[4]);

    srand(24+i); // reproducible


    vector<string> seqs;
    sbwt::SeqIO::Reader<> in(seqfile);
    cerr << "Reading sequences..." << endl;
    while (int64_t len = in.get_next_read_to_buffer()) {
        for (const string& subseq : split_at_non_ACGT(in.read_buf, len)) {
            if (subseq.size() >= sample_read_length) {
                seqs.push_back(subseq);
            }
        }
    }
    cerr << seqs.size() << " sequences read" << endl;

    cerr << "Sampling " << howmany << " subsequences of length " << sample_read_length << endl;

    vector<int64_t> cumul_weight;
    int64_t total_weight = 0;
    for (int64_t seq_id = 0; seq_id < seqs.size(); seq_id++) {
        int64_t weight = seqs[seq_id].size() - sample_read_length + 1;
        total_weight += weight;
        cumul_weight.push_back(total_weight);
    }

    for (int64_t i = 0; i < howmany; i++) {
        int64_t r = rand() % total_weight;
        int64_t seq_id = lower_bound(cumul_weight.begin(), cumul_weight.end(), r) - cumul_weight.begin();
        int64_t start = rand() % (seqs[seq_id].size() - sample_read_length + 1);
        cout << ">" << i << "\n" << seqs[seq_id].substr(start, sample_read_length) << "\n";
    }

    return 0;
}
