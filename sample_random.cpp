#include "./SBWT/include/sbwt/SeqIO.hh"
#include <string>
#include <cstdlib>
#include <iostream>
#include <filesystem>


using namespace std;
using namespace sbwt;


// From Themisto codebase
vector<string> split_at_non_ACGT(const char* S, int64_t S_size){
    vector<string> parts;
    int64_t start = 0;
    for(int64_t end = 0; end <= S_size; end++){ // Exclusive end point
        if(end == S_size || (S[end] != 'A' && S[end] != 'C' && S[end] != 'G' && S[end] != 'T')){
            if(end > start) {
                parts.push_back(string(S + start, S + end));
            }
            start = end + 1;
        }
    }
    return parts;
}

int main(int argc, char** argv){

    srand(time(0));

    string seqfile = argv[1];
    int64_t howmany = stoll(argv[2]);
    int64_t sample_read_length = stoll(argv[3]);

    vector<string> seqs;
    SeqIO::Reader<> in(seqfile);
    cerr << "Reading sequences..." << endl;
    while(int64_t len = in.get_next_read_to_buffer()){
        for(const string& subseq : split_at_non_ACGT(in.read_buf, len)){
            if(subseq.size() >= sample_read_length){
                seqs.push_back(subseq);
            }
        }
    }
    cerr << seqs.size() << " sequences read" << endl;

    cerr << "Sampling " << howmany << " subsequences of length " << sample_read_length << endl;

    // We are going to sample subsequences by first picking the sequence to pick from,
    // and then picking a random starting point. To get an even sample, we must sample
    // the sequences so that the probability of selecting a sequence (the "weight") is
    // proprotional to the number of possible starting points in the sequence. The weight
    // is defined here as the number of substrings of length sample_read_length.
    vector<int64_t> cumul_weight;
    int64_t total_weight = 0;
    for(int64_t seq_id = 0; seq_id < seqs.size(); seq_id++){
        // weight = how many sample_read_length-mers fit into this sequence
        // Is always positive because we discarded sequences that are too short above
        int64_t weight = seqs[seq_id].size() - sample_read_length + 1;
        total_weight += weight;
        cumul_weight.push_back(total_weight);
    }

    for(int64_t i = 0; i < howmany; i++){
        int64_t r = std::rand() % total_weight;

        // Binary search for smallest seq_id such that cumul_weight[seq_id] >= r
        int64_t seq_id = std::lower_bound(cumul_weight.begin(), cumul_weight.end(), r) - cumul_weight.begin();
        int64_t start = std::rand() % (seqs[seq_id].size() - sample_read_length + 1);
        cout << ">" << i << "\n" << seqs[seq_id].substr(start, sample_read_length) << "\n";
    }

}
