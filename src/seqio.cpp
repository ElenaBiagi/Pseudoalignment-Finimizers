#include <string>
#include <algorithm>
#include "SeqIO.hh"

using namespace std;

namespace SeqIO {

void reverse_complement_c_string(char* S, int64_t len) {
    for(int64_t i = 0; i < len; i++) {
        char c = S[i];
        switch(c) {
            case 'A': S[i] = 'T'; break;
            case 'T': S[i] = 'A'; break;
            case 'C': S[i] = 'G'; break;
            case 'G': S[i] = 'C'; break;
            case 'a': S[i] = 't'; break;
            case 't': S[i] = 'a'; break;
            case 'c': S[i] = 'g'; break;
            case 'g': S[i] = 'c'; break;
            default: break; // Keep other characters as-is
        }
    }
    reverse(S, S + len);
}

FileFormat figure_out_file_format(string filename) {
    FileFormat result;
    result.extension = filename;
    result.gzipped = false;
    
    if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz") {
        result.gzipped = true;
        filename = filename.substr(0, filename.size() - 3);
    }
    
    if (filename.size() >= 6 && filename.substr(filename.size() - 6) == ".fasta") {
        result.format = FASTA;
    } else if (filename.size() >= 6 && filename.substr(filename.size() - 6) == ".fastq") {
        result.format = FASTQ;
    } else if (filename.size() >= 2 && filename.substr(filename.size() - 2) == ".f") {
        result.format = FASTA;
    } else if (filename.size() >= 2 && filename.substr(filename.size() - 2) == ".q") {
        result.format = FASTQ;
    } else {
        result.format = FASTA; // Default to FASTA
    }
    
    return result;
}

} // namespace SeqIO
