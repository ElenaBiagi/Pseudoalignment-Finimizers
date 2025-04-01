//g++ -std=c++11 mutate_seq.cpp -o mutate_seq
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <chrono>
#include <iostream>


using namespace std;


#include <fstream>
#include <vector>
#include <string>
#include <iostream>

std::vector<std::string> ReadSequencesFromFile(const std::string& inputFile) {
    std::vector<std::string> sequences;
    std::ifstream inFile(inputFile);
    if (!inFile) {
        std::cerr << "Error: Could not open the file " << inputFile << " for reading.\n";
        return sequences;
    }

    std::string line;
    std::string currentSeq;
    while (std::getline(inFile, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }
        if (line[0] == '>') { // Header line in FASTA format
            if (!currentSeq.empty()) {
                sequences.push_back(currentSeq);
                currentSeq.clear();
            }
        } else {
            currentSeq += line;
        }
    }
    if (!currentSeq.empty()) {
        sequences.push_back(currentSeq);
    }

    inFile.close();
    return sequences;
}

void WriteSequencesToFile(const std::string& outputFile, const std::vector<std::string>& sequences) {
    std::ofstream outFile(outputFile);
    if (!outFile) {
        std::cerr << "Error: Could not open the file " << outputFile << " for writing.\n";
        return;
    }

    for (size_t i = 0; i < sequences.size(); ++i) {
        outFile << "> " << i << "\n";
        outFile << sequences[i] << "\n";
    }

    outFile.close();
}

// Function to introduce random mutations into a DNA sequence
string InsertDNAMutation(const string &seq, double mutation_rate, mt19937 &engine) {
    string alphabet = "ACGT";
    uniform_int_distribution<> distA(0, alphabet.length() - 1);
    uniform_int_distribution<> distM(0, 2); // For random mutation type
    uniform_real_distribution<> distP(0.0, 1.0); // For mutation probability

    string mutated_seq = seq;

    for (size_t pos = 0; pos < mutated_seq.length(); ++pos) {
        if (distP(engine) < mutation_rate) {
            int mutationType = distM(engine);
            if (mutationType == 0) { // INSERTION
                char new_base = alphabet[distA(engine)];
                mutated_seq.insert(pos, 1, new_base);
            } else if (mutationType == 1) { // DELETION
                if (!mutated_seq.empty()) { // Ensure there is something to delete
                    mutated_seq.erase(pos, 1);
                    --pos; // Adjust pos after deletion
                }
            } else { // SUBSTITUTION
                char new_base = alphabet[distA(engine)];
                mutated_seq[pos] = new_base;
            }
        }
    }
    return mutated_seq;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " {inputFile} {outputFile} {percentModifications} {numCopies} {engineNumber}\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = std::string(argv[2]) + argv[3] + ".fa";

    double percent = std::stod(argv[3]);
    int numCopies = std::stoi(argv[4]);
    int i = std::stoi(argv[5]);

    //std::random_device rd;
    //std::mt19937 engine(rd());
    std::mt19937 engine(24+i); // reproducible

    // Read sequences from input file
    std::vector<std::string> sequences = ReadSequencesFromFile(inputFile);
    if (sequences.empty()) {
        std::cerr << "No sequences found in the input file.\n";
        return 1;
    }

    std::vector<std::string> modifiedSequences;
    for (const auto& seq : sequences) {
        for (int i = 0; i < numCopies; ++i) {
            std::string modifiedSeq = InsertDNAMutation(seq, percent, engine);
            modifiedSequences.push_back(modifiedSeq);
        }
    }

    // Write modified sequences to output file
    WriteSequencesToFile(outputFile, modifiedSequences);

    return 0;
}