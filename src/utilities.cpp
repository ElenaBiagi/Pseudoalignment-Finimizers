#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <iostream>

using namespace std;

// Utility functions implementation

vector<string> readlines(string filename) {
    vector<string> lines;
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file: " + filename);
    }
    string line;
    while (getline(file, line)) {
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    file.close();
    return lines;
}

void check_readable(string filename) {
    ifstream file(filename);
    if (!file.good()) {
        throw runtime_error("Cannot read file: " + filename);
    }
}

void check_writable(string filename) {
    ofstream file(filename);
    if (!file.good()) {
        throw runtime_error("Cannot write to file: " + filename);
    }
}
