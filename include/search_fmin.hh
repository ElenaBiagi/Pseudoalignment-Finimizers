#pragma once

#include <string>
#include <cstring>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <variant>

//#include "PackedStrings.hh"
#include "SeqIO.hh"
#include "ColoredFinimizers.hh"

//#include <sdsl/elias_fano_vector.hpp>

using namespace std;

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_queries_streaming(reader_t& reader, out_stream_t& out, const CompressedColoredFinimizers& index, const float& t){

    int64_t total_micros = 0;

    // 1. Read all seq
    vector<string> reads;
    while (true) {
        int64_t len = reader.get_next_read_to_buffer();
        if (len == 0) break;
        reads.push_back(reader.read_buf);  // Copy each read
    }

    const size_t n = reads.size();

    using ResultType = variant<vector<vector<float>>, vector<vector<uint64_t>>>;
    ResultType result;

    if (t > 0) {
        result = vector<vector<float>>{};
    } else {
        result = vector<vector<uint64_t>>{};
    }

    // 2. Parallel search
    int64_t t0 = cur_time_micros();

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        const std::string& seq = reads[i];

        if (auto* res = std::get_if<std::vector<std::vector<float>>>(&result)) {
            index.search(seq, (*res)[i], t);  // Safe because each thread writes to unique index
        } else if (auto* res = std::get_if<std::vector<std::vector<uint64_t>>>(&result)) {
            index.search(seq, (*res)[i]);     // Safe for same reason
        }
    }

    total_micros += cur_time_micros() - t0;

    // 3. Sort output 
    if (auto* res = std::get_if<std::vector<std::vector<uint64_t>>>(&result)) {
        for (size_t j = 0; j < res->size(); ++j) {
            out << j << " ";
            std::vector<std::pair<int, int64_t>> nonzero_entries;

            for (size_t idx = 0; idx < (*res)[j].size(); ++idx) {
                if ((*res)[j][idx] > 0)
                    nonzero_entries.emplace_back(static_cast<int>(idx), (*res)[j][idx]);
            }

            std::sort(nonzero_entries.begin(), nonzero_entries.end(), [](const auto& a, const auto& b) {
                return (a.second > b.second) || (a.second == b.second && a.first < b.first);
            });

            for (const auto& [idx, count] : nonzero_entries) {
                out << idx << ":" << count << " ";
            }
            out << "\n";
        }

    } else if (auto* res = std::get_if<std::vector<std::vector<float>>>(&result)) {
        for (size_t j = 0; j < res->size(); ++j) {
            out << j << " ";
            std::vector<std::pair<int, float>> nonzero_entries;

            for (size_t idx = 0; idx < (*res)[j].size(); ++idx) {
                if ((*res)[j][idx] > 0)
                    nonzero_entries.emplace_back(static_cast<int>(idx), (*res)[j][idx]);
            }

            std::sort(nonzero_entries.begin(), nonzero_entries.end(), [](const auto& a, const auto& b) {
                return (a.second > b.second) || (a.second == b.second && a.first < b.first);
            });

            for (const auto& [idx, score] : nonzero_entries) {
                out << idx << ":" << score << " ";
            }
            out << "\n";
        }
    }

    write_log("us/query (excluding I/O): " + to_string((double)total_micros / n), LogLevel::MAJOR);
    return static_cast<int64_t>(n);
}

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_file(const string& infile, out_stream_t& out, const CompressedColoredFinimizers& index, const float& t){
    reader_t reader(infile);
    //write_log("Running streaming queries from input file " + infile, LogLevel::MAJOR);
    return run_fmin_queries_streaming(reader, out, index, t);
}

// Returns number of queries executed
int64_t run_fmin_queries(const vector<string>& infiles, const optional<vector<string>>& outfiles, const CompressedColoredFinimizers& index, const float& t){

    if(outfiles.has_value()){
        if(infiles.size() != outfiles.value().size()){
            string count1 = to_string(infiles.size());
            string count2 = to_string(outfiles.value().size());
            throw runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
        }
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<ifstream>> in_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input){
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], out, index, t);
            } else { // To stdout
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], cout, index, t);
            }
        }
        else {
            if(outfiles.has_value()){
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], out, index, t);
            } else{ // To stdout
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], cout, index, t);
            }
        }
    }
    return n_queries_run;
}

int search_fmin(int argc, char** argv){

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all Finimizers of all input reads.");

    options.add_options()
        ("o,out-file", "Output filename, or stdout if not given.", cxxopts::value<string>())
        ("i,index-file", "Index filename prefix.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("t", "Threshold", cxxopts::value<float>()->default_value("0"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        cerr << options.help() << endl;
        exit(1);
    }

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> query_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if(multi_file){
        query_files = readlines(queryfile);
    } else{
        query_files = {queryfile};
    }
    for(string file : query_files) check_readable(file);

    // Interpret output file
    optional<vector<string>> output_files;
    try{
        string outfile = opts["out-file"].as<string>();
        if(multi_file){
            output_files = readlines(outfile);
        } else{
            output_files = {outfile};
        }
        for(string file : output_files.value()) check_writable(file);
    } catch(cxxopts::option_has_no_value_exception& e){
        write_log("No output file given, writing to stdout", LogLevel::MAJOR);
    }

    string index_prefix = opts["index-file"].as<string>();

    int64_t number_of_queries = 0;
    float t = opts["t"].as<float>();

    cerr << "Loading index..." << endl;
    CompressedColoredFinimizers index;
    index.load(index_prefix);
    cerr << "Index loaded" << endl;

    number_of_queries += run_fmin_queries(query_files, output_files, index, t); // TODO: Implement this in ColoredFinimizers
    int64_t new_total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)new_total_micros / number_of_queries), LogLevel::MAJOR);
    write_log("total number of queries: " + to_string(number_of_queries), LogLevel::MAJOR);
    
    /* 
    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);
    */

    return 0;

}
