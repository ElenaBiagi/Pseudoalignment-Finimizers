#pragma once

#include <string>
#include <cstring>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <variant>

#include "PackedStrings.hh"
#include "SeqIO.hh"
#include "ColoredFinimizers.hh"

//#include <sdsl/elias_fano_vector.hpp>

using namespace std;

template<typename reader_t, typename out_stream_t>
int64_t run_fmin_queries_streaming(reader_t& reader, out_stream_t& out, const CompressedColoredFinimizers& index, const float& t){

    int k = index.get_k();
    int64_t total_micros = 0;
    int64_t number_of_queries = 1; // TODO remove or fix
    int64_t kmers_count = 0 , kmers_count_rev = 0;
    int64_t total_positive = 0;
    vector<int64_t> out_buffer, out_buffer_rev;

    //vector<vector<float>> result = {};
    using ResultType = variant<vector<vector<float>>, vector<vector<uint64_t>>>;
    
    ResultType result;

    if (t > 0) {
        result = vector<vector<float>>{};
    } else {
        result = vector<vector<uint64_t>>{};
    }
    
    int i=0;
    while(true){
        

        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = cur_time_micros();
        //string seq = remove_N_from_string(reader.read_buf);
        string seq = reader.read_buf;

        if (auto* res = get_if<vector<vector<float>>>(&result)) {
            res->push_back({});
            index.search(seq, (*res)[i], t);
        } else if (auto* res = get_if<vector<vector<uint64_t>>>(&result)) {
            res->push_back({});
            index.search(seq, (*res)[i]);
        }

        i++;
        total_micros += cur_time_micros() - t0;
    }
    write_log("k " + to_string(k), LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    //write_log("Found kmers: " + to_string(kmers_count), LogLevel::MAJOR);
    //write_log("Found kmers reverse : " + to_string(kmers_count_rev), LogLevel::MAJOR);
    //write_log("Total found kmers: " + to_string(total_positive), LogLevel::MAJOR);

    // Compare (genome id, # k-mer matched)
    if (auto* res = std::get_if<std::vector<std::vector<uint64_t>>>(&result)) {
        // Handling case: vector<vector<uint64_t>>
        for (int j = 0; j < i; ++j) {
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

            out << std::endl;
        }
    // Compare (genome id, (% k-mer matched) > t)
    } else if (auto* res = std::get_if<std::vector<std::vector<float>>>(&result)) {
        for (int j = 0; j < i; ++j) {
            out << j << " ";

            std::vector<std::pair<int, float>> nonzero_entries;
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

            out << std::endl;
        }
    }
    return number_of_queries;
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
            throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
        }
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

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
        std::cerr << options.help() << std::endl;
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
