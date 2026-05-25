#pragma once

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <variant>
#include <sstream>

#include "SeqIO.hh"
#include "ColoredFinimizers.hh"
#include <cxxopts.hpp>

using namespace std;



template <typename reader_t, typename out_stream_t>
void run_fmin_queries_batch(reader_t &reader, out_stream_t &out, const CompressedColoredFinimizers &index, const float &t, const uint64_t batch_size)
{
    vector<string> reads;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if (len == 0)
            break;
        const string &seq = reader.read_buf;
        reads.push_back(seq);
    }
    
    // This assumes that all the reads have the same length
    const uint64_t query_len = reads[0].size();

    cerr << "# reads: " << reads.size() <<'\n';
    cerr << "read length: " << query_len <<'\n';
    cerr << "batch size: " << batch_size <<'\n';

    // int i = 0;
    // constexpr size_t flush_t = 8 * 1024 * 1024; // 8 MB //1 << 20; // 1MB

    // std::vector<char> output_buffer;
    // output_buffer.reserve(flush_t * 2);

    // 1. Identify finimizers (Fluke8 and PrefTab, batch)
    // 2. Select the smallest finimizer for each k-mer
    // 3. Deal with colorsets (DONE)
    // 4. Print results (~DONE)
    uint64_t k = index.get_k();
    if (t>0){
        index.search_batch(reads, query_len, batch_size, k, t);
    }
    else{
        index.search_batch(reads, query_len, batch_size, k);
    }
    
    return;
}

template <typename reader_t, typename out_stream_t>
void run_fmin_file(const string &infile, out_stream_t &out, const CompressedColoredFinimizers &index, const float &t, const uint64_t batch_size)
{
    reader_t reader(infile);
    return run_fmin_queries_batch(reader, out, index, t, batch_size);
}

// Returns number of queries executed
void run_fmin_queries(const vector<string> &infiles, const optional<vector<string>> &outfiles, const CompressedColoredFinimizers &index, const float &t, const uint64_t batch_size)
{

    if (outfiles.has_value())
    {
        if (infiles.size() != outfiles.value().size())
        {
            string count1 = to_string(infiles.size());
            string count2 = to_string(outfiles.value().size());
            throw runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
        }
    }

    typedef SeqIO::Reader<Buffered_ifstream<ifstream>> Reader;

    int64_t n_queries_run = 0;
    for (int64_t i = 0; i < infiles.size(); i++)
    {
        if (outfiles.has_value())
        {
            ofstream out(outfiles.value()[i]);
            run_fmin_file<Reader>(infiles[i], out, index, t, batch_size);
        }
        else
        { // To stdout
            run_fmin_file<Reader>(infiles[i], cout, index, t, batch_size);
        }
    }
    return;
}

int search_fmin_batch(int argc, char **argv)
{

    // int64_t micros_start = cur_time_micros();

    // set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all Finimizers of all input reads.");

    options.add_options()
                        ("o,out-file", "Output filename, or stdout if not given.", cxxopts::value<string>())
                        ("i,index-file", "Index filename prefix.", cxxopts::value<string>())
                        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
                        ("t", "Threshold", cxxopts::value<float>()->default_value("0"))
                        ("b", "Batch size", cxxopts::value<uint64_t>()->default_value("1000"))
                        ("h,help", "Print usage");


    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help"))
    {
        cerr << options.help() << endl;
        exit(1);
    }

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> query_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if (multi_file)
    {
        query_files = readlines(queryfile);
    }
    else
    {
        query_files = {queryfile};
    }
    for (string file : query_files)
        check_readable(file);

    // Interpret output file
    optional<vector<string>> output_files;
    try
    {
        string outfile = opts["out-file"].as<string>();
        if (multi_file)
        {
            output_files = readlines(outfile);
        }
        else
        {
            output_files = {outfile};
        }
        for (string file : output_files.value())
            check_writable(file);
    }
    catch (cxxopts::option_has_no_value_exception &e)
    {
        cerr << "No output file given, writing to stdout" << endl;
        // write_log("No output file given, writing to stdout", LogLevel::MAJOR);
    }

    string index_prefix = opts["index-file"].as<string>();

    //int64_t number_of_queries = 0;
    float t = opts["t"].as<float>();
    uint64_t batch_size = opts["b"].as<uint64_t>();

    cerr << "Loading index..." << endl;
    /* int64_t total_micros = 0;
    int64_t t0 = cur_time_micros();
    */
    CompressedColoredFinimizers index;
    {
        // auto start = std::chrono::high_resolution_clock::now();
        index.load(index_prefix);
        // auto end = std::chrono::high_resolution_clock::now();
        // time_index_loading += (end - start);
    }

    /* total_micros += cur_time_micros() - t0;
    cerr << total_micros << endl; */
    cerr << "Index loaded" << endl;

    run_fmin_queries(query_files, output_files, index, t, batch_size);
    // int64_t new_total_micros = cur_time_micros() - micros_start;

    /*
    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);
    */

    return 0;
}
