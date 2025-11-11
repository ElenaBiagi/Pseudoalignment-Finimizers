#pragma once

#include <string>
#include <cstring>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <variant>
#include <sstream>

#include "SeqIO.hh"
#include "ColoredFinimizers.hh"

using namespace std;

template <typename out_stream_t>
void write_out(const char *data, int64_t data_length, out_stream_t &output_writer, vector<char> &buffer, const size_t flush_t = 8 * 1024 * 1024)
{ // TODO find a better way to write threshold
    for (int64_t i = 0; i < data_length; i++)
    {
        buffer.push_back(data[i]);
        if (buffer.size() > flush_t && data[i] == '\n')
        {
            output_writer.write(buffer.data(), buffer.size());
            buffer.clear(); // Let's hope this keeps the reserved capacity of the vector intact
        }
    }
}

uint16_t fast_int_to_string(uint16_t x, char *buffer)
{
    uint16_t i = 0;
    // Write the digits in reverse order (reversed back at the end)
    do
    {
        buffer[i++] = '0' + (x % 10);
        x /= 10;
    } while (x > 0);
    std::reverse(buffer, buffer + i);
    buffer[i] = '\0';
    return i;
}

template <typename reader_t, typename out_stream_t>
int64_t run_fmin_queries_streaming(reader_t &reader, out_stream_t &out, const CompressedColoredFinimizers &index, const float &t)
{
    // int64_t total_micros = 0;

    int i = 0;
    constexpr size_t flush_t = 8 * 1024 * 1024; // 8 MB //1 << 20; // 1MB

    std::vector<char> output_buffer;
    output_buffer.reserve(flush_t * 2);

    // For output writing
    char int_buf[32]; // Enough space for a 64-bit integer in ascii

    vector<int64_t> Finimizers;
    // vector<pair<uint16_t, uint16_t>> ans;
    // ans.resize(index.n_colors)

    vector<int16_t> results(index.n_colors, 0);

    if (t > 0)
    {

        while (true)
        {
            int64_t len = reader.get_next_read_to_buffer();
            if (len == 0)
                break;

            int64_t id_len = fast_int_to_string(i, int_buf);
            write_out(int_buf, id_len, out, output_buffer, flush_t);
            write_out(" ", 1, out, output_buffer, flush_t);

            uint64_t colors_seen_alo = 0;
            uint64_t colors_out = 0;

            const string &seq = reader.read_buf;

            cerr << "n_fmin, n_colorset_ids,  n_colors_seen, n_colors_alo, n_colors_out" << endl;

            const int16_t min_value = index.search(seq, results, t, Finimizers);
            auto start = std::chrono::high_resolution_clock::now();
            for (auto idx = 0; idx < results.size(); idx++)
            {
                const auto &count = results[idx];
                if (count > 0)
                {
                    colors_seen_alo++;
                    if (count >= min_value)
                    {
                        uint16_t idx_len = fast_int_to_string(idx, int_buf);
                        write_out(int_buf, idx_len, out, output_buffer, flush_t);
                        write_out(" ", 1, out, output_buffer, flush_t);
                        colors_out++;
                    }
                }
            }
            cerr << colors_seen_alo << ", " << colors_out << endl;

            /* for (int a = static_cast<int>(ans.size()) - 1; a >= 0; a--)
            {
                const auto &[idx, count] = ans[a];
                if (count >= min_value)
                {
                    uint16_t idx_len = fast_int_to_string(idx, int_buf);
                    write_out(int_buf, idx_len, out, output_buffer, flush_t);
                    write_out(" ", 1, out, output_buffer, flush_t);
                }
            } */

            write_out("\n", 1, out, output_buffer, flush_t);

            i++;
            auto end = std::chrono::high_resolution_clock::now();
            time_output += (end - start);
        }
    }
    else
    { // Print everything
        while (true)
        {
            int64_t len = reader.get_next_read_to_buffer();
            if (len == 0)
            {
                break;
            }
            int64_t id_len = fast_int_to_string(i, int_buf);
            write_out(int_buf, id_len, out, output_buffer, flush_t);
            write_out(" ", 1, out, output_buffer, flush_t);

            const string &seq = reader.read_buf;

            index.search(seq, results, Finimizers);

            auto start = std::chrono::high_resolution_clock::now();
            for (auto idx = 0; idx < results.size(); idx++)
            {
                const auto &count = results[idx];
                if (count >= 0)
                {
                    uint16_t idx_len = fast_int_to_string(idx, int_buf);
                    write_out(int_buf, idx_len, out, output_buffer, flush_t);
                    write_out(" ", 1, out, output_buffer, flush_t);
                }
            }

            /* for (int a = static_cast<int>(ans.size()) - 1; a >= 0; a--)
            {

                const auto &[idx, count] = ans[a];
                if (count >= 0)
                {
                    uint16_t idx_len = fast_int_to_string(idx, int_buf);
                    write_out(int_buf, idx_len, out, output_buffer, flush_t);
                    write_out(" ", 1, out, output_buffer, flush_t);
                }
            } */

            write_out("\n", 1, out, output_buffer, flush_t);

            i++;
            auto end = std::chrono::high_resolution_clock::now();
            time_output += (end - start);
        }
    }
    auto start = std::chrono::high_resolution_clock::now();

    if (!output_buffer.empty())
    {
        out.write(output_buffer.data(), output_buffer.size());
    }

    auto end = std::chrono::high_resolution_clock::now();
    time_output += (end - start);

    print_search_timing_stats();
    return 1;
}

template <typename reader_t, typename out_stream_t>
int64_t run_fmin_file(const string &infile, out_stream_t &out, const CompressedColoredFinimizers &index, const float &t)
{
    reader_t reader(infile);
    return run_fmin_queries_streaming(reader, out, index, t);
}

// Returns number of queries executed
int64_t run_fmin_queries(const vector<string> &infiles, const optional<vector<string>> &outfiles, const CompressedColoredFinimizers &index, const float &t)
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

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<ifstream>> in_no_gzip;

    int64_t n_queries_run = 0;
    for (int64_t i = 0; i < infiles.size(); i++)
    {
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if (gzip_input)
        {
            if (outfiles.has_value())
            {
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], out, index, t);
            }
            else
            { // To stdout
                n_queries_run += run_fmin_file<in_gzip>(infiles[i], cout, index, t);
            }
        }
        else
        {
            if (outfiles.has_value())
            {
                ofstream out(outfiles.value()[i]);
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], out, index, t);
            }
            else
            { // To stdout
                n_queries_run += run_fmin_file<in_no_gzip>(infiles[i], cout, index, t);
            }
        }
    }
    return n_queries_run;
}

int search_fmin(int argc, char **argv)
{

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all Finimizers of all input reads.");

    options.add_options()("o,out-file", "Output filename, or stdout if not given.", cxxopts::value<string>())("i,index-file", "Index filename prefix.", cxxopts::value<string>())("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())("t", "Threshold", cxxopts::value<float>()->default_value("0"))("h,help", "Print usage");

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
        write_log("No output file given, writing to stdout", LogLevel::MAJOR);
    }

    string index_prefix = opts["index-file"].as<string>();

    int64_t number_of_queries = 0;
    float t = opts["t"].as<float>();

    cerr << "Loading index..." << endl;
    /* int64_t total_micros = 0;
    int64_t t0 = cur_time_micros();
    */
    CompressedColoredFinimizers index;
    {
        auto start = std::chrono::high_resolution_clock::now();
        index.load(index_prefix);
        auto end = std::chrono::high_resolution_clock::now();
        time_index_loading += (end - start);
    }

    /* total_micros += cur_time_micros() - t0;
    cerr << total_micros << endl; */
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
