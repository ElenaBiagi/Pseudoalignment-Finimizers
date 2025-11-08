#pragma once

#include "ColoredFinimizers.hh"

int get_stats_buckets(int argc, char **argv)
{
    cxxopts::Options options(argv[0], "Compress colored Finimizers.");

    options.add_options()("i,index-file", "ColoredFinimizers file.", cxxopts::value<string>())("o,out-file", "Output stats.", cxxopts::value<string>()) // as input
                                                                                                                                                        //("p,p_len", "Finimizers prefix length.", cxxopts::value<int64_t>()->default_value(std::to_string(10)))
                                                                                                                                                        //("k", "k-mer length.", cxxopts::value<int64_t>()->default_value(std::to_string(31)))
        ("h,help", "Print usage");

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help"))
    {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    string out_prefix = opts["out-file"].as<string>();

    string index_prefix = opts["index-file"].as<string>();

    cerr << "Loading index..." << endl;
    CompressedColoredFinimizers CCF;
    CCF.load(index_prefix);
    cerr << "Index loaded" << endl;

    CCF.buckets_stats();

    return 0;
}
