#pragma once

#include "ColoredFinimizers.hh"

int build_fmin(int argc, char** argv) {
    cxxopts::Options options(argv[0], "Compress colored Finimizers.");

    options.add_options()
        ("i,index-file", "ColloredFinimizers file.", cxxopts::value<string>())
        ("o,out-file", "Output index filename prefix.", cxxopts::value<string>()) // as input
        ("p,p_len", "Finimizers prefix length.", cxxopts::value<int64_t>()->default_value(std::to_string(10)))
        ("k", "k-mer length.", cxxopts::value<int64_t>()->default_value(std::to_string(31)))
        ("h,help", "Print usage");
    
    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")) {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    string out_prefix = opts["out-file"].as<string>();

    string indexfile = opts["index-file"].as<string>();
    ColoredFinimizers cf;
    ifstream in(indexfile);
    cf.load(in);
    int64_t prefix_len = opts["p"].as<int64_t>();
    int64_t k = opts["k"].as<int64_t>();
    CompressedColoredFinimizers ccf(std::move(cf), prefix_len, k); // Search currently works only with 31
    
    ccf.serialize(out_prefix);
    
    return 0;
}