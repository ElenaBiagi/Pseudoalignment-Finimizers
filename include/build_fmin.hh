#pragma once

#include "ColoredFinimizers.hh"

int build_fmin(int argc, char** argv) {
    cxxopts::Options options(argv[0], "Find all Finimizers of all input reads.");

    options.add_options()
        ("i,index-file", "ColloredFinimizers file.", cxxopts::value<string>())
        //("o,out-file", "Output index filename prefix.", cxxopts::value<string>()) // as input
        ("h,help", "Print usage");
    
    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")) {
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    // string out_prefix = opts["out-file"].as<string>(); // TODO used this

    string indexfile = opts["index-file"].as<string>();
    ColoredFinimizers cf;
    ifstream in(indexfile);
    cf.load(in);
    CompressedColoredFinimizers(std::move(cf), 10, 31);
}