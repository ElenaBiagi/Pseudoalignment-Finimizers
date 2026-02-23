#pragma once

#include "ColoredFinimizers.hh"

int print_fmin(int argc, char** argv) {
    cxxopts::Options options(argv[0], "Compress colored Finimizers.");

    options.add_options()
        ("i,index-file", "ColoredFinimizers file.", cxxopts::value<string>())
        ("o,out-file", "Output fmin as text.", cxxopts::value<string>()) // as input // from finimizer_matrix
        //("p,p_len", "Finimizers prefix length.", cxxopts::value<int64_t>()->default_value(std::to_string(10)))
        //("k", "k-mer length.", cxxopts::value<int64_t>()->default_value(std::to_string(31)))
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
    //int64_t prefix_len = opts["p"].as<int64_t>();
    //int64_t k = opts["k"].as<int64_t>();
    vector<string> fmins = cf.print_fmins();

    // TODO WRITE TO A FILE
    // open output file
    // ofstream out(out_prefix);
    // if (!out) {
    //     cerr << "Error: could not open output file " << out_prefix << endl;
    //     return 1;
    // }

    // for (auto& f: fmins){
    //     out << f << '\n';
    // }

    // cf.print_colors(out_prefix);
    cf.save_colors(out_prefix);
    cf.save_fmins(out_prefix);
    // cf.print_colors_as_list(out_prefix);



    //CompressedColoredFinimizers ccf(std::move(cf), prefix_len, k); // Search currently works only with 31
    
    //ccf.serialize(out_prefix);
    
    return 0;
}