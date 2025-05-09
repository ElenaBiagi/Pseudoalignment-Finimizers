#include "ColoredFinimizers.hh"

int main(int argc, char** argv) {
    ColoredFinimizers cf;
    ifstream in(argv[1]);
    cf.load(in);
    CompressedColoredFinimizers(std::move(cf), 8);
}