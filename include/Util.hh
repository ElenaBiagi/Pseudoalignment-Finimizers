#ifndef _UTIL_H_
#define _UTIL_H_

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <ratio>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

inline char char2bits(int c){
   if(c == 'A') return 0;
   if(c == 'C') return 1;
   if(c == 'G') return 2;
   if(c == 'T') return 3;
   cerr << "Unrecognized symbol: "<<c<<'\n';
   return -1;
}

#endif
