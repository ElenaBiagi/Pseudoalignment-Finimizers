#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

//#include "PackedStrings.hh"
#include "SeqIO.hh"
#include "BoundedDeque.hh"

#include "common.hh"
#include "bitsearch.hh"
#include "Buckets.hh"


// TODO simplify this removing what is not necessary
// Do we want to count the number of found kmers? YES
// set ?
vector<uint64_t> rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint8_t plen, const uint8_t k){ 
    const int64_t str_len = input.size();
    vector<uint64_t> Fmin;// pointer to C

    uint64_t start = 0;
    uint64_t kmer_start = 0;
    bool found;
    uint64_t s_int;

    uint8_t f_len;
    uint64_t f_int;
    uint32_t f_color; 

    tuple<uint8_t, uint64_t, uint32_t, uint64_t> curr_substr;
    tuple<uint8_t, uint64_t, uint32_t, uint64_t> k_fmin = {k+1, 0, 0, str_len};
    BoundedDeque<tuple<uint8_t, uint64_t, uint32_t, uint64_t>> all_fmin(str_len);;
    //fmin = {f_len, f_int, f_color, start}

    /* for(auto x: sB){
        cerr << x.first << ": "<< (int)x.second.first <<endl;
    }
    cerr<< endl; */

    // iterate over input
    cerr << str_len-plen+1<< endl;
    for (start = 0; start < str_len-plen+1; start++ ){ //TODO the last k-plen characters cannot contain a prefix
        //cerr << "START: "<< start<< endl;
        found = false;
        uint64_t int_p = prefix2int(input, start, plen);
        // 1. prefix found
        if (buckets[int_p].has_value()){
            cerr << "prefix found !"<< endl;
            cerr << int_p << endl;
            cerr << input.substr(start, plen) << endl;
            
            // extract the tail
            const Bucket& bucket_p = *buckets[int_p];
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen

            char s_len = (str_len >= start+k) ? k-plen : str_len-start-plen;
            uint64_t int_s = prefix2int(input, start+plen, s_len); // tail
            auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            if (pos > -1){
                cerr << "found! len: " << (int)len << endl;
                found = true;
                f_len = plen + len;
                f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); // TODO ADD PREFIX AND TLEN: shift p_int to the left by len*2, and int_s to the right to remove the unused chars
                f_color = bucket_p.color_set_ids[pos];
            }
        }else{
            cerr << endl;
            cerr << "prefix not found" << endl;
            // 2. Prefix NOT found
            // Start from the longest possible prefix
            // if you find a real match, stop
            uint8_t sp_len = plen-1;
            // shift int_p to the right by (plen-plen-1)*2 characters;
            cerr << "prefix = "<< input.substr(start, plen) << " " << int_p << endl;
            uint64_t int_sp = int_p >> ((plen-sp_len)*2);
            auto it = sB.find(int_sp);
            cerr << input.substr(start, plen-1) << endl;
            while (sp_len > 0){ // sp_len must be < plen as the whole prefix was not found
                    cerr << "int_sp: " << int_sp << endl;

                    cerr << "sp_len: "<< (int)sp_len << endl;
                if (it != sB.end() ){
                    
                    cerr << "found sB" << endl;
                    cerr << "sB len: " << (int)it->second.first << endl;
                    if (sp_len == it->second.first){ // real match
                    found = true;
                    f_len = sp_len;
                    f_int = int_sp;
                    f_color = it->second.second;
                    break;
                }
            }
                sp_len--;
                int_sp = int_p >> ((plen-sp_len)*2); // shift to the right by 2 character less than then prev time
                it = sB.find(int_sp);
            }
            cerr << endl;

        }

        if (found){
            //cerr << "found" << endl;
            //cerr << "START: "<< start<< endl;
            curr_substr = {f_len, f_int, f_color, start};
            // check later if this is a good finimizer
            all_fmin.push_back(curr_substr);
        }
        else{
            //cerr << "not found"<< endl;
        }
        if (start >= k-1){ // we have looked at all the characters of the kmer
            
            //1. check end and start are ok
            
            k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};
            uint64_t i = 1;
            // if all_fmin is not empty and the while loop is not entered, k_fmin is the smallest fmin and is already ok
            // Ends after the kmer? (get<0>(k_fmin) + get<3>(k_fmin) >= kmer_start+k)
            // Starts before the k-mer? get<3>(k_fmin) < kmer_start
            while ( i < all_fmin.size()){
            //loop trough all_fmin
                if ((get<0>(k_fmin) + get<3>(k_fmin) - 1 <= kmer_start+k-1) and get<3>(k_fmin) >= kmer_start){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                    break;
                }
                k_fmin = all_fmin[i]; // if this is ok it will be the smallest
                i++;
            }
            if (get<0>(k_fmin) + get<3>(k_fmin) - 1 <= kmer_start+k-1 and get<3>(k_fmin) >= kmer_start){Fmin.push_back(get<2>(k_fmin));} // Store only that start of the colors set ids in color_set_concat
            else{
                for (int j=0; j<all_fmin.size(); j++){
                    auto ff= all_fmin[j];
                    cerr << (int)get<0>(ff)<< ", " << get<1>(ff) << ", "<< get<2>(ff) << ", "<< get<3>(ff) << endl;
                }
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, std::min((uint64_t)k, str_len - kmer_start)) << endl;}
            
            kmer_start++;
            // THIS IS VERY USELESS NOW
            /* k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};
            while (all_fmin.size()>0 && get<3>(k_fmin) < kmer_start) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                all_fmin.pop_front();
                k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};
            } */

        }

    }
    //TODO check the last plen-1 values
    uint64_t int_p = prefix2int(input, start, plen-1);
        for (uint64_t ss = start; ss < str_len; ss++ ){ //TODO the last k-plen characters cannot contain a prefix
            //cerr << "START: "<< ss<< endl;

            // TODO add a function here
            uint8_t sp_len = 1;
            // shift int_p to the right by plen-1 characters;
            uint64_t int_sp = int_p >> ((plen-1-sp_len)*2);
            auto it = sB.find(int_sp);

            while (sp_len < plen-1){ // sp_len must be < plen as the whole prefix was not found
                if (it != sB.end() && sp_len == it->second.first){ // real match
                    found = true;
                    f_len = sp_len;
                    f_int = int_sp;
                    f_color = it->second.second;
                    break;
                }
                sp_len++;
                int_sp = int_p >> ((plen-1-sp_len)*2); // shift to the right by 2 character less than then prev time
                it = sB.find(int_sp);
            }

            if (found){
                cerr << "found!"<< endl;
                curr_substr = {f_len, f_int, f_color, ss};
                all_fmin.push_back(curr_substr);
            }
            if (ss >= k-1){ // we have looked at all the characters of the kmer
            //1. check end is ok
            // TODO is it not ok to check end only before? no becuase we are modifying k_fmin here
            k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};
            uint64_t i = 1;
            // if all_fmin is not empty this is always true, k+1+kmer_start > kmer_start +k 
            
            while ( i < all_fmin.size()){
            //loop trough all_fmin
                if ((get<0>(k_fmin) + get<3>(k_fmin) - 1 <= kmer_start+k-1) and get<3>(k_fmin) >= kmer_start){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                    break;
                }
                k_fmin = all_fmin[i]; // if this is ok it will be the smallest
                i++;
            }
            if (get<0>(k_fmin) + get<3>(k_fmin) - 1 <= kmer_start+k-1 and get<3>(k_fmin) >= kmer_start){
                cerr << "finimizer : "<< "{ "<< (int)get<0>(k_fmin)<< ", " << get<1>(k_fmin) << ", "<< get<2>(k_fmin) << ", "<< get<3>(k_fmin) << " }"<< endl;
                Fmin.push_back(get<2>(k_fmin));} // Store only that start of the colors set ids in color_set_concat
            else{
                for (int j=0; j<all_fmin.size(); j++){
                    auto ff= all_fmin[j];
                    cerr << (int)get<0>(ff)<< ", " << get<1>(ff) << ", "<< get<2>(ff) << ", "<< get<3>(ff) << endl;
                }
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, std::min((uint64_t)k, str_len - kmer_start)) << endl;}
            
            //2. check that start is ok for the next kmer
            kmer_start++;
            /* k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};

            while (all_fmin.size()>0 && get<3>(k_fmin) < kmer_start) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                all_fmin.pop_front();
                k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};
            } */

        }
    }
    return Fmin;
}


vector<uint64_t> new_rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint8_t plen, const uint8_t k){ 
    const int64_t str_len = input.size();
    vector<uint64_t> Fmin;// pointer to C

    uint64_t start = 0;
    uint64_t kmer_start = 0;
    bool found;
    uint64_t s_int;

    uint8_t f_len;
    uint64_t f_int;
    uint32_t f_color; 

    tuple<uint8_t, uint64_t, uint32_t, uint64_t> curr_substr;
    tuple<uint8_t, uint64_t, uint32_t, uint64_t> k_fmin = {k+1, 0, 0, str_len};
    BoundedDeque<tuple<uint8_t, uint64_t, uint32_t, uint64_t>> all_fmin(input.size()-k+1);;
    //fmin = {f_len, f_int, f_color, start}

    // iterate over input
    for (start = 0; start < str_len-plen+1; start++ ){ //TODO the last k-plen characters cannot contain a prefix
        found = false;
        uint64_t int_p = prefix2int(input, start, plen);
        // 1. prefix found
        if (buckets[int_p].has_value()){
            cerr << "prefix found !"<< endl;
            // extract the tail
            const Bucket& bucket_p = *buckets[int_p];
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen

            char s_len = (str_len >= start+k) ? k-plen : str_len-start;
            uint64_t int_s = prefix2int(input, start+plen, s_len); // tail
            auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            if (pos > -1){
                found = true;
                f_len = plen + len;
                f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); // TODO ADD PREFIX AND TLEN: shift p_int to the left by len*2, and int_s to the right to remove the unused chars
                f_color = bucket_p.color_set_ids[pos];
            }
        }else{
            //cerr << "prefix not found" << endl;
            // 2. Prefix NOT found
            // Start from the shortest possible prefix
            // if you find a real match, stop
            uint8_t sp_len = 1;
            // shift int_p to the right by plen-1 characters;
            uint64_t int_sp = int_p >> ((plen-1)*2);
            auto it = sB.find(int_sp);

            while (sp_len < plen){ // sp_len must be < plen as the whole prefix was not found
                if (it != sB.end() && sp_len == it->second.first){ // real match
                    found = true;
                    f_len = sp_len;
                    f_int = int_sp;
                    f_color = it->second.second;
                    break;
                }
                sp_len++;
                int_sp = int_p >> ((plen-sp_len)*2); // shift to the right by 2 character less than then prev time
                it = sB.find(int_sp);
            }

        }

        if (found){
            curr_substr = {f_len, f_int, f_color, start};
            
            // check the END is within the current k-mer
            if (f_len+start < kmer_start + k){
                if (k_fmin > curr_substr){
                    all_fmin.clear();
                    k_fmin = curr_substr;
                } else{
                   while (all_fmin.size()>0 && all_fmin.back() > curr_substr) { // the start is already ok
                        all_fmin.pop_back();
                    }
                }
            }
            
            all_fmin.push_back(curr_substr);
        }
        if (start > k-1){ // we have looked at all the characters of the kmer
            
            //1. check end is ok
            // TODO is it not ok to check end only before? no becuase we are modifying k_fmin here
            uint64_t i = 0;
            while ( i < all_fmin.size() and  get<0>(k_fmin) + get<3>(k_fmin) >= kmer_start+k) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                //loop trough all_fmin
                k_fmin = all_fmin[i];
                i++;
            }
            if (i<all_fmin.size() and get<0>(k_fmin) + get<3>(k_fmin) < kmer_start+k){Fmin.push_back(get<2>(k_fmin));} // Store only that start of the colors set ids in color_set_concat
            else{cerr << "finimizer not found for kmer " << kmer_start << " "      << input.substr(kmer_start, std::min((size_t)31, input.size() - kmer_start)) << endl;}
            
            //2. check that start is ok for the next kmer
            kmer_start++;
            k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};

            while ( all_fmin.size()>0 && get<3>(k_fmin) < kmer_start) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                all_fmin.pop_front();
                k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};
            }


        }

    }
    //TODO check the last plen-1 values
    uint64_t int_p = prefix2int(input, start, plen-1);
    cerr << "last chars at the end: "<< int_p << endl;
        for (uint64_t ss = start; ss < str_len; ss++ ){ //TODO the last k-plen characters cannot contain a prefix
            // TODO add a function here
            uint8_t sp_len = 1;
            // shift int_p to the right by plen-1 characters;
            uint64_t int_sp = int_p >> ((plen-1-sp_len)*2);
            auto it = sB.find(int_sp);

            while (sp_len < plen-1){ // sp_len must be < plen as the whole prefix was not found
                if (it != sB.end() && sp_len == it->second.first){ // real match
                    found = true;
                    f_len = sp_len;
                    f_int = int_sp;
                    f_color = it->second.second;
                    break;
                }
                sp_len++;
                int_sp = int_p >> ((plen-1-sp_len)*2); // shift to the right by 2 character less than then prev time
                it = sB.find(int_sp);
            }

            if (found){
                cerr << "found!"<< endl;
                curr_substr = {f_len, f_int, f_color, ss};
                
                // check the END is within the current k-mer
                if (f_len+ss < kmer_start + k){
                    if (k_fmin > curr_substr){
                        all_fmin.clear();
                        k_fmin = curr_substr;
                    } else{
                    while (all_fmin.size()>0 && all_fmin.back() > curr_substr) { // the start is already ok
                            all_fmin.pop_back();
                        }
                    }
                }
                
                all_fmin.push_back(curr_substr);
            }
            if (ss > k-1){ // we have looked at all the characters of the kmer
            //1. check end is ok
            // TODO is it not ok to check end only before? no becuase we are modifying k_fmin here
            uint64_t i = 0;
            while ( i < all_fmin.size() and  get<0>(k_fmin) + get<3>(k_fmin) >= kmer_start+k) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                // loop trough all_fmin
                k_fmin = all_fmin[i];
                i++;
            }
            if (i<all_fmin.size() and get<0>(k_fmin) + get<3>(k_fmin) < kmer_start+k){Fmin.push_back(get<2>(k_fmin));} // Store only that start of the colors set ids in color_set_concat
            else{cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, 31) << endl;}
            
            //2. check that start is ok for the next kmer
            kmer_start++;
            k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};

            while (all_fmin.size()>0 && get<3>(k_fmin) < kmer_start) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                all_fmin.pop_front();
                k_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, uint32_t, uint64_t>{k+1,0,0,kmer_start};
            }

        }
    }
    return Fmin;
}

vector<uint64_t> old_rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint8_t plen, const uint8_t k){ 
    const int64_t str_len = input.size();

    vector<uint64_t> Fmin;// pointer to C

    //int64_t last_pos = 0;
    int64_t start = 0;
    int64_t kmer_start = 0; // start of the first k-mer

    bool found = false;
    uint8_t len_fmin = 0;
    uint64_t int_fmin = 0;
    //string str_fmin;
    uint64_t s_int = 0;
    int64_t C_fmin = 0; // Can the result of color index be negative??? If not found??


    BoundedDeque<tuple<uint8_t, uint64_t, int64_t, int64_t>> all_fmin(input.size()-k+1);
    uint64_t m = std::numeric_limits<uint64_t>::max();

    tuple<uint8_t, uint64_t, int64_t, int64_t> curr_substr; // length, fmin, C_offset, start
    tuple<uint8_t, uint64_t, int64_t, int64_t> w_fmin = {k+1,m,0, input.size()}; // start will always be < str_len

    
    // idea: look for prefixes of length p in the hashtable B
    // 1. prefix not found: the finimizer might be smaller
    //     Check sB with a shorter prefix
    //      * found = store finimizer
    //      - not found = continue
    // 2. prefix found
    //    Go to where the pointer takes you in T
    //    Check first the shortest lengths

        
    // TODO select the correct fmin for every k-mer
    int not_found = 0;
    // Look at every possible prefix of length plen
    for (start = 0; start < str_len - plen +1; start++) { // Extract p characters at a time
        found = false;
        // TODO convert 32 values at time = 64 bits
        uint64_t int_p = prefix2int(input, start, plen);
        // Look for the prefix in B
        
        //int64_t tails_so_far = B[int_p];

        if (buckets[int_p].has_value()){
            // 2. Prefix found!
            const Bucket& bucket_p = *buckets[int_p];
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen

            char s_len = (str_len >= start+k) ? k-plen : str_len-start;
            s_int = prefix2int(input, start+plen, s_len); // tail
            auto result = bitMagicSearch(bucket_p.tail_data, s_int, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            
            if (result.first != -1){ 
                // b. Tail Found!
                found = true;
                len_fmin = plen + result.second; // add the correct fmin len
                C_fmin = bucket_p.color_set_ids[result.first]; // TODO NO NEED TO STORE THE COLORS NOW AS LONG AS WE KEEP THE OFFSET 
                // Store the string as a number. OK as only strings of the same length will be compared 
                
                // Extract the number instead of converting again
                int_fmin = (s_int >> ((s_len - len_fmin) * 2)) & ((1ULL << (len_fmin * 2)) - 1); 
                //old //int_fmin =  prefix2int(input, start+plen, len_fmin);                
            }
        } else{
            cerr << "prefix not found" << endl;
            // 1. Prefix NOT found
            // Start from the shortest possible prefix
            // if you find a real match, stop
            uint8_t sp_len = 1;
            uint64_t int_sp = (int_p >> ((plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1); // prefix of length 1 of int_p
            auto it = sB.find(int_sp);

            while (sp_len < plen){ // sp_len must be < plen as the whole prefix was not found
                if (it != sB.end() && sp_len == it->second.first){ // real match
                    found = true;
                    len_fmin = sp_len;
                    int_fmin = int_sp;
                    C_fmin = it->second.second;
                    break;
                }
                sp_len++;
                int_sp = (int_p >> ((plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1);
                it = sB.find(int_sp);
            }
        }
        if (found){
            not_found=0;

            curr_substr = {len_fmin, int_fmin, C_fmin, start}; // {len_fmin, int_fmin, C_fmin, start};
            // still don't know if this is the correct fmin
            // check if the length is ok to be a finimizer of the current kmer
            if ((start+len_fmin)< kmer_start +k){
                if (w_fmin > curr_substr){ // Compare LEXICOGRAPHICALLY Only strings of the same length will be compared so it's ok to use numbers
                    all_fmin.clear();
                    w_fmin = curr_substr;
                } else {
                    while (all_fmin.size()>0 and all_fmin.back() > curr_substr) {
                        all_fmin.pop_back();
                    }
                }
            
            }
            all_fmin.push_back(curr_substr);
        }else{
            not_found++;
            if (not_found == k){ cerr << "not a single finimizer for the whole genome"<< endl;}
        }
        
        // Check if we are in a kmer
        //if (end - kmer_start + 1 == k){
        
        if (start >= k-1){ // Once start reaches k-1 it means that we have looked for all the finimizers in the first k-mer = all possible characters of the kmer
            
            // 1. Check start
            while (all_fmin.size()>0 && get<3>(w_fmin) < kmer_start) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                all_fmin.pop_front();
                w_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, int64_t, int64_t>{k+1,0,0,kmer_start};
            }
            // 2. Check end
            int j = 0;
            // while end > kmer_end
            // TODO: GO THROUGH ALL ELEMENTS OF FMIN UNTIL END OR FINDING ONE THAT IS WITHIN THE K-MER
            while (j < all_fmin.size() and (get<0>(w_fmin)+get<3>(w_fmin)) > kmer_start +k) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                w_fmin = (all_fmin.size()>0) ? all_fmin[j] : tuple<uint8_t, uint64_t, int64_t, int64_t>{k+1,0,0,kmer_start};
                j++;
            }

            if (all_fmin.size()>0 and ((get<0>(w_fmin)+get<3>(w_fmin)) < kmer_start +k)){Fmin.push_back(get<2>(w_fmin));} // Store only that start of the colors set ids in color_set_concat
            else{cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, 31) << endl;}
            
            kmer_start++;
        }
    }
    // TODO check the last characters as well (prefixes shorter then plen)
/*     int max_plen = k - plen;
    for (int ss = start; ss< str_len; ss++){
        found=false;
        //int max_plen = str_len - ss; // k-plen (- 1 at every loop)
        uint64_t int_p = prefix2int(input, ss, max_plen); // 
            
        // Start from the shortest possible prefix
        // if you find a real match, stop
        uint8_t sp_len = 1;
        uint64_t int_sp = (int_p >> ((max_plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1); // prefix of length 1 of int_p
        auto it = sB.find(int_sp);

        while (sp_len < max_plen +1){
            if (it != sB.end() && sp_len == it->second.first){ // real match
                found = true;
                len_fmin = sp_len;
                int_fmin = int_sp;
                C_fmin = it->second.second;
                break;
            }
            sp_len++;
            int_sp = (int_p >> ((plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1);
            it = sB.find(int_sp);
        }
        max_plen--;

        if (found){
            not_found=0;

            curr_substr = {len_fmin, int_fmin, C_fmin, start}; // {len_fmin, int_fmin, C_fmin, start};
            // still don't know if this is the correct fmin
            // check if the length is ok to be a finimizer of the current kmer
            if ((start+len_fmin)< kmer_start +k){
                if (w_fmin > curr_substr){ // Compare LEXICOGRAPHICALLY Only strings of the same length will be compared so it's ok to use numbers
                    all_fmin.clear();
                    w_fmin = curr_substr;
                } else {
                    while (all_fmin.back() > curr_substr) {
                        all_fmin.pop_back();
                    }
                }
            all_fmin.push_back(curr_substr);
            }
        }else{
            not_found++;
            if (not_found == k){ cerr << "not a single finimizer for the whole genome"<< endl;}
        }

        // 1. Check start
            while (get<3>(w_fmin) < kmer_start) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                all_fmin.pop_front();
                w_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, int64_t, int64_t>{k+1,0,0,kmer_start};
            }
            // 2. Check end
            int j = 0;
            // while end > kmer_end
            // TODO: GO THROUGH ALL ELEMENTS OF FMIN UNTIL END OR FINDING ONE THAT IS WITHIN THE K-MER
            while (j < all_fmin.size() and (get<0>(w_fmin)+get<3>(w_fmin)) > kmer_start +k) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                w_fmin = (all_fmin.size()>0) ? all_fmin[j] : tuple<uint8_t, uint64_t, int64_t, int64_t>{k+1,0,0,kmer_start};
                j++;
            }

            if (all_fmin.size()>0 and ((get<0>(w_fmin)+get<3>(w_fmin)) < kmer_start +k)){Fmin.push_back(get<2>(w_fmin));} // Store only that start of the colors set ids in color_set_concat
            else{cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, 31) << endl;}
            
            kmer_start++;
        
    }
     */
    return Fmin;
}
   
   // TODO: FOUND COLORS: vector of bits and flip found and store in a vector, then flip back
   // BITMAPS SETS INSTEAD OF COLORS

    //TODO: int for the number of colors, change if needed
    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<uint64_t>& results){ 
        // count the number of finimizers found
        results.assign(n_colors, 0);
        
        #pragma omp parallel for
        for(const auto& start : Fmin){
            for (uint64_t i=0; i< n_colors; i++){
                #pragma omp atomic
                results[i]+=color_sets_concat[start+i];
            }
        }
        return;
    }

    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<float>& results, const float& t){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size(); // # total finimizers
        //size_t rm_fmin = 0; // finimizers not found in the index

        
        vector<uint64_t> tot_res;
        tot_res.assign(n_colors, 0);

        #pragma omp parallel for
        for(const auto& start : Fmin){
            for (uint64_t i=0; i< n_colors; i++){
                #pragma omp atomic
                tot_res[i]+=color_sets_concat[start+i];
            }
        }
        results.assign(n_colors, 0);
        
        #pragma omp parallel for
        for (uint64_t i=0; i< n_colors; i++){         
            // For every color found, (#finimizers with that color)/(#tot finimizers)
            float fraction = static_cast<float>(tot_res[i]/static_cast<float>(found_fmin));
            
            if (t==1 && fraction ==1){ 
                results[i]=1;
            } else if (fraction > t){
                results[i]=fraction;
            }
        }
        return;
    }

//TODO rewrite
// used in build
/* string print_finimizer_stats(const set<tuple<int64_t, int64_t, int64_t>>& finimizers, int64_t n_kmers, int64_t n_nodes, int64_t t){
    int64_t new_number_of_fmin = finimizers.size();
    int64_t sum_freq = 0;
    int64_t sum_len = 0;
    for (auto x : finimizers){
        sum_freq += get<1>(x);
        sum_len += get<0>(x);
    }

    string result = to_string(new_number_of_fmin) + "," + to_string(sum_freq) + "," + to_string(static_cast<float>(sum_freq) / static_cast<float>(new_number_of_fmin)) + "," + to_string(static_cast<float>(sum_len) / static_cast<float>(new_number_of_fmin)) + "," + to_string(n_kmers);

    write_log(to_string(t) + "," + result, LogLevel::MAJOR);
    write_log("#Distinct finimizers: " + to_string(new_number_of_fmin) , LogLevel::MAJOR);
    write_log("Sum of frequencies: " + to_string(sum_freq) , LogLevel::MAJOR);
    write_log("Avg frequency: " + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    write_log("Avg length: " + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    return result;
} */