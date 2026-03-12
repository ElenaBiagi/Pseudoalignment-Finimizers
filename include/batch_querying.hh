#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

#include "BoundedDeque.hh"
#include "Colors_queries.hh"
#include "Print_output_queries.hh"

//#include "common.hh"
//#include "bitsearch.hh"
//#include "Buckets.hh"

#include "Fluke8.hh"
#include "PrefTab.hh"

#include "sdsl/bit_vectors.hpp"

void PickFinimizer(vector<int64_t> &Fmin, const uint64_t kmer_start, const uint64_t k, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t> &k_fmin){ //, const string& input){
    if (!curr_candidates.empty())
    {
        k_fmin = curr_candidates.front();

        Fmin.push_back(get<2>(k_fmin));
    }

    // 1. Check if this finimizer is good for the next k-mer (still in the window)
    while (!curr_candidates.empty() && get<3>(curr_candidates.front()) <= kmer_start)
    {
        curr_candidates.pop_front();
    }
    k_fmin = (curr_candidates.empty()) ? static_cast<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>(make_tuple(k + 1, 0, 0, kmer_start + 1)) : curr_candidates.front();

    // 2. Check if the NEXT finimizer would be good for the next k-mer
    if (!next_candidates.empty())
    {   
        // ending pos, len, rank, color_set_id[rank] ->  // len, rank, color_set_id[rank], start(i)

        const auto &next_fmin = next_candidates.front(); // tuple<uint64_t, uint64_t, uint64_t, uint64_t>
        tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {get<1>(next_fmin), get<2>(next_fmin), get<3>(next_fmin),  get<0>(next_fmin) +1 - get<1>(next_fmin)};
        if (get<0>(next_fmin) <= kmer_start + k)
        { // end of the fmin is before end of next kmer
            if (new_fmin < k_fmin)
            { // always true if curr_candidates is empty
                curr_candidates.clear();
                k_fmin = new_fmin;
            }
            else
            {
                while (curr_candidates.back() > new_fmin)
                {
                    curr_candidates.pop_back();
                }
            }
            curr_candidates.push_back(new_fmin);
            next_candidates.pop_front();
        }
    }
}

void AddFinimizer(const uint64_t rank, const uint64_t f_len, const uint64_t start, const uint64_t end, const vector<uint64_t> &color_set_ids, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t> &k_fmin){ //, const string& input){
    if ((start + f_len) > end)
    {
        next_candidates.push_back(make_tuple(start + f_len-1, f_len, rank, color_set_ids[rank]));
    } // Sorted based on END
    else
    {
        tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {f_len, rank, color_set_ids[rank], start};
        if (new_fmin < k_fmin)
        {
            curr_candidates.clear();
            k_fmin = new_fmin;
        }
        else
        {
            while (!curr_candidates.empty() && curr_candidates.back() > new_fmin)
            {
                curr_candidates.pop_back();
            }
        }
        curr_candidates.push_back(new_fmin);
    }

}

void FindFinimizers (const vector<std::pair<uint64_t, uint64_t>> &batch, const CompressedColorSets &CCS, const uint64_t n_colors, const vector<uint64_t> &color_set_ids, uint64_t k, const float &t){
    
    // TODO this assumes that all the queries have the same length (1000)
    const vector<uint64_t> query_lens(10000,1000);

    uint64_t q = 0;
    for (auto q_idx = 0; q_idx < query_lens.size(); q_idx++){
        if (query_lens[q_idx]<k){continue;}
        // length, rank, color_set_id[rank], start(i)
        BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> curr_candidates(k); // sort based on len, int (color, start)
        // ending pos, len, rank, color_set_id[rank]
        BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> next_candidates(k); // sort by end (start+len-1)
        tuple<uint64_t, uint64_t, uint64_t, uint64_t> k_fmin = {k + 1, 0, 0, 0};

        
        uint64_t end = k-1;

        // NO decisions in this loop, not enough info
        for (auto i = 0; i < k-1; i++, end++){
            uint64_t rank = (batch[i+q].second >> 32) & 0xFFFF;  // extract upper 32 bits
            uint64_t f_len = batch[i+q].first;
            
            if (rank > 0) // 0 == -1
            {
            rank--; // 0 is a valid result 
            AddFinimizer(rank, f_len, i, end, color_set_ids, curr_candidates, next_candidates, k_fmin);
            }

        // Now we have info on the first k-1 pos and we can make decisions
        vector<int64_t> Fmin;
        vector<int16_t> results;
        for (auto i = k-1; i < query_lens[q_idx]; i++, end++){
            if (rank > 0) // 0 == -1
            {
            rank--; // 0 is a valid result 
            AddFinimizer(rank, f_len, i, end, color_set_ids, curr_candidates, next_candidates, k_fmin);
            }
            PickFinimizer(Fmin, i, k, curr_candidates, next_candidates, k_fmin);
        }
        // TODO modeve this branch before. avoid having it for every query
        if (t>0){
            int64_t T = pseudoalignment_stats(Fmin, CCS, n_colors, results, t);
            print_cout_queries(results,q_idx, T);
        }
        else{
            pseudoalignment_stats(Fmin, CCS, n_colors, results);
            print_cout_queries(results,q_idx);
        }
        
        q+=query_lens[q_idx];
        }
    }
}


void batch_querying(const vector<std::string> &reads, const Fluke8 &f8, const PrefTab &ptab, const uint64_t batch_size, const CompressedColorSets &CCS, const uint64_t n_colors, const vector<uint64_t> &color_set_ids, const uint64_t k, const float &t){
    //do the querying
    auto start = std::chrono::system_clock::now(); 
        
    uint64_t n_f8 = f8.getn();

    uint64_t checksum = 0;
    uint64_t checksum2 = 0;
    std::vector<std::pair<uint64_t, uint64_t>> batch;
    std::vector<std::pair<uint64_t, uint64_t>> batch_sorted;
    batch.reserve(batch_size*1000); //1000 is the read length, need to do this better
    uint64_t f8failedsearches = 0;
    for(uint64_t bi=0; bi<reads.size(); bi+=batch_size){
        uint64_t batchend = std::min((uint64_t)(reads.size()),bi+batch_size);
        batch.clear();
        //prepare the batch
        uint64_t bptr = 0;
        for(uint64_t i=bi;i<batchend;i++){
            std::string r = reads[i];
            //cerr << r << '\n';
            uint64_t packedr = 0;
            for(uint64_t j=0; j<30; j++){
                //pack a binary representation of the 31-mer into packedr
                packedr |= ((uint64_t)char2bits(r[j])) << ((uint64_t)(2*(31-j)));
            } 
            for(uint64_t j=0; j<r.length()-30; j++){
                //packedr |= (uint64_t)char2bits(r[j+30]);
                packedr |= ((uint64_t)char2bits(r[j+30])<<(uint64_t)2);
                batch.push_back({packedr,bptr+j});
                //cerr << "X bptr: " << bptr+j <<'\n';
                //int64_t ret = ptab.finiLookup(packedr);
                //checksum += ret;
                //if(ret != -1){
                //   cerr << i << ": " << ret << '\n';
                //}
                //pair<int64_t,bool> res = ptab.getPred(packedr);
                //cerr << packedr << ' ' << res.first << ' ' << res.second << '\n';
                packedr = packedr << 2;
                //cerr << "----------------------------\n";
                //if (j > 3) exit(1);
            }
            //TODO: next loop makes k-mers with trailing A's overhanging the end of the read
            for(uint64_t j=0; j<30; j++){
                //cerr << "Y bptr: " << bptr+j+r.length()-30 <<'\n';
                //cerr << "'K-mer starts at j: " << j+r.length()-30 << "\n";
                //packedr |= (uint64_t)char2bits(r[j+30]);
                packedr |= ((uint64_t)char2bits('A')<<(uint64_t)2);
                batch.push_back({packedr,bptr+j+r.length()-30});
                //cerr << packedr << '\n';
                //int64_t ret = f8.finiLookup(packedr,30-j);
                packedr = packedr << 2;
                //cerr << "----------------------------\n";
            }

            bptr += r.length(); //-31;
        }
        //cerr << "bptr: "<<bptr<<" batch.size(): "<<batch.size()<<'\n';
        //cerr<<"About tosort\n";
        //sort the batch
        //std::sort(batch.begin(),batch.end());

        //partially sort the batch (bucket the elements in it based on prefix into 2^16 groups) --- much faster than sorting
        batch_sorted.reserve(batch.size());
        uint64_t C[65536];
        for(uint64_t i=0;i<65536;i++) C[i] = 0;
        for(uint64_t i=0;i<batch.size();i++){
            C[batch[i].first>>48]++;
        }
        //cerr<<"Counting done\n";
        uint64_t psum = 0;
        for(uint64_t i=0;i<65536;i++){
            uint64_t count = C[i];
            C[i] = psum;
            psum += count;
            //cerr<<"psum: "<<psum<<"\n";
        }
        //cerr<<"Summing done\n";
        for(uint64_t i=0;i<batch.size();i++){
            uint64_t pos = C[batch[i].first>>48];
            //cerr<<"pos: "<<pos<<"\n";
            //cerr<<"batch[i].first>>48: "<<(batch[i].first>>48)<<"\n";
            batch_sorted[pos]=batch[i];
            C[batch[i].first>>48]++;
        }
        //cerr<<"Rearrangement done\n";

        for(uint64_t i=0;i<batch.size();i++){
            uint64_t len = 1000 - (batch_sorted[i].second%1000);
            len = ((len >= 31) ? 31 : len);
            //cerr << "pos: " << (batch_sorted[i].second%1000) << " len: "<<len<<'\n';
            pair<int64_t,uint64_t> ret = ptab.finiLookup(batch_sorted[i].first,len);
            batch_sorted[i].second = (((uint64_t)(ret.first+1+n_f8)) << 32) | batch_sorted[i].second;
            batch_sorted[i].first = ret.second;
            checksum += (ret.first+1);
        }
        for(uint64_t i=0;i<batch.size();i++){
            if(batch_sorted[i].second>>32 == 0){
                uint64_t len = 1000 - (batch_sorted[i].second%1000);
                //cerr << "pos: " << (batch_sorted[i].second%1000) << " len: "<<len<<'\n';
                len = ((len >= 31) ? 31 : len);
                // ret is the rank of the finimap found
                // todo add the number of short finimizers stored in f8  n_f8
                
                pair<int64_t,uint64_t> ret = f8.finiLookup(batch_sorted[i].first, len);
                //cerr << "ret: "<<ret.first<<'\n';
                batch_sorted[i].second = (((uint64_t)(ret.first+1)) << 32) | batch_sorted[i].second;
                f8failedsearches += (ret.first == -1);
                batch_sorted[i].first = ret.second;
                checksum += (ret.first+1);
                //cerr << "ret: "<<ret<<'\n';
            
            }
        }
        for(uint64_t i=0;i<batch.size();i++){
            batch[batch_sorted[i].second & 0xFFFFFFFF] = batch_sorted[i];
            //int64_t ret = batch_sorted[i].first;
        }
        // for(uint64_t i=0;i<batch.size();i++){
        //     checksum2 += batch[i].second>>32;
        //     //if(batch[i].second>>32){
        //     //   cerr << "Found at: " << i << '\n';
        //     //}
        // }

        // Identify correct finimizers for every pos
        // TODO keep a vector of read lengths or assume they all have the same length
        FindFinimizers (batch, CCS, n_colors, color_set_ids, k, t);
        // print results at the end of each read
    }
    return;
}


void counting_sort(const vector<int16_t> &results, vector<pair<uint16_t, uint16_t>> &ans, const size_t found_fmin, const uint16_t n_colors)
{
    vector<uint16_t> counts(found_fmin + 1);

    for (size_t idx = 0; idx < n_colors; idx++)
    {
        counts[results[idx]]++;
    }

    // Cumulative Sums
    for (size_t c = 1; c < counts.size(); c++)
    {
        counts[c] += counts[c - 1];
    }

    ans.resize(n_colors);
    for (size_t idx = 0; idx < n_colors; idx++)
    {
        ans[counts[results[idx]] - 1] = {static_cast<uint16_t>(idx), results[idx]};
        counts[results[idx]]--;
    }
}