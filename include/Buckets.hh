#pragma once

class Bucket {
public:
    Bucket() = default;
    struct Compact_tails{
        int tlen;
        uint64_t int_tail;
        uint32_t color_set_id;

        bool operator<(const Compact_tails& other) const {
            return tlen < other.tlen;  // sort by tlen // No need for them to be in lexicographic order
        }
    };
    
    // Compressed tail data
    // Bit layout: [tail length: u8][#tails: vbyte][concat of bitpacked tails]
    //             [tail length: u8][#tails: vbyte][concat of bitpacked tails] 
    //             [tail length: u8][#tails: vbyte][concat of bitpacked tails] 
    //             ...
    sdsl::int_vector<1> tail_data;

    // Color set ids for each tail
    vector<uint32_t> color_set_ids;

    Bucket(vector<std::string_view>& tails, vector<uint32_t>& unsorted_color_set_ids) {
        
        // Convert tails to Compact_tails
        vector<Compact_tails> B_tails;
        B_tails.reserve(tails.size());
        for (size_t i = 0; i < tails.size(); i++) {
            Compact_tails ct;
            ct.tlen = tails[i].size();
            ct.int_tail = prefix2int(tails[i], 0, ct.tlen);
            ct.color_set_id = unsorted_color_set_ids[i];
            B_tails.push_back(ct);
        }
        std::sort(B_tails.begin(), B_tails.end()); 
        /* for (size_t i = 0; i < B_tails; i++)
            color_set_ids[i]=B_tails[i].color_set_id;
        } */

        //color_set_ids.reserve(unsorted_color_set_ids.size());
        
        // This permutes the color_set_ids
        WriteTailsVector(B_tails);
    }

    void WriteTailsVector(const vector<Compact_tails>& B_tails){
        // tlen == 0 is a special case -> B_tails.size() == 1
        // Compressed tails
        //sdsl::int_vector<1> tail_data;
        
        // Color set ids for every tail
        //vector<uint32_t> color_set_ids;
        color_set_ids.reserve(B_tails.size()); 

        uint64_t total_bits = 0;
        uint32_t ntails = 0;
        vector <int> tlens;
        unordered_map<int, uint8_t> m_ntails;
        
        // B_tails must be of length at least 0;
        int tlen = B_tails[0].tlen;
        int cur_tlen = tlen;

        if (tlen == 0){
            total_bits += 5; // tlen
        }
        else{
            for (size_t i = 0; i < B_tails.size(); ++i) {
                tlen = B_tails[i].tlen;
                
                if (tlen != cur_tlen){ 
                    m_ntails[cur_tlen]=ntails; // only at this point we know how many tails for the previous tlen
                    total_bits += 5; // cur_tlen (previous tlen)   
		            total_bits += vbyte_encode(ntails).size() * 8; // vbyte encoded tail count
                    total_bits += ntails * cur_tlen * 2; // each tail uses cur_tlen*2 bits
                    ntails = 0;
                    cur_tlen = tlen;

                }
                ntails++;
            }     
            // last tails
            m_ntails[tlen]=ntails;
            total_bits += 5; // tlen
            total_bits += vbyte_encode(ntails).size() * 8;
 
           total_bits += ntails * tlen * 2;
        }
        for (auto t:m_ntails){
            cerr << t.first << ", "<< t.second << endl;
        }
//        cerr << "total_bits = "<< (int)total_bits << endl;
        tail_data.resize(total_bits+128); 
        uint64_t* data = tail_data.data();

  //      cerr << "tail_data.size() = "<< (int)tail_data.size() << endl;
        
        int64_t offset = 0;
        uint64_t word_index = 0;
        uint8_t w_offset = 0;
        
        // if tlen changed it was never 0        
        if (tlen == 0){
            color_set_ids.push_back(B_tails[0].color_set_id); // only one color id so nothing changed

   //         cerr << "Writing tlen... "<< endl;
   //             cerr << "offset = " << (int)offset << endl;
                word_index = offset/64;
                w_offset = offset % 64;
                //cerr << "w_offset = " << (int)w_offset << endl;
                //cerr << "w_offset = " << (int)w_offset << endl;
            sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
            offset += 5;
        } 
        else {
            cur_tlen = 0; // tlen cannot be 0

            for (size_t i=0; i < B_tails.size(); i++){

                color_set_ids.push_back(B_tails[i].color_set_id); // permute the vector of colors

                tlen = B_tails[i].tlen;

                if (tlen != cur_tlen){
                    cur_tlen = tlen;
                    word_index = offset/64;
                    w_offset = offset % 64;
                    sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
                    offset += 5;
//                    cerr << "TLEN = "<< tlen << endl;
                        
                    const uint32_t tnumber = m_ntails[tlen];
                    auto vb = vbyte_encode(tnumber);
                    for (uint8_t b : vb) {
                        cerr << "Writing number of tails... "<< endl;
                        cerr << "offset = " << (int)offset << endl;
                        word_index = offset/64;
                        w_offset = offset % 64;
                        //cerr << "w_offset = " << (int)w_offset << endl;
                        //cerr << "w_offset = " << (int)w_offset << endl;
  //                      cerr << "# TAILS = "<< (int)b << endl;
                        sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
                        offset += 8;
                    }
                } 
                // tails
    //            cerr << "Writing tail... "<< endl;
    //            cerr << "offset = " << (int)offset << endl;
                word_index = offset/64;
                w_offset = offset % 64;
                //cerr << "w_offset = " << (int)w_offset << endl;
                //cerr << "w_offset = " << (int)w_offset << endl;
     //           cerr << "TAIL = "<< B_tails[i].int_tail << endl;
                sdsl::bits::write_int(&data[word_index], B_tails[i].int_tail, w_offset, tlen*2);
                offset += (2 * tlen);
            }
        }
       // cerr << tail_data << endl;
        return;
    }

    void serialize(std::ostream& out) const {
        // tail_data
        sdsl::serialize(tail_data, out);

        // color_set_ids
        size_t num_ids = color_set_ids.size();
        out.write(reinterpret_cast<const char*>(&num_ids), sizeof(num_ids));
        out.write(reinterpret_cast<const char*>(color_set_ids.data()), num_ids * sizeof(uint32_t));
    }

    void load(std::istream& in) {
        // tail_data
        sdsl::load(tail_data, in);

        //color_set_ids
        size_t num_ids;
        in.read(reinterpret_cast<char*>(&num_ids), sizeof(num_ids));
        color_set_ids.resize(num_ids);
        in.read(reinterpret_cast<char*>(color_set_ids.data()), num_ids * sizeof(uint32_t));
    }



};
