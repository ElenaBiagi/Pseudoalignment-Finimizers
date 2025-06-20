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
    // tail length = tlen
    // Bit layout for evry prefix:  [tlen: 5bits][#tails: vbyte][concat of bitpacked tails: #tails * tlen * 2]
    //                              [tlen: 5bits][#tails: vbyte][concat of bitpacked tails: #tails * tlen * 2] 
    //                              [tlen: 5bits][#tails: vbyte][concat of bitpacked tails: #tails * tlen * 2] 
    //                                  ...
    sdsl::int_vector<1> tail_data;

    // Color set ids for each tail
    vector<uint32_t> color_set_ids;

    Bucket(vector<std::string_view>& tails, vector<uint32_t>& unsorted_color_set_ids) {
        
        if (tails.size() == 1){
            int tlen = tails[0].size();
            if (tlen == 0){ // if tlen is not 0, just act normally
                tail_data.resize(5);
                uint64_t* data = tail_data.data();
                color_set_ids.push_back(unsorted_color_set_ids[0]); // only one color id so nothing changed
                sdsl::bits::write_int(&data[0], 0, 0, 5);
                return;
            }
        }
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
        // This permutes also the color_set_ids (later)
        std::sort(B_tails.begin(), B_tails.end()); 
        
        WriteTailsVector(B_tails);
    }

    void WriteTailsVector(const vector<Compact_tails>& B_tails){        
        // Compressed tails
        //sdsl::int_vector<1> tail_data;
        
        // Color set ids for every tail
        //vector<uint32_t> color_set_ids;
        color_set_ids.reserve(B_tails.size()); 

        uint64_t total_bits = 0;
        uint32_t ntails = 0;
        vector <int> tlens;
        vector<uint64_t> v_ntails;


        // tlen == 0 is a special case -> B_tails.size() == 1
        // B_tails must be of size at least 1;
        int tlen = B_tails[0].tlen;
        int cur_tlen = tlen;

        for (size_t i = 0; i < B_tails.size(); ++i) {
            tlen = B_tails[i].tlen;
            if (tlen != cur_tlen){ 
                v_ntails.push_back(ntails);
                total_bits += 5; // cur_tlen (previous tlen)   
                total_bits += vbyte_encode(ntails).size() * 8; // vbyte encoded tail count
                total_bits += ntails * cur_tlen * 2; // each tail uses cur_tlen*2 bits
                ntails = 0;
                cur_tlen = tlen;

            }
            ntails++;
        }     
        // last tails
        v_ntails.push_back(ntails);

        total_bits += 5; // tlen
        total_bits += vbyte_encode(ntails).size() * 8;
 
        total_bits += ntails * tlen * 2;
        tail_data.resize(total_bits+128);// same a bit_resize

        uint64_t* data = tail_data.data();        
        int64_t offset = 0;
        uint64_t word_index = 0;
        uint8_t w_offset = 0;

        cur_tlen = 0; // tlen cannot be 0

        int t=0; // tlen index
        for (size_t i=0; i < B_tails.size(); i++){

            color_set_ids.push_back(B_tails[i].color_set_id); // permute the vector of colors

            tlen = B_tails[i].tlen;

            if (tlen != cur_tlen){
                cur_tlen = tlen;
                word_index = offset/64;
                w_offset = offset % 64;
                
                sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
                offset += 5;
                        
                const uint32_t tnumber = v_ntails[t];
                t++;

                auto vb = vbyte_encode(tnumber);
                for (uint8_t b : vb) {
                    word_index = offset/64;
                    w_offset = offset % 64;

                    sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
                    offset += 8;
                }
            } 
            // tails
            word_index = offset/64;
            w_offset = offset % 64;

            sdsl::bits::write_int(&data[word_index], B_tails[i].int_tail, w_offset, tlen*2);
            offset += (2 * tlen);
        }
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
