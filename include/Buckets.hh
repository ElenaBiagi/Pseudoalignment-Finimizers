#pragma once
struct Compact_tails
{
    int tlen;
    uint64_t int_tail;
    uint64_t color_set_id;
};

class Bucket
{
public:
    Bucket() = default; // TODO store nullptr for missing values instead of optional

    std::unordered_map<int, int> tlen_rank_map;

    void set_tlen_order(const std::vector<uint8_t> &ordered_tlens)
    {
        tlen_rank_map.clear();
        for (size_t i = 0; i < ordered_tlens.size(); ++i)
        {
            tlen_rank_map[ordered_tlens[i]] = static_cast<int>(i);
        }
    }

    // Compressed tail data
    // tail length = tlen
    // Bit layout for evry prefix:  [tlen: 5bits][#tails: vbyte][concat of bitpacked tails: #tails * tlen * 2]
    //                              [tlen: 5bits][#tails: vbyte][concat of bitpacked tails: #tails * tlen * 2]
    //                              [tlen: 5bits][#tails: vbyte][concat of bitpacked tails: #tails * tlen * 2]
    //                                  ...
    sdsl::int_vector<1> tail_data;

    // Color set ids for each tail
    vector<uint64_t> color_set_ids;

    Bucket(vector<std::string_view> &tails, vector<uint64_t> &unsorted_color_set_ids, vector<uint8_t> &tlen_freq)
    {

        if (tails.size() == 1)
        {
            int tlen = tails[0].size();
            if (tlen == 0)
            { // if tlen is not 0, just act normally
                tail_data.resize(5);
                uint64_t *data = tail_data.data();
                color_set_ids.push_back(unsorted_color_set_ids[0]); // only one color id so nothing changed
                sdsl::bits::write_int(&data[0], 0, 0, 5);
                return;
            }
        }
        // Convert tails to Compact_tails
        vector<Compact_tails> B_tails;
        B_tails.reserve(tails.size());
        for (size_t i = 0; i < tails.size(); i++)
        {
            Compact_tails ct;
            ct.tlen = tails[i].size();
            ct.int_tail = prefix2int(tails[i], 0, ct.tlen);
            ct.color_set_id = unsorted_color_set_ids[i];
            B_tails.push_back(ct);
        }
        set_tlen_order(tlen_freq);

        // This permutes also the color_set_ids (later)
        // std::sort(B_tails.begin(), B_tails.end());
        // Sort using lambda comparator with tlen_rank_map
        std::sort(B_tails.begin(), B_tails.end(), [this](const Compact_tails &a, const Compact_tails &b)
                  {
            int rank_a = tlen_rank_map.count(a.tlen) ? tlen_rank_map[a.tlen] : INT_MAX;
            int rank_b = tlen_rank_map.count(b.tlen) ? tlen_rank_map[b.tlen] : INT_MAX;

            if (rank_a != rank_b) return rank_a < rank_b;
            return a.color_set_id < b.color_set_id; });

        WriteTailsVector(B_tails);
    }

    void WriteTailsVector(const vector<Compact_tails> &B_tails)
    {
        // Compressed tails
        // sdsl::int_vector<1> tail_data;

        // Color set ids for every tail
        // vector<uint64_t> color_set_ids;
        color_set_ids.resize(B_tails.size());

        uint64_t total_bits = 0;
        uint64_t ntails = 0;
        vector<int> tlens;
        vector<uint64_t> v_ntails;

        // tlen == 0 is a special case -> B_tails.size() == 1
        // B_tails must be of size at least 1;
        int tlen = B_tails[0].tlen;
        int cur_tlen = tlen;

        for (size_t i = 0; i < B_tails.size(); ++i)
        {
            tlen = B_tails[i].tlen;
            if (tlen != cur_tlen)
            {
                v_ntails.push_back(ntails);
                total_bits += 5;                               // cur_tlen (previous tlen)
                total_bits += vbyte_encode(ntails).size() * 8; // vbyte encoded tail count
                total_bits += ntails * cur_tlen * 2;           // each tail uses cur_tlen*2 bits
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
        tail_data.resize(total_bits + 128); // same a bit_resize

        uint64_t *data = tail_data.data();
        int64_t offset = 0;
        uint64_t word_index = 0;
        uint8_t w_offset = 0;

        cur_tlen = 0; // tlen cannot be 0

        int t = 0; // tlen index
        for (size_t i = 0; i < B_tails.size(); i++)
        {

            color_set_ids[i] = (B_tails[i].color_set_id); // permute the vector of colors

            tlen = B_tails[i].tlen;

            if (tlen != cur_tlen)
            {
                cur_tlen = tlen;
                word_index = offset / 64;
                w_offset = offset % 64;

                sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
                offset += 5;

                const uint64_t tnumber = v_ntails[t];
                t++;

                auto vb = vbyte_encode(tnumber);
                for (uint8_t b : vb)
                {
                    word_index = offset / 64;
                    w_offset = offset % 64;

                    sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
                    offset += 8;
                }
            }
            // tails
            word_index = offset / 64;
            w_offset = offset % 64;

            sdsl::bits::write_int(&data[word_index], B_tails[i].int_tail, w_offset, tlen * 2);
            offset += (2 * tlen);
        }
        return;
    }

    void serialize(std::ostream &out) const
    {
        // tail_data
        sdsl::serialize(tail_data, out);

        // color_set_ids
        size_t num_ids = color_set_ids.size();
        out.write(reinterpret_cast<const char *>(&num_ids), sizeof(num_ids));
        out.write(reinterpret_cast<const char *>(color_set_ids.data()), num_ids * sizeof(uint64_t));
    }

    void load(std::istream &in)
    {
        // tail_data
        sdsl::load(tail_data, in);
        if (tail_data.size() == 0)
        {
            std::cerr << "[ERROR] Bucket::load loaded tail_data.size()==0; file likely corrupted\n";
        }

        // color_set_ids
        size_t num_ids;
        in.read(reinterpret_cast<char *>(&num_ids), sizeof(num_ids));
        color_set_ids.resize(num_ids);
        in.read(reinterpret_cast<char *>(color_set_ids.data()), num_ids * sizeof(uint64_t));
    }
};
