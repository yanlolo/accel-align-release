#include <cstdint>
#include <tuple>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <chrono>
#include <zlib.h>
#include <cmath>
#include <string>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <sys/time.h>
#include <sys/mman.h>
#include <sched.h>
#include <fcntl.h>
#include <omp.h>
#include "hash.hpp"

using namespace std;

const unsigned step = 1;
unsigned kmer_size = 32;
uint64_t mod;
uint32_t xxh_type = 0;
XXHash xxh;
bool rev_comp = true;

// Reference string + some methods
class Reference {
 private:
    string ref;
    uint32_t *keyv, *posv;
    uint32_t nposv;
    uint64_t nkeyv;
    char code[256];

    static bool ichar_equals(char a, char b) {
        return tolower(static_cast<unsigned char>(a)) ==
            tolower(static_cast<unsigned char>(b));
    }

    static bool iequals(const std::string& a, const std::string& b) {
        return a.size() == b.size() &&
            equal(a.begin(), a.end(), b.begin(), ichar_equals);
    }

    static char get_complement(char base) {
        char ubase = toupper(static_cast<unsigned char>(base));
        switch(ubase) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            default:
                // this should never happen in practice
                return ubase;
        }
    }

 public:
    // Load the reference string
    // code --> just a dictionary matching A-0, C-1, G-2, T-3, everything else-4
    bool load_ref(const char *F) {
        char buf[65536];
        for (size_t i = 0; i < 256; i++)
            code[i] = 4;
        code['A'] = code['a'] = 0;
        code['C'] = code['c'] = 1;
        code['G'] = code['g'] = 2;
        code['T'] = code['t'] = 3;

        cerr << "Loading ref\n";
        FILE *f = fopen(F, "rb");
        if (f == NULL)
            return false;
        fseek(f, 0, SEEK_END);
        ref.reserve(ftell(f) + 1);
        fclose(f);
        f = fopen(F, "rt");
        if (f == NULL)
            return false;
        while (fgets(buf, 65536, f) != NULL) {
            if (buf[0] == '>')
            continue;
            for (char *p = buf; *p; p++)
            if (*p >= 33)
                ref.push_back(*p);
        }
        fclose(f);
        cerr << "genome\t" << ref.size() << '\n';
        return true;
    }

    void load_index(const char *F) {
        cerr << "loading hashtable" << endl;
        ifstream fi;
        fi.open(F, ios::binary);
        if (!fi) {
            cerr << "Unable to open index file" << endl;
            exit(0);
        }
        fi.read((char *) &mod, 4);
        fi.read((char *) &xxh_type, 4);
        fi.read((char *) &nposv, 4);
        fi.close();
        // try to detect whether the index support the mod or not
        if ((mod < 128 && mod != 0) || (xxh_type!=0 && xxh_type!=32 && xxh_type!=64)) {
            cerr << "It looks like you are using an old index.\n";
            cerr << "Please, run again ./accindex, and come back later!\n";
            exit(1);
        }
        if (mod == 0)
            mod = MOD_32;
        nkeyv = mod + 1;
        bind_xxhash(xxh_type, xxh);

        cerr << "Mapping keyv of size: " << nkeyv * 4 <<
            " and posv of size " << (size_t) nposv * 4 <<
            " from index file" << endl;
        cerr << "using MOD = " << mod << " and XXH = " << xxh_type << endl;

        size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
        uint64_t keyv_sz = (uint64_t) nkeyv * sizeof(uint32_t);
        int fd = open(F, O_RDONLY);

        #if __linux__
        #include <linux/version.h>
        #if LINUX_VERSION_CODE > KERNEL_VERSION(2,6,22)
        #define _MAP_POPULATE_AVAILABLE
        #endif
        #endif

        #ifdef _MAP_POPULATE_AVAILABLE
        #define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
        #else
        #define MMAP_FLAGS MAP_PRIVATE
        #endif

        char *base = reinterpret_cast<char *>(mmap(NULL, 12 + posv_sz + keyv_sz, PROT_READ, MMAP_FLAGS, fd, 0));
        assert(base != MAP_FAILED);
        posv = (uint32_t * )(base + 12);
        keyv = posv + nposv;
        cerr << "Mapping done" << endl;
        cerr << "done loading hashtable\n";
    }

    uint64_t to_hash(string const& kmer) {
        uint64_t h = 0;
        for (unsigned j = 0; j < kmer_size; j++) {
            h = (h << 2) + code[(int)kmer[j]];
        }
        return uint64_t(xxh(&h) % mod);
    }

    // first element is number of predicted positions, second element is number of actual positions
    pair<size_t,size_t> stats(string const& kmer) {
        // compute key
        uint64_t k = to_hash(kmer);
        uint32_t first_pos_idx, last_pos_idx;
        // get first and last positions
        first_pos_idx = keyv[k];
        last_pos_idx = keyv[k+1];

        pair<size_t,size_t> output(last_pos_idx-first_pos_idx,0);

        size_t actual_pos = 0;

        #pragma omp parallel for reduction(+:actual_pos)
        for (auto i=first_pos_idx; i<last_pos_idx; i++) {
            // get position
            uint32_t pos = posv[i];
            if (iequals(kmer, ref.substr(pos, kmer_size)))
                actual_pos++;
        }
        output.second = actual_pos;
        return output;
    }

    static string reverse_complement(const string& dna) {
        string reversed = dna;
        string complemented;
        // Reverse
        reverse(reversed.begin(), reversed.end());
        // Complement
        transform(reversed.begin(), reversed.end(), back_inserter(complemented), get_complement);
        return complemented;
    }

    ~Reference() {
        size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
        uint64_t keyv_sz = (uint64_t) nkeyv * sizeof(uint32_t);
        char *base = (char *) posv - 12;
        int r = munmap(base, 12 + posv_sz + keyv_sz);
        assert(r == 0);
    }
};

int main(int ac, char **av) {
    if (ac < 3) {
        cerr << "./stats [-l LEN] <ref.fa> <read.fastq>\n";
        return 0;
    }

    unsigned kmer_temp = 0;
    for (int it = 1; it < (ac-2); it++) {
        if (strcmp(av[it], "-l") == 0) {
            kmer_temp = atoi(av[it + 1]);
            it++;
        }
    }

    if (kmer_temp != 0)
        kmer_size = kmer_temp;

    string ref_fn = av[ac-2];            // reference file name
    string idx_fn = ref_fn + ".hash" + to_string(kmer_size); // index file name
    string read_fn = av[ac-1];           // read file name
    string stats_fn = read_fn + ".stats" + to_string(kmer_size);   // stats file name

    Reference ref;
    // Load the reference
    if (!ref.load_ref(ref_fn.c_str()))
        return 1;
    // Load the index
    ref.load_index(idx_fn.c_str());

    // Prepare the stats file
    ofstream stats_f(stats_fn.c_str());
    // Print header
    if (rev_comp)
        stats_f << "read_name,indexed_pos_fwd,actual_pos_fwd,indexed_pos_rev,actual_pos_rev\n";
    else
        stats_f << "read_name,indexed_pos_fwd,actual_pos_fwd\n";

    // Consume the reads
    ifstream read_f(read_fn.c_str());
    string read, name, tmp;
    size_t read_len, limit;
    while(true) {
        // name
        getline(read_f, name);
        if (read_f.eof())
            break;
        // read
        getline(read_f, read);
        // other useless lines
        getline(read_f, tmp);
        getline(read_f, tmp);
        // get read length
        read_len = read.size();
        // get the number of kmers
        limit = read_len - kmer_size + 1;
        size_t pred_pos_fwd=0, actual_pos_fwd=0, pred_pos_rev=0, actual_pos_rev=0;

        #pragma omp parallel for reduction(+:actual_pos_fwd, pred_pos_fwd, actual_pos_rev, pred_pos_rev)
        for (size_t i=0; i<limit; i+=kmer_size) {
            string kmer = read.substr(i,kmer_size);
            auto positions = ref.stats(kmer);
            pred_pos_fwd += positions.first;
            actual_pos_fwd += positions.second;
            if (!rev_comp)
                continue;
            /////////////// if reverse complement ///////////////
            kmer = Reference::reverse_complement(kmer);
            positions = ref.stats(kmer);
            pred_pos_rev += positions.first;
            actual_pos_rev += positions.second;
        }

        if (rev_comp)
            stats_f << name + "," << pred_pos_fwd << "," << actual_pos_fwd << "," << pred_pos_rev << "," << actual_pos_rev << "\n";
        else
            stats_f << name + "," << pred_pos_fwd << "," << actual_pos_fwd << "\n";
        cout << name + "\n";
    }
    return 0;
}