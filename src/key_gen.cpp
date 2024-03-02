#include "header.h"
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;
using namespace std;

const unsigned step = 1;
unsigned kmer;
string dir;

string remove_extension(const string& fn) {
    size_t last_dot = fn.find_last_of(".");
    if (last_dot != string::npos) {
        return fn.substr(0, last_dot);
    }
    return fn;
}

// a structure of pairs (key, position)
// T is supposed to be either uint32_t or uint64_t
template<typename T>
struct Data {
  T key;
  uint32_t pos;
  Data() : key(-1), pos(-1) {}
  Data(T k, uint32_t p) : key(k), pos(p) {}

  bool operator()(const Data &X, const Data &Y) const {
    return X.key == Y.key ? X.pos < Y.pos : X.key < Y.key;
  }

};

class Index {
 private:
  string ref;
 public:
  // function to load the reference string
  // each nucleotide basis is mapped to a 8-bit value
  // A -> 0
  // C -> 1
  // G -> 2
  // T -> 3
  // else (N) -> 4
  bool load_ref(const char *F) {
    char code[256], buf[65536];
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
          ref.push_back(*(code + *p));
    }
    fclose(f);
    cerr << "genome\t" << ref.size() << '\n';
    return true;
  }
  // compute the key value by using 2 bits for each basis (N is discarded)
  // Example: ACCT -> 00-01-01-11
  template <typename T>
  void cal_key(size_t i, vector<Data<T>> &data) {
    T h = 0;
    bool hasn = false;
    for (unsigned j = 0; j < kmer; j++) {
      if (ref[i + j] == 4) {
        hasn = true;
      }
      h = (h << 2) + ref[i + j];
    }
    if (!hasn) {
      data[i / step].key = h;
      data[i / step].pos = i;
    }
  }
  // generates keys
  bool key_gen32();
  bool key_gen64();
};

template<typename T>
class Tbb_cal_key {
  vector<Data<T>> &data;
  Index *index_obj;

 public:
  Tbb_cal_key(vector<Data<T>> &_data, Index *_index_obj) :
      data(_data), index_obj(_index_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      index_obj->cal_key(i, data);
    }
  }
};

bool Index::key_gen32() {
  size_t limit = ref.size() - kmer + 1;
  size_t vsz;
  if (step == 1)
    vsz = limit;
  else
    vsz = ref.size() / step + 1;

  vector<Data<uint32_t>> data(vsz, Data<uint32_t>());
  cerr << "limit = " << limit << ", vsz = " << vsz << endl;

  // get the hash for each key
  tbb::parallel_for(tbb::blocked_range<size_t>(0, limit), Tbb_cal_key<uint32_t>(data, this));

  // sort keys
  //XXX: Parallel sort uses lots of memory. Need to fix this. In general, we
  //use 8 bytes per item. Its a waste.
  try {
    cerr << "Attempting parallel sorting\n";
    tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
    tbb::parallel_sort(data.begin(), data.end(), Data<uint32_t>());
  } catch (std::bad_alloc& e) {
    cerr << "Fall back to serial sorting (low mem)\n";
    sort(data.begin(), data.end(), Data<uint32_t>());
  }
  
  string fn = dir + "keys_uint32";

  ofstream fo_key(fn.c_str(), ios::binary);

  // determine the number of valid and unique entries
  uint32_t prec = uint32_t(-1);
  uint64_t eof = 0;
  size_t valid;   // the number of entries different than -1

  size_t i;

  for (i = 0; i < data.size() && data[i].pos != uint32_t(-1); i++) {
    if (data[i].key != prec) {
        prec = data[i].key;
        eof ++;
      }
  }
  valid = i;
  cerr << "Found " << eof << " valid keys and " << valid << " valid positions out of " <<
      data.size() << " total\n\n";
  // The number of entries is required to be a 64-bit value
  fo_key.write((char *) &eof, 8);

  // write out keys
  try {
    cerr << "Fast writing uint32 keys (" << eof << ")\n";
    size_t elements = eof*2+2;
    uint32_t *buf = new uint32_t[elements];
    // the previous value
    prec = uint32_t(-1); 
    size_t i_buf;

    for (i = 0, i_buf = 0; i < valid && i_buf < (elements-2); i++) {
      if (data[i].key != prec) {
        buf[i_buf++] = data[i].key;
        buf[i_buf++] = i;
        prec = data[i].key;
      }
    }
    ////////// add fake last element //////////
    buf[elements-2] = 0;
    buf[elements-1] = valid;
    ///////////////////////////////////////////
    fo_key.write((char *) buf, elements*sizeof(uint32_t));
    delete[] buf;

  } catch (std::bad_alloc& e) {
    cerr << "Fall back to slow writing keys due to low mem.\n";
    // the previous value
    prec = uint32_t(-1);
    uint32_t buf[2];

    for (size_t i = 0; i < valid; i++) {
      if (data[i].key != prec) {
        buf[0] = data[i].key;
        buf[1] = i;
        fo_key.write((char *) buf, 8);
        prec = data[i].key;
      }
    }
    ////////// add fake last element //////////
    buf[0] = 0;
    buf[1] = valid;
    fo_key.write((char *) buf, 8);
    ///////////////////////////////////////////
  }
  cerr << "Key generation complete!\n\n";
  fo_key.close();

  // now, write positions
  fn = dir + "pos_uint32";

  ofstream fo_pos(fn.c_str(), ios::binary);

  eof = (uint64_t) valid;
  fo_pos.write((char *) &eof, 8);

  try {
    cerr << "Fast writing posv (" << eof << ")\n";
    uint32_t *buf = new uint32_t[eof];
    for (i = 0; i < eof; i++) {
      buf[i] = data[i].pos;
    }
    fo_pos.write((char *) buf, eof * sizeof(uint32_t));
    delete[] buf;
  } catch (std::bad_alloc& e) {
    cerr << "Fall back to slow writing posv due to low mem.\n";
    for (i = 0; i < eof; i++) {
      fo_pos.write((char *) &data[i].pos, 4);
    }
  }

  cerr << "Position generation complete!\n\n";
  fo_pos.close();

  return true;
}

bool Index::key_gen64() {
  size_t limit = ref.size() - kmer + 1;
  size_t vsz;
  if (step == 1)
    vsz = limit;
  else
    vsz = ref.size() / step + 1;

  vector<Data<uint64_t>> data(vsz, Data<uint64_t>());
  cerr << "limit = " << limit << ", vsz = " << vsz << endl;

  // get the hash for each key
  tbb::parallel_for(tbb::blocked_range<size_t>(0, limit), Tbb_cal_key<uint64_t>(data, this));

  // sort keys
  //XXX: Parallel sort uses lots of memory. Need to fix this. In general, we
  //use 8 bytes per item. Its a waste.
  try {
    cerr << "Attempting parallel sorting\n";
    tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
    tbb::parallel_sort(data.begin(), data.end(), Data<uint64_t>());
  } catch (std::bad_alloc& e) {
    cerr << "Fall back to serial sorting (low mem)\n";
    sort(data.begin(), data.end(), Data<uint64_t>());
  }
  
  string fn = dir + "keys_uint64";
  ofstream fo_key(fn.c_str(), ios::binary);

  // determine the number of valid and unique entries
  uint64_t prec = uint64_t(-1);
  uint64_t eof = 0;
  size_t valid;   // the number of entries different than -1

  size_t i;

  for (i = 0; i < data.size() && data[i].pos != uint32_t(-1); i++) {
    if (data[i].key != prec) {
        prec = data[i].key;
        eof ++;
      }
  }
  valid = i;
  cerr << "Found " << eof << " valid keys and " << valid << " valid positions out of " <<
      data.size() << " total\n\n";
  // The number of entries is required to be a 64-bit value
  fo_key.write((char *) &eof, 8);

  // write out keys
  try {
    cerr << "Fast writing uint64 keys (" << eof << ")\n";
    size_t elements = eof*3+3;
    uint32_t *buf = new uint32_t[elements];
    // the previous value
    prec = uint64_t(-1);
    size_t i_buf;
    uint64_t *point;

    for (i = 0, i_buf = 0; i < valid && i_buf < (elements-3); i++) {
      if (data[i].key != prec) {
        //this is what we have to change
        point = reinterpret_cast<uint64_t*>(buf+i_buf);
        *point = data[i].key;
        i_buf += 2;
        buf[i_buf++] = i;
        prec = data[i].key;
      }
    }
    ////////// add fake last element //////////
    buf[elements-3] = 0;
    buf[elements-2] = 0;
    buf[elements-1] = valid;
    ///////////////////////////////////////////
    fo_key.write((char *) buf, elements*sizeof(uint32_t));
    delete[] buf;

  } catch (std::bad_alloc& e) {
    cerr << "Fall back to slow writing keys due to low mem.\n";
    // the previous value
    prec = uint64_t(-1);
    uint32_t buf[3];
    uint64_t *point = reinterpret_cast<uint64_t*>(buf);

    for (size_t i = 0; i < valid; i++) {
      if (data[i].key != prec) {
        *point = data[i].key;
        buf[2] = i;
        fo_key.write((char *) buf, 12);
        prec = data[i].key;
      }
    }
    ////////// add fake last element //////////
    *point = 0;
    buf[2] = valid;
    fo_key.write((char *) buf, 12);
    ///////////////////////////////////////////
  }
  cerr << "Key generation complete!\n\n";
  fo_key.close();

  // now, write positions
  fn = dir + "pos_uint32";

  ofstream fo_pos(fn.c_str(), ios::binary);

  eof = (uint64_t) valid;
  fo_pos.write((char *) &eof, 8);

  try {
    cerr << "Fast writing posv (" << eof << ")\n";
    uint32_t *buf = new uint32_t[eof];
    for (i = 0; i < eof; i++) {
      buf[i] = data[i].pos;
    }
    fo_pos.write((char *) buf, eof * sizeof(uint32_t));
    delete[] buf;
  } catch (std::bad_alloc& e) {
    cerr << "Fall back to slow writing posv due to low mem.\n";
    for (i = 0; i < eof; i++) {
      fo_pos.write((char *) &data[i].pos, 4);
    }
  }

  cerr << "Position generation complete!\n\n";
  fo_pos.close();

  return true;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "key_gen -l LEN <ref.fna>\n";
    return 0;
  }
  unsigned kmer_temp = 0;
  for (int it = 1; it < argc; it++) {
    if (strcmp(argv[it], "-l") == 0)
      kmer_temp = atoi(argv[it + 1]);
  }
  // change default to 32
  kmer = 32;
  if (kmer_temp != 0)
    kmer = kmer_temp;

  cerr << "Using kmer length " << kmer << endl;

  Index i;
  bool good;
  if (!i.load_ref(argv[argc - 1]))
    return 0;

  //make directory
  string prefix = remove_extension(string(argv[argc - 1]));
  dir = prefix + "_index" + to_string(kmer);
  if (!fs::exists(dir)) {
      fs::create_directory(dir);
  }
  dir += "/";

  if (kmer <= 16)
    good = i.key_gen32();
  else
    good = i.key_gen64();
  
  if (!good) return 1;
  // everything is fine
  return 0;
}
