#include "header.h"

// An utility function to perform /path1/path2/path3 --> path_3
std::string get_last_directory(const std::string& directory_path) {
  size_t last_separator = directory_path.find_last_of("/\\");
  if (last_separator != std::string::npos) {
    return directory_path.substr(last_separator + 1);
  }
  return directory_path;
}

class RefParser {
  string &ref;
  char mode;

 public:
  RefParser(string &_ref, char _mode) : ref(_ref), mode(_mode) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    uint8_t code[256];
    for (size_t i = 0; i < 256; i++)
      code[i] = 4;

    code['A'] = code['a'] = 0;
    code['C'] = code['c'] = 1;
    code['G'] = code['g'] = 2;
    code['T'] = code['t'] = 3;
    code['N'] = code['n'] = 0;

    if (mode == 'c')
      code['C'] = code['c'] = 3;
    else if (mode == 'g')
      code['G'] = code['g'] = 0;

    for (size_t i = r.begin(); i != r.end(); ++i)
      ref[i] = *(code + ref[i]);
  }
};

// F should be the directory, that is <reference_prefix>_index32
void Reference::load_index32(const char *F) {

  string keys_f = string(F) + "/keys_uint32";
  string pos_f = string(F) + "/pos_uint32";

  ifstream fi;
  fi.open(keys_f.c_str(), ios::binary);
  if (!fi) {
    cerr << "Unable to open key file" << endl;
    exit(0);
  }
  fi.read((char *) &nkeyv_true, 8);
  nkeyv = (nkeyv_true+1)*2;
  fi.close();

  fi.open(pos_f.c_str(), ios::binary);
  if (!fi) {
    cerr << "Unable to open pos file" << endl;
    exit(0);
  }
  fi.read((char *) &nposv, 8);
  fi.close();

  cerr << "Mapping keyv of size: " << nkeyv * 4 <<
       " and posv of size " << nposv * 4 << endl;
  size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
  size_t keyv_sz = (size_t) nkeyv * sizeof(uint32_t);

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

  int fd = open(keys_f.c_str(), O_RDONLY);
  char *base = reinterpret_cast<char *>(mmap(NULL, 8 + keyv_sz, PROT_READ, MMAP_FLAGS, fd, 0));
  assert(base != MAP_FAILED);
  keyv = (uint32_t * )(base + 8);

  // cerr << "Printing first 4 entries" << endl;
  // cerr << "------ keyv ------" << endl;
  // cerr << keyv[0] << " [" << keyv[1] << " pos]" << endl;
  // cerr << keyv[2] << " [" << keyv[3] << " pos]" << endl;
  // cerr << keyv[4] << " [" << keyv[5] << " pos]" << endl;
  // cerr << keyv[6] << " [" << keyv[7] << " pos]" << endl;

  fd = open(pos_f.c_str(), O_RDONLY);
  base = reinterpret_cast<char *>(mmap(NULL, 8 + posv_sz, PROT_READ, MMAP_FLAGS, fd, 0));
  assert(base != MAP_FAILED);
  posv = (uint32_t * )(base + 8);

  // print first 4 entries
  // cerr << "------ posv ------" << endl;
  // cerr << posv[0] << endl;
  // cerr << posv[1] << endl;
  // cerr << posv[2] << endl;
  // cerr << posv[3] << endl;
  // cerr << "------------------" << endl;

  cerr << "Mapping done" << endl;
  cerr << "done loading hashtable\n";
}

// F should be the directory, that is <reference_prefix>_index64
void Reference::load_index64(const char *F) {

  string keys_f = string(F) + "/keys_uint64";
  string pos_f = string(F) + "/pos_uint32";

  ifstream fi;
  fi.open(keys_f.c_str(), ios::binary);
  if (!fi) {
    cerr << "Unable to open key file" << endl;
    exit(0);
  }
  fi.read((char *) &nkeyv_true, 8);
  nkeyv = (nkeyv_true+1) * 3;
  fi.close();

  fi.open(pos_f.c_str(), ios::binary);
  if (!fi) {
    cerr << "Unable to open pos file" << endl;
    exit(0);
  }
  fi.read((char *) &nposv, 8);
  fi.close();

  cerr << "Mapping keyv of size: " << nkeyv * 4 <<
       " and posv of size " << nposv * 4 << endl;
  size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
  size_t keyv_sz = (size_t) nkeyv * sizeof(uint32_t);

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

  int fd = open(keys_f.c_str(), O_RDONLY);
  char *base = reinterpret_cast<char *>(mmap(NULL, 8 + keyv_sz, PROT_READ, MMAP_FLAGS, fd, 0));
  assert(base != MAP_FAILED);
  keyv = (uint32_t * )(base + 8);

  // cerr << "Printing first 4 entries" << endl;
  // cerr << "------ keyv ------" << endl;
  // cerr << *((uint64_t*)(keyv+0)) << " [" << keyv[2] << " pos]" << endl;
  // cerr << *((uint64_t*)(keyv+3)) << " [" << keyv[5] << " pos]" << endl;
  // cerr << *((uint64_t*)(keyv+6)) << " [" << keyv[8] << " pos]" << endl;
  // cerr << *((uint64_t*)(keyv+9)) << " [" << keyv[11] << " pos]" << endl;

  fd = open(pos_f.c_str(), O_RDONLY);
  base = reinterpret_cast<char *>(mmap(NULL, 8 + posv_sz, PROT_READ, MMAP_FLAGS, fd, 0));
  assert(base != MAP_FAILED);
  posv = (uint32_t * )(base + 8);

  // print first 4 entries
  // cerr << "------ posv ------" << endl;
  // cerr << posv[0] << endl;
  // cerr << posv[1] << endl;
  // cerr << posv[2] << endl;
  // cerr << posv[3] << endl;
  // cerr << "------------------" << endl;

  cerr << "Mapping done" << endl;
  cerr << "done loading hashtable\n";
}

void Reference::load_index_classic(const char *F) {
  string fn;
  if (mode == ' ')
    fn = string(F) + ".hash" + to_string(kmer_len);
  else if (mode == 'c')
    fn = string(F) + ".hash" + to_string(kmer_len) + ".part1";
  else if (mode == 'g')
    fn = string(F) + ".hash" + to_string(kmer_len) + ".part2";

  cerr << "loading hashtable from " << fn << endl;
  ifstream fi;
  fi.open(fn.c_str(), ios::binary);
  if (!fi) {
    cerr << "Unable to open index file " << fn << endl;
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
  nkeyv_true = mod + 1;
  nkeyv = nkeyv_true;
  bind_xxhash(xxh_type, xxh);

  cerr << "Mapping keyv of size: " << nkeyv * 4 <<
       " and posv of size " << (size_t) nposv * 4 <<
       " from index file " << fn << endl;
  cerr << "using MOD = " << mod << " and XXH = " << xxh_type << endl;

  size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
  uint64_t keyv_sz = (uint64_t) nkeyv * sizeof(uint32_t);
  int fd = open(fn.c_str(), O_RDONLY);

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

void Reference::load_reference(const char *F){
  size_t ref_size = 0;
  string bref(F);
  bref += ".bref";
  ifstream brefin(bref.c_str(), ios::in);
  if (brefin) {
    /* Experimentation revealed that using prebuild binary ref leads to
     * embeding procedure taking more latency on memory stalls. Might have
     * something to do with ref now being loaded in cache. So while support
     * is there so that we can test it on other servers, it is disabled.
     */
    cout << "WARNING: USING PREBUILD BINARY REF.\n";
    cout << "-----------------------------------\n";
    cout << "Loading index and binary ref from " << bref << endl;
    size_t count;
    brefin >> count;
    cout << "Done loading " << count << " names " << endl;
    for (size_t i = 0; i < count; ++i) {
      string tmp;
      brefin >> tmp;
      name.push_back(tmp);
    }
    brefin >> count;
    for (size_t i = 0; i < count; ++i) {
      uint32_t tmp;
      brefin >> tmp;
      offset.push_back(tmp);
    }
    cout << "Done loading " << count << " offsets " << endl;
    brefin >> ref_size;
    ref.resize(ref_size);
    brefin.read(&ref[0], ref_size);
    cout << "Done loading binary ref of sz " << ref_size << endl;
    brefin.close();
  } else {
    FILE *f = fopen(F, "rb");
    fseek(f, 0, SEEK_END);
    ref.reserve(ftell(f) + 1);
    ref = "";
    fclose(f);

    ifstream fi(F);
    if (!fi) {
      cerr << "fail to open " << F << '\n';
      return;
    }
    cerr << "Loading reference " << F << "\n";
    string buf = "";
    while (!fi.eof()) {
      getline(fi, buf);
      if (!buf.size())
        continue;

      if (buf[0] == '>') {
        string tmp(buf.c_str() + 1);
        istringstream iss(tmp);
        vector<string> parsed((istream_iterator<string>(iss)),
                              istream_iterator<string>());
        //cerr<<"Loading chromosome " << parsed[0] << "\n";
        name.push_back(parsed[0]);
        offset.push_back(ref_size);
      } else {
        if (isspace(buf[buf.size()-1])) { // isspace(): \t, \n, \v, \f, \r
          ref += buf.substr(0, buf.size()-1);
          ref_size += buf.size() - 1;
        } else {
          ref += buf;
          ref_size += buf.size();
        }
      }
    }

    // process last chromosome
    cerr << "Parsing reference.\n";
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, ref_size),
        RefParser(ref, mode));
    fi.close();
    offset.push_back(ref_size);
    cerr << "Loaded reference size: " << ref_size << endl;
  }
}

void Reference::index_rmi_lookup32(uint64_t key, size_t* b, size_t* e) {
  size_t err;
  uint32_t guess_key;
  uint64_t guess_pos;
  uint64_t l, r;
  uint32_t key32 = (uint32_t) key;
  // call the lookup function of the index
  guess_pos = rmi.lookup(key, &err);

  // set up l and r for the bounded binary search
  l = guess_pos < err? 0 : (guess_pos-err);
  r = (guess_pos+err) < (nkeyv_true-1)? (guess_pos+err) : (nkeyv_true-1);

  // check in the keyv array
  while (l <= r) {
    guess_key = keyv[guess_pos*2];
    // if it's the same, done
    if (guess_key == key32) {
      *b = get_keyv_val(guess_pos);
      *e = get_keyv_val(guess_pos+1);
      return;
    }
    // else, do binary search
    if (guess_key < key32) {
      l = guess_pos + 1;
    } else {
      r = guess_pos - 1;
    }
    // update guess_pos
    guess_pos = l + (r-l)/2;
  }
  // not found :(
  *b = 0;
  *e = 0;
}

void Reference::index_rmi_lookup64(uint64_t key, size_t* b, size_t* e) {
  size_t err;
  uint64_t* guess_key;
  uint64_t guess_pos, l, r;
  // call the lookup function of the index
  guess_pos = rmi.lookup(key, &err);

  // set up l and r for the bounded binary search
  l = guess_pos < err? 0 : (guess_pos-err);
  r = (guess_pos+err) < (nkeyv_true-1)? (guess_pos+err) : (nkeyv_true-1);

  // check in the keyv array
  while (l <= r) {
    guess_key = reinterpret_cast<uint64_t*>(keyv+(guess_pos*3));
    // if it's the same, done
    if (*guess_key == key) {
      *b = get_keyv_val(guess_pos);
      *e = get_keyv_val(guess_pos+1);
      return;
    }
    // else, do binary search
    if (*guess_key < key) {
      l = guess_pos + 1;
    } else {
      r = guess_pos - 1;
    }
    // update guess_pos
    guess_pos = l + (r-l)/2;
  }
  // not found :(
  *b = 0;
  *e = 0;
}

void Reference::index_bin_lookup32(uint64_t key, size_t* b, size_t* e) {
  uint32_t guess_key;
  uint64_t guess_pos;
  uint64_t l, r;
  uint32_t key32 = (uint32_t) key;

  // set up l and r for the bounded binary search
  l = 0;
  r = nkeyv_true-1;
  guess_pos = l + (r-l)/2;

  // check in the keyv array
  while (l <= r) {
    guess_key = keyv[guess_pos*2];
    // if it's the same, done
    if (guess_key == key32) {
      *b = get_keyv_val(guess_pos);
      *e = get_keyv_val(guess_pos+1);
      return;
    }
    // else, do binary search
    if (guess_key < key32) {
      l = guess_pos + 1;
    } else {
      r = guess_pos - 1;
    }
    // update guess_pos
    guess_pos = l + (r-l)/2;
  }
  // not found :(
  *b = 0;
  *e = 0;
}

void Reference::index_bin_lookup64(uint64_t key, size_t* b, size_t* e) {
  uint64_t* guess_key;
  uint64_t guess_pos, l, r;

  l = 0;
  r = nkeyv_true-1;
  guess_pos = l + (r-l)/2;

  // check in the keyv array
  while (l <= r) {
    guess_key = reinterpret_cast<uint64_t*>(keyv+(guess_pos*3));
    // if it's the same, done
    if (*guess_key == key) {
      *b = get_keyv_val(guess_pos);
      *e = get_keyv_val(guess_pos+1);
      return;
    }
    // else, do binary search
    if (*guess_key < key) {
      l = guess_pos + 1;
    } else {
      r = guess_pos - 1;
    }
    // update guess_pos
    guess_pos = l + (r-l)/2;
  }
  // not found :(
  *b = 0;
  *e = 0;
}

void Reference::index_lookup_classic(uint64_t key, size_t* b, size_t* e) {
  // FIXME -- do we need the mask?
  uint64_t hash = uint64_t(xxh(&key) % mod);
  *b = keyv[hash];
  *e = keyv[hash+1];
}

uint32_t Reference::get_keyv_val32(uint64_t idx) {
  return keyv[idx*2+1];
}
uint32_t Reference::get_keyv_val64(uint64_t idx) {
  return keyv[idx*3+2];
}

Reference::Reference(const char *F, unsigned _kmer_len, SType _g_stype, IndexType _index_type, char _mode):
    kmer_len(_kmer_len),
    g_stype(_g_stype),
    index_type(_index_type),
    mode(_mode) {
  auto start = std::chrono::system_clock::now();

  if (g_stype == SType::Minimizer){
    int n_threads = 3;

    mm_idxopt_t ipt;
    mm_idxopt_init(&ipt);
    string fn = string(F) + ".hash";
    const char* fnw = fn.c_str();

    mm_idx_reader_t *idx_rdr = mm_idx_reader_open(fnw, &ipt, nullptr);
    mi = mm_idx_reader_read(idx_rdr, n_threads);

    load_reference(F);
  } else{
    string F_index;
    ////// case HASH //////
    if (index_type == IndexType::HASH_IDX) {
      load_index = std::bind(&Reference::load_index_classic, this, std::placeholders::_1);
      index_lookup = std::bind(&Reference::index_lookup_classic, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
      F_index = F;
    } else {
      // common
      unsigned bit_len = kmer_len > 16? 64 : 32;
      if (bit_len==32) {
        load_index = std::bind(&Reference::load_index32, this, std::placeholders::_1);
        get_keyv_val = std::bind(&Reference::get_keyv_val32, this, std::placeholders::_1);
      } else {
        load_index = std::bind(&Reference::load_index64, this, std::placeholders::_1);
        get_keyv_val = std::bind(&Reference::get_keyv_val64, this, std::placeholders::_1);
      }
      // F          ./data/hg37.fna
      // F_prefix   ./data/hg37
      // F_index    ./data/hg37_index32
      // F_library  ./data/hg37_index32/hg37_index
      string F_prefix = string(F).substr(0, string(F).find_last_of("."));
      F_index = F_prefix  + "_index" + to_string(kmer_len);
      ////// case RMI //////
      if (index_type == IndexType::RMI_IDX) {
        // bind index lookup
        if (bit_len==32)
          index_lookup = std::bind(&Reference::index_rmi_lookup32, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        else
          index_lookup = std::bind(&Reference::index_rmi_lookup64, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

        string F_library = F_index + "/" + get_last_directory(F_index);
        rmi.init(F_library.c_str());
      } else {
        ////// case BINARY //////
        // bind index lookup
        if (bit_len==32)
          index_lookup = std::bind(&Reference::index_bin_lookup32, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        else
          index_lookup = std::bind(&Reference::index_bin_lookup64, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
      }
    }
    thread t([this, F_index]() {load_index(F_index.c_str());});
    load_reference(F);

    t.join(); // wait for index load to finish
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  cerr << "Setup reference in " << elapsed.count() / 1000000 << " secs\n";
}

Reference::~Reference() {
  if (g_stype == SType::Minimizer){
    mm_idx_destroy(mi);
  } else {
    size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
    uint64_t keyv_sz = (uint64_t) nkeyv * sizeof(uint32_t);
    int r;
    if (index_type == IndexType::HASH_IDX) {
      char *base = (char *) posv - 12;
      r = munmap(base, posv_sz + keyv_sz + 12);
      assert(r == 0);
    } else {
      char *base = (char *) keyv - 8;
      r = munmap(base, keyv_sz + 8);
      assert(r == 0);
      base = (char *) posv - 8;
      r = munmap(base, posv_sz + 8);
      assert(r == 0);
    }
  }
}


