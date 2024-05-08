#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <thread>
#include <cassert>
#include <iomanip>
#include <chrono>
#include <experimental/filesystem>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "../strobealign/refs.hpp"
#include "../strobealign/exceptions.hpp"
#include "../strobealign/strobe-index.hpp"
#include "../strobealign/pc.hpp"
#include "../strobealign/aln.hpp"
#include "../strobealign/logger.hpp"
#include "../strobealign/timer.hpp"
#include "../strobealign/readlen.hpp"
//#include "strobealign/version.hpp"
//#include "strobealign/buildconfig.hpp"


#include "../include/header.h"
using namespace std;
namespace fs = std::experimental::filesystem;

uint64_t mod = MOD_29;    // default value is 2^29 - 1
uint32_t mod_tmp;
const unsigned step = 1;
uint32_t xxh_type = 0;
XXHash xxh;
unsigned kmer;
bool enable_idx_minimizer = false, enable_bs = false; //short for bisulfite reads
struct Data {
  uint32_t key, pos;
  Data() : key(-1), pos(-1) {}
  Data(uint32_t k, uint32_t p) : key(k), pos(p) {}
  bool operator()(const Data &X, const Data &Y) const {
    return X.key == Y.key ? X.pos < Y.pos : X.key < Y.key;
  }
};
class Index {
 private:
  string ref;
 public:
  bool load_ref(const char *F, char mode);
  bool make_index(const char *F, int id);
  void cal_key(size_t i, vector<Data> &data);
};

bool Index::load_ref(const char *F, char mode) {
  char code[256], buf[65536];
  for (size_t i = 0; i < 256; i++)
    code[i] = 4;
  code['A'] = code['a'] = 0;
  code['C'] = code['c'] = 1;
  code['G'] = code['g'] = 2;
  code['T'] = code['t'] = 3;
  if (mode == 'c')
    code['C'] = code['c'] = 3;
  else if (mode == 'g')
    code['G'] = code['g'] = 0;
  cerr << "Loading ref\n";
  FILE *f = fopen(F, "rb");
  if (f == NULL){
    printf("Error: File '%s' not found.\n", F);
    return false;
  }
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

void Index::cal_key(size_t i, vector<Data> &data) {
  uint64_t h = 0;
  bool hasn = false;
  for (unsigned j = 0; j < kmer; j++) {
    if (ref[i + j] == 4) {
      hasn = true;
    }
    h = (h << 2) + ref[i + j];
  }
  if (!hasn) {
    data[i / step].key = uint32_t(xxh(&h) % mod);
    data[i / step].pos = i;
  }
}
class Tbb_cal_key {
  vector<Data> &data;
  Index *index_obj;
 public:
  Tbb_cal_key(vector<Data> &_data, Index *_index_obj) :
      data(_data), index_obj(_index_obj) {}
  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      index_obj->cal_key(i, data);
    }
  }
};
bool Index::make_index(const char *F, int id) {
  size_t limit = ref.size() - kmer + 1;
  size_t vsz;
  if (step == 1)
    vsz = limit;
  else
    vsz = ref.size() / step + 1;
  vector<Data> data(vsz, Data());
  cerr << "hashing :limit = " << limit << ", vsz = " << vsz << endl;
  cerr << "using MOD = " << mod << " and XXH = " << xxh_type << endl;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, limit), Tbb_cal_key(data, this));
  cerr << "hash\t" << data.size() << endl;
  //XXX: Parallel sort uses lots of memory. Need to fix this. In general, we
  //use 8 bytes per item. Its a waste.
  try {
    cerr << "Attempting parallel sorting\n";
    tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
    tbb::parallel_sort(data.begin(), data.end(), Data());
  } catch (std::bad_alloc &e) {
    cerr << "Fall back to serial sorting (low mem)\n";
    sort(data.begin(), data.end(), Data());
  }
  cerr << "writing\n";
  string fn = F;
  if (id)
    fn += ".hash" + to_string(kmer) + ".part" + to_string(id);
  else
    fn += ".hash" + to_string(kmer);
  ofstream fo(fn.c_str(), ios::binary);
  // determine the number of valid entries based on first junk entry
  auto joff = std::lower_bound(data.begin(), data.end(), Data(-1, -1), Data());
  size_t eof = joff - data.begin();
  cerr << "Found " << eof << " valid entries out of " <<
       data.size() << " total\n";
  // first, write mod
  fo.write((char *) &mod, 4);
  // then, write xxh_type
  fo.write((char *) &xxh_type, 4);
  // then, write the number of positions
  fo.write((char *) &eof, 4);
  // write out positions
  for (size_t i = eof; i < data.size(); i++)
    assert(data[i].key == (uint32_t) -1);
  try {
    cerr << "Fast writing posv (" << eof << ")\n";
    uint32_t *buf = new uint32_t[eof];
    for (size_t i = 0; i < eof; i++) {
      buf[i] = data[i].pos;
    }
    fo.write((char *) buf, eof * sizeof(uint32_t));
    delete[] buf;
  } catch (std::bad_alloc &e) {
    cerr << "Fall back to slow writing posv due to low mem.\n";
    for (size_t i = 0; i < eof; i++) {
      fo.write((char *) &data[i].pos, 4);
    }
  }

  // write out keys
  size_t last_key = 0, offset;
  try {
    cerr << "Fast writing keyv\n";
    uint64_t buf_idx = 0;
    uint32_t *buf = new uint32_t[mod + 1];
    // for each position
    for (size_t i = 0; i < eof;) {
      assert (data[i].pos != (uint32_t) -1);
      size_t h = data[i].key, n;
      offset = i;
      for (size_t j = last_key; j <= h; j++) {
        buf[buf_idx] = offset;
        ++buf_idx;
      }
      last_key = h + 1;
      for (n = i + 1; n < eof && data[n].key == h; n++);
      i = n;
    }
    offset = eof;
    for (uint64_t j = (uint64_t) last_key; j <= mod; j++) {
      buf[buf_idx] = offset;
      ++buf_idx;
    }
    assert(buf_idx == (mod + 1));
    fo.write((char *) buf, buf_idx * sizeof(uint32_t));
    delete[] buf;
  } catch (std::bad_alloc &e) {
    cerr << "Fall back to slow writing keyv (low mem)\n";
    for (size_t i = 0; i < eof;) {
      assert (data[i].pos != (uint32_t) -1);
      size_t h = data[i].key, n;
      offset = i;
      for (size_t j = last_key; j <= h; j++) {
        fo.write((char *) &offset, 4);
      }
      last_key = h + 1;
      for (n = i + 1; n < eof && data[n].key == h; n++);
      i = n;
    }
    offset = eof;
    for (uint64_t j = (uint64_t) last_key; j <= mod; j++) {
      fo.write((char *) &offset, 4);
    }
  }
  cerr << "Indexing complete\n";
  fo.close();
  return true;
}

static Logger& logger = Logger::get();

void log_parameters(const IndexParameters& index_parameters, const MappingParameters& map_param, const AlignmentParameters& aln_params) {
    logger.debug() << "Using" << std::endl
        << "k: " << index_parameters.syncmer.k << std::endl
        << "s: " << index_parameters.syncmer.s << std::endl
        << "w_min: " << index_parameters.randstrobe.w_min << std::endl
        << "w_max: " << index_parameters.randstrobe.w_max << std::endl
        << "Read length (r): " << map_param.r << std::endl
        << "Maximum seed length: " << index_parameters.randstrobe.max_dist + index_parameters.syncmer.k << std::endl
        << "R: " << map_param.rescue_level << std::endl
        << "Expected [w_min, w_max] in #syncmers: [" << index_parameters.randstrobe.w_min << ", " << index_parameters.randstrobe.w_max << "]" << std::endl
        << "Expected [w_min, w_max] in #nucleotides: [" << (index_parameters.syncmer.k - index_parameters.syncmer.s + 1) * index_parameters.randstrobe.w_min << ", " << (index_parameters.syncmer.k - index_parameters.syncmer.s + 1) * index_parameters.randstrobe.w_max << "]" << std::endl
        << "A: " << aln_params.match << std::endl
        << "B: " << aln_params.mismatch << std::endl
        << "O: " << aln_params.gap_open << std::endl
        << "E: " << aln_params.gap_extend << std::endl
        << "end bonus: " << aln_params.end_bonus << '\n';
}

bool avx2_enabled() {
#ifdef __AVX2__
    return true;
#else
    return false;
#endif
}

string remove_extension(const string& fn) {
  size_t last_dot = fn.find_last_of(".");
  if (last_dot != string::npos) {
    return fn.substr(0, last_dot);
  }
  return fn;
}

int run_strobealign(int argc, char **argv) {
    auto opt = parse_command_line_arguments(argc, argv);

    logger.set_level(opt.verbose ? LOG_DEBUG : LOG_INFO);
    logger.info() << std::setprecision(2) << std::fixed;
    logger.debug() << "AVX2 enabled: " << (avx2_enabled() ? "yes" : "no") << '\n';

    if (opt.c >= 64 || opt.c <= 0) {
        throw BadParameter("c must be greater than 0 and less than 64");
    }
    InputBuffer input_buffer = get_input_buffer(opt);
    if (!opt.r_set && !opt.reads_filename1.empty()) {
        opt.r = estimate_read_length(input_buffer);
        logger.info() << "Estimated read length: " << opt.r << " bp\n";
    }
    input_buffer.rewind_reset();
    IndexParameters index_parameters = IndexParameters::from_read_length(
        opt.r,
        opt.k_set ? opt.k : IndexParameters::DEFAULT,
        opt.s_set ? opt.s : IndexParameters::DEFAULT,
        opt.l_set ? opt.l : IndexParameters::DEFAULT,
        opt.u_set ? opt.u : IndexParameters::DEFAULT,
        opt.c_set ? opt.c : IndexParameters::DEFAULT,
        opt.max_seed_len_set ? opt.max_seed_len : IndexParameters::DEFAULT
    );
    logger.debug() << index_parameters << '\n';
    AlignmentParameters aln_params;
    aln_params.match = opt.A;
    aln_params.mismatch = opt.B;
    aln_params.gap_open = opt.O;
    aln_params.gap_extend = opt.E;
    aln_params.end_bonus = opt.end_bonus;

    MappingParameters map_param;
    map_param.r = opt.r;
    map_param.max_secondary = opt.max_secondary;
    map_param.dropoff_threshold = opt.dropoff_threshold;
    map_param.rescue_level = opt.rescue_level;
    map_param.max_tries = opt.max_tries;
    map_param.is_sam_out = opt.is_sam_out;
    map_param.cigar_ops = opt.cigar_eqx ? CigarOps::EQX : CigarOps::M;
    map_param.output_unmapped = opt.output_unmapped;
    map_param.details = opt.details;
    map_param.verify();

    log_parameters(index_parameters, map_param, aln_params);
    logger.debug() << "Threads: " << opt.n_threads << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");

    // Create index
    References references;
    Timer read_refs_timer;
    references = References::from_fasta(opt.ref_filename);
    logger.info() << "Time reading reference: " << read_refs_timer.elapsed() << " s\n";

    logger.info() << "Reference size: " << references.total_length() / 1E6 << " Mbp ("
        << references.size() << " contig" << (references.size() == 1 ? "" : "s")
        << "; largest: "
        << (*std::max_element(references.lengths.begin(), references.lengths.end()) / 1E6) << " Mbp)\n";
    if (references.total_length() == 0) {
        throw InvalidFasta("No reference sequences found");
    }

    StrobemerIndex index(references, index_parameters, opt.bits);
    logger.debug() << "Bits used to index buckets: " << index.get_bits() << "\n";
    logger.info() << "Indexing ...\n";
    index.populate(opt.f, opt.n_threads);

    logger.debug()
        << "Index statistics\n"
        << "  Total strobemers:    " << std::setw(14) << index.stats.tot_strobemer_count << '\n'
        << "  Distinct strobemers: " << std::setw(14) << index.stats.distinct_strobemers << " (100.00%)\n"
        << "    1 occurrence:      " << std::setw(14) << index.stats.tot_occur_once
            << " (" << std::setw(6) << (100.0 * index.stats.tot_occur_once / index.stats.distinct_strobemers) << "%)\n"
        << "    2..100 occurrences:" << std::setw(14) << index.stats.tot_mid_ab
            << " (" << std::setw(6) << (100.0 * index.stats.tot_mid_ab / index.stats.distinct_strobemers) << "%)\n"
        << "    >100 occurrences:  " << std::setw(14) << index.stats.tot_high_ab
            << " (" << std::setw(6) << (100.0 * index.stats.tot_high_ab / index.stats.distinct_strobemers) << "%)\n"
        ;
    if (opt.only_gen_index) {
        Timer index_writing_timer;
        std::string sti_path = opt.ref_filename + index_parameters.filename_extension();
        logger.info() << "Writing index to " << sti_path << '\n';
        index.write(opt.ref_filename + index_parameters.filename_extension());
        logger.info() << "Total time writing index: " << index_writing_timer.elapsed() << " s\n";
//        return EXIT_SUCCESS;
    }

  //make directory
  string prefix = remove_extension(opt.ref_filename);
  string dir = prefix + "_index" + to_string(32);
  if (!fs::exists(dir)) {
    fs::create_directory(dir);
  }
  dir += "/";

  vector<RefRandstrobe> data = index.randstrobes;
  std::sort(data.begin(), data.end());
  string fn =  dir + "/keys_uint64";
  ofstream fo_key(fn.c_str(), ios::binary);

  // determine the number of valid and unique entries
  uint64_t prec = uint64_t(-1);
  uint64_t eof = 0;
  size_t valid;   // the number of entries different than -1

  size_t i;
  for (i = 0; i < data.size() && data[i].position != uint32_t(-1); i++) {
    if (data[i].hash != prec) {
      prec = data[i].hash;
      eof ++;
    }
  }
  valid = i;
  cerr << "Found " << eof << " valid keys and " << valid << " valid positions out of " << data.size() << " total\n\n";

//  fo_key.write((char *) &index.filter_cutoff, sizeof(int));
//  fo_key.write((char *) &index.bits, sizeof(int));
  fo_key.write((char *) &eof, 8);   // The number of entries is required to be a 64-bit value

  // write out keys
  try {
    cerr << "Fast writing uint64 keys (" << eof << ")\n";
    size_t elements = eof*3+3;
    uint32_t *buf = new uint32_t[elements];
    prec = uint64_t(-1);     // the previous value
    size_t i_buf;
    uint64_t *point;

    for (i = 0, i_buf = 0; i < valid && i_buf < (elements-3); i++) {
      if (data[i].hash != prec) {
        //this is what we have to change
        point = reinterpret_cast<uint64_t*>(buf+i_buf);
        *point = data[i].hash;
        i_buf += 2;
        buf[i_buf++] = i;
        prec = data[i].hash;
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
      if (data[i].hash != prec) {
        *point = data[i].hash;
        buf[2] = i;
        fo_key.write((char *) buf, 12);
        prec = data[i].hash;
      }
    }
    ////////// add fake last element //////////
    *point = 0;
    buf[2] = valid;
    fo_key.write((char *) buf, 12);
    ///////////////////////////////////////////
  }
  fo_key.close();
  cerr << "Key generation complete!\n\n";

  // now, write positions
  fn = dir  + "/pos_uint32";

  ofstream fo_pos(fn.c_str(), ios::binary);

  eof = (uint64_t) valid;
  fo_pos.write((char *) &eof, 8);

  try {
    cerr << "Fast writing posv (" << eof << ")\n";
    uint32_t *buf = new uint32_t[eof];
    for (i = 0; i < eof; i++) {
      buf[i] = data[i].position;
    }
    fo_pos.write((char *) buf, eof * sizeof(uint32_t));
    delete[] buf;
  } catch (std::bad_alloc& e) {
    cerr << "Fall back to slow writing posv due to low mem.\n";
    for (i = 0; i < eof; i++) {
      fo_pos.write((char *) &data[i].position, 4);
    }
  }
  fo_pos.close();
  cerr << "Position generation complete!\n\n";

  return EXIT_SUCCESS;
}


int main(int ac, char **av) {
  int opn = 1;

  if (av[opn] != nullptr && std::string(av[opn]) == "--strobe-mode") {
      try {
        return run_strobealign(ac, av);
      } catch (BadParameter& e) {
        logger.error() << "A parameter is invalid: " << e.what() << std::endl;
      } catch (const std::runtime_error& e) {
        logger.error() << "strobealign: " << e.what() << std::endl;
      }
      return EXIT_FAILURE;
  } else {
    if (ac < 2) {
      cerr << "index [options] <ref.fa>\n";
      cerr << "options:\n";
      cerr << "\t-l INT length of seed [32]\n";
      cerr << "\t-h INT value of hash MOD [2^29-1]\n";
      cerr << "\t   Special string values = 2^29-1, 2^32, prime, lprime\n";
      cerr << "\t-x INT size of xxhash [0]\n";
      cerr << "\t   Values = 0 (xxh not used), 32, 64\n";
      cerr << "\t-m enable minimizer\n";
      cerr << "\t-k minimizer: k, kmer size \n";
      cerr << "\t-w minimizer: w, window size \n";
      cerr << "\t-s bisulfite sequencing read alignment mode \n";
      return 0;
    }

    unsigned kmer_temp = 0, mm_k_tmp = 0, mm_w_tmp = 0;

    for (int it = 1; it < ac; it++) {
      if (strcmp(av[it], "-l") == 0)
        kmer_temp = atoi(av[it + 1]);
      else if (strcmp(av[it], "-m") == 0)
        enable_idx_minimizer = true;
      else if (strcmp(av[it], "-k") == 0)
        mm_k_tmp = atoi(av[it + 1]);
      else if (strcmp(av[it], "-w") == 0)
        mm_w_tmp = atoi(av[it + 1]);
      else if (strcmp(av[it], "-s") == 0)
        enable_bs = true;
      else if (strcmp(av[it], "-h") == 0)
        mod = atoll(av[it + 1]);
      else if (strcmp(av[it], "-x") == 0) {
        xxh_type = atoi(av[it + 1]);
        if (xxh_type!=0 && xxh_type!=32 && xxh_type!=64) {
          cerr << "Unknown value for xxhash. \nSupported values: 0 (xxh not used), 32, 64.\n";
          exit(1);
        }
      }
    }
    string fn = av[ac - 1]; //input ref file name
    bind_xxhash(xxh_type, xxh);

    if (enable_idx_minimizer) {
      int n_threads = 3;
      mm_idxopt_t ipt;
      mm_idxopt_init(&ipt);
      if (mm_k_tmp)
        ipt.k = mm_k_tmp;
      if (mm_w_tmp)
        ipt.w = mm_w_tmp;

      fn += ".hash"; // output hash: xxx.hash

      mm_idx_reader_t *idx_rdr = mm_idx_reader_open(av[ac - 1], &ipt, fn.c_str());
      mm_idx_reader_read(idx_rdr, n_threads, true);
    } else if (enable_bs) {
      kmer = kmer_temp ? kmer_temp: 32;
      cerr << "Using kmer length " << kmer << " and step size " << step << endl;

      cerr << "==== convert reference C to T ===="  << endl;
      Index ic;
      if (!ic.load_ref(fn.c_str(), 'c'))
        return 0;
      if (!ic.make_index(fn.c_str(), 1))
        return 0;

      cerr << "==== convert reference G to A ===="  << endl;
      Index ig;
      if (!ig.load_ref(fn.c_str(), 'g'))
        return 0;
      if (!ig.make_index(fn.c_str(), 2))
        return 0;

      cerr << "hash1 and hash2 have been generated"  << endl;
    } else {
      kmer = kmer_temp ? kmer_temp: 32;

      cerr << "Using kmer length " << kmer << " and step size " << step << endl;

      Index i;
      if (!i.load_ref(fn.c_str(), ' '))
        return 0;
      if (!i.make_index(fn.c_str(), 0))
        return 0;
    }

  }

  return 0;
}

