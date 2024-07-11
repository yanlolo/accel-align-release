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


typedef struct mm_idx_bucket_s {
  mm128_v a;   // (minimizer, position) array
  int32_t n;   // size of the _p_ array
  uint64_t *p; // position array for minimizers appearing >1 times
  void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

using namespace std;
uint64_t mod = MOD_29;    // default value is 2^29 - 1
uint32_t mod_tmp;
const unsigned step = 1;
uint32_t xxh_type = 0;
XXHash xxh;
unsigned kmer;
bool enable_idx_minimizer = false, enable_bs = false; //short for bisulfite reads
struct Data {
  uint32_t key, pos;
  bool is_fwd; // fwd: 1, rev: 0
  Data() : key(-1), pos(-1), is_fwd(true) {}
  Data(uint32_t k, uint32_t p, bool dir) : key(k), pos(p), is_fwd(dir) {}
  bool operator()(const Data &X, const Data &Y) const {
    if (X.key != Y.key ){
      return X.key < Y.key;
    } else if (X.pos != Y.pos){
      return X.pos < Y.pos;
    } else {
      return X.is_fwd > Y.is_fwd;  //fwd before rev
    }
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
  uint64_t h0 = 0, h1 = 0;
  bool hasn = false;
  for (unsigned j = 0; j < kmer; j++) {
    if (ref[i + j] == 4) {
      hasn = true;
    }
    h0 = (h0 << 2) + ref[i + j];
    h1 = (h1 << 2) + (3ULL^ref[i + kmer - 1 - j]);
  }
  uint64_t h = h0 < h1 ? h0 : h1;
  bool is_fwd = h0 < h1 ? 1 : 0;
  if (!hasn) {
    data[i / step].key = uint32_t(xxh(&h) % mod);
    data[i / step].pos = i;
    data[i / step].is_fwd = is_fwd;
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
  size_t vsz = step == 1 ? limit: ref.size() / step + 1;
  vector<Data> data(vsz, Data());
  tbb::parallel_for(tbb::blocked_range<size_t>(0, limit), Tbb_cal_key(data, this));

  cerr << "hashing :limit = " << limit << ", vsz = " << vsz << endl;
  cerr << "using MOD = " << mod << " and XXH = " << xxh_type << endl;
  cerr << "hash\t" << data.size() << endl;

  //XXX: Parallel sort uses lots of memory. Need to fix this. In general, we use 8 bytes per item. Its a waste.
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

//  ofstream fo(fn.c_str(), ios::binary);
  std::ofstream outputFile(fn.c_str());

  // determine the number of valid entries based on first junk entry
  auto joff = std::lower_bound(data.begin(), data.end(), Data(-1, -1, 1), Data());
  size_t eof = joff - data.begin();
  cerr << "Found " << eof << " valid entries out of " << data.size() << " total\n";

  stringstream ss;
  ss << to_string(mod) << ',' << to_string(xxh_type) << ',' << to_string(eof) << endl;

  // write out keys
  cerr << "Number of total keys (" << eof << ")\n";
  size_t last_key = 0, offset, last_offset = 0;
  cerr << "Fall back to slow writing keyv (low mem)\n";
  for (size_t i = 0; i < eof;) {
    assert (data[i].pos != (uint32_t) -1);
    size_t h = data[i].key, n;
    offset = i;
    for (size_t j = last_key; j <= h; j++) {
      size_t tmp = offset - last_offset;
      ss << to_string(j) << "," <<to_string(tmp) << endl;
      last_offset = offset;
    }
    last_key = h + 1;
    for (n = i + 1; n < eof && data[n].key == h; n++);
    i = n;
  }

  offset = eof;
  for (uint64_t j = (uint64_t) last_key; j <= mod; j++) {
    size_t tmp = offset - last_offset;
    ss << to_string(j) << "," << to_string(tmp) << endl;
  }
  outputFile << ss.str();
  outputFile.flush();
  outputFile.close();

//  }

  cerr << "Indexing complete\n";
//  fo.close();
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

  vector<RefRandstrobe> strobe = index.randstrobes;
  std::sort(strobe.begin(), strobe.end());

  // write out keys
  stringstream ss;
  string output_file = opt.ref_filename + ".hash";
  std::ofstream outputFile(output_file);
  cerr << "Number of total strobe (" << strobe.size() << ")\n";
  cerr << "The statistic is output at " << output_file << endl;
  size_t last_key = 0, offset, last_offset = 0;
  for (size_t i = 0; i < strobe.size();) {
    size_t h = strobe[i].hash, n;
    for (n = i + 1; n < strobe.size() && strobe[n].hash == h; n++);
    cerr << to_string(n - i) << endl;
    i = n;
  }

  outputFile << ss.str();
  outputFile.flush();
  outputFile.close();
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
      idx_rdr->opt.bucket_bits = 30; // 2 * k
      mm_idx_t *mi = mm_idx_reader_read(idx_rdr, n_threads, false);
      int max_bin = (1<<mi->b) - 1;
      cerr << "writing\n";
      std::ofstream outputFile(fn.c_str());
      stringstream ss;
      size_t nb_keys = 0, nb_pos = 0;
      for(int i = 0; i < max_bin; ++i){
        mm_idx_bucket_t *b = &mi->B[i];
        size_t nb_pos_per_key = b->a.n;
        if (nb_pos_per_key){
          ++nb_keys;
          nb_pos += nb_pos_per_key;
          ss << to_string(nb_pos_per_key) << endl;
        }
      }
      outputFile << ss.str();
      outputFile.flush();
      outputFile.close();
      cerr << "writing finished \n";
      cerr << "Number of total keys (" << nb_keys << ")\n";
      cerr << "Number of total positions (" << nb_pos << ")\n";


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

