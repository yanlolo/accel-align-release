#pragma once

#include "const.h"
#include "minimap.h"
#include "rmi.h"
#include "hash.hpp"
#include "../strobealign/strobe-index.hpp"
#include "../strobealign/indexparameters.hpp"
#include "../strobealign/aln.hpp"
#include "../strobealign/cmdline.hpp"

//seedtype: how to build the seed
enum class SeedType {
  __NONE__,
  Kmer,
  Strobemer,
  Minimizer
};

//Index type: how to look up the seed
enum IndexType {
  __NONE__,
  RMI_IDX,
  BINARY_IDX,
  HASH_IDX
};

struct Alignment {
  std::string cigar_string;
  int ref_begin;
  int mismatches;
};

struct Interval {
  uint32_t s, e;  // start, end position of the reference match to the whole read
};

struct Region {
  uint32_t rs;  // start position of the reference match to the whole read
  uint32_t qs, qe;  // start, end position of matched seed in the query (read)
  uint16_t cov;
  uint16_t embed_dist;
  bool is_fwd;
  int as;
  std::vector<Interval> matched_intervals;  // list of start pos of matched seeds in read that indicate to this region

  bool operator()(Region &X, Region &Y) {
    if (X.rs == Y.rs)
      return X.qs < Y.qs;
    else
      return X.rs < Y.rs;
  }

  /*
   * if pos is in range of match_interval, merge them togeter; kmer_len is the length of the interval
   */
  void add_match_interval(SeedType g_stype, uint32_t pos, int32_t kmer_len) {
    if (g_stype==SeedType::Minimizer)
      pos = (pos >> 1) - kmer_len + 1; //q_pos format: pos << 1 | z

    if (!matched_intervals.size()){  // no interval yet
      matched_intervals.push_back(Interval{pos, pos + kmer_len});
      return;
    }

    Interval &interval = matched_intervals.back();
    if (pos >= interval.s && pos <= interval.e) { //within range
      interval.e = pos + kmer_len;
    } else {  //out of range
      matched_intervals.push_back(Interval{pos, pos + kmer_len});
    }

    return;
  }

  void extend_interval(const char* ref, char* Q, int rlen) {
    for (size_t i = 0; i < matched_intervals.size(); ++i) {
      Interval &interval = matched_intervals[i];

      size_t left = i == 0 ? 0 : matched_intervals[i - 1].e;
      for(size_t j = 0; j + left < interval.s; ++j){
        if (Q[interval.s - j] != ref[rs + interval.s - j]){
          interval.s = interval.s - j + 1;
          break;
        }
      }

      size_t right = i == matched_intervals.size() - 1 ? rlen : matched_intervals[i + 1].s;
      for(size_t j = 0; interval.e + j < right; ++j){
        if (Q[interval.e + j] != ref[rs + interval.e + j]){
          interval.e = interval.e + j - 1;
          break;
        }
      }
    }
  }

};

struct Read {
  char name[MAX_LEN], qua[MAX_LEN], seq[MAX_LEN], fwd[MAX_LEN], rev[MAX_LEN], rev_str[MAX_LEN], cigar[MAX_LEN];
  int tid, as, nm, best, secBest, best_optional, secBest_optional, rlen, ref_id;
  uint32_t pos;
  short mapq, kmer_step; //kmer_step used that find the seed
  char strand;
  Region best_region;

  char strand_optional;
  Region best_region_optional;

  bool force_align = false;

  friend gzFile &operator>>(gzFile &in, Read &r);

  Read() {
    best = INT_MAX;
    secBest = INT_MAX;
  }

};

class Reference {
 public:
  void load_index32(const char *F);
  void load_index64(const char *F);
  void load_index_classic(const char *F);
  std::function<void(const char*)> load_index;
  void index_rmi_lookup32(uint64_t key, size_t* b, size_t* e);
  void index_rmi_lookup64(uint64_t key, size_t* b, size_t* e);
  void index_bin_lookup32(uint64_t key, size_t* b, size_t* e);
  void index_bin_lookup64(uint64_t key, size_t* b, size_t* e);
  void index_lookup_classic(uint64_t key, size_t* b, size_t* e);
  std::function<void(uint64_t,size_t*,size_t*)> index_lookup;
  uint32_t get_keyv_val32(uint64_t idx);
  uint32_t get_keyv_val64(uint64_t idx);
  std::function<uint32_t(uint64_t)> get_keyv_val;

  void load_reference(const char *F);

  std::string ref;
  std::vector<std::string> name;
  std::vector<uint32_t> offset;
  uint32_t *keyv, *posv;
  bool *is_fwdv;
  uint64_t nposv, nkeyv, nkeyv_true;
  mm_idx_t *mi;
  RMI rmi;  // this is fine
  unsigned kmer_len;
  SeedType g_stype;
  IndexType index_type;
  bool load_accalign_index;
  CommandLineOptions opt;
  StrobemerIndex *strobe_index;
  IndexParameters *index_parameters_reference;

  char mode; // 'c' c-> t; 'g' g->a; ' ' original
  // for classic index
  uint64_t mod;
  uint32_t xxh_type;
  XXHash xxh;
//  Reference(const char *F, SeedType g_stype, char mode, bool load_accalign_index);
  Reference(const char *F, unsigned _kmer_len, SeedType g_stype, IndexType _index_type, char _mode,  IndexParameters *index_parameters_reference);
  ~Reference();
};

typedef std::tuple<Read *, Read *, int> ReadCnt;

typedef std::tuple<Read *, Read *> ReadPair;

typedef struct {
  uint32_t capacity;
  uint32_t n_cigar;
  int32_t dp_score;
  uint32_t cigar[];
} Extension;

// add RMI class
// function typedefs
typedef bool (*load_type)(char const*);
typedef uint64_t (*lookup_type)(uint64_t, size_t*);
typedef void (*cleanup_type)();