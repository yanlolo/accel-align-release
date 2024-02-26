#pragma once

#include "const.h"
#include "minimap.h"
#include "strobealign/strobe-index.hpp"
#include "strobealign/indexparameters.hpp"
#include "strobealign/aln.hpp"

//seedtype
enum class SType {
  Hash,
  Strobemer,
  Minimizer,
  LearntIndex // TODO
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
  int as;
  std::vector<Interval> matched_intervals;  // list of start pos of matched seeds in read that indicate to this region

  bool operator()(Region &X, Region &Y) {
    if (X.rs == Y.rs)
      return X.qs < Y.qs;
    else
      return X.rs < Y.rs;
  }

  // if last_q_pos is in range of match_interval, merge them togeter; k is the length of the interval (kmer_len)
  void merge_interval(SType g_stype, uint32_t last_q_pos, int32_t k) {
    if (g_stype==SType::Minimizer)
      last_q_pos = (last_q_pos >> 1) - k + 1; //q_pos format: pos << 1 | z

    if (!matched_intervals.size())
      matched_intervals.push_back(Interval{last_q_pos, last_q_pos + k});

    Interval &interval = matched_intervals.back();
    if (last_q_pos >= interval.s && last_q_pos <= interval.e) {
      interval.e = last_q_pos + k;
      return;
    }

    //no overlap, new interval
    matched_intervals.push_back(Interval{last_q_pos, last_q_pos + k});
    return;
  }

  void extend_interval(const char*ref, char*Q, int rlen) {
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
  void load_index(const char *F);
  void load_reference(const char *F);

  std::string ref;
  std::vector<std::string> name;
  std::vector<uint32_t> offset;
  uint32_t *keyv, *posv;
  uint32_t nposv, nkeyv;
  mm_idx_t *mi;
  SType g_stype;
  bool load_accalign_index;
  char mode; // 'c' c-> t; 'g' g->a; ' ' original

  Reference(const char *F, SType g_stype, char mode, bool load_accalign_index);

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