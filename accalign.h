#pragma once

#include "mmpriv.h"
#include "strobealign/strobe-index.hpp"
#include "strobealign/aln.hpp"

class AccAlign {
 private:
//  std::string &ref;
//  std::vector<std::string> &name;
//  std::vector<uint32_t> &offset;
  Embedding *embedding;


  float input_io_time, parse_time;
  float seeding_time, hit_count_time, vpair_build_time;
  float swap_time;
  float sw_time, sam_time, sam_pre_time, sam_out_time;
  long vpair_sort_count;

  // output fields
  std::string sam_name;
  std::ofstream sam_stream;
  std::mutex sam_mutex;

  void cpu_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *inputQ,
                   tbb::concurrent_bounded_queue<ReadCnt> *outputQ);
  void output_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *outputQ,
                      tbb::concurrent_bounded_queue<ReadPair> *dataQ);
  void align_wrapper(int tid, int soff, int eoff,
                     Read *ptlread, Read *ptlread2,
                     tbb::concurrent_bounded_queue<ReadPair> *dataQ);
  void embed_wrapper(Read &R, bool ispe, vector<Region> &fregion,
                     vector<Region> &rregion, unsigned &fbest, unsigned &fnext, unsigned &rbest,
                     unsigned &rnext, int &best_threshold, int &next_threshold, int ref_id);
  void pghole_wrapper(Read &R, vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
                      unsigned &fbest, unsigned &rbest, int ref_id);
  void pigeonhole_query_topcov(char *Q, size_t rlen, vector<Region> &candidate_regions, char S, int err_threshold,
                               unsigned kmer_step, unsigned max_occ, unsigned &best, unsigned ori_slide, int ref_id);
  void pghole_wrapper_mates(Read &R, vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
                            unsigned &fbest, unsigned &rbest, unsigned ori_slide, unsigned kmer_step, unsigned max_occ,
                            bool &high_freq, int ref_id);
  void pigeonhole_query(char *Q, size_t rlen, vector<Region> &candidate_regions, char S,
                        unsigned &best, unsigned ori_slide, int err_threshold, unsigned kmer_step,
                        unsigned max_occ, bool &high_freq, int ref_id);
  void pghole_wrapper_pair(Read &mate1, Read &mate2,
                           vector<Region> &region_f1, vector<Region> &region_r1,
                           vector<Region> &region_f2, vector<Region> &region_r2,
                           unsigned &best_f1, unsigned &best_r1, unsigned &best_f2, unsigned &best_r2,
                           unsigned &next_f1, unsigned &next_r1, unsigned &next_f2, unsigned &next_r2,
                           bool *&flag_f1, bool *&flag_r1, bool *&flag_f2, bool *&flag_r2,
                           bool &has_f1r2, bool &has_r1f2, int ref_id);
  void pigeonhole_query_sort(char *Q, size_t rlen, vector<Region> &candidate_regions, char S,
                             unsigned err_threshold, unsigned kmer_step, unsigned max_occ,
                             unsigned &best, unsigned ori_slide, int ref_id);
  bool pairdis_filter(vector<Region> &in_regions1, vector<Region> &in_regions2,
                      bool flag1[], bool flag2[],
                      unsigned &best1, unsigned &next1, unsigned &best2, unsigned &next2);
  void mark_for_extension(Read &read, char S, Region &cregion, int ref_id);
  void save_region(Read &R, size_t rlen, Region &region,
                   Alignment &a);
  void score_region(Read &r, char *qseq, Region &region,
                    Alignment &a);
  void sam_header(void);
  void extend_pair(Read &mate1, Read &mate2,
                   vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                   bool flag_f1[], bool flag_r2[], unsigned &best_f1, unsigned &best_r2,
                   int &best_threshold, int &next_threshold, char strand);
  void embed_wrapper_pair(Read &R1, Read &R2,
                          vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                          bool flag_f1[], bool flag_r2[], unsigned &best_f1, unsigned &best_r2,
                          int &best_threshold, int &next_threshold, char strand, int ref_id);
  void mm(char *Q, size_t rlen, int err_threshold, vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
                    unsigned &fbest, unsigned &rbest);
  inline uint32_t get_global_pos(uint64_t cr, int ref_id);
  inline uint64_t normalize_pos(uint64_t cr, uint32_t q_pos, int k, int rlen);
  void find_candidate_positions_using_strobealign(char *seq, vector<Region> &fcandidate_regions, bool direction, int ref_id);

 public:
  Reference **refs;
  StrobemerIndex *index_reference;
  IndexParameters *index_parameters_reference;
  MappingParameters map_params;
//  uint32_t *keyv, *posv;
//  mm_idx_t *mi;

  std::string& get_ref(int ref_id){
    return refs[ref_id]->ref;
  }

  std::vector<std::string>& get_name(int ref_id){
    return refs[ref_id]->name;
  };

  std::vector<uint32_t>& get_offset(int ref_id){
      return refs[ref_id]->offset;
  };

  uint32_t* get_keyv(int ref_id){
    return refs[ref_id]->keyv;
  };

  uint32_t* get_posv(int ref_id){
    return refs[ref_id]->posv;
  };

  mm_idx_t* get_mi(int ref_id){
    return refs[ref_id]->mi;
  };


  void open_output(std::string &out_file);
  void close_output();
  bool fastq(const char *F1, const char *F2, bool enable_gpu);
  void print_stats();
  void map_read(Read &R, int ref_id);
  void align_read(Read &R);
  void out_sam(string *sam);
  void snprintf_sam(Read &R, string *s);
  void snprintf_pair_sam(Read &R, string *s, Read &R2, string *s2);
  void map_paired_read(Read &mate1, Read &mate2, int ref_id);
  void map_read_wrapper(Read &R);
  void map_paired_read_wrapper(Read &mate1, Read &mate2);
  void wfa_align_read(Read &R);
  void rectify_start_pos(char *strand, Region &region, unsigned rlen, int ref_id);
  void print_paired_sam(Read &R, Read &R2);
  void print_sam(Read &R);
  bool tbb_fastq(const char *F1, const char *F2);
  int get_mapq(int best, int secbest);
  int get_tid(Read &R);
//  void merge_interval(Region &r, uint32_t last_q_pos, int32_t k);
  void extend_interval(Region &r, char*Q, int rlen, int ref_id);
  void collect_seed_hits_priorityqueue(int n_m0, int64_t n_a, size_t rlen, int err_threshold, mm_seed_t* m, vector<Region> &candidate_regions,
                                                 vector<Region> &rcandidate_regions, unsigned &best, unsigned &rbest, int ref_id);
  void fetch_candidates(mm128_v &mv, int32_t mid_occ, size_t rlen, int err_threshold,
                        vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
                        unsigned &fbest, unsigned &rbest, int ref_id);
  AccAlign(Reference **r, StrobemerIndex *_index_reference, IndexParameters *_index_parameters_reference,MappingParameters _map_params);
  ~AccAlign();
};


