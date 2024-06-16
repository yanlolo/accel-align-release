class Tbb_aligner_paired {
  Read *all_reads;
  Read *all_reads2;
  string *sams;
  AccAlign *acc_obj;

 public:
  Tbb_aligner_paired(Read *_all_reads, Read *_all_reads2, string *_sams, AccAlign *_acc_obj) :
      all_reads(_all_reads), all_reads2(_all_reads2), sams(_sams), acc_obj(_acc_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      if ((all_reads + i)->strand != '*') {
        if (enable_wfa_extension)
          acc_obj->wfa_align_read(*(all_reads + i));
        else
          acc_obj->align_read(*(all_reads + i));
      }

      if ((all_reads2 + i)->strand != '*') {
        if (enable_wfa_extension)
          acc_obj->wfa_align_read(*(all_reads2 + i));
        else
          acc_obj->align_read(*(all_reads2 + i));
      }

      acc_obj->snprintf_pair_sam(*(all_reads + i), sams + 2 * i, *(all_reads2 + i), sams + 2 * i + 1);
    }
  }
};
//// This is not updated --> need to check if want to use
//void AccAlign::snprintf_pair_sam(Read &R, string *s, Read &R2, string *s2) {
//  auto start = std::chrono::system_clock::now();
//
//  // 60 is the approximate length for all int
//  int size;
//  if (!enable_extension) {
//    size = 60;
//  } else {
//    size = 60 + strlen(R.seq); //assume the length of cigar will not longer than the read
//  }
//  char strand1 = R.strand;
//  char strand2 = R2.strand;
//
//  //mate 1
//  string rname = R.name;
//  string nn = rname.substr(0, rname.find_last_of("/"));
//
//  uint16_t flag = 0x1;
//  if (strand1 == '*')
//    flag |= 0x4;
//  if (strand2 == '*')
//    flag |= 0x8;
//  if (!(flag & 0x4) && !(flag & 0x8))
//    flag |= 0x2;
//  if (strand1 == '-')
//    flag |= 0x10;
//  if (strand2 == '-')
//    flag |= 0x20;
//  flag |= 0x40;
//
//  int isize = 0;
//  if (R.strand != '*' && R2.strand != '*') {
//    if (R.pos > R2.pos)
//      isize = R2.pos - R.pos - strlen(R.seq);
//    else
//      isize = R2.pos - R.pos + strlen(R2.seq);
//  }
//
//  string format = "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\n";
//  if (R.strand == '+' && !unfill(R2)) {
//    size += strlen(R.name) + get_name(R.ref_id)[R.tid].length() + get_name(R.ref_id)[R2.tid].length() + 2 * strlen(R.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(), nn.c_str(), flag, get_name(R.ref_id)[R.tid].c_str(), R.pos,
//             (int) R.mapq, R.cigar, get_name(R.ref_id)[R.tid] == get_name(R.ref_id)[R2.tid] ? "=" : get_name(R.ref_id)[R2.tid].c_str(),
//             R2.pos, isize, R.seq, R.qua, R.nm, R.as);
//    *s = buf;
//  } else if (R.strand == '-' && !unfill(R2)) {
//    size += strlen(R.name) + get_name(R.ref_id)[R.tid].length() + get_name(R.ref_id)[R2.tid].length() + 2 * strlen(R.seq);
//    char buf[size];
//    std::reverse(R.qua, R.qua + strlen(R.qua));
//    snprintf(buf, size, format.c_str(), nn.c_str(), flag, get_name(R.ref_id)[R.tid].c_str(), R.pos,
//             (int) R.mapq, R.cigar, get_name(R.ref_id)[R.tid] == get_name(R.ref_id)[R2.tid] ? "=" : get_name(R.ref_id)[R2.tid].c_str(),
//             R2.pos, isize, R.rev_str, R.qua, R.nm, R.as);
//    *s = buf;
//  } else if (R.strand == '+' && unfill(R2)) {
//    size += strlen(R.name) + get_name(R.ref_id)[R.tid].length() + 2 * strlen(R.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(),
//             nn.c_str(), flag, get_name(R.ref_id)[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar, "*", 0,
//             isize, R.seq, R.qua, R.nm, R.as);
//    *s = buf;
//  } else if (R.strand == '-' && unfill(R2)) {
//    size += strlen(R.name) + get_name(R.ref_id)[R2.tid].length() + 2 * strlen(R.seq);
//    char buf[size];
//    std::reverse(R.qua, R.qua + strlen(R.qua));
//    snprintf(buf, size, format.c_str(),
//             nn.c_str(), flag, get_name(R.ref_id)[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar, "*", 0,
//             isize, R.rev_str, R.qua, R.nm, R.as);
//    *s = buf;
//  } else if (unfill(R) && R2.strand != '*') {
//    size += strlen(R.name) + get_name(R.ref_id)[R2.tid].length() + 2 * strlen(R.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(),
//             nn.c_str(), flag, "*", 0, 0, "*", get_name(R.ref_id)[R2.tid].c_str(), R2.pos,
//             isize, R.seq, R.qua, R.nm, R.as);
//    *s = buf;
//  } else if (unfill(R) && unfill(R)) {
//    size += strlen(R.name) + 2 * strlen(R.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(),
//             nn.c_str(), flag, "*", 0, 0, "*", "*", 0,
//             isize, R.seq, R.qua, R.nm, R.as);
//    *s = buf;
//  }
//
//  //for mate2
//  string rname2 = R2.name;
//  string nn2 = rname2.substr(0, rname2.find_last_of("/"));
//
//  flag = 0x1;
//  if (strand2 == '*')
//    flag |= 0x4;
//  if (strand1 == '*')
//    flag |= 0x8;
//  if (!(flag & 0x4) && !(flag & 0x8))
//    flag |= 0x2;
//  if (strand2 == '-')
//    flag |= 0x10;
//  if (strand1 == '-')
//    flag |= 0x20;
//  flag |= 0x80;
//
//  if (R2.strand == '+' && R.strand != '*') {
//    size += strlen(R2.name) + get_name(R2.ref_id)[R2.tid].length() + get_name(R2.ref_id)[R.tid].length() + 2 * strlen(R2.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(), nn2.c_str(), flag, get_name(R2.ref_id)[R2.tid].c_str(), R2.pos,
//             (int) R2.mapq, R2.cigar, get_name(R2.ref_id)[R.tid] == get_name(R2.ref_id)[R2.tid] ? "=" : get_name(R2.ref_id)[R.tid].c_str(),
//             R.pos, -isize, R2.seq, R2.qua, R2.nm, R2.as);
//    *s2 = buf;
//  } else if (R2.strand == '-' && R.strand != '*') {
//    size += strlen(R2.name) + get_name(R2.ref_id)[R2.tid].length() + get_name(R2.ref_id)[R.tid].length() + 2 * strlen(R2.seq);
//    char buf[size];
//    std::reverse(R2.qua, R2.qua + strlen(R2.qua));
//    snprintf(buf, size, format.c_str(), nn2.c_str(), flag, get_name(R2.ref_id)[R2.tid].c_str(), R2.pos,
//             (int) R2.mapq, R2.cigar, get_name(R2.ref_id)[R.tid] == get_name(R2.ref_id)[R2.tid] ? "=" : get_name(R2.ref_id)[R.tid].c_str(),
//             R.pos, -isize, R2.rev_str, R2.qua, R2.nm, R2.as);
//    *s2 = buf;
//  } else if (R2.strand == '+' && R.strand == '*') {
//    size += strlen(R2.name) + get_name(R2.ref_id)[R2.tid].length() + 2 * strlen(R2.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(),
//             nn2.c_str(), flag, get_name(R2.ref_id)[R2.tid].c_str(), R2.pos, (int) R2.mapq, R2.cigar, "*", 0,
//             -isize, R2.seq, R2.qua, R2.nm, R2.as);
//    *s2 = buf;
//  } else if (R2.strand == '-' && R.strand == '*') {
//    size += strlen(R2.name) + get_name(R2.ref_id)[R.tid].length() + 2 * strlen(R2.seq);
//    char buf[size];
//    std::reverse(R2.qua, R2.qua + strlen(R2.qua));
//    snprintf(buf, size, format.c_str(),
//             nn2.c_str(), flag, get_name(R2.ref_id)[R2.tid].c_str(), R2.pos, (int) R2.mapq, R2.cigar, "*", 0,
//             -isize, R2.rev_str, R2.qua, R2.nm, R2.as);
//    *s2 = buf;
//  } else if (R2.strand == '*' && R.strand != '*') {
//    size += strlen(R2.name) + get_name(R2.ref_id)[R.tid].length() + 2 * strlen(R2.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(),
//             nn2.c_str(), flag, "*", 0, 0, "*", get_name(R2.ref_id)[R.tid].c_str(), R.pos,
//             -isize, R2.seq, R2.qua, R2.nm, R2.as);
//    *s2 = buf;
//  } else if (R2.strand == '*' && R.strand == '*') {
//    size += strlen(R2.name) + 2 * strlen(R2.seq);
//    char buf[size];
//    snprintf(buf, size, format.c_str(),
//             nn2.c_str(), flag, "*", 0, 0, "*", "*", 0,
//             -isize, R2.seq, R2.qua, R2.nm, R2.as);
//    *s2 = buf;
//  }
//
//  auto end = std::chrono::system_clock::now();
//  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//  sam_pre_time += elapsed.count();
//}

void AccAlign::pigeonhole_query_topcov(char *Q,
                                       size_t rlen,
                                       vector<Region> &candidate_regions,
                                       char S,
                                       int err_threshold,
                                       unsigned kmer_step,
                                       unsigned max_occ,
                                       unsigned &best,
                                       unsigned ori_slide,
                                       int ref_id) {
  int max_cov = 0;
  unsigned nkmers = (rlen - ori_slide - kmer_len) / kmer_step + 1;
  unsigned nkmers_slct = (rlen - ori_slide - kmer_len) / kmer_len + 1;
  size_t ntotal_hits = 0;
  size_t b[nkmers], e[nkmers];
  unsigned kmer_idx = 0;
  unsigned ori_slide_bk = ori_slide;
  unsigned nseed_freq = 0;
  bool high_freq = false;
  vector<size_t> cnt;
  cnt.reserve(nkmers);

  // Take non-overlapping seeds and find all hits
  auto start = std::chrono::system_clock::now();
  for (size_t i = ori_slide; i + kmer_len <= rlen; i += kmer_step) {
    uint64_t k = 0;
    for (size_t j = i; j < i + kmer_len; j++)
//      k = (k << 2) + *(Q + j);
  }
  // lookup to get the position of the hash
  get_lookup(ref_id, k, b+kmer_idx, e+kmer_idx);
  cnt.push_back(e[kmer_idx] - b[kmer_idx]);

  if (e[kmer_idx] - b[kmer_idx] >= max_occ)
    nseed_freq++;
  kmer_idx++;
}
assert(kmer_idx == nkmers);
auto end = std::chrono::system_clock::now();
auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
keyvTime += elapsed.count();

vector<int> slct_seed_idx = findBestIndices(cnt, nkmers_slct);

if (nseed_freq > nkmers / 2)
high_freq = true;

for (size_t i = 0; i < nkmers; i++) {
if (is_selected_seed(high_freq, max_occ, cnt[i], i, slct_seed_idx))
ntotal_hits += (e[i] - b[i]);
}

// if we have no hits, we are done
if (!ntotal_hits)
return;

uint32_t top_pos[nkmers];
int rel_off[nkmers];
uint32_t MAX_POS = numeric_limits<uint32_t>::max();

start = std::chrono::system_clock::now();
// initialize top values with first values for each kmer.
for (unsigned i = 0; i < nkmers; i++) {
if (b[i] < e[i] && is_selected_seed(high_freq, max_occ, cnt[i], i, slct_seed_idx)) {
top_pos[i] = get_posv(ref_id)[b[i]];
rel_off[i] = i * kmer_step;
uint32_t shift_pos = rel_off[i] + ori_slide_bk;
//TODO: for each chrome, happen to < the start pos
if (top_pos[i] < shift_pos)
top_pos[i] = 0; // there is insertion before this kmer
else
top_pos[i] -= shift_pos;
} else {
top_pos[i] = MAX_POS;
}
}
end = std::chrono::system_clock::now();
elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
posvTime += elapsed.count();

size_t nprocessed = 0;
uint32_t last_pos = MAX_POS, last_qs = ori_slide_bk; //last query start pos
int last_cov = 0;

start = std::chrono::system_clock::now();

vector<Region> unique_regions;
unique_regions.reserve(ntotal_hits);
size_t idx = 0;
Region r;
r.matched_intervals.reserve(nkmers);
//  Region unique_regions[ntotal_hits];
//  size_t idx = 0;
//  Region *r = unique_regions + idx;

while (nprocessed < ntotal_hits) {
//find min
uint32_t *min_item = min_element(top_pos, top_pos + nkmers);
uint32_t min_pos = *min_item;
int min_kmer = min_item - top_pos;

if (is_selected_seed(high_freq, max_occ, cnt[min_kmer], min_kmer, slct_seed_idx)) {
// kick off prefetch for next round
__builtin_prefetch(get_posv(ref_id) + b[min_kmer] + 1);

// if previous min element was same as current one, increment coverage.
// otherwise, check if last min element's coverage was high enough to make it a candidate region

if (min_pos == last_pos) {
add_match_interval(r, last_qs, kmer_len);
//        r.matched_intervals.push_back(Interval{last_qs, last_qs + kmer_len});
last_cov++;
} else {
if (nprocessed != 0) {
r.cov = last_cov;
r.rs = last_pos;
add_match_interval(r, last_qs, kmer_len);
extend_interval( r, Q,  rlen,  ref_id);
//          r.matched_intervals.push_back(Interval{last_qs, last_qs + kmer_len});
r.qs = r.matched_intervals[0].s; //the first match seed, so left extension could be accurate
r.qe = r.matched_intervals[0].e;

if (last_cov > max_cov)
max_cov = last_cov;

assert(r.rs != MAX_POS && r.rs < MAX_POS);
unique_regions.push_back(move(r));
++idx;
//          r =  unique_regions + idx;
}

last_cov = 1;
}
last_qs = min_kmer * kmer_step + ori_slide_bk;
last_pos = min_pos;
}

// add next element
b[min_kmer]++;
uint32_t next_pos = b[min_kmer] < e[min_kmer] ? get_posv(ref_id)[b[min_kmer]] : MAX_POS;
if (next_pos != MAX_POS) {
uint32_t shift_pos = rel_off[min_kmer] + ori_slide_bk;
//TODO: for each chrome, happen to < the start pos
if (next_pos < shift_pos)
*min_item = 0; // there is insertion before this kmer
else
*min_item = next_pos - shift_pos;
} else
*min_item = MAX_POS;

++nprocessed;
}

// we will have the last few positions not processed. check here.
if (last_pos != MAX_POS) {
r.cov = last_cov;
r.rs = last_pos;
add_match_interval(r, last_qs, kmer_len);
extend_interval( r, Q,  rlen,  ref_id);
//    r.matched_intervals.push_back(Interval{last_qs, last_qs + kmer_len});
r.qs = r.matched_intervals[0].s; //the first match seed, so left extension could be accurate
r.qe = r.matched_intervals[0].e;

if (last_cov > max_cov)
max_cov = last_cov;

assert(r.rs != MAX_POS && r.rs < MAX_POS);
unique_regions.push_back(move(r));
++idx;
}

err_threshold = max(err_threshold, max_cov - 1);
assert(idx <= ntotal_hits);
for (size_t i = 0; i < idx; i++) {
Region &r = unique_regions[i];
if (r.cov >= err_threshold) {
if (r.cov == max_cov)
best = candidate_regions.size();

candidate_regions.push_back(move(r));
}
}

end = std::chrono::system_clock::now();
elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
hit_count_time += elapsed.count();
}

void AccAlign::pigeonhole_query_sort(char *Q,
                                     size_t rlen,
                                     vector<Region> &candidate_regions,
                                     char S,
                                     unsigned err_threshold,
                                     unsigned kmer_step,
                                     unsigned max_occ,
                                     unsigned &best,
                                     unsigned ori_slide,
                                     int ref_id) {
  unsigned max_cov = 0;
  unsigned nkmers = (rlen - ori_slide - kmer_len) / kmer_step + 1;
  size_t ntotal_hits = 0;
  size_t b[nkmers], e[nkmers];
  unsigned kmer_idx = 0;
  unsigned nseed_freq = 0;
  bool high_freq = false;

  // Take non-overlapping seeds and find all hits
  auto start = std::chrono::system_clock::now();
  for (size_t i = ori_slide; i + kmer_len <= rlen; i += kmer_step) {
    uint64_t k = 0;
    for (size_t j = i; j < i + kmer_len; j++)
      k = (k << 2) + *(Q + j);
    size_t hash = (k & mask) % MOD;
    b[kmer_idx] = get_keyv(ref_id)[hash];
    e[kmer_idx] = get_keyv(ref_id)[hash + 1];
    if (e[kmer_idx] - b[kmer_idx] >= max_occ)
      nseed_freq++;
//    if (e[kmer_idx] - b[kmer_idx] < max_occ) {
//      ntotal_hits += (e[kmer_idx] - b[kmer_idx]);
//    }
    kmer_idx++;
  }
  assert(kmer_idx == nkmers);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  keyvTime += elapsed.count();

  if (nseed_freq > nkmers / 2)
    high_freq = true;

  for (size_t i = 0; i < nkmers; i++) {
    if ((!high_freq && e[i] - b[i] < max_occ) || high_freq)
      ntotal_hits += (e[i] - b[i]);
  }

  // if we have no hits, we are done
  if (!ntotal_hits)
    return;

  start = std::chrono::system_clock::now();
  // initialize top values with first values for each kmer.
  uint32_t MAX_POS = numeric_limits<uint32_t>::max();
  vector<Region> regions;
  regions.reserve(ntotal_hits);
  for (unsigned i = 0; i < nkmers; i++) {
    if (b[i] < e[i] && ((!high_freq && e[i] - b[i] < max_occ) || high_freq)) {
//    if (b[i] < e[i] && e[i] - b[i] < max_occ) {
      for (uint32_t j = b[i]; j < e[i]; j++) {
        Region r;
        r.rs = get_posv(ref_id)[j];
        r.qs = i * kmer_step + ori_slide;
        r.rs -= min(r.rs, r.qs);
        regions.push_back(r);
        // rs can't be samller than 0, if insertion before this kmer, set rs to 0 instead of -1
      }
    }
  }
  assert(regions.size() == ntotal_hits);
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  posvTime += elapsed.count();

  start = std::chrono::system_clock::now();

  sort(regions.begin(), regions.end(), Region());

  size_t nprocessed = 0, last_cov = 0;
  uint32_t last_pos = MAX_POS;

  while (nprocessed < ntotal_hits) {

    if (regions[nprocessed].rs == last_pos) {
      last_cov++;
    } else {
      if (last_cov >= err_threshold) {
        Region r;
        r.cov = last_cov;
        r.rs = last_pos;
        for (unsigned i = nprocessed - last_cov; i < nprocessed; i++)
          r.matched_intervals.push_back(Interval{regions[i].qs, regions[i].qs + kmer_len});
        r.qs = r.matched_intervals[0].s; //the first match seed, so left extension could be accurate
        r.qe = r.qs + kmer_len;

        assert(r.rs < MAX_POS);

        if (last_cov >= max_cov) {
          max_cov = last_cov;
          best = candidate_regions.size();
        }
        candidate_regions.push_back(r);
      }
      last_cov = 1;
    }
    last_pos = regions[nprocessed].rs;

    ++nprocessed;
  }

  // we will have the last few positions not processed. check here.
  if (last_cov >= err_threshold && last_pos != MAX_POS) {
    Region r;
    r.cov = last_cov;
    r.rs = last_pos;
    for (unsigned i = nprocessed - last_cov; i < nprocessed; i++)
      r.matched_intervals.push_back(Interval{regions[i].qs, regions[i].qs + kmer_len});
    r.qs = r.matched_intervals[0].s; //the first match seed, so left extension could be accurate
    r.qe = r.qs + kmer_len;
    assert(r.rs < MAX_POS);

    if (last_cov >= max_cov) {
      max_cov = last_cov;
      best = candidate_regions.size();
    }
    candidate_regions.push_back(r);
  }

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  posvTime += elapsed.count();
}

//#include "ksort.h"
//#define heap_lt(a, b) ((a).x > (b).x)
//KSORT_INIT(heap, mm128_t, heap_lt)
//static mm128_t *collect_seed_hits_heap(void *km,
//                                       const mm_mapopt_t *opt,
//                                       int max_occ,
//                                       const mm_idx_t *mi,
//                                       const char *qname,
//                                       const mm128_v *mv,
//                                       int qlen,
//                                       int64_t *n_a,
//                                       int *rep_len,
//                                       int *n_mini_pos,
//                                       uint64_t **mini_pos) {
//  int i, n_m, heap_size = 0;
//  int64_t j, n_for = 0, n_rev = 0;
//  mm_seed_t *m;
//  mm128_t *a, *heap;
//
//  m = mm_collect_matches(km,
//                         &n_m,
//                         qlen,
//                         max_occ,
//                         opt->max_max_occ,
//                         opt->occ_dist,
//                         mi,
//                         mv,
//                         n_a,
//                         rep_len,
//                         n_mini_pos,
//                         mini_pos);
//
//  fprintf(stderr, "\n hahaha: %d, %d", n_m, *n_a);
//
//  heap = (mm128_t *) kmalloc(km, n_m * sizeof(mm128_t));
//  a = (mm128_t *) kmalloc(km, *n_a * sizeof(mm128_t));
//
//  for (i = 0, heap_size = 0; i < n_m; ++i) {
//    if (m[i].n > 0) {
//      heap[heap_size].x = m[i].cr[0];
//      heap[heap_size].y = (uint64_t) i << 32;
//      ++heap_size;
//    }
//  }
//  ks_heapmake_heap(heap_size, heap);
//  while (heap_size > 0) {
//    mm_seed_t *q = &m[heap->y >> 32];
//    mm128_t *p;
//    uint64_t r = heap->x;
//    int32_t is_self, rpos = (uint32_t) r >> 1;
////    if (!skip_seed(opt->flag, r, q, qname, qlen, mi, &is_self)) {
//    if ((r & 1) == (q->q_pos & 1)) { // forward strand
//      p = &a[n_for++];
//      p->x = (r & 0xffffffff00000000ULL) | rpos;
//      p->y = (uint64_t) q->q_span << 32 | q->q_pos >> 1;
//    } else { // reverse strand
//      p = &a[(*n_a) - (++n_rev)];
//      p->x = 1ULL << 63 | (r & 0xffffffff00000000ULL) | rpos;
//      p->y = (uint64_t) q->q_span << 32 | (qlen - ((q->q_pos >> 1) + 1 - q->q_span) - 1);
//    }
//    p->y |= (uint64_t) q->seg_id << MM_SEED_SEG_SHIFT;
//    if (q->is_tandem) p->y |= MM_SEED_TANDEM;
//    if (is_self) p->y |= MM_SEED_SELF;
////    }
//    // update the heap
//    if ((uint32_t) heap->y < q->n - 1) {
//      ++heap[0].y;
//      heap[0].x = m[heap[0].y >> 32].cr[(uint32_t) heap[0].y];
//    } else {
//      heap[0] = heap[heap_size - 1];
//      --heap_size;
//    }
//    ks_heapdown_heap(0, heap_size, heap);
//  }
//  kfree(km, m);
//  kfree(km, heap);
//
//  // reverse anchors on the reverse strand, as they are in the descending order
//  for (j = 0; j < n_rev >> 1; ++j) {
//    mm128_t t = a[(*n_a) - 1 - j];
//    a[(*n_a) - 1 - j] = a[(*n_a) - (n_rev - j)];
//    a[(*n_a) - (n_rev - j)] = t;
//  }
//  if (*n_a > n_for + n_rev) {
//    memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mm128_t));
//    *n_a = n_for + n_rev;
//  }
//  return a;
//}

//void AccAlign::mm(char *Q, size_t rlen, int err_threshold,
//                  vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
//                  unsigned &fbest, unsigned &rbest){
//  mm128_v mv = {0,0,0};
//  void *km = nullptr;
//
//  // cal minimizer
//  auto start = std::chrono::system_clock::now();
//  mm_sketch(km, Q, rlen, mi->w, mi->k, 0, mi->flag&MM_I_HPC, &mv);
//  auto end = std::chrono::system_clock::now();
//  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//  mm_cal += elapsed.count();
//
//  fetch_candidates(mv, int32_t mid_occ, size_t rlen, int err_threshold,
//      vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
//      unsigned &fbest, unsigned &rbest)
//
//  if (!fcandidate_regions.size() && !rcandidate_regions.size()){
//    mid_occ = 5000;
//    start = std::chrono::system_clock::now();
//    m = mm_collect_matches(km, &n_m0, rlen, mid_occ, max_max_occ, occ_dist,
//                           mi, &mv, &n_a, &rep_len, &n_mini_pos, &mini_pos);
//    end = std::chrono::system_clock::now();
//    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    mm_fetch += elapsed.count();
//
//    start = std::chrono::system_clock::now();
//    collect_seed_hits_priorityqueue(n_m0, n_a, rlen, err_threshold,  m, fcandidate_regions, rcandidate_regions, fbest, rbest);
//    end = std::chrono::system_clock::now();
//    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    mm_hit_cnt+= elapsed.count();
//  }
//}