bool AccAlign::fastq(const char *F1, const char *F2, bool enable_gpu) {

  bool is_paired = false;

  gzFile in1 = gzopen(F1, "rt");
  if (in1 == Z_NULL)
    return false;

  gzFile in2 = Z_NULL;
  if (strlen(F2) > 0) {
    is_paired = true;
    in2 = gzopen(F2, "rt");
    if (in2 == Z_NULL)
      return false;
  }

  cerr << "Reading fastq file " << F1 << ", " << F2 << "\n";

  // start CPU and GPU master threads, they consume reads from inputQ
  // dataQ is to re-use
  tbb::concurrent_bounded_queue<ReadCnt> inputQ;
  tbb::concurrent_bounded_queue<ReadCnt> outputQ;
  tbb::concurrent_bounded_queue<ReadPair> dataQ;

  thread cpu_thread = thread(&AccAlign::cpu_root_fn, this, &inputQ, &outputQ);
  thread out_thread = thread(&AccAlign::output_root_fn, this, &outputQ, &dataQ);

  auto start = std::chrono::system_clock::now();

  int total_nreads = 0, nreads_per_vec = 0, vec_index = 0, vec_size = 50;

  int batch_size = BATCH_SIZE;
  if (is_paired)
    batch_size /= 2;
  Read *reads[vec_size];
  reads[vec_index] = new Read[batch_size];
  Read *reads2[vec_size];
  if (is_paired) {
    reads2[vec_index] = new Read[batch_size];
  }

  bool neof1 = (!gzeof(in1) && gzgetc(in1) != EOF);
  bool neof2 = (!is_paired || (!gzeof(in2) && gzgetc(in2) != EOF));
  while (vec_index < vec_size && neof1 && neof2) {
    Read &r = *(reads[vec_index] + nreads_per_vec);
    in1 >> r;

    if (!strlen(r.seq)) {
      break;
    }

    if (is_paired) {
      Read &r2 = *(reads2[vec_index] + nreads_per_vec);
      in2 >> r2;

      if (!strlen(r2.seq)) {
        break;
      }
    }
    neof1 = (!gzeof(in1) && gzgetc(in1) != EOF);
    neof2 = (!is_paired || (!gzeof(in2) && gzgetc(in2) != EOF));

    ++nreads_per_vec;

    if (nreads_per_vec == batch_size) {
      if (is_paired)
        inputQ.push(make_tuple(reads[vec_index], reads2[vec_index], batch_size));
      else
        inputQ.push(make_tuple(reads[vec_index], (Read *) NULL, batch_size));

      vec_index++;

      if (vec_index < vec_size) {
        reads[vec_index] = new Read[batch_size];
        if (is_paired)
          reads2[vec_index] = new Read[batch_size];
      }

      total_nreads += nreads_per_vec;
      nreads_per_vec = 0;
    }
  }

  ReadPair cur_vec = make_tuple((Read *) NULL, (Read *) NULL);

  // the nb of reads is less than vec_size *BATCH_SIZE, and there are some reads not pushed to inputQ
  if (nreads_per_vec && vec_index < vec_size) {
    // the remaining reads
    if (is_paired)
      inputQ.push(make_tuple(reads[vec_index], reads2[vec_index], nreads_per_vec));
    else
      inputQ.push(make_tuple(reads[vec_index], (Read *) NULL, nreads_per_vec));

    total_nreads += nreads_per_vec;
  } else {
    // still have reads not loaded
    dataQ.pop(cur_vec);

    while (neof1 && neof2) {
      Read &r = *(std::get<0>(cur_vec) + nreads_per_vec);
      in1 >> r;

      if (!strlen(r.seq)) {
        break;
      }

      if (is_paired) {
        Read &r2 = *(std::get<1>(cur_vec) + nreads_per_vec);
        in2 >> r2;

        if (!strlen(r2.seq)) {
          break;
        }
      }
      neof1 = (!gzeof(in1) && gzgetc(in1) != EOF);
      neof2 = (!is_paired || (!gzeof(in2) && gzgetc(in2) != EOF));

      ++nreads_per_vec;

      if (nreads_per_vec == batch_size) {
        inputQ.push(make_tuple(std::get<0>(cur_vec), std::get<1>(cur_vec), batch_size));
        dataQ.pop(cur_vec);
        total_nreads += nreads_per_vec;
        nreads_per_vec = 0;
      }
    }

    // the remaining reads
    if (nreads_per_vec) {
      total_nreads += nreads_per_vec;
      inputQ.push(make_tuple(std::get<0>(cur_vec), std::get<1>(cur_vec), nreads_per_vec));
    }
  }

  gzclose(in1);
  if (is_paired) {
    gzclose(in2);
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  input_io_time += elapsed.count();

  cerr << "done reading " << total_nreads << " reads from fastq file " << F1 << ", " << F2 << " in " <<
       input_io_time / 1000000.0 << " secs\n";

  ReadCnt sentinel = make_tuple((Read *) NULL, (Read *) NULL, 0);
  inputQ.push(sentinel);

  int size = vec_index < vec_size ? vec_index : vec_size;
  start = std::chrono::system_clock::now();
  if (total_nreads % batch_size == 0) {
    //because the last popped cur_vec has not been pushed back
    size -= 1;
    delete[] std::get<0>(cur_vec);

    if (is_paired)
      delete[] std::get<1>(cur_vec);
  }
  for (int i = 0; i < size; i++) {
    dataQ.pop(cur_vec);
    delete[] std::get<0>(cur_vec);

    if (is_paired)
      delete[] std::get<1>(cur_vec);
  }
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  delTime += elapsed.count();

  cpu_thread.join();

  outputQ.push(sentinel);
  out_thread.join();

  cerr << "Processed " << total_nreads << " in total \n";
  return true;
}


void AccAlign::output_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *outputQ,
                              tbb::concurrent_bounded_queue<ReadPair> *dataQ) {
  cerr << "Extension and output function starting.." << endl;

  unsigned nreads = 0;
  tbb::concurrent_bounded_queue<ReadCnt> *targetQ = outputQ;
  do {
    ReadCnt gpu_reads;
    targetQ->pop(gpu_reads);
    nreads = std::get<2>(gpu_reads);
    if (!nreads) {
      targetQ->push(gpu_reads);   //put sentinel back
      break;
    }
    align_wrapper(0, 0, nreads, std::get<0>(gpu_reads), std::get<1>(gpu_reads), dataQ);
  } while (1);

  cerr << "Extension and output function quitting...\n";
}

class Parallel_mapper {
  Read *all_reads1;
  Read *all_reads2;
  AccAlign *acc_obj;

 public:
  Parallel_mapper(Read *_all_reads1, Read *_all_reads2, AccAlign *_acc_obj) :
      all_reads1(_all_reads1), all_reads2(_all_reads2), acc_obj(_acc_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    if (!all_reads2) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        acc_obj->map_read_wrapper(*(all_reads1 + i));
      }
    } else {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        acc_obj->map_paired_read_wrapper(*(all_reads1 + i), *(all_reads2 + i));
      }
    }
  }
};

void AccAlign::cpu_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *inputQ,
                           tbb::concurrent_bounded_queue<ReadCnt> *outputQ) {
  cerr << "CPU Root function starting.." << endl;

  tbb::concurrent_bounded_queue<ReadCnt> *targetQ = inputQ;
  int nreads = 0, total = 0;
  do {
    ReadCnt cpu_readcnt;
    targetQ->pop(cpu_readcnt);
    nreads = std::get<2>(cpu_readcnt);
    total += nreads;
    if (nreads == 0) {
      inputQ->push(cpu_readcnt);    // push sentinel back
      break;
    }

    tbb::task_scheduler_init init(g_ncpus);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nreads),
                      Parallel_mapper(std::get<0>(cpu_readcnt), std::get<1>(cpu_readcnt), this)
    );

    outputQ->push(cpu_readcnt);
  } while (1);

  cerr << "Processed " << total << " reads in cpu \n";
  cerr << "CPU Root function quitting.." << endl;
}

void AccAlign::align_wrapper(int tid, int soff, int eoff, Read *ptlread, Read *ptlread2,
                             tbb::concurrent_bounded_queue<ReadPair> *dataQ) {

  if (!ptlread2) {
    // single-end read alignment
    string sams[eoff];
    tbb::task_scheduler_init init(g_ncpus);
    tbb::parallel_for(tbb::blocked_range<size_t>(soff, eoff), Tbb_aligner(ptlread, sams, this));

    auto start = std::chrono::system_clock::now();
    for (int i = soff; i < eoff; i++) {
      out_sam(sams + i);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    sam_time += elapsed.count();

    dataQ->push(make_tuple(ptlread, (Read *) NULL));
  }
//  else {
//    string sams[2 * eoff];
//    tbb::task_scheduler_init init(g_ncpus);
//    tbb::parallel_for(tbb::blocked_range<size_t>(soff, eoff), Tbb_aligner_paired(ptlread, ptlread2, sams, this));
//
//    auto start = std::chrono::system_clock::now();
//    for (int i = soff; i < 2 * eoff; i++) {
//      out_sam(sams + i);
//    }
//    auto end = std::chrono::system_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    sam_time += elapsed.count();
//
//    dataQ->push(make_tuple(ptlread, ptlread2));
//  }
}

class Tbb_aligner {
  Read *all_reads;
  string *sams;
  AccAlign *acc_obj;

 public:
  Tbb_aligner(Read *_all_reads, string *_sams, AccAlign *_acc_obj) :
      all_reads(_all_reads), sams(_sams), acc_obj(_acc_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      if ((all_reads + i)->strand != '*') {
        if (enable_wfa_extension)
          acc_obj->wfa_align_read(*(all_reads + i));
        else
          acc_obj->align_read(*(all_reads + i));
      }
      acc_obj->snprintf_sam(*(all_reads + i), sams + i);
    }
  }
};

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

void AccAlign::snprintf_sam(Read &R, string *s) {
  auto start = std::chrono::system_clock::now();

  // 50 is the approximate length for all int
  int size;
  if (!enable_extension) {
    size = 50;
  } else {
    size = 50 + strlen(R.seq); //assume the length of cigar will not longer than the read
  }

  string format = "%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\n";
  if (R.strand == '*') {
    size += strlen(R.name) + 2 * strlen(R.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             R.name, 0, "*", 0, 0, "*", R.seq, R.qua, 0, 0);
    *s = buf;
  } else {
    size += strlen(R.name) + get_name(R.ref_id)[R.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    if (R.strand == '+') {
      snprintf(buf, size, format.c_str(),
               R.name, 0, get_name(R.ref_id)[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar,
               R.seq, R.qua, R.nm, R.as);
    } else {
      std::reverse(R.qua, R.qua + strlen(R.qua));
      snprintf(buf, size, format.c_str(),
               R.name, 16, get_name(R.ref_id)[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar,
               R.rev_str, R.qua, R.nm, R.as);
    }
    *s = buf;
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sam_pre_time += elapsed.count();
}

void AccAlign::out_sam(string *s) {
  auto start = std::chrono::system_clock::now();
  {
    if (sam_name.length()) {
      sam_stream << *s;
    } else {
      cout << *s;
    }
  }
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sam_out_time += elapsed.count();
}