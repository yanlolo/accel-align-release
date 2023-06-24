#include "cmdline.hpp"

#include "ext/args.hxx"
#include "arguments.hpp"

class Version {};


CommandLineOptions parse_command_line_arguments(int argc, char **argv, bool use_strobealign) {

  args::ArgumentParser parser("Accel-Align with Support for Strobealign Extension");
  parser.helpParams.showTerminator = false;
  parser.helpParams.helpindent = 20;
  parser.helpParams.width = 90;
  parser.helpParams.programName = "strobeacc";
  parser.helpParams.shortSeparator = " ";

  args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});
  args::ActionFlag version(parser, "version", "Print version and exit", {"version"}, []() { throw Version(); });


  if(use_strobealign) {
    return do_strobealign_setup(parser, argc, argv);
  }
  return do_accalign_setup(parser, argc, argv);
}



CommandLineOptions do_strobealign_setup(args::ArgumentParser parser, int argc, char **argv) {
  // Threading
  args::ValueFlag<int> threads(parser, "INT", "Number of threads [3]", {'t', "threads"});
  args::ValueFlag<int> chunk_size(parser, "INT", "Number of reads processed by a worker thread at once [10000]", {"chunk-size"}, args::Options::Hidden);

  args::Group io(parser, "Input/output:");
  args::ValueFlag<std::string> o(parser, "PATH", "redirect output to file [stdout]", {'o'});
  args::Flag v(parser, "v", "Verbose output", {'v'});
  args::Flag no_progress(parser, "no-progress", "Disable progress report (enabled by default if output is a terminal)", {"no-progress"});
  args::Flag eqx(parser, "eqx", "Emit =/X instead of M CIGAR operations", {"eqx"});
  args::Flag x(parser, "x", "Only map reads, no base level alignment (produces PAF file)", {'x'});
  args::Flag U(parser, "U", "Suppress output of unmapped reads", {'U'});
  args::Flag interleaved(parser, "interleaved", "Interleaved input", {"interleaved"});
  args::ValueFlag<std::string> rgid(parser, "ID", "Read group ID", {"rg-id"});
  args::ValueFlagList<std::string> rg(parser, "TAG:VALUE", "Add read group metadata to SAM header (can be specified multiple times). Example: SM:samplename", {"rg"});

  args::ValueFlag<int> N(parser, "INT", "Retain at most INT secondary alignments (is upper bounded by -M and depends on -S) [0]", {'N'});
  args::ValueFlag<std::string> index_statistics(parser, "PATH", "Print statistics of indexing to PATH", {"index-statistics"});
  args::Flag i(parser, "index", "Do not map reads; only generate the strobemer index and write it to disk. If read files are provided, they are used to estimate read length", {"create-index", 'i'});
  args::Flag use_index(parser, "use_index", "Use a pre-generated index previously written with --create-index.", { "use-index" });
  args::Flag use_strobealign(parser, "use_strobealign", "Use Strobealign mode", { "strobe-mode" });


  args::Group seeding_group(parser, "Seeding:");
  SeedingArguments *seedingArgumentsPointer = NULL;


  auto seeding = SeedingArguments{parser};
  seedingArgumentsPointer = &seeding;


  args::Group alignment(parser, "Alignment:");
  args::ValueFlag<int> A(parser, "INT", "Matching score [2]", {'A'});
  args::ValueFlag<int> B(parser, "INT", "Mismatch penalty [8]", {'B'});
  args::ValueFlag<int> O(parser, "INT", "Gap open penalty [12]", {'O'});
  args::ValueFlag<int> E(parser, "INT", "Gap extension penalty [1]", {'E'});
  args::ValueFlag<int> end_bonus(parser, "INT", "Soft clipping penalty [10]", {'L'});

  args::Group search(parser, "Search parameters:");
  args::ValueFlag<float> f(parser, "FLOAT", "Top fraction of repetitive strobemers to filter out from sampling [0.0002]", {'f'});
  args::ValueFlag<float> S(parser, "FLOAT", "Try candidate sites with mapping score at least S of maximum mapping score [0.5]", {'S'});
  args::ValueFlag<int> M(parser, "INT", "Maximum number of mapping sites to try [20]", {'M'});
  args::ValueFlag<int> R(parser, "INT", "Rescue level. Perform additional search for reads with many repetitive seeds filtered out. This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than default makes strobealign significantly slower but more accurate. R <= 1 deactivates rescue and is the fastest.", {'R'});

  args::Positional<std::string> ref_filename(parser, "reference", "Reference in FASTA format", args::Options::Required);
  args::Positional<std::string> reads1_filename(parser, "reads1", "Reads 1 in FASTA or FASTQ format, optionally gzip compressed");
  args::Positional<std::string> reads2_filename(parser, "reads2", "Reads 2 in FASTA or FASTQ format, optionally gzip compressed");

  try {
    parser.ParseCLI(argc, argv);
  }
  catch (const args::Completion& e) {
    std::cout << e.what();
    exit(EXIT_SUCCESS);
  }
  catch (const args::Help&) {
    std::cout << parser;
    exit(EXIT_SUCCESS);
  }
  catch (const Version& e) {
    std::cout << "Version 0.1" << std::endl;
    exit(EXIT_SUCCESS);
  }
  catch (const args::Error& e) {
    std::cerr << parser;
    std::cerr << "Error: " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  CommandLineOptions opt;

  // Threading
  if (threads) { opt.n_threads = args::get(threads); }
  if (chunk_size) { opt.chunk_size = args::get(chunk_size); }

  // Input/output
  if (o) { opt.output_file_name = args::get(o); opt.write_to_stdout = false; }
  if (v) { opt.verbose = true; }
  if (no_progress) { opt.show_progress = false; }
  if (eqx) { opt.cigar_eqx = true; }
  if (x) { opt.is_sam_out = false; }
  if (U) { opt.output_unmapped = false; }
  if (rgid) { opt.read_group_id = args::get(rgid); }
  if (rg) { opt.read_group_fields = args::get(rg); }
  if (N) { opt.max_secondary = args::get(N); }
  if (index_statistics) { opt.logfile_name = args::get(index_statistics); }
  if (i) { opt.only_gen_index = true; }
  if (use_index) { opt.use_index = true; }
  if (use_strobealign) {opt.use_strobealign = true;}

  // Seeding
  if(use_strobealign) {
    if (seedingArgumentsPointer->r) { opt.r = args::get(seedingArgumentsPointer->r); opt.r_set = true; }
    if (seedingArgumentsPointer->m) { opt.max_seed_len = args::get(seedingArgumentsPointer->m); opt.max_seed_len_set = true; }
    if (seedingArgumentsPointer->k) { opt.k = args::get(seedingArgumentsPointer->k); opt.k_set = true; }
    if (seedingArgumentsPointer->l) { opt.l = args::get(seedingArgumentsPointer->l); opt.l_set = true; }
    if (seedingArgumentsPointer->u) { opt.u = args::get(seedingArgumentsPointer->u); opt.u_set = true; }
    if (seedingArgumentsPointer->s) { opt.s = args::get(seedingArgumentsPointer->s); opt.s_set = true; }
    if (seedingArgumentsPointer->c) { opt.c = args::get(seedingArgumentsPointer->c); opt.c_set = true; }
  }

  // Alignment
  // if (n) { n = args::get(n); }
  if (A) { opt.A = args::get(A); }
  if (B) { opt.B = args::get(B); }
  if (O) { opt.O = args::get(O); }
  if (E) { opt.E = args::get(E); }
  if (end_bonus) { opt.end_bonus = args::get(end_bonus); }

  // Search parameters
  if (f) { opt.f = args::get(f); }
  if (S) { opt.dropoff_threshold = args::get(S); }
  if (M) { opt.maxTries = args::get(M); }
  if (R) { opt.R = args::get(R); }

  // Reference and read files
  opt.ref_filename = args::get(ref_filename);
  opt.reads_filename1 = args::get(reads1_filename);
  opt.is_interleaved = bool(interleaved);

  if (reads2_filename) {
    opt.reads_filename2 = args::get(reads2_filename);
    opt.is_SE = false;
  } else if (interleaved) {
    opt.is_SE = false;
  } else {
    opt.reads_filename2 = std::string();
    opt.is_SE = true;
  }

  if (opt.use_index && opt.only_gen_index) {
    std::cerr << "Error: Options -i and --use-index cannot be used at the same time" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (opt.reads_filename1.empty() && !opt.only_gen_index) {
    std::cerr << "Error: At least one file with reads must be specified." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (opt.only_gen_index && !(opt.r_set || !opt.reads_filename1.empty())) {
    std::cerr << "Error: The target read length needs to be known when generating an index.\n"
                 "Use -r to set it explicitly or let the program estimate it by providing at least one read file.\n";
    exit(EXIT_FAILURE);
  }

  return opt;
}



CommandLineOptions do_accalign_setup(args::ArgumentParser parser, int argc, char **argv){

  args::Flag use_strobealign(parser, "use_strobealign", "Use Strobealign mode", { "strobe-mode" });

  // t
  args::ValueFlag<int> threads(parser, "INT", "Number of threads [3]", {'t', "threads"});

  // l
  args::ValueFlag<int> l(parser, "INT", "k-mer size", {'l'});

  // o
  args::ValueFlag<std::string> o(parser, "STRING", "Name of the output file", {'o'});

  // e
  args::ValueFlag<std::string> e(parser, "STRING", "Name of embed file", {'e'});

  // b
  args::ValueFlag<std::string> b(parser, "STRING", "Name of batch file", {'b'});

  // p
  args::ValueFlag<int> p(parser, "INT", "Paired end distance", {'p'});

  // x
  args::Flag x(parser, "x", "Alignment-free mode", {'x'});

  // w
  args::Flag w(parser, "w", "Use WFA for extension. KSW used by default", {'w'});

  // d
  args::Flag d(parser, "d", "Disable embedding, extend all candidates from seeding (this mode is super slow, only for benchmark)", {'d'});

  // m
  args::Flag m(parser, "m", "Enable Minimizer", {'m'});

  // s
  args::Flag s(parser, "s", "Use bisulfite sequencing alignment mode", {'s'});

  args::Positional<std::string> ref_filename(parser, "reference", "Reference in FASTA format", args::Options::Required);
  args::Positional<std::string> reads1_filename(parser, "reads1", "Reads 1 in FASTA or FASTQ format, optionally gzip compressed");
  args::Positional<std::string> reads2_filename(parser, "reads2", "Reads 2 in FASTA or FASTQ format, optionally gzip compressed");

  try {
    parser.ParseCLI(argc, argv);
  }
  catch (const args::Completion& e) {
    std::cout << e.what();
    exit(EXIT_SUCCESS);
  }
  catch (const args::Help&) {
    std::cout << parser;
    exit(EXIT_SUCCESS);
  }
  catch (const Version& e) {
    std::cout << "Version 0.1" << std::endl;
    exit(EXIT_SUCCESS);
  }
  catch (const args::Error& e) {
    std::cerr << parser;
    std::cerr << "Error: " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  CommandLineOptions opt;

  // Threading
  if (threads) { opt.n_threads = args::get(threads); }

  if(l) {opt.l = args::get(l);}
  if(o) {opt.o = args::get(o);}
  if(e) {opt.e = args::get(e);}
  if(b) {opt.b = args::get(b);}
  if(p) {opt.p = args::get(p);}
  if(x) {opt.x = false;}
  if(w) {opt.w = true;}
  if(d) {opt.d = true;}
  if(m) {opt.m = true;}
  if(s) {opt.bs = true;}


  // Reference and read files
  opt.ref_filename = args::get(ref_filename);
  opt.reads_filename1 = args::get(reads1_filename);

  if (reads2_filename) {
    opt.reads_filename2 = args::get(reads2_filename);
    opt.is_SE = false;
  } else {
    opt.reads_filename2 = std::string();
    opt.is_SE = true;
  }

  if (opt.reads_filename1.empty() && !opt.only_gen_index) {
    std::cerr << "Error: At least one file with reads must be specified." << std::endl;
    exit(EXIT_FAILURE);
  }

  return opt;
}
