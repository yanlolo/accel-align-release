#!/bin/bash
set -e  # Stop if there is a failure
_redo_=0

# Default values
ref_name=""
thread_number=$(nproc --all)
#clean=false
kmer_len=32

# Function to display usage instructions
usage() {
    echo -e "\n\033[1;96mbash rmi.sh [OPTIONS] <reference.fna>\033[0m"
    echo -e "Builds a learned index for the <reference.fna> reference string."
    echo "Options:"
    echo "  -t, --threads  THREADS  The number of threads to be used [all]"
    echo "  -l, --len      LEN      The length of the kmer [32]"
    #echo "      --clean             Delete the index related to <reference.fna>"
    echo -e "  -h, --help              Display this help message\n"
    exit 1
}

################################### INTRO ###################################

# Parse command line options
while [[ $# -gt 1 ]]; do
    key="$1"
    case $key in
        -t|--threads)
            thread_number=$2
            shift 2
            ;;
        -l|--len)
            kmer_len=$2
            shift 2
            ;;
        # --clean)
        #     clean=true
        #     shift 1
        #     ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $key"
            usage
            ;;
    esac
done
ref_name="$1"
if [ -z "$ref_name" ] || [ ! -e "$ref_name" ]; then
    echo "Error: Please provide a valid reference genome file."
    usage
fi
ref_name=$(realpath -se $ref_name)              # ./data/hg37.fna

INITIAL_DIR=$(pwd)
_source_dir_=$(dirname "$0")
BASE_DIR=$(readlink -f "$_source_dir_")     # /home/ilaria/Documents/uny/project/accel-align-rmi-dev
cd $BASE_DIR
                             
# The output file will be
dir_name=$(dirname $ref_name)               # ./data
base_name=$(basename $ref_name .fna)        # hg37

# Make output directory
OUTPUT_DIR="${dir_name}/${base_name}_index${kmer_len}"   # ./data/hg37_index32
OUTPUT_DIR=$(realpath -se $OUTPUT_DIR)

echo -e "\n\033[1;96m [rmi.sh] \033[0mBuilding index on file $base_name.fna"
echo -e "            --- Using $thread_number threads."
echo -e "            --- kmer length is $kmer_len."

# Compute the number of bit in the key
bit_len=$((kmer_len * 2))
echo -e "            --- expected key length is $bit_len bits."
if ((bit_len > 64)); then
    echo -e "\n\033[1;31m  [error!]  \033[0mSorry, this kmer length is too big!\n            The maximum supported length is 32.\n"
    exit 1
fi
if ((bit_len > 32)); then
    bit_len=64
fi
if ((bit_len < 32)); then
    bit_len=32
fi
keys_name="${OUTPUT_DIR}/keys_uint${bit_len}"       # ./data/hg37_index/keys_unit32 OR ./data/hg37_index/keys_unit64
pos_name="${OUTPUT_DIR}/pos_uint32"       # ./data/hg37_index/pos_unit32

echo -e "            --- generating files 'keys_uint${bit_len}' and 'pos_uint32'."

################################### KEY_GEN ###################################

echo -e "\n\033[1;96m [rmi.sh] \033[0mCompiling the key_gen program..."
make key_gen

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR                             # ----> NOW WE ARE IN hg37_index/keys_unit32

echo -e "\n\033[1;96m [rmi.sh] \033[0mRunning key_gen..."
if [ ! -e $keys_name ] || [ ! -e $pos_name ]; then
  # The file does not exist, so execute the command
  "${BASE_DIR}/key_gen" -l $kmer_len $ref_name
else
  # The file exists, so ask the user before executing
  read -ep $'\033[1;33m [rmi.sh] \033[0mkey_gen output already exists. Do you want to execute the command anyway? [y/N] ' choice
  case "$choice" in 
    y|Y )
      _redo_=1 
      "${BASE_DIR}/key_gen" -l $kmer_len $ref_name 
      ;;
    * ) 
      echo -e "\033[1;33m [rmi.sh] \033[0mcommand not executed" ;;
  esac
fi

################################### INDEX ###################################

echo -e "\n\033[1;96m [rmi.sh] \033[0mCompiling the index generation program..."
# Build the index - if it has not being compiled yet
if ! [ -e "${BASE_DIR}/rmi/target/release/rmi" ] && ! [ -e ./rmi ]; then
  cd "${BASE_DIR}/rmi" && cargo build --release
  cd $OUTPUT_DIR
fi
if [[ ! -e ./rmi ]]; then
  # Copy it
  cp "${BASE_DIR}/rmi/target/release/rmi" .
fi

echo -e "\n\033[1;96m [rmi.sh] \033[0mRunning the optimizer..."
# Run RMI optimization - if not present
if [ ! -e optimizer.out ] || [ "$_redo_" == "1" ]; then
  # The file does not exist, so execute the command
  ./rmi --threads $thread_number --optimize optimizer.json "./keys_uint${bit_len}" > optimizer.out
else
  # The file exists, so ask the user before executing
  read -ep $'\033[1;33m [rmi.sh] \033[0moptimizer output already exists. Do you want to execute the command anyway? [y/N] ' choice
  case "$choice" in 
    y|Y )
      _redo_=1  
      ./rmi --threads $thread_number --optimize optimizer.json "./keys_uint${bit_len}" > optimizer.out
      ;;
    * ) 
      echo -e "\033[1;33m [rmi.sh] \033[0mcommand not executed" ;;
  esac
fi

# Print optimization result
echo -n -e "\t" && head -1 optimizer.out && tail -10 optimizer.out | cat -n

# Chose the best model and train it
# The chosen parameters in rmi_type.txt as type, branching_factor, size (KB), avg_log2_err
if [ ! -e rmi_type.txt ] || [ "$_redo_" == "1" ]; then
  # The file does not exist, so execute the command
  bash "${BASE_DIR}/utilities/model_select.sh"
else
  # The file exists, so ask the user before executing
  read type branching size avg_err max_err b_time < rmi_type.txt
  echo -e "\033[1;33m [rmi.sh] \033[0mmodel has already been chosen:"
  echo -e "               | MODEL\t\t\t$type"
  echo -e "               | BRANCHING FACTOR\t$branching"
  echo -e "               | SIZE (B)\t\t$size"
  echo -e "               | AVG LOG2 ERROR\t\t$avg_err"
  echo -e "               | MAX LOG2 ERROR\t\t$max_err"
  echo -e "               | BUILD TIME (ns)\t$b_time"
  read -ep $'            do you want to train a new one? [y/N] ' choice
  case "$choice" in 
    y|Y ) 
      _redo_=1 
      bash "${BASE_DIR}/utilities/model_select.sh"
      ;;
    * ) 
      echo -e "\033[1;33m [rmi.sh] \033[0mcommand not executed" ;;
  esac
fi

echo -e "\n\033[1;96m [rmi.sh] \033[0mBuilding the index..."
if [ ! -e rmi.h ] || [ "$_redo_" == "1" ]; then
  # Get the parameters
  read type branching _ _ _ _ < rmi_type.txt
  echo -e "\n\033[1;96m [rmi.sh] \033[0mTraining the $type ($branching) index"
  # Train the model
  ./rmi "./keys_uint${bit_len}" rmi $type $branching
else
  echo -e "\n\033[1;96m [rmi.sh] \033[0mIndex already exists!\n"
fi

################################### SHARED OBJECT ###################################
# Copy Makefile in the current directory
if [ ! -e Makefile ]; then
  cp "${BASE_DIR}/rmi/Makefile" .
fi

# make ! -> create the shared library
make

cd $INITIAL_DIR

echo -e "\n\033[1;32m [rmi.sh] \033[0mDone!\n"