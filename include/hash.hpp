#pragma once

// ------------------- xxhash ------------------- //
#define XXH_STATIC_LINKING_ONLY /* access advanced declarations */
#define XXH_IMPLEMENTATION      /* access definitions */
#define XXH_INLINE_ALL          /* make all functions inline -- improvements of 200%! */

#include "xxhash.h"
#include <functional>

// --------------- hard-coded mods --------------- //
#define MOD_29      ((1UL<<29)-1)
#define MOD_PRIME   1073741651
#define MOD_LPRIME  2861333663
#define MOD_32      (1UL<<32)

#define SEED        426942
#define IN_LEN      8

// --------------- xxhash wrapper --------------- //
using XXHash = std::function<uint64_t(const uint64_t*)>;

inline void bind_xxhash(const unsigned type, XXHash& hash_fn) {
  if (type==32) {
    hash_fn = [](const uint64_t* input) -> uint64_t {
      return static_cast<uint64_t>(XXH32(input, IN_LEN, SEED));
    };
  } else if (type==64) {
    hash_fn = [](const uint64_t* input) -> uint64_t {
      return XXH64(input, IN_LEN, SEED);
    };
  } else {
    hash_fn = [](const uint64_t* input) -> uint64_t {
      return *input;
    };
  }
}

// XXH32_hash_t XXH32 (const void* input, size_t length, XXH32_hash_t seed);
// XXH64_hash_t XXH64 (const void* input, size_t length, XXH64_hash_t seed);