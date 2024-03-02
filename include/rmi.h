#ifndef __RMI__
#define __RMI__

#include <cstdint>
#include <cstddef>

typedef bool (*load_type)(char const*);
typedef uint64_t (*lookup_type)(uint64_t, size_t*);
typedef void (*cleanup_type)();

class RMI {
  // private
  void* library_handle;
  bool is_init;
  load_type rmi_load;
  lookup_type rmi_lookup;
  cleanup_type rmi_cleanup;

  public:
    RMI();
    ~RMI();
    void init(const char *library_prefix);
    uint64_t lookup(uint64_t key, size_t* err);
};

#endif