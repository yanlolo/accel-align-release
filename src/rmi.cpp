#include "header.h"

// A wrapper for the RMI library 

// An utility function to perform /path1/path2/file --> /path1/path2
std::string get_parent_directory(const std::string& path) {
    // Find the last occurrence of the directory separator
    size_t last_separator = path.find_last_of("/\\");

    if (last_separator != std::string::npos) {
        // Extract the substring up to the last separator
        return path.substr(0, last_separator);
    }

    // No separator found, return the original path
    return path;
}

/**
 * Initializes the wrapper of the RMI library.
 *
 * @param library_prefix The prefix of the library we want to use, that is, the path and the name WITHOUT the .so extension.
 * @return A wrapper object.
 */
void RMI::init(const char *library_prefix) {                  // ./data/hg37_index/hg37_index
    auto start = std::chrono::system_clock::now();
    auto lib_str = std::string(library_prefix);
    std::string library_name = lib_str + ".so";         // ./data/hg37_index/hg37_index.so
    std::string library_sym = lib_str + ".sym";         // ./data/hg37_index/hg37_index.sym

    std::cerr << "Loading library " << library_name << std::endl;
    // first, open the library
    library_handle = dlopen(library_name.c_str(), RTLD_LAZY);
    assert(library_handle);
    // now, read the symbols from library_sym and load them 
    std::ifstream f(library_sym);
    std::string line;
    /*      -------- LOAD --------      */
    std::getline(f, line);
    rmi_load = reinterpret_cast<load_type>(dlsym(library_handle, line.c_str()));
    assert(rmi_load);
    /*      ------- LOOKUP -------      */
    std::getline(f, line);
    rmi_lookup = reinterpret_cast<lookup_type>(dlsym(library_handle, line.c_str()));
    assert(rmi_lookup);
    /*      ------- CLEANUP ------      */
    std::getline(f, line);
    rmi_cleanup = reinterpret_cast<cleanup_type>(dlsym(library_handle, line.c_str()));
    assert(rmi_cleanup);
    f.close();
    
    // load the index
    std::string param_path = get_parent_directory(lib_str) + "/rmi_data";
    int done = rmi_load(param_path.c_str());
    assert(done);
    
    // print elapsed time
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cerr << "Setup RMI library in " << elapsed.count() << " us" << std::endl;
    is_init = true;
}
// Constructor -> simply initializes is_init
RMI::RMI() {
    is_init = false;
}

// Destructor -> removes the index and closes the dynamic library object
RMI::~RMI() {
    if (is_init) {
        rmi_cleanup();
        dlclose(library_handle);
    }
}

// lookup
uint64_t RMI::lookup(uint64_t key, size_t* err) {
    assert(is_init);
    return rmi_lookup(key, err);
}