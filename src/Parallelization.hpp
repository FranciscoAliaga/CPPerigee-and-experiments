#ifndef PARALLELIZATION_UTILS
#define PARALLELIZATION_UTILS

#include <iostream>
#include <thread>
#include <utility> // std::pair
#include <vector>


// parallelization utils
const unsigned int DEFAULT_PARALLELIZATION_THREADS = 4;

using thread_work_range = std::pair<unsigned int,unsigned int>;
using thread_work_ranges = std::vector<thread_work_range>;

// splits the range [0...V-1] into P ranges for parallelized tasks
thread_work_ranges getWorkRanges(unsigned int V, unsigned int P){
    thread_work_ranges ranges{};
    for(unsigned int part = 0; part < P; part++){
        auto lower = ( part*V   ) / P;
        auto upper = ((part+1)*V) / P;
        if (part==P-1) upper = V;
        ranges.push_back({lower,upper});
    }
    return ranges;
    // P tasks should be performed on ranges [lower,upper-1]
}

void GET_HARDWARE_THREADS(unsigned int* threads){
    *threads = std::thread::hardware_concurrency();
    if (*threads < DEFAULT_PARALLELIZATION_THREADS){
        *threads = DEFAULT_PARALLELIZATION_THREADS;
    }
    std::cout << "\n\navailable hardware threads : " << *threads << "\n\n";
}


#endif