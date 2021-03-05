//
// Created by wyx on 19-6-20.
//

#ifndef CONFIGURABLE_DEDUP_GLOBAL_H
#define CONFIGURABLE_DEDUP_GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <ctime>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <list>
#include <sstream>
#include <time.h>
#include <chrono>
#include <deque>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <bitset>
#include <utility>
#include <math.h>
#include <inttypes.h>
#include <set>
#include <unordered_set>
#include <queue>
using namespace std;

#define PRIMER_CAPACITY 736

extern string g_blast_Result;
extern long g_chunk_size;
extern int g_primer_capacity;
extern int g_tube_capacity;
extern long long int g_total_strand_number;
extern bool g_if_decomposition_on_primer_graph;
extern bool g_if_baseline;
extern bool g_if_control_payload_totalsize;
extern double g_threshold_of_totalsize;


int Parse(string cfgfile);

typedef std::uint64_t hash_t;
constexpr hash_t prime = 0x100000001B3ull;
constexpr hash_t basis = 0xCBF29CE484222325ull;
constexpr hash_t hash_(char const* str, hash_t last_value = basis)
{
    return *str ? hash_(str+1, (*str ^ last_value) * prime) : last_value;
}
#endif //CONFIGURABLE_DEDUP_GLOBAL_H