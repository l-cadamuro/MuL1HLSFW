#ifndef MINI_SORTER_H
#define MINI_SORTER_H

#include "dataformats.h"

// #define MAX_SPAR 2147483647
// #define NULL_DATA 2147483647

#define N_MU 12

void find_clusters(muon_t in_muons[N_MU], hw_minv_t out_masses[N_MU]);

// void tau_3mu(muon_t in_muons[N_MU], hw_minv_t out_masses[N_MU]);

#endif