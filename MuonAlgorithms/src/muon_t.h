#ifndef MUON_T_H
#define MUON_T_H

#include "ap_int.h"

// FIXME : only test values
struct muon_t {
    ap_uint<15> pt;
    // ap_uint<12> theta;
    ap_int<16> eta; // cfr the track class
    ap_int<12> phi;
};

struct iso_muon_t : public muon_t {
    ap_uint<1> iso;
};

#endif // MUON_T_H
