#ifndef TRACK_T_H
#define TRACK_T_H

#include "ap_int.h"

#define TRACK_QOVERR_W       15
#define TRACK_PHI_W          12
#define TRACK_TANLAMBDA_W    16
#define TRACK_Z0_W           12
#define TRACK_D0_W           13
#define TRACK_CHI2_W         4
#define TRACK_BENDCHI2_W     3
#define TRACK_HITMASK_W      7
#define TRACK_TKQUALMVA_W    3
#define TRACK_OTHERQUALMVA_W 6
#define TRACK_SPARE_W        5

// from the TDR version 3
struct track_t {
    ap_int<TRACK_QOVERR_W>        qOverR;
    ap_int<TRACK_PHI_W>           phi;
    ap_uint<TRACK_TANLAMBDA_W>    tanLambda;
    ap_uint<TRACK_Z0_W>           z0;
    ap_uint<TRACK_D0_W>           d0;
    ap_uint<TRACK_CHI2_W>         chi2;
    ap_uint<TRACK_BENDCHI2_W>     bendChi2;
    ap_uint<TRACK_HITMASK_W>      hitMask;
    ap_uint<TRACK_TKQUALMVA_W>    tkQualMVA;
    ap_uint<TRACK_OTHERQUALMVA_W> otherQualMVA;
    ap_uint<TRACK_SPARE_W>        spare;
};


// my DF that will be propagated to the algorithm

#define TRACKCONV_PT_W         15
#define TRACKCONV_PHI_W        12
#define TRACKCONV_ETA_W        16
#define TRACKCONV_Z0_W         12
#define TRACKCONV_D0_W         13
#define TRACKCONV_CHI2_W       4

struct track_conv_t
{
    ap_uint<TRACKCONV_PT_W>    pt;
    ap_uint<1>                 q; // 1 : +, 0 : -
    ap_int<TRACKCONV_PHI_W>    phi; 
    ap_int<TRACKCONV_ETA_W>    eta; 
    ap_uint<TRACKCONV_Z0_W>    z0;
    ap_uint<TRACKCONV_D0_W>    d0;
    ap_uint<TRACKCONV_CHI2_W>  chi2;  
};

#endif // TRACK_T_H
