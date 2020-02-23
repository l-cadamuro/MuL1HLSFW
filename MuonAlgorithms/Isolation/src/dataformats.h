#ifndef DATAFORMATS_H
#define DATAFORMATS_H

#include "ap_int.h"
#include "ap_fixed.h"

///////////////////////////////////////////////////////////////
////// format of the individual data

// // definition of the dataformat here
// #define PT_W        18
// #define PT_INT_W    8 //11

// // #define THETA_W     16
// // #define THETA_INT_W 3

// #define ETA_W       16
// #define ETA_INT_W   3

// #define PHI_W       12
// #define PHI_INT_W   3

// // z ranges approx from -40 to +40 [cm] -> 7 bits to express it in -64, 64
// // resolution on z is >= 0.1 cm -> 4 bits for the decimal part are e ough
// #define Z0_W        11
// #define Z0_INT_W    7  

// typedef ap_ufixed<PT_W, PT_INT_W>                    hw_pt_t;
// // typedef ap_ufixed<THETA_W, THETA_INT_W>              hw_theta_t;
// typedef ap_fixed<THETA_W, THETA_INT_W>               hw_eta_t;
// // typedef ap_fixed<PHI_W, PHI_INT_W, AP_TRN, AP_SAT>   hw_phi_t;
// typedef ap_fixed<PHI_W, PHI_INT_W>                   hw_phi_t; // no truncation for phi makes dphi roll over the pi boundary
// typedef ap_fixed<Z0_W,  Z0_INT_W,  AP_TRN, AP_SAT>   hw_z0_t;

// DEBUG : back to 15 bits for pt
// typedef ap_ufixed <15, 8, AP_TRN, AP_SAT>  hw_pt_t;
// typedef ap_ufixed <15, 8>  hw_pt_t;
// typedef ap_fixed  <16, 3>  hw_eta_t;
// typedef ap_fixed  <12, 3>  hw_phi_t;
// typedef ap_fixed  <11, 7>  hw_z0_t;

typedef ap_ufixed <18, 8, AP_TRN, AP_SAT>  hw_pt_t;
typedef ap_fixed  <16, 3, AP_TRN, AP_SAT>  hw_eta_t;
typedef ap_fixed  <12, 3>                  hw_phi_t; // no truncation for phi makes dphi roll over the pi boundary
typedef ap_fixed  <11, 7, AP_TRN, AP_SAT>  hw_z0_t;

// ap_fixed <total_bits, integer_part_bits, rounding mode, overflow mode, >
struct muon_t {

    hw_pt_t     pt;
    // hw_theta_t  theta;
    hw_eta_t    eta;
    hw_phi_t    phi;
    hw_z0_t     z0;
};

struct track_t {

    hw_pt_t     pt;
    // hw_theta_t  theta;
    hw_eta_t    eta;
    hw_phi_t    phi;
    hw_z0_t     z0;
};

///////////////////////////////////////////////////////////////
////// format of the transmission (N inputs, TMT, ...)

#define N_TRK_LINKS    18 // total number of input links (each gives 1 trk / clk)
// #define N_TRK_PER_LINK 18 // total number of track serially trasmitted on each link
#define N_MUON 12

struct track_data_t {
    track_t trk_0;
    track_t trk_1;
    track_t trk_2;
    track_t trk_3;
    track_t trk_4;
    track_t trk_5;
    track_t trk_6;
    track_t trk_7;
    track_t trk_8;
    track_t trk_9;
    track_t trk_10;
    track_t trk_11;
    track_t trk_12;
    track_t trk_13;
    track_t trk_14;
    track_t trk_15;
    track_t trk_16;
    track_t trk_17;
};

struct muon_data_t {
    muon_t mu_0;
    muon_t mu_1;
    muon_t mu_2;
    muon_t mu_3;
    muon_t mu_4;
    muon_t mu_5;
    muon_t mu_6;
    muon_t mu_7;
    muon_t mu_8;
    muon_t mu_9;
    muon_t mu_10;
    muon_t mu_11;
};

struct muon_isodata_t {
    ap_uint<1> isomu_0;
    ap_uint<1> isomu_1;
    ap_uint<1> isomu_2;
    ap_uint<1> isomu_3;
    ap_uint<1> isomu_4;
    ap_uint<1> isomu_5;
    ap_uint<1> isomu_6;
    ap_uint<1> isomu_7;
    ap_uint<1> isomu_8;
    ap_uint<1> isomu_9;
    ap_uint<1> isomu_10;
    ap_uint<1> isomu_11;
};
///////////////////////////////////////////////////////////////
////// formats and constants specific to the isolation algorithms
////// FIXME : move to another "constants.h" header?

typedef ap_ufixed<5, 5, AP_TRN, AP_SAT> iso_accum_t; // threshold is typically integer in GeV and of a few GeV -> 5 bits, all integers

// #define c_iso_sumpt_thr  12  // < 
#define c_iso_sumpt_thr  1  // purely for debug
#define c_iso_dangle_max 0.2 // < 
#define c_iso_dz_max     1.0 // < 
#define c_iso_pt_min     3.0 // >= 

typedef ap_ufixed<hw_eta_t::width-1, hw_eta_t::iwidth-1> deta_t ;
typedef ap_ufixed<hw_phi_t::width-1, hw_phi_t::iwidth-1> dphi_t ;
typedef  ap_ufixed<hw_z0_t::width-1,  hw_z0_t::iwidth-1> dz0_t ;

#endif