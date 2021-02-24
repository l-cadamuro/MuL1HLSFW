#ifndef DATAFORMATS_V2_H
#define DATAFORMATS_V2_H

#include "ap_int.h"
#include "ap_fixed.h"

// #define USE_REL_ISO 1 // 0 : absolute iso, 1 : relative iso

///////////////////////////////////////////////////////////////
////// basic data types definitions

typedef ap_uint<14> hw_pt_t;   // LSB : 1/2^5 GeV
typedef ap_int<13>  hw_eta_t;  // LSB : 2pi/2^13 rad
typedef ap_int<13>  hw_phi_t;  // LSB : 2pi/2^13 rad
typedef ap_int<12>  hw_z0_t;   // LSB : 60/2^12 cm
typedef ap_uint<1>  hw_q_t;    // LSB : 1

// unused in algo, but they build a "full" dataformat
typedef ap_int<13>   hw_d0_t;     // d0 max: 15.4 cm
typedef ap_uint<62>  hw_spare_t;  // full muon track has 128 bits : 128 - 66 = 62 bits

typedef ap_uint<2>  hw_iso_t;  // two bits for isolation WPs

// ap_fixed <total_bits, integer_part_bits, rounding mode, overflow mode, >
struct muon_t {
    hw_pt_t     pt;
    hw_eta_t    eta;
    hw_phi_t    phi;
    hw_z0_t     z0;
    hw_q_t      q;
    hw_d0_t     d0;
    hw_spare_t  spare;
};

struct isomuon_t {
    muon_t    mu;
    hw_iso_t  isoflags;
};

struct track_t {

    hw_pt_t     pt;
    hw_eta_t    eta;
    hw_phi_t    phi;
    hw_z0_t     z0;
    hw_q_t      q;
};

///////////////////////////////////////////////////////////////
////// formats specific to the isolation algorithms

// threshold is typically integer in GeV and of a few GeV -> minimise the number of bits needed
// 14 bits -> 512 GeV ==> 14 - 6  = 8 bits can express 512 / 2^6 = 8 GeV

// NOTE : less resources can be used by summing reduced precision pT values, i.e. accum += (pt >> N)
// typedef ap_uint<8, AP_TRN, AP_SAT> iso_accum_t;
typedef ap_ufixed<8, 8, AP_TRN, AP_SAT> iso_accum_t;

typedef ap_uint<hw_eta_t::width>   deta_t ;
typedef ap_uint<hw_phi_t::width-1> dphi_t ; // as dphi rolls over the boundary, need 1/2 the range to express it -> 1 bit less
typedef ap_uint<hw_z0_t::width>    dz0_t ;

///////////////////////////////////////////////////////////////
////// format of the transmission (N inputs, TMT, ...)

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

struct track_data_9_t {
    track_t trk_0;
    track_t trk_1;
    track_t trk_2;
    track_t trk_3;
    track_t trk_4;
    track_t trk_5;
    track_t trk_6;
    track_t trk_7;
    track_t trk_8;
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

struct isomuon_data_t {
    isomuon_t mu_0;
    isomuon_t mu_1;
    isomuon_t mu_2;
    isomuon_t mu_3;
    isomuon_t mu_4;
    isomuon_t mu_5;
    isomuon_t mu_6;
    isomuon_t mu_7;
    isomuon_t mu_8;
    isomuon_t mu_9;
    isomuon_t mu_10;
    isomuon_t mu_11;
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
////// constants used to drive the isolation algorithm
////// Note : must be aligned with the df (LSB)

#define c_iso_dangle_max 260 // <  , 260 x 2pi/2^13 = 0.2 rad
#define c_iso_dz_max     68  // <  , 68  x 60/2^12  = 1 cm
#define c_iso_pt_min     96  // >= , 96  x 1/2^5    = 3 GeV

#endif