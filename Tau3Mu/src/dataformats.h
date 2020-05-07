#ifndef DATAFORMATS_H
#define DATAFORMATS_H

#include "ap_int.h"
#include "ap_fixed.h"

// definition of the dataformat here
// #define PT_W        18
// #define PT_INT_W    8 //11
// #define THETA_W     16
// #define THETA_INT_W 3
// #define ETA_W       16
// #define ETA_INT_W   3
// #define PHI_W       12
// #define PHI_INT_W   3

typedef ap_ufixed<18, 8> hw_pt_t;
typedef ap_ufixed<16, 3> hw_theta_t;
typedef ap_fixed<16, 3>  hw_eta_t;
typedef ap_fixed<12, 3, AP_TRN, AP_SAT> hw_phi_t;
typedef ap_uint<50> hw_spare_t;

struct muon_t
{
    hw_pt_t  pt;
    hw_eta_t eta;
    hw_phi_t phi;
    hw_spare_t spares; // extra bits, unused, but give a reasonable size to the input word
};

typedef hw_phi_t hw_minv_t; // FIXME

// struct iso_muon_t : public muon_t
// {
//     ap_uint<1> iso;
// };

// an input/output set of muons - all independent
// struct muon_gr_t
// {
//     muon_t mu0;
//     muon_t mu1;
//     muon_t mu2;
//     muon_t mu3;
//     muon_t mu4;
//     muon_t mu5;
//     muon_t mu6;
//     muon_t mu7;
//     muon_t mu8;
//     muon_t mu9;
//     muon_t mu10;
//     muon_t mu11;
//     muon_t mu12;
// };


#endif