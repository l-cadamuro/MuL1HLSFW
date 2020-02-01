#include "invariant_mass.h"

ap_uint<c_minv_w> invariant_mass(ap_uint<c_minv_pt_w> pt1, ap_int<c_minv_phi_w> phi1, ap_int<c_minv_eta_w> eta1, ap_uint<c_minv_pt_w> pt2, ap_int<c_minv_phi_w> phi2, ap_int<c_minv_eta_w> eta2)
{
    #pragma HLS pipeline II=1

    // compute the pt product
    ap_uint<2*c_minv_pt_w + 1> ptprod = 2*pt1*pt2; // is vivado intelligent enough to do 
    
    // compute abs(dphi)
    ap_uint < c_minv_phi_w - 1 > absdphi; // using parity to use 1 bit less 
    ap_uint < c_minv_phi_w >     dphi;
    dphi = phi2-phi1;
    if (dphi < 0)
        absdphi = -1*dphi;
    else
        absdphi = dphi;
    
    // compute the theta angular part
    

}