#ifndef INVARIANT_MASS_H
#define INVARIANT_MASS_H

#include "constants.h"

ap_uint<c_minv_w> invariant_mass(ap_uint<c_minv_pt_w> pt1, ap_int<c_minv_phi_w> phi1, ap_int<c_minv_eta_w> eta1, ap_uint<c_minv_pt_w> pt2, ap_int<c_minv_phi_w> phi2, ap_int<c_minv_eta_w> eta2);

#endif