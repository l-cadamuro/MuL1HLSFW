#ifndef INV_MASS_DF_H
#define INV_MASS_DF_H

#include "ap_int.h"
#include "ap_fixed.h"

// definition of the dataformat here
#define PT_W        18
#define PT_INT_W    8 //11
#define THETA_W     16
#define THETA_INT_W 3
#define ETA_W       16
#define ETA_INT_W   3
#define PHI_W       12
#define PHI_INT_W   3

typedef ap_ufixed<PT_W, PT_INT_W>                    hw_pt_t;
typedef ap_ufixed<THETA_W, THETA_INT_W>              hw_theta_t;
typedef ap_fixed<THETA_W, THETA_INT_W>               hw_eta_t;
typedef ap_fixed<PHI_W, PHI_INT_W, AP_TRN, AP_SAT>   hw_phi_t;

/////////////////
// LUT and object definition
#define COS_PHI_ADDR_W     11 
#define COS_PHI_ADDR_W_INT 3 // max dphi is 2pi < 2**3 = 8
#define COS_PHI_OUT_W      18

#define COSH_ETA_ADDR_W      11
#define COSH_ETA_ADDR_W_INT  3  // max eta is 2.8 -> max deta is 5.6 < 2**3 = 8
#define COSH_ETA_OUT_W       18
#define COSH_ETA_OUT_W_INT   8  // because cosh (5.6) = 135


// ap_fixed <total_bits, integer_part_bits, rounding mode, overflow mode, >
struct p3_polar_tr_digi {
    
    hw_pt_t     pt;   
    hw_theta_t  theta;
    hw_eta_t    eta;
    hw_phi_t    phi;
};

typedef p3_polar_tr_digi hw_p3_t  ;

// and precision of the output of the inv mass calculation here
// that actually returns m^2/2
typedef ap_ufixed<18, 9, AP_TRN, AP_SAT>   hw_minv2over2_t;

#endif //INV_MASS_DF_H