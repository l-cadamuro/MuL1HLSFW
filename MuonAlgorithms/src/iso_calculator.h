#ifndef ISO_CALCULATOR_H
#define ISO_CALCULATOR_H

#include "ap_fixed.h"
#include "ap_int.h"

#include "muon_t.h"
#include "track_t.h"
#include "track_input.h"
#include "muon_input.h"

#include "constants.h"

// compute the isolation of one muon (mu) with N tracks streaming in
// returns true if isolated, else false

class iso_calculator
{
    public:
        iso_calculator();
        ~iso_calculator();
        ap_uint<1> isolation (muon_t mu_in, track_conv_input trks_in, ap_uint<1> last_in);

        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_0  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_1  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_2  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_3  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_4  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_5  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_6  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_7  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_8  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_9  ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_10 ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_11 ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_12 ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_13 ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_14 ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_15 ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_16 ; // = 0;
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_17 ; // = 0;        

        // void reset_accums();
};

// does abs(v2 - v1)
template <class tin, class tout>
tout abs_delta_isocalc (tin v1, tin v2)
{
    tin  dbuf = v2-v1;
    tout dout;
    if (dbuf < 0)
        dout = -1*dbuf;
    else
        dout = dbuf;
    
    return dout;
}

#endif