#include "isolation.h"
#include "constants.h"

#ifndef __SYNTHESIS__
#include <iostream>
#endif

ap_uint<1> isolation (muon_t mu_in, track_conv_input trks_in, ap_uint<1> last_in)
{
    #pragma HLS pipeline II=1

    //test_out = trk_in.trk0.qOverR + mu_in.mu0.pt;
    ap_uint<1> out;

    // the accumlators - there is one for every track incoming
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_0  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_1  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_2  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_3  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_4  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_5  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_6  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_7  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_8  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_9  = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_10 = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_11 = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_12 = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_13 = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_14 = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_15 = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_16 = 0;
    static ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> accum_17 = 0;

    #ifndef __SYNTHESIS__
    std::cout << " ~~~~ DEBUG : 0 : accum0 = " << accum_0.to_uint() << " += " << trks_in.trk0.pt.to_uint() << std::endl;
    #endif

    // the conditions to add a track

    ap_uint<TRACKCONV_ETA_W> deta0 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk0.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi0 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk0.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum0;
    if (deta0 < c_deta_iso && dphi0 < c_dphi_iso)
        psum0 = trks_in.trk0.pt;
    else
        psum0 = 0;

    ap_uint<TRACKCONV_ETA_W> deta1 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk1.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi1 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk1.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum1;
    if (deta1 < c_deta_iso && dphi1 < c_dphi_iso)
        psum1 = trks_in.trk1.pt;
    else
        psum1 = 0;

    ap_uint<TRACKCONV_ETA_W> deta2 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk2.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi2 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk2.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum2;
    if (deta2 < c_deta_iso && dphi2 < c_dphi_iso)
        psum2 = trks_in.trk2.pt;
    else
        psum2 = 0;

    ap_uint<TRACKCONV_ETA_W> deta3 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk3.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi3 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk3.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum3;
    if (deta3 < c_deta_iso && dphi3 < c_dphi_iso)
        psum3 = trks_in.trk3.pt;
    else
        psum3 = 0;

    ap_uint<TRACKCONV_ETA_W> deta4 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk4.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi4 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk4.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum4;
    if (deta4 < c_deta_iso && dphi4 < c_dphi_iso)
        psum4 = trks_in.trk4.pt;
    else
        psum4 = 0;

    ap_uint<TRACKCONV_ETA_W> deta5 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk5.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi5 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk5.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum5;
    if (deta5 < c_deta_iso && dphi5 < c_dphi_iso)
        psum5 = trks_in.trk5.pt;
    else
        psum5 = 0;

    ap_uint<TRACKCONV_ETA_W> deta6 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk6.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi6 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk6.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum6;
    if (deta6 < c_deta_iso && dphi6 < c_dphi_iso)
        psum6 = trks_in.trk6.pt;
    else
        psum6 = 0;

    ap_uint<TRACKCONV_ETA_W> deta7 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk7.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi7 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk7.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum7;
    if (deta7 < c_deta_iso && dphi7 < c_dphi_iso)
        psum7 = trks_in.trk7.pt;
    else
        psum7 = 0;

    ap_uint<TRACKCONV_ETA_W> deta8 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk8.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi8 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk8.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum8;
    if (deta8 < c_deta_iso && dphi8 < c_dphi_iso)
        psum8 = trks_in.trk8.pt;
    else
        psum8 = 0;

    ap_uint<TRACKCONV_ETA_W> deta9 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk9.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi9 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk9.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum9;
    if (deta9 < c_deta_iso && dphi9 < c_dphi_iso)
        psum9 = trks_in.trk9.pt;
    else
        psum9 = 0;

    ap_uint<TRACKCONV_ETA_W> deta10 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk10.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi10 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk10.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum10;
    if (deta10 < c_deta_iso && dphi10 < c_dphi_iso)
        psum10 = trks_in.trk10.pt;
    else
        psum10 = 0;

    ap_uint<TRACKCONV_ETA_W> deta11 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk11.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi11 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk11.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum11;
    if (deta11 < c_deta_iso && dphi11 < c_dphi_iso)
        psum11 = trks_in.trk11.pt;
    else
        psum11 = 0;

    ap_uint<TRACKCONV_ETA_W> deta12 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk12.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi12 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk12.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum12;
    if (deta12 < c_deta_iso && dphi12 < c_dphi_iso)
        psum12 = trks_in.trk12.pt;
    else
        psum12 = 0;

    ap_uint<TRACKCONV_ETA_W> deta13 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk13.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi13 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk13.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum13;
    if (deta13 < c_deta_iso && dphi13 < c_dphi_iso)
        psum13 = trks_in.trk13.pt;
    else
        psum13 = 0;

    ap_uint<TRACKCONV_ETA_W> deta14 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk14.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi14 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk14.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum14;
    if (deta14 < c_deta_iso && dphi14 < c_dphi_iso)
        psum14 = trks_in.trk14.pt;
    else
        psum14 = 0;

    ap_uint<TRACKCONV_ETA_W> deta15 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk15.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi15 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk15.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum15;
    if (deta15 < c_deta_iso && dphi15 < c_dphi_iso)
        psum15 = trks_in.trk15.pt;
    else
        psum15 = 0;

    ap_uint<TRACKCONV_ETA_W> deta16 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk16.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi16 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk16.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum16;
    if (deta16 < c_deta_iso && dphi16 < c_dphi_iso)
        psum16 = trks_in.trk16.pt;
    else
        psum16 = 0;

    ap_uint<TRACKCONV_ETA_W> deta17 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk17.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi17 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk17.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum17;
    if (deta17 < c_deta_iso && dphi17 < c_dphi_iso)
        psum17 = trks_in.trk17.pt;
    else
        psum17 = 0;


    accum_0 += psum0;
    accum_1 += psum1;
    accum_2 += psum2;
    accum_3 += psum3;
    accum_4 += psum4;
    accum_5 += psum5;
    accum_6 += psum6;
    accum_7 += psum7;
    accum_8 += psum8;
    accum_9 += psum9;
    accum_10 += psum10;
    accum_11 += psum11;
    accum_12 += psum12;
    accum_13 += psum13;
    accum_14 += psum14;
    accum_15 += psum15;
    accum_16 += psum16;
    accum_17 += psum17;

    #ifndef __SYNTHESIS__
    std::cout << " ~~~~~~~~> accum0 now is = " << accum_0.to_uint() << std::endl;
    #endif

    if (last_in == 0)
        return ap_uint<1>(0); // an empty result

    else{
        // compute the final energy sum
        ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> ptsum;
        ptsum = accum_0 +
                accum_1 +
                accum_2 +
                accum_3 +
                accum_4 +
                accum_5 +
                accum_6 +
                accum_7 +
                accum_8 +
                accum_9 +
                accum_10 +
                accum_11 +
                accum_12 +
                accum_13 +
                accum_14 +
                accum_15 +
                accum_16 +
                accum_17 ;

        if (ptsum >= c_iso_thr)
            out = 0;
        else
            out = 1;

        return out;
    }
}