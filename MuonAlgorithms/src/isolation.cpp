#include "isolation.h"
#include "constants.h"

#ifndef __SYNTHESIS__
#include <iostream>
#endif

#define DEBUG_ISOLATION false

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
    if (DEBUG_ISOLATION) std::cout << " ~~~~ DEBUG : 0 : accum14 = " << accum_14.to_uint() << " += " << trks_in.trk14.pt.to_uint() << " (if matched)" << std::endl;
    #endif

    // the conditions to add a track

    ap_uint<TRACKCONV_ETA_W> deta0 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk0.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi0 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk0.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum0;
    if (deta0 < c_deta_iso && dphi0 < c_dphi_iso && deta0 > 1 && dphi0 > 1 && trks_in.trk0.pt > c_trk_iso_thresh)
        psum0 = trks_in.trk0.pt;
    else
        psum0 = 0;

    ap_uint<TRACKCONV_ETA_W> deta1 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk1.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi1 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk1.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum1;
    if (deta1 < c_deta_iso && dphi1 < c_dphi_iso && deta1 > 1 && dphi1 > 1 && trks_in.trk1.pt > c_trk_iso_thresh)
        psum1 = trks_in.trk1.pt;
    else
        psum1 = 0;

    ap_uint<TRACKCONV_ETA_W> deta2 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk2.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi2 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk2.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum2;
    if (deta2 < c_deta_iso && dphi2 < c_dphi_iso && deta2 > 1 && dphi2 > 1 && trks_in.trk2.pt > c_trk_iso_thresh)
        psum2 = trks_in.trk2.pt;
    else
        psum2 = 0;

    ap_uint<TRACKCONV_ETA_W> deta3 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk3.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi3 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk3.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum3;
    if (deta3 < c_deta_iso && dphi3 < c_dphi_iso && deta3 > 1 && dphi3 > 1 && trks_in.trk3.pt > c_trk_iso_thresh)
        psum3 = trks_in.trk3.pt;
    else
        psum3 = 0;

    ap_uint<TRACKCONV_ETA_W> deta4 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk4.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi4 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk4.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum4;
    if (deta4 < c_deta_iso && dphi4 < c_dphi_iso && deta4 > 1 && dphi4 > 1 && trks_in.trk4.pt > c_trk_iso_thresh)
        psum4 = trks_in.trk4.pt;
    else
        psum4 = 0;

    ap_uint<TRACKCONV_ETA_W> deta5 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk5.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi5 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk5.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum5;
    if (deta5 < c_deta_iso && dphi5 < c_dphi_iso && deta5 > 1 && dphi5 > 1 && trks_in.trk5.pt > c_trk_iso_thresh)
        psum5 = trks_in.trk5.pt;
    else
        psum5 = 0;

    ap_uint<TRACKCONV_ETA_W> deta6 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk6.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi6 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk6.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum6;
    if (deta6 < c_deta_iso && dphi6 < c_dphi_iso && deta6 > 1 && dphi6 > 1 && trks_in.trk6.pt > c_trk_iso_thresh)
        psum6 = trks_in.trk6.pt;
    else
        psum6 = 0;

    ap_uint<TRACKCONV_ETA_W> deta7 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk7.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi7 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk7.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum7;
    if (deta7 < c_deta_iso && dphi7 < c_dphi_iso && deta7 > 1 && dphi7 > 1 && trks_in.trk7.pt > c_trk_iso_thresh)
        psum7 = trks_in.trk7.pt;
    else
        psum7 = 0;

    ap_uint<TRACKCONV_ETA_W> deta8 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk8.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi8 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk8.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum8;
    if (deta8 < c_deta_iso && dphi8 < c_dphi_iso && deta8 > 1 && dphi8 > 1 && trks_in.trk8.pt > c_trk_iso_thresh)
        psum8 = trks_in.trk8.pt;
    else
        psum8 = 0;

    ap_uint<TRACKCONV_ETA_W> deta9 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk9.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi9 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk9.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum9;
    if (deta9 < c_deta_iso && dphi9 < c_dphi_iso && deta9 > 1 && dphi9 > 1 && trks_in.trk9.pt > c_trk_iso_thresh)
        psum9 = trks_in.trk9.pt;
    else
        psum9 = 0;

    ap_uint<TRACKCONV_ETA_W> deta10 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk10.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi10 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk10.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum10;
    if (deta10 < c_deta_iso && dphi10 < c_dphi_iso && deta10 > 1 && dphi10 > 1 && trks_in.trk10.pt > c_trk_iso_thresh)
        psum10 = trks_in.trk10.pt;
    else
        psum10 = 0;

    ap_uint<TRACKCONV_ETA_W> deta11 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk11.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi11 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk11.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum11;
    if (deta11 < c_deta_iso && dphi11 < c_dphi_iso && deta11 > 1 && dphi11 > 1 && trks_in.trk11.pt > c_trk_iso_thresh)
        psum11 = trks_in.trk11.pt;
    else
        psum11 = 0;

    ap_uint<TRACKCONV_ETA_W> deta12 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk12.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi12 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk12.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum12;
    if (deta12 < c_deta_iso && dphi12 < c_dphi_iso && deta12 > 1 && dphi12 > 1 && trks_in.trk12.pt > c_trk_iso_thresh)
        psum12 = trks_in.trk12.pt;
    else
        psum12 = 0;

    ap_uint<TRACKCONV_ETA_W> deta13 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk13.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi13 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk13.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum13;
    if (deta13 < c_deta_iso && dphi13 < c_dphi_iso && deta13 > 1 && dphi13 > 1 && trks_in.trk13.pt > c_trk_iso_thresh)
        psum13 = trks_in.trk13.pt;
    else
        psum13 = 0;

    ap_uint<TRACKCONV_ETA_W> deta14 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk14.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi14 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk14.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum14;
    if (deta14 < c_deta_iso && dphi14 < c_dphi_iso && deta14 > 1 && dphi14 > 1 && trks_in.trk14.pt > c_trk_iso_thresh)
        psum14 = trks_in.trk14.pt;
    else
        psum14 = 0;

    ap_uint<TRACKCONV_ETA_W> deta15 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk15.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi15 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk15.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum15;
    if (deta15 < c_deta_iso && dphi15 < c_dphi_iso && deta15 > 1 && dphi15 > 1 && trks_in.trk15.pt > c_trk_iso_thresh)
        psum15 = trks_in.trk15.pt;
    else
        psum15 = 0;

    ap_uint<TRACKCONV_ETA_W> deta16 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk16.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi16 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk16.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum16;
    if (deta16 < c_deta_iso && dphi16 < c_dphi_iso && deta16 > 1 && dphi16 > 1 && trks_in.trk16.pt > c_trk_iso_thresh)
        psum16 = trks_in.trk16.pt;
    else
        psum16 = 0;

    ap_uint<TRACKCONV_ETA_W> deta17 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk17.eta, mu_in.eta);
    ap_uint<TRACKCONV_PHI_W> dphi17 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk17.phi, mu_in.phi);
    ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum17;
    if (deta17 < c_deta_iso && dphi17 < c_dphi_iso && deta17 > 1 && dphi17 > 1 && trks_in.trk17.pt > c_trk_iso_thresh)
        psum17 = trks_in.trk17.pt;
    else
        psum17 = 0;

    #ifndef __SYNTHESIS__
    if (DEBUG_ISOLATION) std::cout << " ~~~~~~~~> partial sum 14 is = " << psum14.to_uint() << std::endl;
    #endif

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
    if (DEBUG_ISOLATION) std::cout << " ~~~~~~~~> accum14 now is = " << accum_14.to_uint() << std::endl;
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

        #ifndef __SYNTHESIS__
        if (DEBUG_ISOLATION){
            std::cout << " ~~~~~~~~|||> the 18 partial sums are = " << std::endl;
            std::cout << "     " << std::endl;
            std::cout << "     " << accum_0.to_uint() << std::endl;
            std::cout << "     " << accum_1.to_uint() << std::endl;
            std::cout << "     " << accum_2.to_uint() << std::endl;
            std::cout << "     " << accum_3.to_uint() << std::endl;
            std::cout << "     " << accum_4.to_uint() << std::endl;
            std::cout << "     " << accum_5.to_uint() << std::endl;
            std::cout << "     " << accum_6.to_uint() << std::endl;
            std::cout << "     " << accum_7.to_uint() << std::endl;
            std::cout << "     " << accum_8.to_uint() << std::endl;
            std::cout << "     " << accum_9.to_uint() << std::endl;
            std::cout << "     " << accum_10.to_uint() << std::endl;
            std::cout << "     " << accum_11.to_uint() << std::endl;
            std::cout << "     " << accum_12.to_uint() << std::endl;
            std::cout << "     " << accum_13.to_uint() << std::endl;
            std::cout << "     " << accum_14.to_uint() << std::endl;
            std::cout << "     " << accum_15.to_uint() << std::endl;
            std::cout << "     " << accum_16.to_uint() << std::endl;
            std::cout << "     " << accum_17.to_uint() << std::endl;
            std::cout << " ~~~~~||||||||> the total sum is = " << ptsum.to_uint() << std::endl;
            std::cout << " ~~~~~||||||||> compared to a mu pt of = " << mu_in.pt.to_uint() << std::endl;
        }
        #endif

        if (ptsum >= c_iso_thr)
            out = 0;
        else
            out = 1;

        #ifndef __SYNTHESIS__
        if (DEBUG_ISOLATION)  std::cout << " ~~~~~~~~> isolation flag is = " << out.to_uint() << std::endl;
        #endif

        // reset the accumulators
        accum_0  = 0;
        accum_1  = 0;
        accum_2  = 0;
        accum_3  = 0;
        accum_4  = 0;
        accum_5  = 0;
        accum_6  = 0;
        accum_7  = 0;
        accum_8  = 0;
        accum_9  = 0;
        accum_10 = 0;
        accum_11 = 0;
        accum_12 = 0;
        accum_13 = 0;
        accum_14 = 0;
        accum_15 = 0;
        accum_16 = 0;
        accum_17 = 0;

        return out;
    }
}

ap_uint<1> isolation_class (muon_t mu_in, track_conv_input trks_in, ap_uint<1> last_in)
{
    #pragma HLS pipeline II=1

    static iso_calculator ic;
    return ic.isolation(mu_in, trks_in, last_in);
}

iso_muon_input isolation_class_allmu (muon_input mus_in, track_conv_input trks_in, ap_uint<1> last_in)
{
    #pragma HLS pipeline II=1

    iso_muon_input isomus;

    static iso_calculator ic0;
    static iso_calculator ic1;
    static iso_calculator ic2;
    static iso_calculator ic3;
    static iso_calculator ic4;
    static iso_calculator ic5;
    static iso_calculator ic6;
    static iso_calculator ic7;
    static iso_calculator ic8;
    static iso_calculator ic9;
    static iso_calculator ic10;
    static iso_calculator ic11;
    static iso_calculator ic12;
    static iso_calculator ic13;
    static iso_calculator ic14;
    static iso_calculator ic15;
    static iso_calculator ic16;
    static iso_calculator ic17;

    ap_uint<1> iso0  = ic0.isolation(mus_in.mu0,  trks_in, last_in);
    ap_uint<1> iso1  = ic1.isolation(mus_in.mu1,  trks_in, last_in);
    ap_uint<1> iso2  = ic2.isolation(mus_in.mu2,  trks_in, last_in);
    ap_uint<1> iso3  = ic3.isolation(mus_in.mu3,  trks_in, last_in);
    ap_uint<1> iso4  = ic4.isolation(mus_in.mu4,  trks_in, last_in);
    ap_uint<1> iso5  = ic5.isolation(mus_in.mu5,  trks_in, last_in);
    ap_uint<1> iso6  = ic6.isolation(mus_in.mu6,  trks_in, last_in);
    ap_uint<1> iso7  = ic7.isolation(mus_in.mu7,  trks_in, last_in);
    ap_uint<1> iso8  = ic8.isolation(mus_in.mu8,  trks_in, last_in);
    ap_uint<1> iso9  = ic9.isolation(mus_in.mu9,  trks_in, last_in);
    ap_uint<1> iso10 = ic10.isolation(mus_in.mu10, trks_in, last_in);
    ap_uint<1> iso11 = ic11.isolation(mus_in.mu11, trks_in, last_in);
    ap_uint<1> iso12 = ic12.isolation(mus_in.mu12, trks_in, last_in);
    ap_uint<1> iso13 = ic13.isolation(mus_in.mu13, trks_in, last_in);
    ap_uint<1> iso14 = ic14.isolation(mus_in.mu14, trks_in, last_in);
    ap_uint<1> iso15 = ic15.isolation(mus_in.mu15, trks_in, last_in);
    ap_uint<1> iso16 = ic16.isolation(mus_in.mu16, trks_in, last_in);
    ap_uint<1> iso17 = ic17.isolation(mus_in.mu17, trks_in, last_in);

    isomus.mu0.iso = iso0;
    isomus.mu1.iso = iso1;
    isomus.mu2.iso = iso2;
    isomus.mu3.iso = iso3;
    isomus.mu4.iso = iso4;
    isomus.mu5.iso = iso5;
    isomus.mu6.iso = iso6;
    isomus.mu7.iso = iso7;
    isomus.mu8.iso = iso8;
    isomus.mu9.iso = iso9;
    isomus.mu10.iso = iso10;
    isomus.mu11.iso = iso11;
    isomus.mu12.iso = iso12;
    isomus.mu13.iso = iso13;
    isomus.mu14.iso = iso14;
    isomus.mu15.iso = iso15;
    isomus.mu16.iso = iso16;
    isomus.mu17.iso = iso17;

    return isomus;
}
