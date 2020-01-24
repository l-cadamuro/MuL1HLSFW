#include "test_algo.h"
#include "df_conversions.h"
#include "constants.h"

#ifndef __SYNTHESIS__
#include <iostream>
#endif

// some macros to reduce verbosity
#define COPY_MU(FROM, TO) \
    TO .  pt     = FROM .  pt ; \
    TO .  theta  = FROM .  theta ; \
    TO .  phi    = FROM .  phi ; \


iso_muons_out test_algo(muon_input mu_in, track_input trk_in, ap_uint<1> last_in)
{
    #pragma HLS pipeline II=1

    //test_out = trk_in.trk0.qOverR + mu_in.mu0.pt;
    iso_muons_out out;

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

    // // every track has pt in 15 bits, and I expect up to 18 tracks for each fiber -> the sum can be expressed on 20 bits
    // static ap_uint<20> accum_0 = 0;
    // static ap_uint<20> accum_1 = 0;
    // static ap_uint<20> accum_2 = 0;
    // static ap_uint<20> accum_3 = 0;
    // static ap_uint<20> accum_4 = 0;
    // static ap_uint<20> accum_5 = 0;
    // static ap_uint<20> accum_6 = 0;
    // static ap_uint<20> accum_7 = 0;
    // static ap_uint<20> accum_8 = 0;
    // static ap_uint<20> accum_9 = 0;
    // static ap_uint<20> accum_10 = 0;
    // static ap_uint<20> accum_11 = 0;
    // static ap_uint<20> accum_12 = 0;
    // static ap_uint<20> accum_13 = 0;
    // static ap_uint<20> accum_14 = 0;
    // static ap_uint<20> accum_15 = 0;
    // static ap_uint<20> accum_16 = 0;
    // static ap_uint<20> accum_17 = 0;


    #ifndef __SYNTHESIS__
    std::cout << " ~~~~ DEBUG : 0 : accum0 = " << accum_0.to_uint() << " += " << trk_in.trk0.qOverR.to_uint() << std::endl;
    #endif

    // first convert the muons to the new df with pt etc..
    track_conv_t trkconv0;
    track_conv_t trkconv1;
    track_conv_t trkconv2;
    track_conv_t trkconv3;
    track_conv_t trkconv4;
    track_conv_t trkconv5;
    track_conv_t trkconv6;
    track_conv_t trkconv7;
    track_conv_t trkconv8;
    track_conv_t trkconv9;
    track_conv_t trkconv10;
    track_conv_t trkconv11;
    track_conv_t trkconv12;
    track_conv_t trkconv13;
    track_conv_t trkconv14;
    track_conv_t trkconv15;
    track_conv_t trkconv16;
    track_conv_t trkconv17;

    track_to_trackconv (trk_in.trk0, trkconv0);
    track_to_trackconv (trk_in.trk1, trkconv1);
    track_to_trackconv (trk_in.trk2, trkconv2);
    track_to_trackconv (trk_in.trk3, trkconv3);
    track_to_trackconv (trk_in.trk4, trkconv4);
    track_to_trackconv (trk_in.trk5, trkconv5);
    track_to_trackconv (trk_in.trk6, trkconv6);
    track_to_trackconv (trk_in.trk7, trkconv7);
    track_to_trackconv (trk_in.trk8, trkconv8);
    track_to_trackconv (trk_in.trk9, trkconv9);
    track_to_trackconv (trk_in.trk10, trkconv10);
    track_to_trackconv (trk_in.trk11, trkconv11);
    track_to_trackconv (trk_in.trk12, trkconv12);
    track_to_trackconv (trk_in.trk13, trkconv13);
    track_to_trackconv (trk_in.trk14, trkconv14);
    track_to_trackconv (trk_in.trk15, trkconv15);
    track_to_trackconv (trk_in.trk16, trkconv16);
    track_to_trackconv (trk_in.trk17, trkconv17);

    // accum_0 += trk_in.trk0.qOverR;
    // accum_1 += trk_in.trk1.qOverR;
    // accum_2 += trk_in.trk2.qOverR;
    // accum_3 += trk_in.trk3.qOverR;
    // accum_4 += trk_in.trk4.qOverR;
    // accum_5 += trk_in.trk5.qOverR;
    // accum_6 += trk_in.trk6.qOverR;
    // accum_7 += trk_in.trk7.qOverR;
    // accum_8 += trk_in.trk8.qOverR;
    // accum_9 += trk_in.trk9.qOverR;
    // accum_10 += trk_in.trk10.qOverR;
    // accum_11 += trk_in.trk11.qOverR;
    // accum_12 += trk_in.trk12.qOverR;
    // accum_13 += trk_in.trk13.qOverR;
    // accum_14 += trk_in.trk14.qOverR;
    // accum_15 += trk_in.trk15.qOverR;
    // accum_16 += trk_in.trk16.qOverR;
    // accum_17 += trk_in.trk17.qOverR;

    // accum_0 += trk_in.trk0.qOverR;
    // accum_1 += trk_in.trk1.qOverR;
    // accum_2 += trk_in.trk2.qOverR;
    // accum_3 += trk_in.trk3.qOverR;
    // accum_4 += trk_in.trk4.qOverR;
    // accum_5 += trk_in.trk5.qOverR;
    // accum_6 += trk_in.trk6.qOverR;
    // accum_7 += trk_in.trk7.qOverR;
    // accum_8 += trk_in.trk8.qOverR;
    // accum_9 += trk_in.trk9.qOverR;
    // accum_10 += trk_in.trk10.qOverR;
    // accum_11 += trk_in.trk11.qOverR;
    // accum_12 += trk_in.trk12.qOverR;
    // accum_13 += trk_in.trk13.qOverR;
    // accum_14 += trk_in.trk14.qOverR;
    // accum_15 += trk_in.trk15.qOverR;
    // accum_16 += trk_in.trk16.qOverR;
    // accum_17 += trk_in.trk17.qOverR;

    // // compute delta theta and delta phi
    // ap_uint<TRACKCONV_ETA_W> dthehta0 = delta<ap_uint<TRACKCONV_ETA_W>> (trkconv0.eta, trkconv1.eta);
    // ap_uint<TRACKCONV_ETA_W> dthehta0 = delta<ap_uint<TRACKCONV_ETA_W>> (trkconv0.eta, trkconv1.eta);

    accum_0 += trkconv0.pt;
    accum_1 += trkconv1.pt;
    accum_2 += trkconv2.pt;
    accum_3 += trkconv3.pt;
    accum_4 += trkconv4.pt;
    accum_5 += trkconv5.pt;
    accum_6 += trkconv6.pt;
    accum_7 += trkconv7.pt;
    accum_8 += trkconv8.pt;
    accum_9 += trkconv9.pt;
    accum_10 += trkconv10.pt;
    accum_11 += trkconv11.pt;
    accum_12 += trkconv12.pt;
    accum_13 += trkconv13.pt;
    accum_14 += trkconv14.pt;
    accum_15 += trkconv15.pt;
    accum_16 += trkconv16.pt;
    accum_17 += trkconv17.pt;

    #ifndef __SYNTHESIS__
    std::cout << " ~~~~~~~~> accum0 now is = " << accum_0.to_uint() << std::endl;
    #endif

    if (last_in == 0)
        return iso_muons_out(); // an empty result

    else{

        // copy the properties

        COPY_MU(mu_in.mu0, out.mu0) ;
        COPY_MU(mu_in.mu1, out.mu1) ;
        COPY_MU(mu_in.mu2, out.mu2) ;
        COPY_MU(mu_in.mu3, out.mu3) ;
        COPY_MU(mu_in.mu4, out.mu4) ;
        COPY_MU(mu_in.mu5, out.mu5) ;
        COPY_MU(mu_in.mu6, out.mu6) ;
        COPY_MU(mu_in.mu7, out.mu7) ;
        COPY_MU(mu_in.mu8, out.mu8) ;
        COPY_MU(mu_in.mu9, out.mu9) ;
        COPY_MU(mu_in.mu10, out.mu10) ;
        COPY_MU(mu_in.mu11, out.mu11) ;
        COPY_MU(mu_in.mu12, out.mu12) ;
        COPY_MU(mu_in.mu13, out.mu13) ;
        COPY_MU(mu_in.mu14, out.mu14) ;
        COPY_MU(mu_in.mu15, out.mu15) ;
        COPY_MU(mu_in.mu16, out.mu16) ;
        COPY_MU(mu_in.mu17, out.mu17) ;

        // set the computed isolation
        // dummy threshold of 100, to be implemented in a LUT maybe? (vs pt mu)

        // was 1000 earlier

        out.mu0.iso  = (accum_0  > c_iso_thr ? 0 : 1);
        out.mu1.iso  = (accum_1  > c_iso_thr ? 0 : 1);
        out.mu2.iso  = (accum_2  > c_iso_thr ? 0 : 1);
        out.mu3.iso  = (accum_3  > c_iso_thr ? 0 : 1);
        out.mu4.iso  = (accum_4  > c_iso_thr ? 0 : 1);
        out.mu5.iso  = (accum_5  > c_iso_thr ? 0 : 1);
        out.mu6.iso  = (accum_6  > c_iso_thr ? 0 : 1);
        out.mu7.iso  = (accum_7  > c_iso_thr ? 0 : 1);
        out.mu8.iso  = (accum_8  > c_iso_thr ? 0 : 1);
        out.mu9.iso  = (accum_9  > c_iso_thr ? 0 : 1);
        out.mu10.iso = (accum_10 > c_iso_thr ? 0 : 1);
        out.mu11.iso = (accum_11 > c_iso_thr ? 0 : 1);
        out.mu12.iso = (accum_12 > c_iso_thr ? 0 : 1);
        out.mu13.iso = (accum_13 > c_iso_thr ? 0 : 1);
        out.mu14.iso = (accum_14 > c_iso_thr ? 0 : 1);
        out.mu15.iso = (accum_15 > c_iso_thr ? 0 : 1);
        out.mu16.iso = (accum_16 > c_iso_thr ? 0 : 1);
        out.mu17.iso = (accum_17 > c_iso_thr ? 0 : 1);

        // reset the accumulators
        accum_0 = 0;
        accum_1 = 0;
        accum_2 = 0;
        accum_3 = 0;
        accum_4 = 0;
        accum_5 = 0;
        accum_6 = 0;
        accum_7 = 0;
        accum_8 = 0;
        accum_9 = 0;
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