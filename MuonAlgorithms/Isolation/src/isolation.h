#ifndef ISOLATION_H
#define ISOLATION_H

// #include "dataformats.h"
#include "dataformats_v2.h"

#include "external_constants.h" // FIXME!


#ifndef __SYNTHESIS__
#define ISODEBUG false
#include <iostream>
using namespace std;
#endif

// checks the iso cone and returns partial sum
iso_accum_t compute_trk_iso (muon_t in_mu, track_t in_trk);

// template <int i_module>
// ap_uint<1> isolation_single(muon_t in_mu, track_t in_tracks[N_TRK_LINKS], ap_uint<1> is_last);

// void isolation(muon_t in_muons[N_MUON], track_t in_tracks[N_TRK_LINKS], ap_uint<1> is_last, ap_uint<1> iso_flags[N_MUON]);

// void isolation_allmu_flags      (muon_data_t in_muons, track_data_t  in_tracks,  ap_uint<1> is_last, iso_accum_t iso_threshold, muon_isodata_t& iso_flags);
// void isolation_allmu_9trk_flags (muon_data_t in_muons, track_data_9_t in_tracks, ap_uint<1> is_last, iso_accum_t iso_threshold, muon_isodata_t& iso_flags);

void isolation_allmu      (muon_data_t in_muons, track_data_t  in_tracks,  ap_uint<1> is_last, iso_accum_t iso_threshold_1, iso_accum_t iso_threshold_2, isomuon_data_t& iso_muons);
void isolation_allmu_9trk (muon_data_t in_muons, track_data_9_t in_tracks, ap_uint<1> is_last, iso_accum_t iso_threshold_1, iso_accum_t iso_threshold_2, isomuon_data_t& iso_muons);


// ap_uint<1> isolation_single_wrap(muon_t in_mu, track_t in_tracks[N_TRK_LINKS], ap_uint<1> is_last);
// ap_uint<1> isolation_single_muon_wrap(muon_t in_mu, track_data_t in_tracks, ap_uint<1> is_last);

// // the delta rolls at overflow (dphi)
// template <typename tret, typename tin>
// tret abs_delta_roll (tin v1, tin v2);

// // the delta saturates at overflow (deta, dz)
// template <typename tret, typename tin>
// tret abs_delta_sat (tin v1, tin v2);

// template <class tin, class tout>
// tout abs_delta (tin v1, tin v2)
// {
//     tin  dbuf = v2-v1;
//     tout dout;
//     if (dbuf < 0)
//         dout = -1*dbuf;
//     else
//         dout = dbuf;
    
//     return dout;
// }

///////////////////////////////////////////////////////////////////////////////
//// implementation of the templated function

// template <typename tret, typename tin>
// tret abs_delta_roll (tin v1, tin v2)
// {
//     #pragma HLS INLINE

//     tret result;
//     ap_fixed<tin::width+1, tin::iwidth+1> delta;    
    
//     delta = v1 - v2;
    
//     if (delta < 0)
//         delta = -1*delta;
//     else
//         delta = delta;
    
//     result = delta;
//     return result;
// }

// template <typename tret, typename tin>
// tret abs_delta_sat (tin v1, tin v2)
// {
//     #pragma HLS INLINE

//     tret result;
//     ap_fixed<tin::width+1, tin::iwidth+1, AP_TRN, AP_SAT> delta;    
    
//     delta = v1 - v2;
    
//     if (delta < 0)
//         delta = -1*delta;
//     else
//         delta = delta;
    
//     result = delta;
//     return result;
// }

// handles deta that does not rollover (deta, dz0)
template <typename tret, typename tin>
tret abs_delta_noroll (tin v1, tin v2)
{
    #pragma HLS inline

    tret result;
    if (v1 >= v2)
        result = v1 - v2;
    else
        result = v2 - v1;

    return result;
}

template <typename tret, typename tin>
tret abs_delta_roll (tin v1, tin v2)
{
    #pragma HLS inline

    tret result;
    tin  delta = v1 - v2;
    if (delta >= 0)
        result = delta;
    else
        result = -1 * delta;

    return result;
}

template <int i_module>
hw_iso_t isolation_single_muon(muon_t in_mu, track_data_t in_tracks, ap_uint<1> is_last, iso_accum_t iso_threshold_1, iso_accum_t iso_threshold_2)
{
    #pragma HLS pipeline II=1
    #pragma HLS interface ap_stable port=iso_threshold_1
    #pragma HLS interface ap_stable port=iso_threshold_2
    #pragma HLS inline
    
    // the accumulators of the energy sums
    static iso_accum_t accum_0  = 0;
    static iso_accum_t accum_1  = 0;
    static iso_accum_t accum_2  = 0;
    static iso_accum_t accum_3  = 0;
    static iso_accum_t accum_4  = 0;
    static iso_accum_t accum_5  = 0;
    static iso_accum_t accum_6  = 0;
    static iso_accum_t accum_7  = 0;
    static iso_accum_t accum_8  = 0;
    static iso_accum_t accum_9  = 0;
    static iso_accum_t accum_10 = 0;
    static iso_accum_t accum_11 = 0;
    static iso_accum_t accum_12 = 0;
    static iso_accum_t accum_13 = 0;
    static iso_accum_t accum_14 = 0;
    static iso_accum_t accum_15 = 0;
    static iso_accum_t accum_16 = 0;
    static iso_accum_t accum_17 = 0;

    #ifndef __SYNTHESIS__
    #if ISODEBUG
    cout << "--- before incrementing " << endl;
    cout << "... accum 0 : " << accum_0.to_double() << endl;
    cout << "... accum 1 : " << accum_1.to_double() << endl;
    cout << "... accum 2 : " << accum_2.to_double() << endl;
    cout << "... accum 3 : " << accum_3.to_double() << endl;
    cout << "... accum 4 : " << accum_4.to_double() << endl;
    cout << "... accum 5 : " << accum_5.to_double() << endl;
    cout << "... accum 6 : " << accum_6.to_double() << endl;
    cout << "... accum 7 : " << accum_7.to_double() << endl;
    cout << "... accum 8 : " << accum_8.to_double() << endl;
    cout << "... accum 9 : " << accum_9.to_double() << endl;
    cout << "... accum 10 : " << accum_10.to_double() << endl;
    cout << "... accum 11 : " << accum_11.to_double() << endl;
    cout << "... accum 12 : " << accum_12.to_double() << endl;
    cout << "... accum 13 : " << accum_13.to_double() << endl;
    cout << "... accum 14 : " << accum_14.to_double() << endl;
    cout << "... accum 15 : " << accum_15.to_double() << endl;
    cout << "... accum 16 : " << accum_16.to_double() << endl;
    cout << "... accum 17 : " << accum_17.to_double() << endl;
    cout << endl;
    #endif
    #endif

    hw_iso_t result;

    iso_accum_t psum_0 = compute_trk_iso (in_mu, in_tracks.trk_0);
    iso_accum_t psum_1 = compute_trk_iso (in_mu, in_tracks.trk_1);
    iso_accum_t psum_2 = compute_trk_iso (in_mu, in_tracks.trk_2);
    iso_accum_t psum_3 = compute_trk_iso (in_mu, in_tracks.trk_3);
    iso_accum_t psum_4 = compute_trk_iso (in_mu, in_tracks.trk_4);
    iso_accum_t psum_5 = compute_trk_iso (in_mu, in_tracks.trk_5);
    iso_accum_t psum_6 = compute_trk_iso (in_mu, in_tracks.trk_6);
    iso_accum_t psum_7 = compute_trk_iso (in_mu, in_tracks.trk_7);
    iso_accum_t psum_8 = compute_trk_iso (in_mu, in_tracks.trk_8);
    iso_accum_t psum_9 = compute_trk_iso (in_mu, in_tracks.trk_9);
    iso_accum_t psum_10 = compute_trk_iso (in_mu, in_tracks.trk_10);
    iso_accum_t psum_11 = compute_trk_iso (in_mu, in_tracks.trk_11);
    iso_accum_t psum_12 = compute_trk_iso (in_mu, in_tracks.trk_12);
    iso_accum_t psum_13 = compute_trk_iso (in_mu, in_tracks.trk_13);
    iso_accum_t psum_14 = compute_trk_iso (in_mu, in_tracks.trk_14);
    iso_accum_t psum_15 = compute_trk_iso (in_mu, in_tracks.trk_15);
    iso_accum_t psum_16 = compute_trk_iso (in_mu, in_tracks.trk_16);
    iso_accum_t psum_17 = compute_trk_iso (in_mu, in_tracks.trk_17);

    #ifndef __SYNTHESIS__
    #if ISODEBUG
    cout << "... psum 0 : " << psum_0.to_double() << endl;
    cout << "... psum 1 : " << psum_1.to_double() << endl;
    cout << "... psum 2 : " << psum_2.to_double() << endl;
    cout << "... psum 3 : " << psum_3.to_double() << endl;
    cout << "... psum 4 : " << psum_4.to_double() << endl;
    cout << "... psum 5 : " << psum_5.to_double() << endl;
    cout << "... psum 6 : " << psum_6.to_double() << endl;
    cout << "... psum 7 : " << psum_7.to_double() << endl;
    cout << "... psum 8 : " << psum_8.to_double() << endl;
    cout << "... psum 9 : " << psum_9.to_double() << endl;
    cout << "... psum 10 : " << psum_10.to_double() << endl;
    cout << "... psum 11 : " << psum_11.to_double() << endl;
    cout << "... psum 12 : " << psum_12.to_double() << endl;
    cout << "... psum 13 : " << psum_13.to_double() << endl;
    cout << "... psum 14 : " << psum_14.to_double() << endl;
    cout << "... psum 15 : " << psum_15.to_double() << endl;
    cout << "... psum 16 : " << psum_16.to_double() << endl;
    cout << "... psum 17 : " << psum_17.to_double() << endl;
    cout << endl;
    #endif
    #endif

    accum_0 += psum_0;
    accum_1 += psum_1;
    accum_2 += psum_2;
    accum_3 += psum_3;
    accum_4 += psum_4;
    accum_5 += psum_5;
    accum_6 += psum_6;
    accum_7 += psum_7;
    accum_8 += psum_8;
    accum_9 += psum_9;
    accum_10 += psum_10;
    accum_11 += psum_11;
    accum_12 += psum_12;
    accum_13 += psum_13;
    accum_14 += psum_14;
    accum_15 += psum_15;
    accum_16 += psum_16;
    accum_17 += psum_17;

    #ifndef __SYNTHESIS__
    #if ISODEBUG
    cout << "--- after incrementing " << endl;
    cout << "... accum 0 : " << accum_0.to_double() << endl;
    cout << "... accum 1 : " << accum_1.to_double() << endl;
    cout << "... accum 2 : " << accum_2.to_double() << endl;
    cout << "... accum 3 : " << accum_3.to_double() << endl;
    cout << "... accum 4 : " << accum_4.to_double() << endl;
    cout << "... accum 5 : " << accum_5.to_double() << endl;
    cout << "... accum 6 : " << accum_6.to_double() << endl;
    cout << "... accum 7 : " << accum_7.to_double() << endl;
    cout << "... accum 8 : " << accum_8.to_double() << endl;
    cout << "... accum 9 : " << accum_9.to_double() << endl;
    cout << "... accum 10 : " << accum_10.to_double() << endl;
    cout << "... accum 11 : " << accum_11.to_double() << endl;
    cout << "... accum 12 : " << accum_12.to_double() << endl;
    cout << "... accum 13 : " << accum_13.to_double() << endl;
    cout << "... accum 14 : " << accum_14.to_double() << endl;
    cout << "... accum 15 : " << accum_15.to_double() << endl;
    cout << "... accum 16 : " << accum_16.to_double() << endl;
    cout << "... accum 17 : " << accum_17.to_double() << endl;
    cout << endl;
    #endif
    #endif

    if (is_last == 0)
        result = 0; // an empty result

    else
    {
        // compute the final energy sum
        iso_accum_t tot_sum;
        tot_sum = accum_0
            + accum_1
            + accum_2
            + accum_3
            + accum_4
            + accum_5
            + accum_6
            + accum_7
            + accum_8
            + accum_9
            + accum_10
            + accum_11
            + accum_12
            + accum_13
            + accum_14
            + accum_15
            + accum_16
            + accum_17;

        // if (tot_sum < iso_threshold_1)
        //     result.range(0,0) = 1;
        // else
        //     result.range(0,0) = 0;

        // if (tot_sum < iso_threshold_2)
        //     result.range(1,1) = 1;
        // else
        //     result.range(1,1) = 0;

        result.range(0,0) = (tot_sum < iso_threshold_1 ? 1 : 0);
        result.range(1,1) = (tot_sum < iso_threshold_2 ? 1 : 0);

        #ifndef __SYNTHESIS__
        #if ISODEBUG
        cout << ".......  is last, tot sum is " << tot_sum.to_double() << endl;
        cout << ".......  is last, iso result is " << result.to_double() << endl;
        #endif
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
    }

    return result;
}


template <int i_module>
hw_iso_t isolation_single_muon_9trk(muon_t in_mu, track_data_9_t in_tracks, ap_uint<1> is_last, iso_accum_t iso_threshold_1, iso_accum_t iso_threshold_2)
{
    #pragma HLS pipeline II=1
    #pragma HLS interface ap_stable port=iso_threshold_1
    #pragma HLS interface ap_stable port=iso_threshold_2
    #pragma HLS inline

    // the accumulators of the energy sums
    static iso_accum_t accum_0  = 0;
    static iso_accum_t accum_1  = 0;
    static iso_accum_t accum_2  = 0;
    static iso_accum_t accum_3  = 0;
    static iso_accum_t accum_4  = 0;
    static iso_accum_t accum_5  = 0;
    static iso_accum_t accum_6  = 0;
    static iso_accum_t accum_7  = 0;
    static iso_accum_t accum_8  = 0;

    #ifndef __SYNTHESIS__
    #if ISODEBUG
    cout << "--- before incrementing " << endl;
    cout << "... accum 0 : " << accum_0.to_double() << endl;
    cout << "... accum 1 : " << accum_1.to_double() << endl;
    cout << "... accum 2 : " << accum_2.to_double() << endl;
    cout << "... accum 3 : " << accum_3.to_double() << endl;
    cout << "... accum 4 : " << accum_4.to_double() << endl;
    cout << "... accum 5 : " << accum_5.to_double() << endl;
    cout << "... accum 6 : " << accum_6.to_double() << endl;
    cout << "... accum 7 : " << accum_7.to_double() << endl;
    cout << "... accum 8 : " << accum_8.to_double() << endl;
    cout << endl;
    #endif
    #endif

    hw_iso_t result;

    iso_accum_t psum_0 = compute_trk_iso (in_mu, in_tracks.trk_0);
    iso_accum_t psum_1 = compute_trk_iso (in_mu, in_tracks.trk_1);
    iso_accum_t psum_2 = compute_trk_iso (in_mu, in_tracks.trk_2);
    iso_accum_t psum_3 = compute_trk_iso (in_mu, in_tracks.trk_3);
    iso_accum_t psum_4 = compute_trk_iso (in_mu, in_tracks.trk_4);
    iso_accum_t psum_5 = compute_trk_iso (in_mu, in_tracks.trk_5);
    iso_accum_t psum_6 = compute_trk_iso (in_mu, in_tracks.trk_6);
    iso_accum_t psum_7 = compute_trk_iso (in_mu, in_tracks.trk_7);
    iso_accum_t psum_8 = compute_trk_iso (in_mu, in_tracks.trk_8);

    #ifndef __SYNTHESIS__
    #if ISODEBUG
    cout << "... psum 0 : " << psum_0.to_double() << endl;
    cout << "... psum 1 : " << psum_1.to_double() << endl;
    cout << "... psum 2 : " << psum_2.to_double() << endl;
    cout << "... psum 3 : " << psum_3.to_double() << endl;
    cout << "... psum 4 : " << psum_4.to_double() << endl;
    cout << "... psum 5 : " << psum_5.to_double() << endl;
    cout << "... psum 6 : " << psum_6.to_double() << endl;
    cout << "... psum 7 : " << psum_7.to_double() << endl;
    cout << "... psum 8 : " << psum_8.to_double() << endl;
    cout << endl;
    #endif
    #endif

    accum_0 += psum_0;
    accum_1 += psum_1;
    accum_2 += psum_2;
    accum_3 += psum_3;
    accum_4 += psum_4;
    accum_5 += psum_5;
    accum_6 += psum_6;
    accum_7 += psum_7;
    accum_8 += psum_8;

    #ifndef __SYNTHESIS__
    #if ISODEBUG
    cout << "--- after incrementing " << endl;
    cout << "... accum 0 : " << accum_0.to_double() << endl;
    cout << "... accum 1 : " << accum_1.to_double() << endl;
    cout << "... accum 2 : " << accum_2.to_double() << endl;
    cout << "... accum 3 : " << accum_3.to_double() << endl;
    cout << "... accum 4 : " << accum_4.to_double() << endl;
    cout << "... accum 5 : " << accum_5.to_double() << endl;
    cout << "... accum 6 : " << accum_6.to_double() << endl;
    cout << "... accum 7 : " << accum_7.to_double() << endl;
    cout << "... accum 8 : " << accum_8.to_double() << endl;
    cout << endl;
    #endif
    #endif

    if (is_last == 0)
        result = 0; // an empty result

    else
    {
        // compute the final energy sum
        iso_accum_t tot_sum;
        tot_sum = accum_0
            + accum_1
            + accum_2
            + accum_3
            + accum_4
            + accum_5
            + accum_6
            + accum_7
            + accum_8;

        // if (tot_sum < iso_threshold_1)
        //     result.range(0,0) = 1;
        // else
        //     result.range(0,0) = 0;

        // if (tot_sum < iso_threshold_2)
        //     result.range(1,1) = 1;
        // else
        //     result.range(1,1) = 0;

        result.range(0,0) = (tot_sum < iso_threshold_1 ? 1 : 0);
        result.range(1,1) = (tot_sum < iso_threshold_2 ? 1 : 0);

        #ifndef __SYNTHESIS__
        #if ISODEBUG
        cout << ".......  is last, tot sum is " << tot_sum.to_double() << endl;
        cout << ".......  is last, iso result is " << result.to_double() << endl;
        #endif
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
    }

    return result;
}






// template <int i_module>
// ap_uint<1> isolation_single(muon_t in_mu, track_t in_tracks[N_TRK_LINKS], ap_uint<1> is_last)
// {
//     #pragma HLS pipeline II=1
//     #pragma HLS ARRAY_PARTITION variable=in_tracks complete

//     static iso_accum_t accums[N_TRK_LINKS]; // for every trk, compute in parallel accumulated pt sums
//     // FIXME : should I init them explicitly at the 1st clk? -> needs list generation at compile time

//     #pragma HLS ARRAY_PARTITION variable=accums complete

//     ap_uint<1> result;

//     for (size_t itrk = 0; itrk < N_TRK_LINKS; ++itrk)
//     {
//         #pragma HLS unroll

//         track_t trk = in_tracks[itrk];

//         iso_accum_t part_sum;

//         // the eta df shoudl be able to handle the max deta accross the tracker width
//         deta_t deta = abs_delta_sat<deta_t>   (in_mu.eta, trk.eta);
//         dphi_t dphi = abs_delta_roll<dphi_t>  (in_mu.phi, trk.phi);
//         dz0_t  dz0  = abs_delta_sat<dz0_t>    (in_mu.z0, trk.z0);

//         // debug - simpler deltas conditions
//         // deta_t deta = in_mu.eta -  trk.eta;
//         // dphi_t dphi = in_mu.phi -  trk.phi;
//         // dz0_t  dz0  = in_mu.z0  -  trk.z0;

//         bool pass_deta   = (deta < deta_t(c_iso_dangle_max) ? true : false);
//         bool pass_dphi   = (dphi < dphi_t(c_iso_dangle_max) ? true : false);
//         bool pass_dz0    = (dz0  < dz0_t(c_iso_dz_max)      ? true : false);
//         bool pass_trkpt  = (trk.pt >= hw_pt_t(c_iso_pt_min) ? true : false);

//         bool pass_ovrl_deta   = (deta > 0 ? true : false);
//         bool pass_ovrl_dphi   = (dphi > 0 ? true : false);


//         // // debug - remove deta / dphi calculation
//         // bool pass_deta   = (in_mu.eta > trk.eta  ? true : false);
//         // bool pass_dphi   = (in_mu.phi > trk.phi  ? true : false);
//         // bool pass_dz0    = (in_mu.z0  > trk.z0   ? true : false);
//         // bool pass_trkpt  = (trk.pt > 0           ? true : false);

//         // bool pass_ovrl_deta   = (in_mu.eta > 0 ? true : false);
//         // bool pass_ovrl_dphi   = (in_mu.phi > 0 ? true : false);

//         if (
//             // match conditions
//             pass_deta  &&
//             pass_dphi  &&
//             pass_dz0   &&
//             pass_trkpt &&
//             pass_ovrl_deta &&
//             pass_ovrl_dphi
//         )
//             part_sum = trk.pt;
//         else
//             part_sum = 0;

//         // // debug : always accept
//         // part_sum = trk.pt;        
//         // LARGELY reduces resources

//         accums[itrk] += part_sum;
//     }

//     if (is_last == 0)
//         result = 0; // an empty result

//     else
//     {
//         // compute the final energy sum
//         iso_accum_t tot_sum = 0;
//         for (size_t itrk = 0; itrk < N_TRK_LINKS; ++itrk)
//         {
//             #pragma HLS unroll
//             tot_sum += accums[itrk];
//         }

//         if (tot_sum < c_iso_sumpt_thr)
//             result = 1;
//         else
//             result = 0;

//         // reset the accumulators
//         for (size_t itrk = 0; itrk < N_TRK_LINKS; ++itrk)
//         {
//             #pragma HLS unroll
//             accums[itrk] = 0;
//         }
//     }

//     return result;
// }

// ///////////// DEBUG copy below (no loops)
// template <int i_module>
// ap_uint<1> isolation_single_copy(muon_t in_mu, track_t in_tracks[N_TRK_LINKS], ap_uint<1> is_last)
// {
//     #pragma HLS pipeline II=1
//     #pragma HLS ARRAY_PARTITION variable=in_tracks complete

//     static iso_accum_t accums[N_TRK_LINKS]; // for every trk, compute in parallel accumulated pt sums
//     // FIXME : should I init them explicitly at the 1st clk? -> needs list generation at compile time

//     #pragma HLS ARRAY_PARTITION variable=accums complete

//     ap_uint<1> result;


//     track_t trk_0 = in_tracks[0];
//     iso_accum_t part_sum_0;
    
//     deta_t deta_0 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[0].eta);
//     dphi_t dphi_0 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[0].phi);
//     dz0_t  dz0_0  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[0].z0);

//     bool pass_deta_0        = (deta_0      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_0        = (dphi_0      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_0         = (dz0_0       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_0       = (in_tracks[0].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_0   = (deta_0 > 0 ? true : false);
//     bool pass_ovrl_dphi_0   = (dphi_0 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_0  &&
//         pass_dphi_0  &&
//         pass_dz0_0   &&
//         pass_trkpt_0 &&
//         pass_ovrl_deta_0 &&
//         pass_ovrl_dphi_0
//     )
//         part_sum_0 = in_tracks[0].pt;
//     else
//         part_sum_0 = 0;
//     accums[0] += part_sum_0;



//     track_t trk_1 = in_tracks[1];
//     iso_accum_t part_sum_1;
    
//     deta_t deta_1 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[1].eta);
//     dphi_t dphi_1 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[1].phi);
//     dz0_t  dz0_1  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[1].z0);

//     bool pass_deta_1        = (deta_1      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_1        = (dphi_1      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_1         = (dz0_1       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_1       = (in_tracks[1].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_1   = (deta_1 > 0 ? true : false);
//     bool pass_ovrl_dphi_1   = (dphi_1 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_1  &&
//         pass_dphi_1  &&
//         pass_dz0_1   &&
//         pass_trkpt_1 &&
//         pass_ovrl_deta_1 &&
//         pass_ovrl_dphi_1
//     )
//         part_sum_1 = in_tracks[1].pt;
//     else
//         part_sum_1 = 0;
//     accums[1] += part_sum_1;



//     track_t trk_2 = in_tracks[2];
//     iso_accum_t part_sum_2;
    
//     deta_t deta_2 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[2].eta);
//     dphi_t dphi_2 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[2].phi);
//     dz0_t  dz0_2  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[2].z0);

//     bool pass_deta_2        = (deta_2      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_2        = (dphi_2      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_2         = (dz0_2       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_2       = (in_tracks[2].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_2   = (deta_2 > 0 ? true : false);
//     bool pass_ovrl_dphi_2   = (dphi_2 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_2  &&
//         pass_dphi_2  &&
//         pass_dz0_2   &&
//         pass_trkpt_2 &&
//         pass_ovrl_deta_2 &&
//         pass_ovrl_dphi_2
//     )
//         part_sum_2 = in_tracks[2].pt;
//     else
//         part_sum_2 = 0;
//     accums[2] += part_sum_2;



//     track_t trk_3 = in_tracks[3];
//     iso_accum_t part_sum_3;
    
//     deta_t deta_3 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[3].eta);
//     dphi_t dphi_3 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[3].phi);
//     dz0_t  dz0_3  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[3].z0);

//     bool pass_deta_3        = (deta_3      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_3        = (dphi_3      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_3         = (dz0_3       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_3       = (in_tracks[3].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_3   = (deta_3 > 0 ? true : false);
//     bool pass_ovrl_dphi_3   = (dphi_3 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_3  &&
//         pass_dphi_3  &&
//         pass_dz0_3   &&
//         pass_trkpt_3 &&
//         pass_ovrl_deta_3 &&
//         pass_ovrl_dphi_3
//     )
//         part_sum_3 = in_tracks[3].pt;
//     else
//         part_sum_3 = 0;
//     accums[3] += part_sum_3;



//     track_t trk_4 = in_tracks[4];
//     iso_accum_t part_sum_4;
    
//     deta_t deta_4 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[4].eta);
//     dphi_t dphi_4 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[4].phi);
//     dz0_t  dz0_4  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[4].z0);

//     bool pass_deta_4        = (deta_4      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_4        = (dphi_4      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_4         = (dz0_4       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_4       = (in_tracks[4].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_4   = (deta_4 > 0 ? true : false);
//     bool pass_ovrl_dphi_4   = (dphi_4 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_4  &&
//         pass_dphi_4  &&
//         pass_dz0_4   &&
//         pass_trkpt_4 &&
//         pass_ovrl_deta_4 &&
//         pass_ovrl_dphi_4
//     )
//         part_sum_4 = in_tracks[4].pt;
//     else
//         part_sum_4 = 0;
//     accums[4] += part_sum_4;



//     track_t trk_5 = in_tracks[5];
//     iso_accum_t part_sum_5;
    
//     deta_t deta_5 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[5].eta);
//     dphi_t dphi_5 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[5].phi);
//     dz0_t  dz0_5  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[5].z0);

//     bool pass_deta_5        = (deta_5      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_5        = (dphi_5      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_5         = (dz0_5       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_5       = (in_tracks[5].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_5   = (deta_5 > 0 ? true : false);
//     bool pass_ovrl_dphi_5   = (dphi_5 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_5  &&
//         pass_dphi_5  &&
//         pass_dz0_5   &&
//         pass_trkpt_5 &&
//         pass_ovrl_deta_5 &&
//         pass_ovrl_dphi_5
//     )
//         part_sum_5 = in_tracks[5].pt;
//     else
//         part_sum_5 = 0;
//     accums[5] += part_sum_5;



//     track_t trk_6 = in_tracks[6];
//     iso_accum_t part_sum_6;
    
//     deta_t deta_6 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[6].eta);
//     dphi_t dphi_6 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[6].phi);
//     dz0_t  dz0_6  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[6].z0);

//     bool pass_deta_6        = (deta_6      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_6        = (dphi_6      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_6         = (dz0_6       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_6       = (in_tracks[6].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_6   = (deta_6 > 0 ? true : false);
//     bool pass_ovrl_dphi_6   = (dphi_6 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_6  &&
//         pass_dphi_6  &&
//         pass_dz0_6   &&
//         pass_trkpt_6 &&
//         pass_ovrl_deta_6 &&
//         pass_ovrl_dphi_6
//     )
//         part_sum_6 = in_tracks[6].pt;
//     else
//         part_sum_6 = 0;
//     accums[6] += part_sum_6;



//     track_t trk_7 = in_tracks[7];
//     iso_accum_t part_sum_7;
    
//     deta_t deta_7 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[7].eta);
//     dphi_t dphi_7 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[7].phi);
//     dz0_t  dz0_7  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[7].z0);

//     bool pass_deta_7        = (deta_7      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_7        = (dphi_7      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_7         = (dz0_7       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_7       = (in_tracks[7].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_7   = (deta_7 > 0 ? true : false);
//     bool pass_ovrl_dphi_7   = (dphi_7 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_7  &&
//         pass_dphi_7  &&
//         pass_dz0_7   &&
//         pass_trkpt_7 &&
//         pass_ovrl_deta_7 &&
//         pass_ovrl_dphi_7
//     )
//         part_sum_7 = in_tracks[7].pt;
//     else
//         part_sum_7 = 0;
//     accums[7] += part_sum_7;



//     track_t trk_8 = in_tracks[8];
//     iso_accum_t part_sum_8;
    
//     deta_t deta_8 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[8].eta);
//     dphi_t dphi_8 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[8].phi);
//     dz0_t  dz0_8  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[8].z0);

//     bool pass_deta_8        = (deta_8      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_8        = (dphi_8      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_8         = (dz0_8       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_8       = (in_tracks[8].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_8   = (deta_8 > 0 ? true : false);
//     bool pass_ovrl_dphi_8   = (dphi_8 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_8  &&
//         pass_dphi_8  &&
//         pass_dz0_8   &&
//         pass_trkpt_8 &&
//         pass_ovrl_deta_8 &&
//         pass_ovrl_dphi_8
//     )
//         part_sum_8 = in_tracks[8].pt;
//     else
//         part_sum_8 = 0;
//     accums[8] += part_sum_8;



//     track_t trk_9 = in_tracks[9];
//     iso_accum_t part_sum_9;
    
//     deta_t deta_9 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[9].eta);
//     dphi_t dphi_9 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[9].phi);
//     dz0_t  dz0_9  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[9].z0);

//     bool pass_deta_9        = (deta_9      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_9        = (dphi_9      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_9         = (dz0_9       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_9       = (in_tracks[9].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_9   = (deta_9 > 0 ? true : false);
//     bool pass_ovrl_dphi_9   = (dphi_9 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_9  &&
//         pass_dphi_9  &&
//         pass_dz0_9   &&
//         pass_trkpt_9 &&
//         pass_ovrl_deta_9 &&
//         pass_ovrl_dphi_9
//     )
//         part_sum_9 = in_tracks[9].pt;
//     else
//         part_sum_9 = 0;
//     accums[9] += part_sum_9;



//     track_t trk_10 = in_tracks[10];
//     iso_accum_t part_sum_10;
    
//     deta_t deta_10 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[10].eta);
//     dphi_t dphi_10 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[10].phi);
//     dz0_t  dz0_10  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[10].z0);

//     bool pass_deta_10        = (deta_10      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_10        = (dphi_10      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_10         = (dz0_10       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_10       = (in_tracks[10].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_10   = (deta_10 > 0 ? true : false);
//     bool pass_ovrl_dphi_10   = (dphi_10 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_10  &&
//         pass_dphi_10  &&
//         pass_dz0_10   &&
//         pass_trkpt_10 &&
//         pass_ovrl_deta_10 &&
//         pass_ovrl_dphi_10
//     )
//         part_sum_10 = in_tracks[10].pt;
//     else
//         part_sum_10 = 0;
//     accums[10] += part_sum_10;



//     track_t trk_11 = in_tracks[11];
//     iso_accum_t part_sum_11;
    
//     deta_t deta_11 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[11].eta);
//     dphi_t dphi_11 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[11].phi);
//     dz0_t  dz0_11  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[11].z0);

//     bool pass_deta_11        = (deta_11      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_11        = (dphi_11      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_11         = (dz0_11       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_11       = (in_tracks[11].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_11   = (deta_11 > 0 ? true : false);
//     bool pass_ovrl_dphi_11   = (dphi_11 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_11  &&
//         pass_dphi_11  &&
//         pass_dz0_11   &&
//         pass_trkpt_11 &&
//         pass_ovrl_deta_11 &&
//         pass_ovrl_dphi_11
//     )
//         part_sum_11 = in_tracks[11].pt;
//     else
//         part_sum_11 = 0;
//     accums[11] += part_sum_11;



//     track_t trk_12 = in_tracks[12];
//     iso_accum_t part_sum_12;
    
//     deta_t deta_12 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[12].eta);
//     dphi_t dphi_12 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[12].phi);
//     dz0_t  dz0_12  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[12].z0);

//     bool pass_deta_12        = (deta_12      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_12        = (dphi_12      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_12         = (dz0_12       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_12       = (in_tracks[12].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_12   = (deta_12 > 0 ? true : false);
//     bool pass_ovrl_dphi_12   = (dphi_12 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_12  &&
//         pass_dphi_12  &&
//         pass_dz0_12   &&
//         pass_trkpt_12 &&
//         pass_ovrl_deta_12 &&
//         pass_ovrl_dphi_12
//     )
//         part_sum_12 = in_tracks[12].pt;
//     else
//         part_sum_12 = 0;
//     accums[12] += part_sum_12;



//     track_t trk_13 = in_tracks[13];
//     iso_accum_t part_sum_13;
    
//     deta_t deta_13 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[13].eta);
//     dphi_t dphi_13 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[13].phi);
//     dz0_t  dz0_13  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[13].z0);

//     bool pass_deta_13        = (deta_13      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_13        = (dphi_13      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_13         = (dz0_13       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_13       = (in_tracks[13].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_13   = (deta_13 > 0 ? true : false);
//     bool pass_ovrl_dphi_13   = (dphi_13 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_13  &&
//         pass_dphi_13  &&
//         pass_dz0_13   &&
//         pass_trkpt_13 &&
//         pass_ovrl_deta_13 &&
//         pass_ovrl_dphi_13
//     )
//         part_sum_13 = in_tracks[13].pt;
//     else
//         part_sum_13 = 0;
//     accums[13] += part_sum_13;



//     track_t trk_14 = in_tracks[14];
//     iso_accum_t part_sum_14;
    
//     deta_t deta_14 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[14].eta);
//     dphi_t dphi_14 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[14].phi);
//     dz0_t  dz0_14  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[14].z0);

//     bool pass_deta_14        = (deta_14      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_14        = (dphi_14      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_14         = (dz0_14       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_14       = (in_tracks[14].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_14   = (deta_14 > 0 ? true : false);
//     bool pass_ovrl_dphi_14   = (dphi_14 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_14  &&
//         pass_dphi_14  &&
//         pass_dz0_14   &&
//         pass_trkpt_14 &&
//         pass_ovrl_deta_14 &&
//         pass_ovrl_dphi_14
//     )
//         part_sum_14 = in_tracks[14].pt;
//     else
//         part_sum_14 = 0;
//     accums[14] += part_sum_14;



//     track_t trk_15 = in_tracks[15];
//     iso_accum_t part_sum_15;
    
//     deta_t deta_15 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[15].eta);
//     dphi_t dphi_15 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[15].phi);
//     dz0_t  dz0_15  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[15].z0);

//     bool pass_deta_15        = (deta_15      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_15        = (dphi_15      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_15         = (dz0_15       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_15       = (in_tracks[15].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_15   = (deta_15 > 0 ? true : false);
//     bool pass_ovrl_dphi_15   = (dphi_15 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_15  &&
//         pass_dphi_15  &&
//         pass_dz0_15   &&
//         pass_trkpt_15 &&
//         pass_ovrl_deta_15 &&
//         pass_ovrl_dphi_15
//     )
//         part_sum_15 = in_tracks[15].pt;
//     else
//         part_sum_15 = 0;
//     accums[15] += part_sum_15;



//     track_t trk_16 = in_tracks[16];
//     iso_accum_t part_sum_16;
    
//     deta_t deta_16 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[16].eta);
//     dphi_t dphi_16 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[16].phi);
//     dz0_t  dz0_16  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[16].z0);

//     bool pass_deta_16        = (deta_16      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_16        = (dphi_16      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_16         = (dz0_16       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_16       = (in_tracks[16].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_16   = (deta_16 > 0 ? true : false);
//     bool pass_ovrl_dphi_16   = (dphi_16 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_16  &&
//         pass_dphi_16  &&
//         pass_dz0_16   &&
//         pass_trkpt_16 &&
//         pass_ovrl_deta_16 &&
//         pass_ovrl_dphi_16
//     )
//         part_sum_16 = in_tracks[16].pt;
//     else
//         part_sum_16 = 0;
//     accums[16] += part_sum_16;



//     track_t trk_17 = in_tracks[17];
//     iso_accum_t part_sum_17;
    
//     deta_t deta_17 = abs_delta_sat<deta_t>   (in_mu.eta, in_tracks[17].eta);
//     dphi_t dphi_17 = abs_delta_roll<dphi_t>  (in_mu.phi, in_tracks[17].phi);
//     dz0_t  dz0_17  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_tracks[17].z0);

//     bool pass_deta_17        = (deta_17      < deta_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dphi_17        = (dphi_17      < dphi_t(c_iso_dangle_max)    ? true : false);
//     bool pass_dz0_17         = (dz0_17       < dz0_t(c_iso_dz_max)         ? true : false);
//     bool pass_trkpt_17       = (in_tracks[17].pt >= hw_pt_t(c_iso_pt_min)  ? true : false);
//     bool pass_ovrl_deta_17   = (deta_17 > 0 ? true : false);
//     bool pass_ovrl_dphi_17   = (dphi_17 > 0 ? true : false);

//     if (
//         // match conditions
//         pass_deta_17  &&
//         pass_dphi_17  &&
//         pass_dz0_17   &&
//         pass_trkpt_17 &&
//         pass_ovrl_deta_17 &&
//         pass_ovrl_dphi_17
//     )
//         part_sum_17 = in_tracks[17].pt;
//     else
//         part_sum_17 = 0;
//     accums[17] += part_sum_17;


    
//     if (is_last == 0)
//         result = 0; // an empty result

//     else
//     {
//         // compute the final energy sum
//         iso_accum_t tot_sum;
//         tot_sum = accums[0] +
//             accums[1] +  
//             accums[2] +  
//             accums[3] +  
//             accums[4] +  
//             accums[5] +  
//             accums[6] +  
//             accums[7] +  
//             accums[8] +  
//             accums[9] +  
//             accums[10] +  
//             accums[11] +  
//             accums[12] +  
//             accums[13] +  
//             accums[14] +  
//             accums[15] +  
//             accums[16] +  
//             accums[17] ;

//         if (tot_sum < c_iso_sumpt_thr)
//             result = 1;
//         else
//             result = 0;

//         // reset the accumulators
//         accums[0] = 0;
//         accums[1] = 0;
//         accums[2] = 0;
//         accums[3] = 0;
//         accums[4] = 0;
//         accums[5] = 0;
//         accums[6] = 0;
//         accums[7] = 0;
//         accums[8] = 0;
//         accums[9] = 0;
//         accums[10] = 0;
//         accums[11] = 0;
//         accums[12] = 0;
//         accums[13] = 0;
//         accums[14] = 0;
//         accums[15] = 0;
//         accums[16] = 0;
//         accums[17] = 0;        
//     }

//     return result;
// }


#endif