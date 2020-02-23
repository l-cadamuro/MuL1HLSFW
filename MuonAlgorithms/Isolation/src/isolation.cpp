#include "isolation.h"

iso_accum_t compute_trk_iso (muon_t in_mu, track_t in_trk)
{

    // ap_uint<TRACKCONV_ETA_W> deta16 = abs_delta<ap_uint<TRACKCONV_ETA_W>, ap_int<TRACKCONV_ETA_W> > (trks_in.trk16.eta, mu_in.eta);
    // ap_uint<TRACKCONV_PHI_W> dphi16 = abs_delta<ap_uint<TRACKCONV_PHI_W>, ap_int<TRACKCONV_PHI_W> > (trks_in.trk16.phi, mu_in.phi);
    // ap_ufixed<c_iso_thr_w, c_iso_thr_w, AP_TRN, AP_SAT> psum16;
    // if (deta16 < c_deta_iso && dphi16 < c_dphi_iso && deta16 > 1 && dphi16 > 1 && trks_in.trk16.pt > c_trk_iso_thresh)
    //     psum16 = trks_in.trk16.pt;
    // else
    //     psum16 = 0;

    #pragma HLS INLINE

    deta_t deta = abs_delta_sat<deta_t>   (in_mu.eta, in_trk.eta);
    dphi_t dphi = abs_delta_roll<dphi_t>  (in_mu.phi, in_trk.phi);
    dz0_t  dz0  = abs_delta_sat<dz0_t>    (in_mu.z0,  in_trk.z0);

    // ap_uint<deta_t::width> deta = abs_delta<ap_uint<deta_t::width>, deta_t > (in_trk.eta, in_mu.eta);
    // ap_uint<dphi_t::width> dphi = abs_delta<ap_uint<dphi_t::width>, dphi_t > (in_trk.phi, in_mu.phi);

    iso_accum_t part_sum;

    bool pass_deta        = (deta      < deta_t(c_iso_dangle_max)    ? true : false);
    bool pass_dphi        = (dphi      < dphi_t(c_iso_dangle_max)    ? true : false);
    bool pass_dz0         = (dz0       < dz0_t(c_iso_dz_max)        ? true : false);
    bool pass_trkpt       = (in_trk.pt >= hw_pt_t(c_iso_pt_min)       ? true : false);
    bool pass_ovrl_deta   = (deta > 0 ? true : false);
    bool pass_ovrl_dphi   = (dphi > 0 ? true : false);

    if (
        // match conditions
        pass_deta  &&
        pass_dphi  &&
        pass_dz0   &&
        pass_trkpt &&
        pass_ovrl_deta &&
        pass_ovrl_dphi
    )
        part_sum = in_trk.pt;
    else
        part_sum = 0;

    return part_sum;
}

void isolation(muon_t in_muons[N_MUON], track_t in_tracks[N_TRK_LINKS], ap_uint<1> is_last, ap_uint<1> iso_flags[N_MUON])
{
    #pragma HLS pipeline II=1
    
    #pragma HLS ARRAY_PARTITION variable=in_muons   complete
    #pragma HLS ARRAY_PARTITION variable=in_tracks  complete
    #pragma HLS ARRAY_PARTITION variable=iso_flags  complete

    // for (size_t imu = 0; imu < N_MUON; ++imu)
    // {
    //     #pragma HLS unroll
    //     iso_flags[imu] = isolation_single <imu> (in_muons[imu], in_tracks, is_last);
    // }

    // FIXME : to be unrolled automatically (loop idx cannot be used in template - need preprocessor loop expansion)

    iso_flags[0] = isolation_single <0> (in_muons[0], in_tracks, is_last);
    iso_flags[1] = isolation_single <1> (in_muons[1], in_tracks, is_last);
    iso_flags[2] = isolation_single <2> (in_muons[2], in_tracks, is_last);
    iso_flags[3] = isolation_single <3> (in_muons[3], in_tracks, is_last);
    iso_flags[4] = isolation_single <4> (in_muons[4], in_tracks, is_last);
    iso_flags[5] = isolation_single <5> (in_muons[5], in_tracks, is_last);
    iso_flags[6] = isolation_single <6> (in_muons[6], in_tracks, is_last);
    iso_flags[7] = isolation_single <7> (in_muons[7], in_tracks, is_last);
    iso_flags[8] = isolation_single <8> (in_muons[8], in_tracks, is_last);
    iso_flags[9] = isolation_single <9> (in_muons[9], in_tracks, is_last);
    iso_flags[10] = isolation_single <10> (in_muons[10], in_tracks, is_last);
    iso_flags[11] = isolation_single <11> (in_muons[11], in_tracks, is_last);

}


void isolation_allmu(muon_data_t in_muons, track_data_t in_tracks, ap_uint<1> is_last, muon_isodata_t& iso_flags)
{
    #pragma HLS pipeline II=1

    iso_flags.isomu_0  = isolation_single_muon <0>  (in_muons.mu_0,  in_tracks, is_last);
    iso_flags.isomu_1  = isolation_single_muon <1>  (in_muons.mu_1,  in_tracks, is_last);
    iso_flags.isomu_2  = isolation_single_muon <2>  (in_muons.mu_2,  in_tracks, is_last);
    iso_flags.isomu_3  = isolation_single_muon <3>  (in_muons.mu_3,  in_tracks, is_last);
    iso_flags.isomu_4  = isolation_single_muon <4>  (in_muons.mu_4,  in_tracks, is_last);
    iso_flags.isomu_5  = isolation_single_muon <5>  (in_muons.mu_5,  in_tracks, is_last);
    iso_flags.isomu_6  = isolation_single_muon <6>  (in_muons.mu_6,  in_tracks, is_last);
    iso_flags.isomu_7  = isolation_single_muon <7>  (in_muons.mu_7,  in_tracks, is_last);
    iso_flags.isomu_8  = isolation_single_muon <8>  (in_muons.mu_8,  in_tracks, is_last);
    iso_flags.isomu_9  = isolation_single_muon <9>  (in_muons.mu_9,  in_tracks, is_last);
    iso_flags.isomu_10 = isolation_single_muon <10> (in_muons.mu_10, in_tracks, is_last);
    iso_flags.isomu_11 = isolation_single_muon <11> (in_muons.mu_11, in_tracks, is_last);
    
}

ap_uint<1> isolation_single_muon_wrap(muon_t in_mu, track_data_t in_tracks, ap_uint<1> is_last)
{
    #pragma HLS pipeline II=1
    return isolation_single_muon <999> (in_mu, in_tracks, is_last);
}



ap_uint<1> isolation_single_wrap(muon_t in_mu, track_t in_tracks[N_TRK_LINKS], ap_uint<1> is_last)
{
    #pragma HLS pipeline II=1
    #pragma HLS ARRAY_PARTITION variable=in_tracks  complete
    return isolation_single <999> (in_mu, in_tracks, is_last);
    // return isolation_single_copy <999> (in_mu, in_tracks, is_last);
}

