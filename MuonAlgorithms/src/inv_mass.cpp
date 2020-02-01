//////////// IMPLEMENTATION
// Same algo as in the current uGT : m^2 = 2 pT(1) pT(2)  [ cosh (delta_eta) - cos (delta_phi) ]
// and the algorithm returns m^2/2


//////////////////
/// -------- NB : THE TEXT BELOW MAY BE OUTDATED --------
//////////////////
///
// the invariant of two massless particles can be expressed as:
// ptprod * angprod
// with 
// ptprod  = 2. * v1.Pt()*v2.Pt()
// angprod = 1./(sin(th1)*sin(th2)) - cos(phi1-phi2) - 1/(tan(th1)*tan(th2))
// where th1, phi1, th2, phi2 are the angles of the two spatial momenta, in the ROOT TLorentzVector coordinate system, i.e.
// x = p  sin(theta) cos(phi)
// y = p  sin(theta) sin(phi)
// z = p  cos(theta)
//
// therefore the implementation is:
// compute ptprod with a DSP
// compute the 3 partial sums with 3 LUT (for the cos one, input dphi)
// sum the 3 partial sums
// multiplicate all

#include "inv_mass.h"
#include "inv_mass_df.h"
#include "cos_dphi_LUT.h"
#include "cosh_deta_LUT.h"
#include "utils/x_hls_traits.h"

hw_minv2over2_t inv_mass(hw_p3_t part1, hw_p3_t part2)
{
    #pragma HLS pipeline II=1

    // pt part
    ap_fixed<24, 13, AP_TRN, AP_SAT> hw_ptprod  = part1.pt * part2.pt; // sum the output bits in the multiplication - for a DSP with 18*24
    
    // phi part
    ap_fixed<COS_PHI_ADDR_W+1, COS_PHI_ADDR_W_INT+1> dphi_sign =  part1.phi - part2.phi;
    if (dphi_sign < 0)
        dphi_sign = -1*dphi_sign;
    ap_ufixed<COS_PHI_ADDR_W, COS_PHI_ADDR_W_INT> dphi = dphi_sign;
    
    // ap_fixed<COS_PHI_OUT_W, 1> hw_phipart = -1*cos_lut<COS_PHI_OUT_W> (dphi);
    ap_fixed<COS_PHI_OUT_W, 1> hw_phipart = cos_dphi_lut[dphi.range()]; // the LUT is indexed by the integer repr of the number, returned by range()

    // eta part
    ap_fixed<COSH_ETA_ADDR_W+1, COSH_ETA_ADDR_W_INT+1> deta_sign =  part1.eta - part2.eta;
    if (deta_sign < 0)
        deta_sign = -1*deta_sign;
    ap_ufixed<COSH_ETA_ADDR_W, COSH_ETA_ADDR_W_INT> deta = deta_sign;
    // ap_ufixed<COSH_ETA_OUT_W, COSH_ETA_OUT_W_INT> hw_etapart = cosh_lut<COSH_ETA_OUT_W, COSH_ETA_OUT_W_INT> (deta);
    ap_ufixed<COSH_ETA_OUT_W, COSH_ETA_OUT_W_INT> hw_etapart = cosh_deta_lut[deta.range()];

    // double ptprod = hw_ptprod.to_double();
    // double a1     = hw_phipart.to_double();
    // double a2     = hw_etapart.to_double();

    // FIXME : the type of the pt sum
    ap_ufixed<18, 8, AP_TRN, AP_SAT> hw_angsum = hw_etapart - hw_phipart;

    hw_minv2over2_t minv = hw_ptprod * hw_angsum;

    return minv;
}

hw_minv2over2_t inv_mass_3body (hw_p3_t part1, hw_p3_t part2, hw_p3_t part3)
{
    #pragma HLS pipeline II=1

    hw_minv2over2_t m12;
    hw_minv2over2_t m23;
    hw_minv2over2_t m31;

    m12 = inv_mass(part1, part2);
    m23 = inv_mass(part2, part3);
    m31 = inv_mass(part3, part1);

    hw_minv2over2_t mSum;
    mSum = m12 + m23 + m31;

    return mSum;

}
