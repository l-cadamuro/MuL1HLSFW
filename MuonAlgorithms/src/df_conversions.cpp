#include "df_conversions.h"
#include "qOverR_to_pt_lut.h"

void qOverR_to_pt (ap_int<TRACK_QOVERR_W> qor_in, ap_uint<TRACKCONV_PT_W>& pt_out, ap_uint<1>& q_out)
{

    ap_uint<TRACK_QOVERR_W-1> a_qor_in;
    if (qor_in < 0)
        a_qor_in = -1*qor_in;
    else
        a_qor_in = qor_in;

    pt_out = qOverR_to_pt_lut[a_qor_in];
    q_out  = (qor_in >= 0 ? 1 : 0);

    return;
}


void track_to_trackconv (track_t t_in, track_conv_t& tconv_out)
{
    ap_uint<TRACKCONV_PT_W> ptbuf;
    ap_uint<1>              qbuf;
    qOverR_to_pt(t_in.qOverR, ptbuf, qbuf);
    
    tconv_out.pt   = ptbuf; 
    tconv_out.q    = qbuf;
    tconv_out.eta  = t_in.tanLambda; // FIXME!
    tconv_out.z0   = t_in.z0;
    tconv_out.d0   = t_in.d0;
    tconv_out.chi2 = t_in.chi2;
}