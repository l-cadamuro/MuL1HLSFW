#ifndef DF_CONVERSIONS_H
#define DF_CONVERSIONS_H

#include "track_t.h"

void qOverR_to_pt (ap_int<TRACK_QOVERR_W> qor_in, ap_uint<TRACKCONV_PT_W>& pt_out, ap_uint<1>& q_out);
void track_to_trackconv (track_t t_in, track_conv_t& tconv_out);

#endif