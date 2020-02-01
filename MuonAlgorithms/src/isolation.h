#ifndef ISOLATION_H
#define ISOLATION_H

#include "ap_fixed.h"
#include "ap_int.h"

#include "muon_t.h"
#include "track_t.h"
#include "track_input.h"
#include "muon_input.h"

#include "iso_calculator.h"

// compute the isolation of one muon (mu) with N tracks streaming in
// returns true if isolated, else false
ap_uint<1> isolation (muon_t mu_in, track_conv_input trks_in, ap_uint<1> last_in);

ap_uint<1> isolation_class (muon_t mu_in, track_conv_input trks_in, ap_uint<1> last_in);

iso_muon_input isolation_class_allmu (muon_input mus_in, track_conv_input trks_in, ap_uint<1> last_in);

// does abs(v2 - v1)
template <class tin, class tout>
tout abs_delta (tin v1, tin v2)
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