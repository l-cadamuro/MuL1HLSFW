#ifndef TEST_ALGO_H
#define TEST_ALGO_H

// data formats
#include "muon_t.h"
#include "track_t.h"

// organisation of the inputs
#include "muon_input.h"
#include "track_input.h"
#include "df_conversions.h"

// struct iso_muon_t : public muon_t {
//     ap_uint<1> iso;
// };

struct iso_muons_out {
    iso_muon_t mu0;
    iso_muon_t mu1;
    iso_muon_t mu2;
    iso_muon_t mu3;
    iso_muon_t mu4;
    iso_muon_t mu5;
    iso_muon_t mu6;
    iso_muon_t mu7;
    iso_muon_t mu8;
    iso_muon_t mu9;
    iso_muon_t mu10;
    iso_muon_t mu11;
    iso_muon_t mu12;
    iso_muon_t mu13;
    iso_muon_t mu14;
    iso_muon_t mu15;
    iso_muon_t mu16;
    iso_muon_t mu17;
};

iso_muons_out test_algo(muon_input mu_in, track_input trk_in, ap_uint<1> last_in);

// does v2 - v1
template <class t>
t delta (t v1, t v2)
{
    return v2-v1;
}

#endif // TEST_ALGO_H
