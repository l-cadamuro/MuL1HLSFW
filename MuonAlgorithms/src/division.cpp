#include "division.h"

ap_fixed<24, 12> division(ap_fixed<12, 6> num, ap_fixed<12, 6> den)
{
    #pragma HLS pipeline II=1
    ap_fixed<24, 12> res = num/den;
    return res;
}