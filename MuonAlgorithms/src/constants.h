#ifndef CONSTANTS_H
#define CONSTANTS_H

// this header defines the constants that are used thoughout the code
// the scales, thresholds, etc.. have to be adjusted here

// assuming my muon is on a scale expressed in LSB of 0.25 GeV (to make in GeV : number >> 2)

#define c_iso_thr   40 // 10 GeV * 4
#define c_iso_thr_w 6  // c_iso_thr must be expressed in c_iso_thr_w bits

// #define c_iso_thr   1000
// #define c_iso_thr_w 10  

// for eta / phi : express them over 9 bits, i.e. in units of 2*pi / 2^9 = 0.0122718463
// [cfr : eta is in units of 0.010875 in the  GMT ]
// for eta this is slightly suboptimal but let's consider it is fine

#define c_deta_iso  33 // 0.4 radius cone (actually a square, +/ 0.4) => 0.4/0.0122718463 = 32.6 ~ 33
#define c_dphi_iso  33 // 

#endif