// Constants that are used to drive the isolation FW but are meant to be supplied externally
// These include
// - numbers that characterise the data sent (N muons, N tracks)
// - constants that are configured externally and passed as registered data ports
// These values do not enter the FW synthesis

#ifndef EXTERNAL_CONSTANTS_H
#define EXTERNAL_CONSTANTS_H

// NOTE : these numbers must be in line with the dataformat

// the track LSB is 1/2^5 GeV
// the angle LSB is (2pi)/(2^13) rad ~= 0.000767 rad
#define c_iso_sumpt_thr_1  96  // < , 96  x 1/2^5    = 3 GeV
#define c_iso_sumpt_thr_2  96  // < , 96  x 1/2^5    = 3 GeV

#define N_TRK_LINKS    18 // total number of input links (each gives 1 trk / clk)
// #define N_TRK_PER_LINK 18 // total number of track serially trasmitted on each link
#define N_MUON 12

#define RUN_DOUBLE_CLK_SPEED 1 // 1 : N_TRK_LINKS/2 tracks/clk, 0 : N_TRK_LINKS tracks/clk


#endif // EXTERNAL_CONSTANTS_H