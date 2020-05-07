#ifndef MINISORTER_H
#define MINISORTER_H

#include "dataformats.h"

typedef muon_t data_c_t;
typedef hw_phi_t spar_c_t; // the sorting parameter used

class miniSorter {
    public:
        miniSorter();
        void make_null_data(data_c_t *data);
        void copy_data(data_c_t from, data_c_t *to);
        void get_sort_par(data_c_t data_in, spar_c_t *spar_out);
        void abs_delta(spar_c_t p1, spar_c_t p2, spar_c_t *pdelta);
        void mini_closer_elem(data_c_t data_in, data_c_t target, ap_uint<1> reset, data_c_t *data_out_0, data_c_t *data_out_1, data_c_t *data_out_2);

    private:
        spar_c_t shift_reg_spar[3];
        data_c_t shift_reg_data[3];
};



#endif