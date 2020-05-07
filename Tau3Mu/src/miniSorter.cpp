#include "miniSorter.h"

miniSorter::miniSorter() {
    
    #pragma HLS inline
    #pragma HLS array_partition variable = shift_reg_spar complete
    #pragma HLS array_partition variable = shift_reg_data complete

}

void miniSorter::make_null_data(data_c_t *data)
{
    #pragma HLS inline
    (*data).phi = 0;
}

void miniSorter::copy_data(data_c_t from, data_c_t *to)
{
    #pragma HLS inline
    (*to).phi = from.phi;
}
void miniSorter::get_sort_par(data_c_t data_in, spar_c_t *spar_out)
{
    #pragma HLS inline
    *spar_out = data_in.phi;
}

void miniSorter::abs_delta(spar_c_t p1, spar_c_t p2, spar_c_t *pdelta)
{
    #pragma HLS inline

    spar_c_t delta = p1 - p2;
    if (delta > 0)
        *pdelta = delta;
    else
        *pdelta = -1 * delta;
}

void miniSorter::mini_closer_elem(data_c_t data_in, data_c_t target, ap_uint<1> reset, data_c_t *data_out_0, data_c_t *data_out_1, data_c_t *data_out_2)
{
    // #pragma HLS inline
    #pragma HLS pipeline

    if (reset)
    {
    Reset_to_default: for (size_t i = 0; i < 3; ++i)
        {
            #pragma HLS unroll
            shift_reg_spar[i] = 3.14;
            // shift_reg_data[i] = NULL_DATA;
            make_null_data(&shift_reg_data[i]);
        }
    }

    // spar_c_t val_for_spar;
    // get_sort_par(data_in, &val_for_spar);

    // spar_c_t val_for_spar_target;
    // get_sort_par(target, &val_for_spar_target);

    // spar_c_t spar_in;
    // abs_delta(val_for_spar, val_for_spar_target, &spar_in);

    spar_c_t spar_in;
    abs_delta(data_in.phi, target.phi, &spar_in);

    bool gtr_0 = (spar_in < shift_reg_spar[0]);
    bool gtr_1 = (spar_in < shift_reg_spar[1]);
    bool gtr_2 = (spar_in < shift_reg_spar[2]);

    // #ifndef __SYNTHESIS__
    // if (instance == 0)
    //     std::cout << "---------- instance " << instance << " enters " << data_in.phi << " target " << target.phi << " results " << gtr_0 << " " << gtr_1 << " " << gtr_2 << std::endl;
    // #endif

    if (gtr_0)
    {
        // shift_reg_data[2] = shift_reg_data[1];
        // shift_reg_data[1] = shift_reg_data[0];
        // shift_reg_data[0] = data_in;
        copy_data(shift_reg_data[1], &shift_reg_data[2]);
        copy_data(shift_reg_data[0], &shift_reg_data[1]);
        copy_data(data_in, &shift_reg_data[0]);

        shift_reg_spar[2] = shift_reg_spar[1];
        shift_reg_spar[1] = shift_reg_spar[0];
        shift_reg_spar[0] = spar_in;
    }

    else if (gtr_1)
    {
        // shift_reg_data[2] = shift_reg_data[1];
        // shift_reg_data[1] = data_in;
        copy_data(shift_reg_data[1], &shift_reg_data[2]);
        copy_data(data_in, &shift_reg_data[1]);

        shift_reg_spar[2] = shift_reg_spar[1];
        shift_reg_spar[1] = spar_in;
    }

    else if (gtr_2)
    {
        // shift_reg_data[2] = data_in;
        copy_data(data_in, &shift_reg_data[2]);
        shift_reg_spar[2] = spar_in;
    }

    // *data_out_0 = shift_reg_data[0];
    // *data_out_1 = shift_reg_data[1];
    // *data_out_2 = shift_reg_data[2];

    copy_data(shift_reg_data[0], data_out_0);
    copy_data(shift_reg_data[1], data_out_1);
    copy_data(shift_reg_data[2], data_out_2);
}