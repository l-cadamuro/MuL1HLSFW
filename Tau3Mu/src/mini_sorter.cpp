#include "mini_sorter.h" // with the template function
// #include "miniSorter.h"  // by hand

// #ifndef __SYNTHESIS__
// #include <iostream>
// #endif

typedef muon_t data_t;
typedef hw_phi_t spar_t; // the sorting parameter used

void get_sort_par(data_t data_in, spar_t *spar_out)
{
    #pragma HLS inline
    *spar_out = data_in.phi;
}

void abs_delta(spar_t p1, spar_t p2, spar_t *pdelta)
{
    #pragma HLS inline

    spar_t delta = p1 - p2;
    if (delta > 0)
        *pdelta = delta;
    else
        *pdelta = -1 * delta;
}

void make_null_data(data_t* data)
{
    #pragma HLS inline
    (*data).phi = 0;
}

void copy_data (data_t from, data_t* to)
{
    #pragma HLS inline
    (*to).phi = from.phi;
}

// templated to create multiple copies
template <size_t instance>
void mini_closer_elem(data_t data_in, data_t target, ap_uint<1> reset, data_t *data_out_0, data_t *data_out_1, data_t *data_out_2)
{
    #pragma HLS pipeline
    // #pragma HLS inline

    // #pragma HLS interface ap_none port=data_out_0
    // #pragma HLS interface ap_none port=data_out_1
    // #pragma HLS interface ap_none port=data_out_2

    static spar_t shift_reg_spar[3];// = {MAX_SPAR, MAX_SPAR, MAX_SPAR};    // try reset value here
    static data_t shift_reg_data[3];// = {NULL_DATA, NULL_DATA, NULL_DATA}; // try reset value here

    // static spar_t shift_reg_spar[3] = {3.14, 3.14, 3.14};    // try reset value here
    // static data_t shift_reg_data[3];
    // // static data_t shift_reg_data[3] = {NULL_DATA, NULL_DATA, NULL_DATA}; // try reset value here

    #pragma HLS array_partition variable = shift_reg_spar complete
    #pragma HLS array_partition variable = shift_reg_data complete

    spar_t spar_in;
    abs_delta(data_in.phi, target.phi, &spar_in);

    if (reset) {
        Reset_to_default : for (size_t i = 1; i < 3; ++i){
            #pragma HLS unroll
            shift_reg_spar[i] = 3.14;
            // shift_reg_data[i] = NULL_DATA;
            make_null_data(&shift_reg_data[i]);
        }
        
        // the first will always go to the queue since other data are invalid
        shift_reg_spar[0] = spar_in;
        copy_data(data_in, &shift_reg_data[0]);
    }

    // spar_t val_for_spar;
    // get_sort_par(data_in, &val_for_spar);

    // spar_t val_for_spar_target;
    // get_sort_par(target, &val_for_spar_target);

    // spar_t spar_in;
    // abs_delta(val_for_spar, val_for_spar_target, &spar_in);

    else
    {    

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

    }

    copy_data(shift_reg_data[0], data_out_0);
    copy_data(shift_reg_data[1], data_out_1);
    copy_data(shift_reg_data[2], data_out_2);

    // if (reset)
    // {
    //     Reset_to_default: for (size_t i = 0; i < 3; ++i)
    //     {
    //         #pragma HLS unroll
    //         shift_reg_spar[i] = 3.14;
    //         // shift_reg_data[i] = NULL_DATA;
    //         make_null_data(&shift_reg_data[i]);
    //     }
    // }
}


void find_clusters(muon_t in_muons[N_MU], hw_minv_t out_masses[N_MU])
{

    #pragma HLS array_partition variable = in_muons complete
    #pragma HLS array_partition variable = out_masses complete

    // the 3 grouped muons are stored here, all fully partitioned (check resources)
    muon_t groups[N_MU][3];
    #pragma HLS array_partition variable = groups complete dim = 2
    #pragma HLS array_partition variable = groups complete dim = 1

    Process_in_mu : for (size_t i = 0; i < N_MU; ++i)
    {
        #pragma HLS pipeline

        ap_uint<1> reset;
        if (i == 0)
            reset = 1;
        else
            reset = 0;

        // ap_uint<1> reset;
        // if (i == N_MU - 1)
        //     reset = 1;
        // else
        //     reset = 0;

        mini_closer_elem<0> (in_muons[i], in_muons[0],  reset, &groups[0][0],  &groups[0][1],  &groups[0][2]);
        mini_closer_elem<1> (in_muons[i], in_muons[1],  reset, &groups[1][0],  &groups[1][1],  &groups[1][2]);
        mini_closer_elem<2> (in_muons[i], in_muons[2],  reset, &groups[2][0],  &groups[2][1],  &groups[2][2]);
        mini_closer_elem<3> (in_muons[i], in_muons[3],  reset, &groups[3][0],  &groups[3][1],  &groups[3][2]);
        mini_closer_elem<4> (in_muons[i], in_muons[4],  reset, &groups[4][0],  &groups[4][1],  &groups[4][2]);
        mini_closer_elem<5> (in_muons[i], in_muons[5],  reset, &groups[5][0],  &groups[5][1],  &groups[5][2]);
        mini_closer_elem<6> (in_muons[i], in_muons[6],  reset, &groups[6][0],  &groups[6][1],  &groups[6][2]);
        mini_closer_elem<7> (in_muons[i], in_muons[7],  reset, &groups[7][0],  &groups[7][1],  &groups[7][2]);
        mini_closer_elem<8> (in_muons[i], in_muons[8],  reset, &groups[8][0],  &groups[8][1],  &groups[8][2]);
        mini_closer_elem<9> (in_muons[i], in_muons[9],  reset, &groups[9][0],  &groups[9][1],  &groups[9][2]);
        mini_closer_elem<10>(in_muons[i], in_muons[10], reset, &groups[10][0], &groups[10][1], &groups[10][2]);
        mini_closer_elem<11>(in_muons[i], in_muons[11], reset, &groups[11][0], &groups[11][1], &groups[11][2]);
    }

    // once done, send out the result - temporary : just sending out one property of the last muon
    Copy_to_output : for(size_t i = 0; i < N_MU; ++i)
    {
        #pragma HLS unroll
        out_masses[i] = groups[i][2].phi;

        // #ifndef __SYNTHESIS__
        // std::cout << "...... inside algo. Mu " << i << " " << groups[i][0].phi << " " << groups[i][1].phi << " " << groups[i][2].phi << std::endl;
        // #endif
    }
}

/*
void find_clusters_class(muon_t in_muons[N_MU], hw_minv_t out_masses[N_MU])
{
    #pragma HLS inline

    #pragma HLS array_partition variable = in_muons complete
    #pragma HLS array_partition variable = out_masses complete

    // the 3 grouped muons are stored here, all fully partitioned (check resources)
    muon_t groups[N_MU][3];
    #pragma HLS array_partition variable = groups complete dim = 2
    #pragma HLS array_partition variable = groups complete dim = 1

    static miniSorter ms[N_MU];
#pragma HLS array_partition variable = ms complete

Process_in_mu:
    for (size_t i = 0; i < N_MU; ++i)
    {
        #pragma HLS pipeline
        ap_uint<1> reset;
        if (i == 0)
            reset = 1;
        else
            reset = 0;

        ms[0].mini_closer_elem(in_muons[i], in_muons[0], reset, &groups[0][0], &groups[0][1], &groups[0][2]);
        ms[1].mini_closer_elem(in_muons[i], in_muons[1], reset, &groups[1][0], &groups[1][1], &groups[1][2]);
        ms[2].mini_closer_elem(in_muons[i], in_muons[2], reset, &groups[2][0], &groups[2][1], &groups[2][2]);
        ms[3].mini_closer_elem(in_muons[i], in_muons[3], reset, &groups[3][0], &groups[3][1], &groups[3][2]);
        ms[4].mini_closer_elem(in_muons[i], in_muons[4], reset, &groups[4][0], &groups[4][1], &groups[4][2]);
        ms[5].mini_closer_elem(in_muons[i], in_muons[5], reset, &groups[5][0], &groups[5][1], &groups[5][2]);
        ms[6].mini_closer_elem(in_muons[i], in_muons[6], reset, &groups[6][0], &groups[6][1], &groups[6][2]);
        ms[7].mini_closer_elem(in_muons[i], in_muons[7], reset, &groups[7][0], &groups[7][1], &groups[7][2]);
        ms[8].mini_closer_elem(in_muons[i], in_muons[8], reset, &groups[8][0], &groups[8][1], &groups[8][2]);
        ms[9].mini_closer_elem(in_muons[i], in_muons[9], reset, &groups[9][0], &groups[9][1], &groups[9][2]);
        ms[10].mini_closer_elem(in_muons[i], in_muons[10], reset, &groups[10][0], &groups[10][1], &groups[10][2]);
        ms[11].mini_closer_elem(in_muons[i], in_muons[11], reset, &groups[11][0], &groups[11][1], &groups[11][2]);
    }

// once done, send out the result - temporary : just sending out one property of the last muon
    Copy_to_output: for (size_t i = 0; i < N_MU; ++i)
    {
        #pragma HLS unroll
        out_masses[i] = groups[i][2].phi;

        // #ifndef __SYNTHESIS__
        // std::cout << "...... inside algo. Mu " << i << " " << groups[i][0].phi << " " << groups[i][1].phi << " " << groups[i][2].phi << std::endl;
        // #endif
    }
}
*/

// void tau_3mu(muon_t in_muons[N_MU], hw_minv_t out_masses[N_MU])
// {
//     find_clusters(in_muons, out_masses);
//     // find_clusters_class(in_muons, out_masses);
//     // find_clusters_unroll(in_muons, out_masses);
// }
