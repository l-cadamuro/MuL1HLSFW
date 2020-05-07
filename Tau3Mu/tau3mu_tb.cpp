#include <stdlib.h>
#include <iostream>
#include <vector>
#include <utility>
#include "src/mini_sorter.h"

#define NTEST 100

std::vector<hw_phi_t> cluster_tb(std::vector<muon_t> data)
{
    std::vector<hw_phi_t> result;
    for (uint i = 0; i < data.size(); ++i)
    {
        std::vector <std::pair<hw_phi_t, muon_t>> deltas_and_mus; 
        for (uint j = 0; j < data.size(); ++j){
            hw_phi_t delta = data.at(i).phi - data.at(j).phi;
            if (delta < 0) delta = -1*delta;
            deltas_and_mus.push_back(std::make_pair(delta, data.at(j)));
        }
        sort(deltas_and_mus.begin(), deltas_and_mus.end(), [](const std::pair<hw_phi_t, muon_t> &a, const std::pair<hw_phi_t, muon_t> &b) -> bool {return a.first < b.first;});
        result.push_back(deltas_and_mus.at(2).second.phi);
    }
    return result;
}

int main()
{
    srand(123456);

    std::vector<muon_t> test_data [NTEST];
    for (uint i = 0; i < NTEST; ++i) {
        for (uint j = 0; j < N_MU; ++j){
            float r2 = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / 3.14)); // uniform 0 - 3.14
            muon_t mu;
            mu.phi = r2;
            test_data[i].push_back(mu);
        }
    }

    int n_errors = 0;
    for (int itest = 0; itest < NTEST; ++itest)
    {
        muon_t data_in[N_MU];
        for (int i = 0; i < N_MU; ++i){
            data_in[i] = test_data[itest][i];
        }

        // debug the input test data
        std::cout << "------ TEST : " << itest << std::endl;
        for (int i = 0; i < N_MU; ++i){
            std::cout << "- imu : " << i << " phi = " << data_in[i].phi.to_float() << std::endl;
        }

        hw_minv_t data_out[N_MU];
        tau_3mu(data_in, data_out);

        std::vector<muon_t> data_in_v;
        for (int i = 0; i < N_MU; ++i)
            data_in_v.push_back(data_in[i]);
        auto exp_result = cluster_tb(data_in_v);

        // debug the answer test data
        std::cout << "------ TEST : " << itest << std::endl;
        for (int i = 0; i < N_MU; ++i)
        {
            std::cout << "- imu : " << i << " phi = " << data_in[i].phi.to_float() << " ---- gives " << data_out[i] << " .. expect " << exp_result.at(i) << std::endl;
            if (exp_result.at(i) != data_out[i]){
                n_errors += 1;
                std::cout << " !!! WARNING : found an error" << std::endl;
            }
        }
    }

    if (n_errors == 0)
        return 0;
    else
        return 1;
}