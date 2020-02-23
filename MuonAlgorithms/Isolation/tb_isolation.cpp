#include <cstdio>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include "src/isolation.h"
#include "src/dataformats.h"
#include "ap_int.h"
#include "ap_fixed.h"

// #define NTEST 3000 // input pattern has 3035
#define NTEST 100 // input pattern has 3035

int roundDouble (double x){
    int r = (int) (x + 0.5);
    return r;
}

std::vector<std::vector<double> > read_in_data (std::string fin_name)
{
    std::cout << " ... reading from " << fin_name << std::endl;

    std::fstream fin (fin_name.c_str());
    std::string line;
    int nlines = 0;

    std::vector<std::vector<double> > read_lines;
    
    while (std::getline(fin, line))
    {
        read_lines.push_back(std::vector<double>(0));
        std::istringstream iss(line);
        double buf; // all input data are integers
        while (iss >> buf)
        {
            read_lines.back().push_back(buf);
        }
        ++nlines;
    }
    
    std::cout << " ... read " << nlines << " lines" << std::endl;
    return read_lines;
}

int main()
{
    std::vector<std::vector<double> > in_data = read_in_data("/home/lcadamur/MuL1HLSFW/MuonAlgorithms/patterns/mu_track_infolist_PU200.txt");

    std::string ofile = "/home/lcadamur/MuL1HLSFW/MuonAlgorithms/Isolation/isolation_results.txt";
    std::ofstream fout (ofile);
    int nwritten = 0;

    for (uint itest = 0; itest < NTEST; ++itest)
    {

        std::cout << " ........ THIS IS THE TEST NUMBER .......... " << itest << std::endl;
 
        const int nmu         = N_MUON;      // 12; // NB : must be aligned with N_MUON in the dataformats
        const int ntrklinks   = N_TRK_LINKS; // 18; // NB : must be aligned with N_TRK_LINKS in the dataformats
        const int ntrkperlink = 100;         // used to read data input and raise "last" flag -> must be in synch with the input patterns

        const int ntrkperlink_to_send = 18; // must be <= to ntrkperlink, to send a subset of the read tracks (for faster debug, but will not give meaningful performance)

        if (ntrkperlink_to_send > ntrkperlink)
            throw std::runtime_error("too many tracks were set to be sent");

        std::vector<int> pts_mu (nmu);
        std::vector<int> pts_trks (ntrklinks*ntrkperlink);

        std::vector<int> etas_mu (nmu);
        std::vector<int> etas_trks (ntrklinks*ntrkperlink);

        std::vector<int> phis_mu (nmu);
        std::vector<int> phis_trks (ntrklinks*ntrkperlink);

        std::vector<int> z0s_mu (nmu);
        std::vector<int> z0s_trks (ntrklinks*ntrkperlink);

       const int nobj_per_ev = nmu + ntrklinks*ntrkperlink; // all track and muons from a collision

        // trk_obj  = ['L1TT_trk_pt', 'L1TT_trk_eta', 'L1TT_trk_phi', 'L1TT_trk_z',  'L1TT_trk_charge', 'L1TT_trk_chi2']
        // muon_obj = ['L1_TkMu_pt',  'L1_TkMu_eta',  'L1_TkMu_phi',  'L1_TkMu_z',   'L1_TkMu_charge']

       // with a offset of 3 for iEv, itype [mu=0, trk=1], iobj

        for (uint i = 0; i < pts_mu.size(); ++i){
            pts_mu.at(i)  = in_data.at(itest*nobj_per_ev + i).at(3+0);
            etas_mu.at(i) = in_data.at(itest*nobj_per_ev + i).at(3+1);
            phis_mu.at(i) = in_data.at(itest*nobj_per_ev + i).at(3+2);
            z0s_mu.at(i)  = in_data.at(itest*nobj_per_ev + i).at(3+3);
        }
        for (uint i = 0; i < pts_trks.size(); ++i){
            pts_trks.at(i)  = in_data.at(itest*nobj_per_ev + nmu + i).at(3+0);
            etas_trks.at(i) = in_data.at(itest*nobj_per_ev + nmu + i).at(3+1);
            phis_trks.at(i) = in_data.at(itest*nobj_per_ev + nmu + i).at(3+2);
            z0s_trks.at(i)  = in_data.at(itest*nobj_per_ev + nmu + i).at(3+3);
        }
 
        // verify that data are read correctly
        int iEv_mufirst = in_data.at(itest*nobj_per_ev + 0).at(0);
        int iEv_mulast  = in_data.at(itest*nobj_per_ev + pts_mu.size() -1).at(0);

        int iEv_trkfirst = in_data.at(itest*nobj_per_ev + nmu + 0).at(0);
        int iEv_trklast  = in_data.at(itest*nobj_per_ev + nmu + pts_trks.size()-1).at(0);

        if (! (iEv_mufirst == iEv_mulast  && iEv_mulast == iEv_trkfirst && iEv_trkfirst ==  iEv_trklast) )
            throw std::runtime_error(std::string("data misread"));

        // build the fw objects from the read data
        muon_t muons[nmu];
        for (size_t imu = 0; imu < nmu; ++imu)
        {
            muons[imu].pt  = pts_mu.at(imu);
            muons[imu].eta = etas_mu.at(imu);
            muons[imu].phi = phis_mu.at(imu);
            muons[imu].z0  = z0s_mu.at(imu);
        }

        muon_data_t mudata;
        mudata.mu_0 = muons[0];
        mudata.mu_1 = muons[1];
        mudata.mu_2 = muons[2];
        mudata.mu_3 = muons[3];
        mudata.mu_4 = muons[4];
        mudata.mu_5 = muons[5];
        mudata.mu_6 = muons[6];
        mudata.mu_7 = muons[7];
        mudata.mu_8 = muons[8];
        mudata.mu_9 = muons[9];
        mudata.mu_10 = muons[10];
        mudata.mu_11 = muons[11];

        // for (uint iclk = 0; iclk < ntrkperlink; ++iclk) // to send all trks
        for (uint iclk = 0; iclk < ntrkperlink_to_send; ++iclk) // to just send a subset of nclks for faster debug
        {
            track_t tracks[ntrklinks];
            for (size_t itrk = 0; itrk < ntrklinks; ++itrk)
            {
                tracks[itrk].pt  = pts_trks.at(itrk  + ntrklinks*iclk);
                tracks[itrk].eta = etas_trks.at(itrk + ntrklinks*iclk);
                tracks[itrk].phi = phis_trks.at(itrk + ntrklinks*iclk);
                tracks[itrk].z0  = z0s_trks.at(itrk  + ntrklinks*iclk);
            }

            track_data_t tracks_data;

            tracks_data.trk_0 = tracks[0];
            tracks_data.trk_1 = tracks[1];
            tracks_data.trk_2 = tracks[2];
            tracks_data.trk_3 = tracks[3];
            tracks_data.trk_4 = tracks[4];
            tracks_data.trk_5 = tracks[5];
            tracks_data.trk_6 = tracks[6];
            tracks_data.trk_7 = tracks[7];
            tracks_data.trk_8 = tracks[8];
            tracks_data.trk_9 = tracks[9];
            tracks_data.trk_10 = tracks[10];
            tracks_data.trk_11 = tracks[11];
            tracks_data.trk_12 = tracks[12];
            tracks_data.trk_13 = tracks[13];
            tracks_data.trk_14 = tracks[14];
            tracks_data.trk_15 = tracks[15];
            tracks_data.trk_16 = tracks[16];
            tracks_data.trk_17 = tracks[17];


            // ap_uint<1> iso_flags[nmu];
            // ap_uint<1> is_last = (iclk == ntrkperlink - 1 ? 1 : 0);
            // isolation(muons, tracks, is_last, iso_flags);

            // if (is_last)
            // {
            //     std::cout << "... isolation results" << std::endl;
            //     for (size_t imu = 0; imu < nmu; ++imu)
            //     {
            //         std::cout << " imu : " << imu << " " << iso_flags[imu].to_uint() << std::endl;
            //     }
            // }

            ap_uint<1> is_last = (iclk == ntrkperlink_to_send - 1 ? 1 : 0);
            // ap_uint<1> iso_flag = isolation_single_wrap(muons[0], tracks, is_last);
            // ap_uint<1> iso_flag = isolation_single_muon_wrap(muons[0], tracks_data, is_last);

            // if (is_last)
            // {
            //     std::cout << "... isolation results " << iso_flag.to_uint() << std::endl;
            // }

            muon_isodata_t iso_flags;
            isolation_allmu(mudata, tracks_data, is_last, iso_flags);

            if (is_last)
            {
                std::cout << "... isolation results mu 0 " << iso_flags.isomu_0.to_uint() << std::endl;
            }

        }
    }
}