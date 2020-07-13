#include <cstdio>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "src/isolation.h"
// #include "src/dataformats.h"
#include "src/dataformats_v2.h"
#include "src/external_constants.h"
#include "ap_int.h"
#include "ap_fixed.h"

// #define NTEST 3000 // input pattern has 3035
#define NSKIP 0 
#define NTEST 100 // for a quick test

// #define NTEST 1
// #define NSKIP 11 // number of events to skip from the beginning

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

void dump_muon(muon_t mu, std::string prefix = "")
{
    if (prefix.size() > 0)
        cout << prefix << " ";
    cout << "muon pt : " << mu.pt.to_double()
         << "  eta : "   << mu.eta.to_double()
         << "  phi : "   << mu.phi.to_double()
         << "  z0 : "   << mu.z0.to_double()
    << endl;
}

void dump_track(track_t trk, std::string prefix = "")
{
    if (prefix.size() > 0)
        cout << prefix << " ";
    cout << "track pt : " << trk.pt.to_double()
         << "  eta : "   << trk.eta.to_double()
         << "  phi : "   << trk.phi.to_double()
         << "  z0 : "   << trk.z0.to_double()
    << endl;
}

// return integer that expresses number in lsb
// note   : sign is taken into account (e.g., an angle of -pi/2 on a LSB of pi/2^10 will return (-pi/2) / (pi/2^10) = -512)
// note 2 : works up tp 16 bits
long long int quantize (double number, double lsb)
{
    double div = number/lsb;
    long long int result = (long long int) floor(div);
    
    // cross check
    if (std::abs(result*lsb - number) > lsb)
        cout << "Quantization error! " << number << " in units of " << lsb << " yields " << result << endl;

    return result;
}

void dump_results(int ievt, muon_data_t in_mu, muon_isodata_t iso_flags, std::ofstream& fout, bool drop_invalid = true)
{
    std::vector<muon_t> muons = {
        in_mu.mu_0,
        in_mu.mu_1,
        in_mu.mu_2,
        in_mu.mu_3,
        in_mu.mu_4,
        in_mu.mu_5,
        in_mu.mu_6,
        in_mu.mu_7,
        in_mu.mu_8,
        in_mu.mu_9,
        in_mu.mu_10,
        in_mu.mu_11,
    };

    std::vector<hw_iso_t> isoflags = {
        iso_flags.isomu_0,
        iso_flags.isomu_1,
        iso_flags.isomu_2,
        iso_flags.isomu_3,
        iso_flags.isomu_4,
        iso_flags.isomu_5,
        iso_flags.isomu_6,
        iso_flags.isomu_7,
        iso_flags.isomu_8,
        iso_flags.isomu_9,
        iso_flags.isomu_10,
        iso_flags.isomu_11,
    };


    for (size_t i = 0; i < muons.size(); ++i)
    {
        muon_t     thismu  = muons.at(i);
        hw_iso_t   thisiso = isoflags.at(i);
        if (drop_invalid && thismu.pt.to_double() < 1)
            continue;
        
        fout << ievt
             << " " << i << " " << thismu.pt.to_double() << " " << thismu.eta.to_double() << " " << thismu.phi.to_double() << " " << thismu.z0.to_double() << " "
             << thisiso.to_uint() << endl;
    }    
}

void convert_and_dump_results(int ievt, muon_data_t in_mu, muon_isodata_t iso_flags, std::ofstream& fout, 
    double lsb_pt, double lsb_eta, double lsb_phi, double lsb_z0,
    bool drop_invalid = true)
{
    std::vector<muon_t> muons = {
        in_mu.mu_0,
        in_mu.mu_1,
        in_mu.mu_2,
        in_mu.mu_3,
        in_mu.mu_4,
        in_mu.mu_5,
        in_mu.mu_6,
        in_mu.mu_7,
        in_mu.mu_8,
        in_mu.mu_9,
        in_mu.mu_10,
        in_mu.mu_11,
    };

    std::vector<hw_iso_t> isoflags = {
        iso_flags.isomu_0,
        iso_flags.isomu_1,
        iso_flags.isomu_2,
        iso_flags.isomu_3,
        iso_flags.isomu_4,
        iso_flags.isomu_5,
        iso_flags.isomu_6,
        iso_flags.isomu_7,
        iso_flags.isomu_8,
        iso_flags.isomu_9,
        iso_flags.isomu_10,
        iso_flags.isomu_11,
    };


    for (size_t i = 0; i < muons.size(); ++i)
    {
        muon_t     thismu  = muons.at(i);
        hw_iso_t   thisiso = isoflags.at(i);
        if (drop_invalid && thismu.pt.to_double() < 1)
            continue;
        
        fout << ievt
             << " " << i << " " << thismu.pt.to_double()*lsb_pt << " " << thismu.eta.to_double()*lsb_eta << " " << thismu.phi.to_double()*lsb_phi << " " << thismu.z0.to_double()*lsb_z0 << " "
             << thisiso.to_uint() << endl;
    }    
}

void convert_and_dump_results_fulldata(int ievt, isomuon_data_t out_muons, std::ofstream& fout, 
    double lsb_pt, double lsb_eta, double lsb_phi, double lsb_z0,
    bool drop_invalid = true)
{
    std::vector<muon_t> muons = {
        out_muons.mu_0.mu,
        out_muons.mu_1.mu,
        out_muons.mu_2.mu,
        out_muons.mu_3.mu,
        out_muons.mu_4.mu,
        out_muons.mu_5.mu,
        out_muons.mu_6.mu,
        out_muons.mu_7.mu,
        out_muons.mu_8.mu,
        out_muons.mu_9.mu,
        out_muons.mu_10.mu,
        out_muons.mu_11.mu,
    };

    std::vector<hw_iso_t> isoflags = {
        out_muons.mu_0.isoflags,
        out_muons.mu_1.isoflags,
        out_muons.mu_2.isoflags,
        out_muons.mu_3.isoflags,
        out_muons.mu_4.isoflags,
        out_muons.mu_5.isoflags,
        out_muons.mu_6.isoflags,
        out_muons.mu_7.isoflags,
        out_muons.mu_8.isoflags,
        out_muons.mu_9.isoflags,
        out_muons.mu_10.isoflags,
        out_muons.mu_11.isoflags,
    };


    for (size_t i = 0; i < muons.size(); ++i)
    {
        muon_t     thismu  = muons.at(i);
        hw_iso_t   thisiso = isoflags.at(i);
        if (drop_invalid && thismu.pt.to_double() < 1)
            continue;
        
        fout << ievt
             << " " << i << " " << thismu.pt.to_double()*lsb_pt << " " << thismu.eta.to_double()*lsb_eta << " " << thismu.phi.to_double()*lsb_phi << " " << thismu.z0.to_double()*lsb_z0 << " "
             << thisiso.to_uint() << endl;
    }    
}

int main()
{
    std::vector<std::vector<double> > in_data = read_in_data("/home/lcadamur/MuL1HLSFW/MuonAlgorithms/patterns/mu_track_infolist_PU200.txt");

    std::string ofile = "/home/lcadamur/MuL1HLSFW/MuonAlgorithms/Isolation/isolation_results.txt";
    std::ofstream fout (ofile);
    int nwritten = 0;


    const int nmu         = N_MUON;      // 12; // NB : must be aligned with N_MUON in the dataformats
    const int ntrklinks   = N_TRK_LINKS; // 18; // NB : must be aligned with N_TRK_LINKS in the dataformats
    const int ntrkperlink = 100;         // used to read data input and raise "last" flag -> must be in synch with the input patterns

    // const int ntrkperlink_to_send = 5; // must be <= to ntrkperlink, to send a subset of the read tracks (for faster debug, but will not give meaningful performance)
    // const int ntrkperlink_to_send = 18; // must be <= to ntrkperlink, to send a subset of the read tracks (for faster debug, but will not give meaningful performance)
    const int ntrkperlink_to_send = 100; // must be <= to ntrkperlink, to send a subset of the read tracks (for faster debug, but will not give meaningful performance)

    const double lsb_pt  = 512./pow(2, hw_pt_t::width);
    const double lsb_eta = 2.*M_PI/pow(2, hw_eta_t::width);
    const double lsb_phi = 2.*M_PI/pow(2, hw_phi_t::width);
    const double lsb_z0  = 60./pow(2, hw_z0_t::width);

    cout << "---- LSB considered ---" << endl;
    cout << "lsb_pt   = " << lsb_pt   << endl;
    cout << "lsb_eta  = " << lsb_eta  << endl;
    cout << "lsb_phi  = " << lsb_phi  << endl;
    cout << "lsb_z0   = " << lsb_z0   << endl;

    for (uint itest = 0 + NSKIP; itest < NTEST + NSKIP; ++itest)
    {
        std::cout << " ........ THIS IS THE TEST NUMBER .......... " << itest << std::endl;
 
        if (ntrkperlink_to_send > ntrkperlink)
            throw std::runtime_error("too many tracks were set to be sent");

        std::vector<double> pts_mu (nmu);
        std::vector<double> pts_trks (ntrklinks*ntrkperlink);

        std::vector<double> etas_mu (nmu);
        std::vector<double> etas_trks (ntrklinks*ntrkperlink);

        std::vector<double> phis_mu (nmu);
        std::vector<double> phis_trks (ntrklinks*ntrkperlink);

        std::vector<double> z0s_mu (nmu);
        std::vector<double> z0s_trks (ntrklinks*ntrkperlink);

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
            muons[imu].pt  = quantize (pts_mu.at(imu),  lsb_pt);
            muons[imu].eta = quantize (etas_mu.at(imu), lsb_eta);
            muons[imu].phi = quantize (phis_mu.at(imu), lsb_phi);
            muons[imu].z0  = quantize (z0s_mu.at(imu),  lsb_z0);

            // muons[imu].pt  = pts_mu.at(imu);
            // muons[imu].eta = etas_mu.at(imu);
            // muons[imu].phi = phis_mu.at(imu);
            // muons[imu].z0  = z0s_mu.at(imu);
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

        #if ISODEBUG
        cout << "--- tb : muons of this event ---" << endl;
        dump_muon(mudata.mu_0, "imu = 0");
        #endif

        // for (uint iclk = 0; iclk < ntrkperlink; ++iclk) // to send all trks
        for (uint iclk = 0; iclk < ntrkperlink_to_send; ++iclk) // to just send a subset of nclks for faster debug
        {
            track_t tracks[ntrklinks];
            for (size_t itrk = 0; itrk < ntrklinks; ++itrk)
            {
                tracks[itrk].pt  = quantize (pts_trks.at(itrk  + ntrklinks*iclk), lsb_pt);
                tracks[itrk].eta = quantize (etas_trks.at(itrk + ntrklinks*iclk), lsb_eta);
                tracks[itrk].phi = quantize (phis_trks.at(itrk + ntrklinks*iclk), lsb_phi);
                tracks[itrk].z0  = quantize (z0s_trks.at(itrk  + ntrklinks*iclk), lsb_z0);

                // tracks[itrk].pt  = pts_trks.at(itrk  + ntrklinks*iclk);
                // tracks[itrk].eta = etas_trks.at(itrk + ntrklinks*iclk);
                // tracks[itrk].phi = phis_trks.at(itrk + ntrklinks*iclk);
                // tracks[itrk].z0  = z0s_trks.at(itrk  + ntrklinks*iclk);
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

            #if ISODEBUG
            cout << "--- tb : tracks of iclk " << iclk << " ---" << endl;
            dump_track(tracks_data.trk_0, "itrk = 0");
            dump_track(tracks_data.trk_1, "itrk = 1");
            dump_track(tracks_data.trk_2, "itrk = 2");
            dump_track(tracks_data.trk_3, "itrk = 3");
            dump_track(tracks_data.trk_4, "itrk = 4");
            dump_track(tracks_data.trk_5, "itrk = 5");
            dump_track(tracks_data.trk_6, "itrk = 6");
            dump_track(tracks_data.trk_7, "itrk = 7");
            dump_track(tracks_data.trk_8, "itrk = 8");
            dump_track(tracks_data.trk_9, "itrk = 9");
            dump_track(tracks_data.trk_10, "itrk = 10");
            dump_track(tracks_data.trk_11, "itrk = 11");
            dump_track(tracks_data.trk_12, "itrk = 12");
            dump_track(tracks_data.trk_13, "itrk = 13");
            dump_track(tracks_data.trk_14, "itrk = 14");
            dump_track(tracks_data.trk_15, "itrk = 15");
            dump_track(tracks_data.trk_16, "itrk = 16");
            dump_track(tracks_data.trk_17, "itrk = 17");
            #endif

            track_data_9_t tracks_data_9_p1;
            track_data_9_t tracks_data_9_p2;

            tracks_data_9_p1.trk_0 = tracks_data.trk_0;
            tracks_data_9_p1.trk_1 = tracks_data.trk_1;
            tracks_data_9_p1.trk_2 = tracks_data.trk_2;
            tracks_data_9_p1.trk_3 = tracks_data.trk_3;
            tracks_data_9_p1.trk_4 = tracks_data.trk_4;
            tracks_data_9_p1.trk_5 = tracks_data.trk_5;
            tracks_data_9_p1.trk_6 = tracks_data.trk_6;
            tracks_data_9_p1.trk_7 = tracks_data.trk_7;
            tracks_data_9_p1.trk_8 = tracks_data.trk_8;
            
            tracks_data_9_p2.trk_0 = tracks_data.trk_9;
            tracks_data_9_p2.trk_1 = tracks_data.trk_10;
            tracks_data_9_p2.trk_2 = tracks_data.trk_11;
            tracks_data_9_p2.trk_3 = tracks_data.trk_12;
            tracks_data_9_p2.trk_4 = tracks_data.trk_13;
            tracks_data_9_p2.trk_5 = tracks_data.trk_14;
            tracks_data_9_p2.trk_6 = tracks_data.trk_15;
            tracks_data_9_p2.trk_7 = tracks_data.trk_16;
            tracks_data_9_p2.trk_8 = tracks_data.trk_17;


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

            // ap_uint<1> is_last = (iclk == ntrkperlink_to_send - 1 ? 1 : 0);
            // // ap_uint<1> iso_flag = isolation_single_wrap(muons[0], tracks, is_last);
            // ap_uint<1> iso_flag = isolation_single_muon_wrap(muons[0], tracks_data, is_last);

            // if (is_last)
            // {
            //     std::cout << "... isolation results " << iso_flag.to_uint() << std::endl;
            // }

            iso_accum_t iso_threshold_1 = c_iso_sumpt_thr_1;
            iso_accum_t iso_threshold_2 = c_iso_sumpt_thr_2;
            
            // muon_isodata_t iso_flags;
            isomuon_data_t iso_muons;

            ap_uint<1> is_last = (iclk == ntrkperlink_to_send - 1 ? 1 : 0);            
            
            //-------------------- 18 tracks per link
            #if RUN_DOUBLE_CLK_SPEED == 0
            // isolation_allmu(mudata, tracks_data, is_last, iso_threshold, iso_flags);
            isolation_allmu(mudata, tracks_data, is_last, iso_threshold_1, iso_threshold_2, iso_muons);
            //--------------------------------------------------------------


            // //-------------------- 9 tracks per link
            #else
            // isolation_allmu_9trk(mudata, tracks_data_9_p1, 0, iso_threshold, iso_flags);
            // isolation_allmu_9trk(mudata, tracks_data_9_p2, is_last, iso_threshold, iso_flags);
            isolation_allmu_9trk(mudata, tracks_data_9_p1, 0,       iso_threshold_1, iso_threshold_2, iso_muons);
            isolation_allmu_9trk(mudata, tracks_data_9_p2, is_last, iso_threshold_1, iso_threshold_2, iso_muons);
            #endif

            // //--------------------------------------------------------------

            if (is_last)
            {
                // std::cout << "... isolation results mu 0 " << iso_flags.isomu_0.to_uint() << std::endl;
                std::cout << "... isolation results mu 0 " << iso_muons.mu_0.isoflags.to_uint() << std::endl;
                
                // dump_results(itest, mudata, iso_flags, fout);
                // convert_and_dump_results(itest, mudata, iso_flags, fout, lsb_pt, lsb_eta, lsb_phi, lsb_z0);
                convert_and_dump_results_fulldata(itest, iso_muons, fout, lsb_pt, lsb_eta, lsb_phi, lsb_z0);
            }

        }
    }
}