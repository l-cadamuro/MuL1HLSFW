#include <cstdio>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include "src/isolation.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 6500 // input pattern has 6552
#define NTMT 18

#define VERBOSE false

int roundDouble (double x){
    int r = (int) (x + 0.5);
    return r;
}

std::vector<std::vector<int> > read_in_data (std::string fin_name)
{
    std::cout << " ... reading from " << fin_name << std::endl;

    std::fstream fin (fin_name.c_str());
    std::string line;
    int nlines = 0;

    std::vector<std::vector<int> > read_lines;
    
    while (std::getline(fin, line))
    {
        read_lines.push_back(std::vector<int>(0));
        std::istringstream iss(line);
        int buf; // all input data are integers
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

    // srand (123456);
    std::default_random_engine rndmEngine;
    std::uniform_real_distribution<double> distribution(10., 100.);

    std::vector<std::vector<int> > in_data = read_in_data("/home/zynq/luca/MuL1HLSFW/MuonAlgorithms/patterns/mu_track_infolist.txt");

    std::string ofile = "/home/zynq/luca/MuL1HLSFW/MuonAlgorithms/isolation_results.txt";
    std::ofstream fout (ofile);
    int nwritten = 0;

    // for (int il = 0; il < in_data.size(); ++il){
    //     std::cout << " LINE : " << il << " ";
    //     for (int iv = 0; iv < in_data.at(il).size(); ++iv)
    //         std::cout << in_data.at(il).at(iv) << " ";
    //     std::cout << std::endl;
    // }

    for (uint itest = 0; itest < NTEST; ++itest)
    {

        std::cout << " ........ THIS IS THE TEST NUMBER .......... " << itest << std::endl;
 
        // generate the muons for this test
        // I only care for the pT actually, so let's not bother to generate other variables
        std::vector<int> pts_mu (18); // there will be 18 muons
        std::vector<int> pts_trks (18*NTMT); // and 18 tracks sent for every TMT period

        std::vector<int> etas_mu (18); // there will be 18 muons
        std::vector<int> etas_trks (18*NTMT); // and 18 tracks sent for every TMT period

        std::vector<int> phis_mu (18); // there will be 18 muons
        std::vector<int> phis_trks (18*NTMT); // and 18 tracks sent for every TMT period

        // for (uint i = 0; i < pts_mu.size(); ++i){
        //     pts_mu.at(i) = roundDouble(distribution(rndmEngine));
        // }

        // for (uint i = 0; i < pts_trks.size(); ++i){
        //     pts_trks.at(i) = roundDouble(distribution(rndmEngine));
        // }

        // from patterns
        // the first 3 fields are iEv, itype, iobj
        // followed by
        // muon_obj = ['EMTF_mu_pt',  'EMTF_mu_eta', 'EMTF_mu_phi', 'EMTF_mu_charge']
        // trk_obj  = ['L1TT_trk_pt', 'L1TT_trk_eta', 'L1TT_trk_phi', 'L1TT_trk_charge', 'L1TT_trk_chi2']

        const int nobj_per_ev = 18 + 18*NTMT; // all track and muons from a collision

        for (uint i = 0; i < pts_mu.size(); ++i){
            pts_mu.at(i)  = in_data.at(itest*nobj_per_ev + i).at(3+0);
            etas_mu.at(i) = in_data.at(itest*nobj_per_ev + i).at(3+1);
            phis_mu.at(i) = in_data.at(itest*nobj_per_ev + i).at(3+2);
        }
        for (uint i = 0; i < pts_trks.size(); ++i){
            pts_trks.at(i)  = in_data.at(itest*nobj_per_ev + 18 + i).at(3+0);
            etas_trks.at(i) = in_data.at(itest*nobj_per_ev + 18 + i).at(3+1);
            phis_trks.at(i) = in_data.at(itest*nobj_per_ev + 18 + i).at(3+2);
        }

        track_conv_input trk_in;
        muon_input       mu_in;

        mu_in.mu0.pt = pts_mu.at(0);
        mu_in.mu1.pt = pts_mu.at(1);
        mu_in.mu2.pt = pts_mu.at(2);
        mu_in.mu3.pt = pts_mu.at(3);
        mu_in.mu4.pt = pts_mu.at(4);
        mu_in.mu5.pt = pts_mu.at(5);
        mu_in.mu6.pt = pts_mu.at(6);
        mu_in.mu7.pt = pts_mu.at(7);
        mu_in.mu8.pt = pts_mu.at(8);
        mu_in.mu9.pt = pts_mu.at(9);
        mu_in.mu10.pt = pts_mu.at(10);
        mu_in.mu11.pt = pts_mu.at(11);
        mu_in.mu12.pt = pts_mu.at(12);
        mu_in.mu13.pt = pts_mu.at(13);
        mu_in.mu14.pt = pts_mu.at(14);
        mu_in.mu15.pt = pts_mu.at(15);
        mu_in.mu16.pt = pts_mu.at(16);
        mu_in.mu17.pt = pts_mu.at(17);

        mu_in.mu0.eta = etas_mu.at(0);
        mu_in.mu1.eta = etas_mu.at(1);
        mu_in.mu2.eta = etas_mu.at(2);
        mu_in.mu3.eta = etas_mu.at(3);
        mu_in.mu4.eta = etas_mu.at(4);
        mu_in.mu5.eta = etas_mu.at(5);
        mu_in.mu6.eta = etas_mu.at(6);
        mu_in.mu7.eta = etas_mu.at(7);
        mu_in.mu8.eta = etas_mu.at(8);
        mu_in.mu9.eta = etas_mu.at(9);
        mu_in.mu10.eta = etas_mu.at(10);
        mu_in.mu11.eta = etas_mu.at(11);
        mu_in.mu12.eta = etas_mu.at(12);
        mu_in.mu13.eta = etas_mu.at(13);
        mu_in.mu14.eta = etas_mu.at(14);
        mu_in.mu15.eta = etas_mu.at(15);
        mu_in.mu16.eta = etas_mu.at(16);
        mu_in.mu17.eta = etas_mu.at(17);

        mu_in.mu0.phi = phis_mu.at(0);
        mu_in.mu1.phi = phis_mu.at(1);
        mu_in.mu2.phi = phis_mu.at(2);
        mu_in.mu3.phi = phis_mu.at(3);
        mu_in.mu4.phi = phis_mu.at(4);
        mu_in.mu5.phi = phis_mu.at(5);
        mu_in.mu6.phi = phis_mu.at(6);
        mu_in.mu7.phi = phis_mu.at(7);
        mu_in.mu8.phi = phis_mu.at(8);
        mu_in.mu9.phi = phis_mu.at(9);
        mu_in.mu10.phi = phis_mu.at(10);
        mu_in.mu11.phi = phis_mu.at(11);
        mu_in.mu12.phi = phis_mu.at(12);
        mu_in.mu13.phi = phis_mu.at(13);
        mu_in.mu14.phi = phis_mu.at(14);
        mu_in.mu15.phi = phis_mu.at(15);
        mu_in.mu16.phi = phis_mu.at(16);
        mu_in.mu17.phi = phis_mu.at(17);

        for (uint itmt = 0; itmt < NTMT; ++itmt)
        {
            trk_in.trk0.pt = pts_trks.at(0 + NTMT*itmt);
            trk_in.trk1.pt = pts_trks.at(1 + NTMT*itmt);
            trk_in.trk2.pt = pts_trks.at(2 + NTMT*itmt);
            trk_in.trk3.pt = pts_trks.at(3 + NTMT*itmt);
            trk_in.trk4.pt = pts_trks.at(4 + NTMT*itmt);
            trk_in.trk5.pt = pts_trks.at(5 + NTMT*itmt);
            trk_in.trk6.pt = pts_trks.at(6 + NTMT*itmt);
            trk_in.trk7.pt = pts_trks.at(7 + NTMT*itmt);
            trk_in.trk8.pt = pts_trks.at(8 + NTMT*itmt);
            trk_in.trk9.pt = pts_trks.at(9 + NTMT*itmt);
            trk_in.trk10.pt = pts_trks.at(10 + NTMT*itmt);
            trk_in.trk11.pt = pts_trks.at(11 + NTMT*itmt);
            trk_in.trk12.pt = pts_trks.at(12 + NTMT*itmt);
            trk_in.trk13.pt = pts_trks.at(13 + NTMT*itmt);
            trk_in.trk14.pt = pts_trks.at(14 + NTMT*itmt);
            trk_in.trk15.pt = pts_trks.at(15 + NTMT*itmt);
            trk_in.trk16.pt = pts_trks.at(16 + NTMT*itmt);
            trk_in.trk17.pt = pts_trks.at(17 + NTMT*itmt);            

            trk_in.trk0.eta = etas_trks.at(0 + NTMT*itmt);
            trk_in.trk1.eta = etas_trks.at(1 + NTMT*itmt);
            trk_in.trk2.eta = etas_trks.at(2 + NTMT*itmt);
            trk_in.trk3.eta = etas_trks.at(3 + NTMT*itmt);
            trk_in.trk4.eta = etas_trks.at(4 + NTMT*itmt);
            trk_in.trk5.eta = etas_trks.at(5 + NTMT*itmt);
            trk_in.trk6.eta = etas_trks.at(6 + NTMT*itmt);
            trk_in.trk7.eta = etas_trks.at(7 + NTMT*itmt);
            trk_in.trk8.eta = etas_trks.at(8 + NTMT*itmt);
            trk_in.trk9.eta = etas_trks.at(9 + NTMT*itmt);
            trk_in.trk10.eta = etas_trks.at(10 + NTMT*itmt);
            trk_in.trk11.eta = etas_trks.at(11 + NTMT*itmt);
            trk_in.trk12.eta = etas_trks.at(12 + NTMT*itmt);
            trk_in.trk13.eta = etas_trks.at(13 + NTMT*itmt);
            trk_in.trk14.eta = etas_trks.at(14 + NTMT*itmt);
            trk_in.trk15.eta = etas_trks.at(15 + NTMT*itmt);
            trk_in.trk16.eta = etas_trks.at(16 + NTMT*itmt);
            trk_in.trk17.eta = etas_trks.at(17 + NTMT*itmt);    

            trk_in.trk0.phi = phis_trks.at(0 + NTMT*itmt);
            trk_in.trk1.phi = phis_trks.at(1 + NTMT*itmt);
            trk_in.trk2.phi = phis_trks.at(2 + NTMT*itmt);
            trk_in.trk3.phi = phis_trks.at(3 + NTMT*itmt);
            trk_in.trk4.phi = phis_trks.at(4 + NTMT*itmt);
            trk_in.trk5.phi = phis_trks.at(5 + NTMT*itmt);
            trk_in.trk6.phi = phis_trks.at(6 + NTMT*itmt);
            trk_in.trk7.phi = phis_trks.at(7 + NTMT*itmt);
            trk_in.trk8.phi = phis_trks.at(8 + NTMT*itmt);
            trk_in.trk9.phi = phis_trks.at(9 + NTMT*itmt);
            trk_in.trk10.phi = phis_trks.at(10 + NTMT*itmt);
            trk_in.trk11.phi = phis_trks.at(11 + NTMT*itmt);
            trk_in.trk12.phi = phis_trks.at(12 + NTMT*itmt);
            trk_in.trk13.phi = phis_trks.at(13 + NTMT*itmt);
            trk_in.trk14.phi = phis_trks.at(14 + NTMT*itmt);
            trk_in.trk15.phi = phis_trks.at(15 + NTMT*itmt);
            trk_in.trk16.phi = phis_trks.at(16 + NTMT*itmt);
            trk_in.trk17.phi = phis_trks.at(17 + NTMT*itmt);    

            ap_uint<1> is_last = (itmt == NTMT - 1 ? 1 : 0);

            if (VERBOSE)
            {
                std::cout << " @@@ TMT # " << itmt << " / " << NTMT << " ---- ISLAST? " << is_last.to_uint() << std::endl;

                // some debug cout
                std::cout << " +++++ muon : pt " << mu_in.mu0.pt << "   .  eta = " << mu_in.mu0.eta << "   . phi = " << mu_in.mu0.phi << std::endl;
                std::cout << " +++++ trk0 : pt " << trk_in.trk0.pt << "   .  eta = " << trk_in.trk0.eta << "   . phi = " << trk_in.trk0.phi << std::endl;
                std::cout << " +++++ trk1 : pt " << trk_in.trk1.pt << "   .  eta = " << trk_in.trk1.eta << "   . phi = " << trk_in.trk1.phi << std::endl;
                std::cout << " +++++ trk17 : pt " << trk_in.trk17.pt << "   .  eta = " << trk_in.trk17.eta << "   . phi = " << trk_in.trk17.phi << std::endl;
            }

            // iso_muons_out outmus = test_algo(mu_in, trk_in, is_last);
            // ap_uint<1> isomu = isolation(mu_in.mu0, trk_in, is_last);
            // ap_uint<1> isomu = isolation_class(mu_in.mu0, trk_in, is_last);
            // ap_uint<1> isomu = isolation_class(mu_in.mu0, trk_in, is_last);
            
            iso_muon_input isomus;
            isomus = isolation_class_allmu(mu_in, trk_in, is_last);
            ap_uint<1> isomu = isomus.mu0.iso;

            // print the result
            if (is_last.to_uint() == 1)
            {
                if(VERBOSE)
                    std::cout << "... mu0  : " <<  isomu.to_uint()  << std::endl;

                // std::cout << "... mu0  : " <<  outmus.mu0.iso.to_uint()  << std::endl;
                // std::cout << "... mu1  : " <<  outmus.mu1.iso.to_uint()  << std::endl;
                // std::cout << "... mu2  : " <<  outmus.mu2.iso.to_uint()  << std::endl;
                // std::cout << "... mu3  : " <<  outmus.mu3.iso.to_uint()  << std::endl;
                // std::cout << "... mu4  : " <<  outmus.mu4.iso.to_uint()  << std::endl;
                // std::cout << "... mu5  : " <<  outmus.mu5.iso.to_uint()  << std::endl;
                // std::cout << "... mu6  : " <<  outmus.mu6.iso.to_uint()  << std::endl;
                // std::cout << "... mu7  : " <<  outmus.mu7.iso.to_uint()  << std::endl;
                // std::cout << "... mu8  : " <<  outmus.mu8.iso.to_uint()  << std::endl;
                // std::cout << "... mu9  : " <<  outmus.mu9.iso.to_uint()  << std::endl;
                // std::cout << "... mu10 : " <<  outmus.mu10.iso.to_uint() << std::endl;
                // std::cout << "... mu11 : " <<  outmus.mu11.iso.to_uint() << std::endl;
                // std::cout << "... mu12 : " <<  outmus.mu12.iso.to_uint() << std::endl;
                // std::cout << "... mu13 : " <<  outmus.mu13.iso.to_uint() << std::endl;
                // std::cout << "... mu14 : " <<  outmus.mu14.iso.to_uint() << std::endl;
                // std::cout << "... mu15 : " <<  outmus.mu15.iso.to_uint() << std::endl;
                // std::cout << "... mu16 : " <<  outmus.mu16.iso.to_uint() << std::endl;
                // std::cout << "... mu17 : " <<  outmus.mu17.iso.to_uint() << std::endl;
                // write the result to the output
                fout << itest << " " << mu_in.mu0.pt.to_uint() << " " << mu_in.mu0.eta.to_int() << " " << mu_in.mu0.phi.to_int() << " " << isomu.to_uint() << std::endl;
                ++nwritten;
            }
        }

        // ap_uint<18> test_out;

        // trk_in.trk0.qOverR = itest;
        // trk_in.trk1.qOverR = itest;
        // trk_in.trk2.qOverR = itest;
        // trk_in.trk3.qOverR = itest;
        // trk_in.trk4.qOverR = itest;
        // trk_in.trk5.qOverR = itest;
        // trk_in.trk6.qOverR = itest;
        // trk_in.trk7.qOverR = itest;
        // trk_in.trk8.qOverR = itest;
        // trk_in.trk9.qOverR = itest;
        // trk_in.trk10.qOverR = itest;
        // trk_in.trk11.qOverR = itest;
        // trk_in.trk12.qOverR = itest;
        // trk_in.trk13.qOverR = itest;
        // trk_in.trk14.qOverR = itest;
        // trk_in.trk15.qOverR = itest;
        // trk_in.trk16.qOverR = itest;
        // trk_in.trk17.qOverR = itest;

        // mu_in.mu0.pt = 10+itest;
        // mu_in.mu1.pt = 10+itest;
        // mu_in.mu2.pt = 10+itest;
        // mu_in.mu3.pt = 10+itest;
        // mu_in.mu4.pt = 10+itest;
        // mu_in.mu5.pt = 10+itest;
        // mu_in.mu6.pt = 10+itest;
        // mu_in.mu7.pt = 10+itest;
        // mu_in.mu8.pt = 10+itest;
        // mu_in.mu9.pt = 10+itest;
        // mu_in.mu10.pt = 10+itest;
        // mu_in.mu11.pt = 10+itest;
        // mu_in.mu12.pt = 10+itest;
        // mu_in.mu13.pt = 10+itest;
        // mu_in.mu14.pt = 10+itest;
        // mu_in.mu15.pt = 10+itest;
        // mu_in.mu16.pt = 10+itest;
        // mu_in.mu17.pt = 10+itest;

        // test_out = 0;

        // test_algo(mu_in, trk_in, test_out);

        // std::cout << " .. itest : " << itest << " -> " << test_out.to_uint() << std::endl;
    }
    fout.close();
    std::cout << " ... " << nwritten << " lines written to file " << ofile << std::endl;

}
