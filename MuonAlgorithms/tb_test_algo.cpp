#include <cstdio>
#include <iostream>
#include <random>
#include "src/test_algo.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 20
#define NTMT 18

int roundDouble (double x){
    int r = (int) (x + 0.5);
    return r;
}

int main()
{

    // srand (123456);
    std::default_random_engine rndmEngine;
    std::uniform_real_distribution<double> distribution(10., 100.);


    for (uint itest = 0; itest < NTEST; ++itest)
    {

        std::cout << " ........ THIS IS THE TEST NUMBER .......... " << itest << std::endl;
 
        // generate the muons for this test
        // I only care for the pT actually, so let's not bother to generate other variables
        std::vector<int> pts_mu (18); // there will be 18 muons
        std::vector<int> pts_trks (18*NTMT); // and 18 tracks sent for every TMT period

        for (uint i = 0; i < pts_mu.size(); ++i){
            pts_mu.at(i) = roundDouble(distribution(rndmEngine));
        }

        for (uint i = 0; i < pts_trks.size(); ++i){
            pts_trks.at(i) = roundDouble(distribution(rndmEngine));
        }

        track_input trk_in;
        muon_input  mu_in;


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

        for (uint itmt = 0; itmt < NTMT; ++itmt)
        {
            trk_in.trk0.qOverR = pts_trks.at(0 + NTMT*itmt);
            trk_in.trk1.qOverR = pts_trks.at(1 + NTMT*itmt);
            trk_in.trk2.qOverR = pts_trks.at(2 + NTMT*itmt);
            trk_in.trk3.qOverR = pts_trks.at(3 + NTMT*itmt);
            trk_in.trk4.qOverR = pts_trks.at(4 + NTMT*itmt);
            trk_in.trk5.qOverR = pts_trks.at(5 + NTMT*itmt);
            trk_in.trk6.qOverR = pts_trks.at(6 + NTMT*itmt);
            trk_in.trk7.qOverR = pts_trks.at(7 + NTMT*itmt);
            trk_in.trk8.qOverR = pts_trks.at(8 + NTMT*itmt);
            trk_in.trk9.qOverR = pts_trks.at(9 + NTMT*itmt);
            trk_in.trk10.qOverR = pts_trks.at(10 + NTMT*itmt);
            trk_in.trk11.qOverR = pts_trks.at(11 + NTMT*itmt);
            trk_in.trk12.qOverR = pts_trks.at(12 + NTMT*itmt);
            trk_in.trk13.qOverR = pts_trks.at(13 + NTMT*itmt);
            trk_in.trk14.qOverR = pts_trks.at(14 + NTMT*itmt);
            trk_in.trk15.qOverR = pts_trks.at(15 + NTMT*itmt);
            trk_in.trk16.qOverR = pts_trks.at(16 + NTMT*itmt);
            trk_in.trk17.qOverR = pts_trks.at(17 + NTMT*itmt);            

            ap_uint<1> is_last = (itmt == NTMT - 1 ? 1 : 0);

            std::cout << " @@@ TMT # " << itmt << " / " << NTMT << " ---- ISLAST? " << is_last.to_uint() << std::endl;

            iso_muons_out outmus = test_algo(mu_in, trk_in, is_last);

            // print the result
            if (is_last.to_uint() == 1)
            {
                std::cout << "... mu0  : " <<  outmus.mu0.iso.to_uint()  << std::endl;
                std::cout << "... mu1  : " <<  outmus.mu1.iso.to_uint()  << std::endl;
                std::cout << "... mu2  : " <<  outmus.mu2.iso.to_uint()  << std::endl;
                std::cout << "... mu3  : " <<  outmus.mu3.iso.to_uint()  << std::endl;
                std::cout << "... mu4  : " <<  outmus.mu4.iso.to_uint()  << std::endl;
                std::cout << "... mu5  : " <<  outmus.mu5.iso.to_uint()  << std::endl;
                std::cout << "... mu6  : " <<  outmus.mu6.iso.to_uint()  << std::endl;
                std::cout << "... mu7  : " <<  outmus.mu7.iso.to_uint()  << std::endl;
                std::cout << "... mu8  : " <<  outmus.mu8.iso.to_uint()  << std::endl;
                std::cout << "... mu9  : " <<  outmus.mu9.iso.to_uint()  << std::endl;
                std::cout << "... mu10 : " <<  outmus.mu10.iso.to_uint() << std::endl;
                std::cout << "... mu11 : " <<  outmus.mu11.iso.to_uint() << std::endl;
                std::cout << "... mu12 : " <<  outmus.mu12.iso.to_uint() << std::endl;
                std::cout << "... mu13 : " <<  outmus.mu13.iso.to_uint() << std::endl;
                std::cout << "... mu14 : " <<  outmus.mu14.iso.to_uint() << std::endl;
                std::cout << "... mu15 : " <<  outmus.mu15.iso.to_uint() << std::endl;
                std::cout << "... mu16 : " <<  outmus.mu16.iso.to_uint() << std::endl;
                std::cout << "... mu17 : " <<  outmus.mu17.iso.to_uint() << std::endl;
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
}
