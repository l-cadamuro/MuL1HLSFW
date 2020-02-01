#include <cstdio>
#include <iostream>
#include <cmath>
// #include <random>
#include <vector>
#include <utility>
#include <tuple>
#include <fstream>
#include <sstream>
#include <string>
#include "src/inv_mass.h"
#include "src/inv_mass_df.h"
#include "ap_int.h"
#include "ap_fixed.h"

using namespace std;

#define NTEST 100000000 // input pattern has 6552

// #define VERBOSE false

// int roundDouble (double x){
//     int r = (int) (x + 0.5);
//     return r;
// }

// df to interface with my input file text
struct p3_polar_tr {
    double pt;
    double theta;
    double eta;
    double phi;
};


vector<pair<p3_polar_tr, p3_polar_tr> > read_input(string filename, std::vector<double>& v_minv, std::vector<double>& v_minv_massless, bool applyFiducialCuts = true)
{
    std::cout << " ... reading from " << filename << std::endl;

    double pt1, theta1, phi1, eta1, pt2, theta2, phi2, eta2, minv, minv_massless;
    std::fstream fin (filename.c_str());
    std::string line;
    int nlines = 0;
    int nskip  = 0;

    std::vector<pair<p3_polar_tr, p3_polar_tr> > read_values;
    v_minv.clear();
    v_minv_massless.clear();

    while (std::getline(fin, line))
    {
        std::istringstream iss(line);
        double buf; // all input data are integers
        iss >> pt1 >> theta1 >> phi1 >> eta1 >> pt2 >> theta2 >> phi2 >> eta2 >> minv >> minv_massless;

        if (applyFiducialCuts)
        {
            if (abs(eta1) > 3.0 || abs(eta2) > 3.0){
                ++nskip;
                continue;
            }
        }

        p3_polar_tr p1;
        p3_polar_tr p2;
        
        p1.pt    = pt1;
        p1.theta = theta1;
        p1.eta   = eta1;
        p1.phi   = phi1;
        
        p2.pt    = pt2;
        p2.theta = theta2;
        p2.eta   = eta2;
        p2.phi   = phi2;        

        read_values.push_back(make_pair(p1, p2));
        v_minv.push_back(minv);
        v_minv_massless.push_back(minv_massless);

        ++nlines;
    }
    
    std::cout << " ... read " << nlines << " lines, and " << nskip << " skipped because out of fiducial cuts" << std::endl;
    return read_values;
}

vector<tuple<p3_polar_tr, p3_polar_tr, p3_polar_tr> > read_input_3body(string filename, std::vector<double>& v_minv, std::vector<double>& v_minv_massless)
{
    std::cout << " ... reading from " << filename << std::endl;

    double pt1, theta1, phi1, eta1, pt2, theta2, phi2, eta2, pt3, theta3, phi3, eta3, minv, minv_massless;
    std::fstream fin (filename.c_str());
    std::string line;
    int nlines = 0;
    // int nskip  = 0;

    std::vector<tuple<p3_polar_tr, p3_polar_tr, p3_polar_tr> > read_values;
    v_minv.clear();
    v_minv_massless.clear();

    while (std::getline(fin, line))
    {
        std::istringstream iss(line);
        double buf; // all input data are integers
        iss >> pt1 >> theta1 >> phi1 >> eta1 >> pt2 >> theta2 >> phi2 >> eta2 >> pt3 >> theta3 >> phi3 >> eta3 >> minv >> minv_massless;

        // if (applyFiducialCuts)
        // {
        //     if (abs(eta1) > 3.0 || abs(eta2) > 3.0){
        //         ++nskip;
        //         continue;
        //     }
        // }

        p3_polar_tr p1;
        p3_polar_tr p2;
        p3_polar_tr p3;
        
        p1.pt    = pt1;
        p1.theta = theta1;
        p1.eta   = eta1;
        p1.phi   = phi1;
        
        p2.pt    = pt2;
        p2.theta = theta2;
        p2.eta   = eta2;
        p2.phi   = phi2;        

        p3.pt    = pt3;
        p3.theta = theta3;
        p3.eta   = eta3;
        p3.phi   = phi3;        

        read_values.push_back(make_tuple(p1, p2, p3));
        v_minv.push_back(minv);
        v_minv_massless.push_back(minv_massless);

        ++nlines;
    }
    
    // std::cout << " ... read " << nlines << " lines, and " << nskip << " skipped because out of fiducial cuts" << std::endl;
    std::cout << " ... read " << nlines << " lines" << std::endl;
    return read_values;
}

int main()
{

    // srand (123456);
    std::vector<double> v_minv;
    std::vector<double> v_minv_massless;

    // auto in_data        = read_input("/home/zynq/luca/MuL1HLSFW/DevelopmentsCode/JPsiMuMu.txt", v_minv, v_minv_massless);
    // std::string ofile = "/home/zynq/luca/MuL1HLSFW/MuonAlgorithms/inv_mass_hw_results.txt";

    auto in_data_tau3mu = read_input_3body("/home/zynq/luca/MuL1HLSFW/DevelopmentsCode/ta3mu_input_data.txt", v_minv, v_minv_massless);
    std::string ofile = "/home/zynq/luca/MuL1HLSFW/MuonAlgorithms/tau3mu_inv_mass_hw_results.txt";
    
    std::ofstream fout (ofile);
    int nwritten = 0;

    for (uint itest = 0; itest < NTEST; ++itest)
    {
        if (itest % 5000 == 0)
            std::cout << " ........ THIS IS THE TEST NUMBER .......... " << itest << std::endl;
        
        // for the 2 body invariant mass
        /*
        if (itest >= in_data.size())
            break; // patterns finished

        // double m = sqrt(minv2(trp.first, trp.second));
        auto& trp            = in_data.at(itest);
        double minv          = v_minv.at(itest);
        double minv_massless = v_minv_massless.at(itest);

        hw_p3_t part1, part2;

        part1.pt    = trp.first.pt;
        part1.theta = trp.first.theta;
        part1.eta   = trp.first.eta;
        part1.phi   = trp.first.phi;

        part2.pt    = trp.second.pt;
        part2.theta = trp.second.theta;
        part2.eta   = trp.second.eta;
        part2.phi   = trp.second.phi;
        
        // double m = sqrt(minv2_hw_floatprecision (part1, part2));
        // double m = sqrt(minv2_hw (part1, part2));
        // double m = sqrt(minv2_hw_uGTalgo_floatprecision(part1, part2));
        // double m = sqrt(minv2_hw_uGTAlgo(part1, part2));
        hw_minv2over2_t hw_minv2over2 = inv_mass(part1, part2);
        double m = sqrt(2.*hw_minv2over2.to_double());
        */

        // for the 3-body invariant mass
        if (itest >= in_data_tau3mu.size())
            break; // patterns finished

        auto& trp            = in_data_tau3mu.at(itest);
        double minv          = v_minv.at(itest);
        double minv_massless = v_minv_massless.at(itest);

        hw_p3_t part1, part2, part3;

        part1.pt    = get<0>(trp).pt;
        part1.theta = get<0>(trp).theta;
        part1.eta   = get<0>(trp).eta;
        part1.phi   = get<0>(trp).phi;

        part2.pt    = get<1>(trp).pt;
        part2.theta = get<1>(trp).theta;
        part2.eta   = get<1>(trp).eta;
        part2.phi   = get<1>(trp).phi;

        part3.pt    = get<2>(trp).pt;
        part3.theta = get<2>(trp).theta;
        part3.eta   = get<2>(trp).eta;
        part3.phi   = get<2>(trp).phi;

        hw_minv2over2_t hw_minv2over2 = inv_mass_3body(part1, part2, part3);
        double m = sqrt(2.*hw_minv2over2.to_double());
        
        // write to output

        fout << m << endl;
        // fout_minv_true << minv << endl;
        // fout_minv_massless << minv_massless << endl;
        ++nwritten;

    }

    fout.close();
    std::cout << " ... " << nwritten << " lines written to file " << ofile << std::endl;


}
