// the code takes the invariant mass of two tracks and computes them using the ap_* classes as the algo would do
// useful to tune the parameters for resolution
// the code takes in input a txt file with the two tracks to be made the inv mass per line
// pt1 theta1 phi1 pt2 theta2 phi2 
// NB: particles in the final state are supposed to be massless

// c++ -std=c++11 -lm -o inv_mass inv_mass.cpp -I /home/zynq/software/Xilinx/Vivado/2018.3/include/  # `root-config --glibs --cflags`

#include <vector>
#include <utility>
#include <cmath>
#include <iostream>
#include <fstream>
#include "ap_int.h"
#include "ap_fixed.h"

// #include "TLorentzVector.h"

using namespace std;

struct p3_polar_tr {
    double pt;
    double theta;
    double eta;
    double phi;
};

// about the track pt : 10 bits for the track decimal part give a LSB of 1/2**10 = 0.00098 which is close to 0.00123 that one gets from expressing the track 1/R on 14 bits (+ 1 sign bit)
// and for the integer part 8 bits are enough because at most the pT goes to 256 (beyond 200 GeV we do not care too much)

// definition of the dataformat here
#define PT_W        18
#define PT_INT_W    8 //11
#define THETA_W     16
#define THETA_INT_W 3
#define ETA_W       16
#define ETA_INT_W   3
#define PHI_W       12
#define PHI_INT_W   3

///// definitio of the intermediate widths
// implement as a 4k x 9 LUT
// >> better implement as a 2k x 18 LUT
// while 1k x 36 has worse resolution
#define COS_PHI_ADDR_W     11 
#define COS_PHI_ADDR_W_INT 3 // max dphi is 2pi < 2**3 = 8
#define COS_PHI_OUT_W      18

// implement as 2 LUTs of 4k x 9  [6+6 bit addr each]
#define LUT_THETA_HALF_ADDR_W     6 
#define LUT_THETA_HALF_ADDR_W_INT 3 // max theta is pi < 2**2 = 4 + 1 sign bit
#define LUT_THETA_OUT_W     32 // 
#define LUT_THETA_OUT_W_INT 14

/// for the cos(delta)/cos(sum) implementation
#define COS_THETA_ADDR_W     12  // 11 
#define COS_THETA_ADDR_W_INT 3   // 3 // max sum(theta) is 2pi < 2**3 = 8 (for dtheta extra optimisation could be done to just go up to pi)
#define COS_THETA_OUT_W      9  // 18

// implement as a 2k x 18 LUT
#define DIV_LUT_NUM_W 6
#define DIV_LUT_DEN_W 6
#define DIV_LUT_OUT_W 9
#define DIV_LUT_OUT_INT_W 5

// remember
// 1k x 36 : 10 in -> 36 out
// 2k x 18 : 11 in -> 18 out
// 4k x 9  : 12 in -> 9 out
#define COSH_ETA_ADDR_W      11
#define COSH_ETA_ADDR_W_INT  3  // max eta is 2.8 -> max deta is 5.6 < 2**3 = 8
#define COSH_ETA_OUT_W       18
#define COSH_ETA_OUT_W_INT   8  // because cosh (5.6) = 135

// here what is used for inputs to the LUT
typedef ap_ufixed<COS_PHI_ADDR_W, COS_PHI_ADDR_W_INT>   dphi_t;
typedef ap_ufixed<COSH_ETA_ADDR_W, COSH_ETA_ADDR_W_INT> deta_t;

typedef ap_ufixed<PT_W, PT_INT_W>                    hw_pt_t;
typedef ap_ufixed<THETA_W, THETA_INT_W>              hw_theta_t;
typedef ap_fixed<THETA_W, THETA_INT_W>               hw_eta_t;
typedef ap_fixed<PHI_W, PHI_INT_W, AP_TRN, AP_SAT>   hw_phi_t;

// ap_fixed <total_bits, integer_part_bits, rounding mode, overflow mode, >
struct p3_polar_tr_digi {
    
    hw_pt_t     pt;   
    hw_theta_t  theta;
    hw_eta_t    eta;
    hw_phi_t    phi;
};

typedef p3_polar_tr_digi hw_p3_t  ;

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


double minv2 (p3_polar_tr part1, p3_polar_tr part2)
{
    // the actual calculation
    // m2 = 2.*p1*p2*(1.-costh) #massless
    // which becomes

    double ptprod  = 2. * part1.pt*part2.pt;
    double a1 = -1.*cos(part1.phi - part2.phi);
    double a2 = (cos(part1.theta-part2.theta) + cos(part1.theta+part2.theta) - 2.) / (cos(part1.theta+part2.theta) - cos(part1.theta-part2.theta));
    return ptprod * (a1 + a2);

}

//////////////////////////////////////////////////////////////////////////////////////

double minv2_hw_floatprecision (hw_p3_t part1, hw_p3_t part2)
{
    double ptprod  = 2. * part1.pt.to_double() * part2.pt.to_double();    
    double a1      = -1.*cos(part1.phi.to_double() - part2.phi.to_double());
    double a2      = (cos(part1.theta.to_double()-part2.theta.to_double()) + cos(part1.theta.to_double()+part2.theta.to_double()) - 2.) / (cos(part1.theta.to_double()+part2.theta.to_double()) - cos(part1.theta.to_double()-part2.theta.to_double()));
    return ptprod * (a1 + a2);
}

double minv2_hw_uGTalgo_floatprecision (hw_p3_t part1, hw_p3_t part2)
{
    double ptprod  = 2. * part1.pt.to_double() * part2.pt.to_double();    
    double a1      = -1.*cos(part1.phi.to_double() - part2.phi.to_double());
    double a2      = cosh(part1.eta.to_double() - part2.eta.to_double());

    // outliers happen only if 
    // if (sqrt(ptprod * (a1 + a2)) > 10)
    // {
    //     cout << "HIGH MASS? " << sqrt(ptprod * (a1 + a2)) << " " << ptprod << " " << a1 << a2
    //          << " ... ptetaphi1 " << part1.pt << " " << part1.eta << " " << part1.phi
    //          << " ... ptetaphi2 " << part2.pt << " " << part2.eta << " " << part2.phi
    //     << endl;
    // }

    return ptprod * (a1 + a2);
}

//////////////////////////////////////////////////////////////////////////////////////

template <size_t out_w, class type_in>
ap_fixed<out_w, 1> cos_lut (type_in input)
{
    double f_result = cos(input.to_double()); // the quantization of the input comes from the evaluation of ap.to_double()
    ap_fixed<out_w, 1> result = f_result;
    return result;
}

template <size_t out_w, size_t out_w_i, class type_in>
ap_ufixed<out_w, out_w_i> cosh_lut (type_in input)
{
    double f_result = cosh(input.to_double()); // the quantization of the input comes from the evaluation of ap.to_double()
    ap_ufixed<out_w, out_w_i> result = f_result;
    return result;
}

//////////////////////////////////////////////////////////////////////////////////////

// out_w  : the total number of bits in the output of the LUT
// int_w  : how many of those bits are allocated for the integer part
// addr_w : how many bits should be used for addressing *each* theta
// addr_int_w : how many those bits are allocated for the integer part
// NOTE: the final LUT will have an address builts as [..th1..][..th2..] so 2*addr_w in width
template <size_t out_w, size_t int_w, size_t addr_w, size_t addr_int_w> 
ap_fixed<out_w, int_w, AP_TRN, AP_SAT> angular_theta_LUT (hw_theta_t t1, hw_theta_t t2)
{
    // accounts for the reduction of the precision on theta1, theta2
    ap_fixed<addr_w, addr_int_w> t1red = t1;
    ap_fixed<addr_w, addr_int_w> t2red = t2;

    double ft1 = t1red.to_double();
    double ft2 = t2red.to_double();

    ap_fixed<out_w, int_w, AP_TRN, AP_SAT> result = 0.0;

    // here the big expression
    double num  = cos(ft1+ft2) + cos(ft1-ft2) -2.0;
    double den  = cos(ft1+ft2) - cos(ft1-ft2);

    if (den == 0)
        return result; // just to avoid domain error, never happens at the LHC

    double fres = num/den;

    result = fres; // here the digitalisation done by the LUT happens
    return result;
}

// //////////////////////////////////////////////////////////////////////////////////////
// // being a cosine, one bit is enough for the integer part (it's the sign)
// template <size_t out_w, size_t addr_w, size_t addr_int_w> 
// // out_w  : the total number of bits in the output of the LUT
// // addr_w : how many bits should be used for addressing the deltaPhi
// // addr_int_w : how many those bits are allocated for the integer part
// ap_fixed<out_w, 1, AP_TRN, AP_SAT> angular_phi_LUT (hw_phi_t p1, hw_phi_t p2)
// {
//     hw_phi_t dphi = p2-p1;
//     ap_fixed<addr_w, addr_int_w> dphired = dphi; // here the address is compressed in precision
//     ap_fixed<out_w, 1> result = cos_lut<out_w>(dphired);
//     return result;
// }

//////////////////////////////////////////////////////////////////////////////////////
// being a cosine, one bit is enough for the integer part (it's the sign)
template <size_t out_w, size_t int_w> 
ap_ufixed<out_w, int_w> pt_prod (hw_pt_t pt1, hw_pt_t pt2)
{
    ap_ufixed<out_w, int_w> ptprod = pt1*pt2;
    return ptprod;
}

//////////////////////////////////////////////////////////////////////////////////////
// division LUT
// NOTE: return absolute value for best precision
template <size_t out_w, size_t out_int_w, class num_t, class den_t>
ap_ufixed<out_w, out_int_w> divide(num_t num, den_t den)
{
    ap_ufixed<out_w, out_int_w> res = 0;
    double fnum = num.to_double();
    double fden = den.to_double();
    if (fden == 0)
        return res;
    double result = fnum/fden;
    res = result;
    return res;
}

//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


template <size_t out_w, class type_in>
void make_cos_lut (string foutname, bool addFormat = true, string lutname = "cos_dphi_lut")
{
    cout << "... building the cos lut on file " << foutname << endl;
    std::ofstream fout (foutname);
    size_t nbits = type_in::width; // total number of bits of the ap_* input type
    cout << "....... address expressed on " << nbits << " bits" << endl;
    cout << "....... output expressed as ap_fixed with " << out_w << " bits total, 1 for integer part" << endl;
    size_t maxn = (1 << nbits);
    
    if (addFormat)
        fout << "const ap_fixed<" << out_w << ", 1> " << lutname << "[1 << " << nbits << "] = {" << endl;


    for (size_t i = 0; i < maxn; ++i)
    {
        type_in addr;
        addr.range() = i; // will assing to its digital representation the integer value of i (address)
        ap_fixed<out_w, 1> out = cos_lut<out_w> (addr);
        
        string ostr = "";
        if (addFormat){
            ostr += "   ";
            ostr += to_string_with_precision(out.to_double(), out_w); // clearly the precision of out_w on the double is too much, but will be trimmed by vivado synthesis anyway
            ostr += ",";
        }
        else
            ostr += to_string(out.to_double());

        fout << ostr << endl;
    }
    if (addFormat)
        fout << "};" << endl;

}

//////////////////////////////////////////////////////////////////////////////////////
template <size_t out_w, size_t out_w_i, class type_in>
void make_cosh_lut (string foutname, bool addFormat = true, string lutname = "cosh_deta_lut")
{
    cout << "... building the cosh lut on file " << foutname << endl;
    std::ofstream fout (foutname);
    size_t nbits = type_in::width; // total number of bits of the ap_* input type
    cout << "....... address expressed on " << nbits << " bits" << endl;
    cout << "....... output expressed as ap_ufixed with " << out_w << " bits total, " << out_w_i << " for integer part" << endl;
    size_t maxn = (1 << nbits);

    if (addFormat)
        fout << "const ap_ufixed<" << out_w << ", " << out_w_i << " > " << lutname << "[1 << " << nbits << "] = {" << endl;

    for (size_t i = 0; i < maxn; ++i)
    {
        type_in addr;
        addr.range() = i; // will assing to its digital representation the integer value of i (address)
        ap_ufixed<out_w, out_w_i> out = cosh_lut<out_w, out_w_i> (addr);
        string ostr = "";
        if (addFormat){
            ostr += "   ";
            ostr += to_string_with_precision(out.to_double(), out_w); // clearly the precision of out_w on the double is too much, but will be trimmed by vivado synthesis anyway
            ostr += ",";
        }
        else
            ostr += to_string(out.to_double());

        fout << ostr << endl;
    }
    if (addFormat)
        fout << "};" << endl;

}

//////////////////////////////////////////////////////////////////////////////////////

double minv2_hw (hw_p3_t part1, hw_p3_t part2)
{
    // pt part
    ap_fixed<2*PT_W, 2*PT_INT_W> hw_ptprod  = part1.pt * part2.pt; // sum the output bits in the multiplication - for a DSP with 18*24
    
    // phi part
    ap_fixed<COS_PHI_ADDR_W+1, COS_PHI_ADDR_W_INT+1> dphi_sign =  part1.phi - part2.phi;
    if (dphi_sign < 0)
        dphi_sign = -1*dphi_sign;
    dphi_t dphi = dphi_sign;
    ap_fixed<COS_PHI_OUT_W, 1> hw_phipart = -1*cos_lut<COS_PHI_OUT_W> (dphi);

    // theta part
    // NOTE: by indexing two thetas, the precision on the two becomes very reduced
    // ap_fixed<LUT_THETA_OUT_W, LUT_THETA_OUT_W_INT, AP_TRN, AP_SAT> hw_thetapart = angular_theta_LUT<LUT_THETA_OUT_W, LUT_THETA_OUT_W_INT, LUT_THETA_HALF_ADDR_W, LUT_THETA_HALF_ADDR_W_INT>  (part1.theta, part2.theta);

    ap_fixed<COS_THETA_ADDR_W+1, COS_THETA_ADDR_W_INT+1> dtheta_sign =  part1.theta - part2.theta;
    if (dtheta_sign < 0)
        dtheta_sign = -1*dtheta_sign;
    
    ap_ufixed<COS_THETA_ADDR_W, COS_THETA_ADDR_W_INT> dtheta   = dtheta_sign;
    ap_ufixed<COS_THETA_ADDR_W, COS_THETA_ADDR_W_INT> sumtheta = part1.theta + part2.theta;

    ap_fixed<COS_THETA_OUT_W, 1> cos_sumtheta = cos_lut<COS_THETA_OUT_W> (sumtheta);
    ap_fixed<COS_THETA_OUT_W, 1> cos_dtheta   = cos_lut<COS_THETA_OUT_W> (dtheta);

    // now I have to make sums that will make num in range -4, 0 and den in range -2, 2
    // so we need to reexpress the result as:
    // num : unsigned (will multiply by -1 after) , 3 integer bits
    // den : signed   , 2+1 integer bits
    // NOTE: maybe it is enough to keep one bit less (number will never saturate to +4, +/-2, but it's still probably close enough)

    ap_ufixed<COS_THETA_OUT_W, 3> num = 2;
    num = num - cos_sumtheta;
    num = num - cos_dtheta;

    ap_fixed<COS_THETA_OUT_W, 3> den = cos_sumtheta;
    den = den - cos_dtheta;

    // FIXME : float division
    // double a2 = -1.*num.to_double()/den.to_double();

    // compress the inputs for the division LUT
    ap_ufixed<DIV_LUT_NUM_W, 3> lutnum = abs(num.to_double());
    ap_ufixed<DIV_LUT_DEN_W, 3> lutden = abs(den.to_double());
    ap_ufixed<DIV_LUT_OUT_W, DIV_LUT_OUT_INT_W> hw_thetapart_us = divide<DIV_LUT_OUT_W, DIV_LUT_OUT_INT_W>(lutnum, lutden);

    // fix the sign
    ap_fixed<DIV_LUT_OUT_W+1, DIV_LUT_OUT_INT_W+1> hw_thetapart = hw_thetapart_us;
    if (num*den < 0)
        hw_thetapart = -1*hw_thetapart;
    // and the -1 I am still missing from before
    hw_thetapart = -1*hw_thetapart;

    // another test : ap_fixed division

    ap_fixed <32, 16, AP_TRN, AP_SAT> hw_thetapart_hlsdiv;
    if (den.to_double() == 0)
        hw_thetapart_hlsdiv = 0;
    else{
        hw_thetapart_hlsdiv = num/den;
        // cout << "... doing " << num.to_double() << " / " << den.to_double() << " = " << hw_thetapart_hlsdiv.to_double() << endl;
    }
    hw_thetapart_hlsdiv = -1*hw_thetapart_hlsdiv;

    // cout << "target : " << num.to_double() << " / " << den.to_double() << " = " << -1 * num.to_double()/den.to_double() << " my LUT : " << hw_thetapart.to_double() << " " << " HLS div : " << hw_thetapart_hlsdiv.to_double() << endl;

    double ptprod = hw_ptprod.to_double();
    double a1     = hw_phipart.to_double();
    // double a2     = hw_thetapart.to_double();
    double a2     = hw_thetapart_hlsdiv.to_double();
    // double a2      = (cos(part1.theta.to_double()-part2.theta.to_double()) + cos(part1.theta.to_double()+part2.theta.to_double()) - 2.) / (cos(part1.theta.to_double()+part2.theta.to_double()) - cos(part1.theta.to_double()-part2.theta.to_double()));

    // double target_ptprod  = part1.pt.to_double() * part2.pt.to_double();    // the factor 2 is removed for now
    // double target_a1      = -1.*cos(part1.phi.to_double() - part2.phi.to_double());
    // double target_a2      = (cos(part1.theta.to_double()-part2.theta.to_double()) + cos(part1.theta.to_double()+part2.theta.to_double()) - 2.) / (cos(part1.theta.to_double()+part2.theta.to_double()) - cos(part1.theta.to_double()-part2.theta.to_double()));
    // cout << " ptprod : " <<  ptprod << " expect " << target_ptprod << endl;
    // cout << " a1 : "     <<  a1     << " expect " << target_a1     << endl;
    // cout << " a2 : "     <<  a2     << " expect " << target_a2     << endl;


//     double a1      = -1.*cos(part1.phi.to_double() - part2.phi.to_double());
//     double a2      = (cos(part1.theta.to_double()-part2.theta.to_double()) + cos(part1.theta.to_double()+part2.theta.to_double()) - 2.) / (cos(part1.theta.to_double()+part2.theta.to_double()) - cos(part1.theta.to_double()-part2.theta.to_double()));
    return 2. * ptprod * (a1 + a2);
}

//////////////////////////////////////////////////////////////////////////////////////

double minv2_hw_uGTAlgo (hw_p3_t part1, hw_p3_t part2)
{
    // pt part
    ap_fixed<2*PT_W, 2*PT_INT_W> hw_ptprod  = part1.pt * part2.pt; // sum the output bits in the multiplication - for a DSP with 18*24
    
    // phi part
    ap_fixed<COS_PHI_ADDR_W+1, COS_PHI_ADDR_W_INT+1> dphi_sign =  part1.phi - part2.phi;
    if (dphi_sign < 0)
        dphi_sign = -1*dphi_sign;
    ap_ufixed<COS_PHI_ADDR_W, COS_PHI_ADDR_W_INT> dphi = dphi_sign;
    ap_fixed<COS_PHI_OUT_W, 1> hw_phipart = -1*cos_lut<COS_PHI_OUT_W> (dphi);

    // eta part
    ap_fixed<COSH_ETA_ADDR_W+1, COSH_ETA_ADDR_W_INT+1> deta_sign =  part1.eta - part2.eta;
    if (deta_sign < 0)
        deta_sign = -1*deta_sign;
    deta_t deta = deta_sign;
    ap_ufixed<COSH_ETA_OUT_W, COSH_ETA_OUT_W_INT> hw_etapart = cosh_lut<COSH_ETA_OUT_W, COSH_ETA_OUT_W_INT> (deta);

    double ptprod = hw_ptprod.to_double();
    double a1     = hw_phipart.to_double();
    double a2     = hw_etapart.to_double();

    return 2. * ptprod * (a1 + a2);
}


//////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
{
    // string fin_name = "BsMuMu.txt";
    // std::ofstream fout ("minv_Bsmumu.txt");

    string fin_name = "JPsiMuMu.txt";
    string fout_name = "minv_JPsiMuMu_test.txt";
    
    if (argc > 1)
        fout_name = argv[1];
    if (argc > 2)
        fin_name = argv[2];

    cout << " ... in file  : " << fin_name << endl;
    cout << " ... out file : " << fout_name << endl;

    std::ofstream fout (fout_name);
    std::ofstream fout_minv_true ("minv_JPsiMuMu_true.txt");
    std::ofstream fout_minv_massless ("minv_JPsiMuMu_massless.txt");

    std::vector<double> v_minv;
    std::vector<double> v_minv_massless;
    auto trackpairs = read_input(fin_name, v_minv, v_minv_massless);

    // for (auto& trp : trackpairs)
    for (size_t itrp = 0; itrp < trackpairs.size(); ++itrp)
    {
        // double m = sqrt(minv2(trp.first, trp.second));
        auto& trp   = trackpairs.at(itrp);
        double minv = v_minv.at(itrp);
        double minv_massless = v_minv_massless.at(itrp);

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
        double m = sqrt(minv2_hw_uGTAlgo(part1, part2));

        fout << m << endl;
        fout_minv_true << minv << endl;
        fout_minv_massless << minv_massless << endl;
    }

    // ---- build the LUTs for the calculation in HLS
    if (true)
    {
        make_cos_lut<COS_PHI_OUT_W, dphi_t> ("cos_dphi_LUT.h");
        make_cosh_lut<COSH_ETA_OUT_W, COSH_ETA_OUT_W_INT, deta_t>("cosh_deta_LUT.h");
    }

    // checking the traits
    // will infer from the type width and integer width (for fixed/ufixed)
    // cout << " fixed " << ap_fixed  <10, 6> :: width << " " << ap_fixed  <10, 6> :: iwidth << endl;
    // cout << " ufixed " << ap_ufixed <10, 6> :: width << " " << ap_ufixed <10, 6> :: iwidth << endl;
    // cout << " int " << ap_int    <10> :: width << endl;
    // cout << " uint " << ap_uint   <10> :: width << endl;

    // converting a ap_fixed into its integer representation - useful for indexing a LUT
    // ap_fixed<3, 2, AP_TRN, AP_SAT> x;
    // ap_fixed<3, 2, AP_TRN, AP_SAT> y;
    // double step = 0.5;
    // double start = -4;
    // double stop = 4;
    // for (uint istep = 0; true; ++istep)
    // {
    //     double d = start + (istep*step);
    //     if (d > stop)
    //         break;
    //     x = d;
    //     int i = x.range();
    //     y.range() = i; // can also be assigned reversed
    //     cout << d << " ... ) " << x << " -> " << i << " -> " << y << endl;
    // }

    // // checking boundaries of my objects
    // all OK within their range
    // hw_pt_t    pt_min     = 2.0; 
    // hw_theta_t theta_min  = 0.0;
    // hw_eta_t   eta_min    = -2.5;
    // hw_phi_t   phi_min    = -3.1415926536;
    // //
    // hw_pt_t    pt_max     = 200.0; 
    // hw_theta_t theta_max  = 3.1415926536;
    // hw_eta_t   eta_max    = 2.5;
    // hw_phi_t   phi_max    = 3.1415926536;
    // //
    // cout << "pt_min = "     << pt_min.to_double()     << endl;
    // cout << "theta_min = "  << theta_min.to_double()  << endl;
    // cout << "eta_min = "    << eta_min.to_double()    << endl;
    // cout << "phi_min = "    << phi_min.to_double()    << endl;
    // cout << "pt_max = "     << pt_max.to_double()     << endl;
    // cout << "theta_max = "  << theta_max.to_double()  << endl;
    // cout << "eta_max = "    << eta_max.to_double()    << endl;
    // cout << "phi_max = "    << phi_max.to_double()    << endl;

    // // numbers with > 32 bits
    // double d = 260144.921921929139;
    // ap_ufixed<36, 18> bigfloat = d;
    // cout << std::setprecision(10) << bigfloat.to_float() << " " << bigfloat.to_double() << " " << d << endl; // all the 3 seem OK

    // // delta of unsigned quantities
    // ap_ufixed <4, 4> a = 1;
    // ap_ufixed <4, 4> b = 15;
    // ap_ufixed<4,4> udelta = a-b;  // gives 2
    // ap_fixed<5,5> sdelta  = a-b;  // good, gives 14
    // cout << a.to_float() << " " << b.to_float() << " " << udelta.to_float() << " " << sdelta.to_float() << endl; // works fine


    // // fixed to ufixed
    // ap_fixed <7, 4> is = -4.5213; // 3 bits for 
    // ap_ufixed <6, 3> ius = -1*is;
    // cout << is.to_float() << " " << ius.to_float() << endl; // works fine

    // // reducing the precision of an ap_fixed
    // ap_fixed <10, 6> i1 = 3.1022342;
    // ap_fixed <6, 4>  i2 = i1;
    // ap_fixed <10, 6> i3 = i2;
    // cout << i1.to_float() << " " << i2.to_float() << " " << i3.to_float() << endl;


    // handling truncations and saturation the right way
    // ap_fixed<4, 2, AP_TRN, AP_SAT>  x = 9999;
    // cout << x.to_float() << endl; // 1.75
    // x = -9999;
    // cout << x.to_float() << endl; // -2

    // // // example of using the cosine function
    // const int pi = 3.14159265359;
    // std::vector<ap_fixed<6, 3> > angl = {
    //     -2*pi, 
    //     -pi/2,
    //     0,
    //     pi,
    //     pi/2
    // };
    // for (auto a : angl){
    //     auto x = cos_lut<3>(a); // type in the cos_lut template can be auto deduced!
    //     cout << a.to_float() << " " << x << " cfr " << cos(a.to_float()) << endl;
    // }

    // ap_fixed<10, 3> v1 =  3.14;
    // ap_fixed<10, 3> v2 =  -3.14;
    // ap_fixed<10, 4> vd =  v2 - v1; // GOOD LOGIC : the delta accounts for the extra needed bit
    // cout << v1.to_float() << endl;
    // cout << v2.to_float() << endl;
    // cout << vd.to_float() << endl;

}
