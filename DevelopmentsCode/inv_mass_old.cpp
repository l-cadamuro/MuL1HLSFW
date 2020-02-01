// the code takes the invariant mass of two tracks and computes them using the ap_* classes as the algo would do
// useful to tune the parameters for resolution
// the code takes in input a txt file with the two tracks to be made the inv mass per line
// pt1 theta1 phi1 pt2 theta2 phi2 
// NB: particles in the final state are supposed to be massless

// c++ -std=c++11 -lm -o inv_mass inv_mass.cpp -I /home/zynq/software/Xilinx/Vivado/2018.3/include/ `root-config --glibs --cflags`

#include <vector>
#include <utility>
#include <cmath>
#include <iostream>
#include <fstream>
#include "ap_int.h"
#include "ap_fixed.h"

#include "TLorentzVector.h"

using namespace std;

struct p3_polar_tr {
    double pt;
    double theta;
    double phi;
};

// ap_fixed <total_bits, integer_part_bits, rounding mode, overflow mode, >
template <int pt_w, int pt_dec_w, int theta_w, int phi_w>
struct p3_polar_tr_digi {
    
    // to auto-deduce after
    const int t_pt_w     = pt_w;
    const int t_pt_dec_w = pt_dec_w;
    const int t_theta_w  = theta_w;
    const int t_phi_w    = phi_w;

    typedef ap_ufixed<pt_w, pt_w - pt_dec_w>   pt_type;
    typedef ap_fixed<theta_w, 3>               theta_type;
    typedef ap_fixed<phi_w, 3, AP_TRN, AP_SAT> phi_type;

    pt_type     pt; // for the integer part keep it fixed, assuring a max pt can be reached. For a pT scale of 0.25 GeV the decimal part can be expressed in 2 bits (1/4 GeV)
    theta_type  theta; // instead for the angular part 2 bits are enough to go up to pi + a bit sign. Note that for the integer part one bit is considered for the sign (so 2 bits = 2 bits value + 1 sign)
    phi_type    phi;

    // ap_ufixed<pt_w, pt_w - pt_dec_w>   pt; // for the integer part keep it fixed, assuring a max pt can be reached. For a pT scale of 0.25 GeV the decimal part can be expressed in 2 bits (1/4 GeV)
    // ap_fixed<theta_w, 3> theta; // instead for the angular part 2 bits are enough to go up to pi + a bit sign. Note that for the integer part one bit is considered for the sign (so 2 bits = 2 bits value + 1 sign)
    // ap_fixed<phi_w, 3, AP_TRN, AP_SAT>   phi;
};

typedef p3_polar_tr_digi<15, 2, 16, 12> fw_type  ;

vector<pair<p3_polar_tr, p3_polar_tr> > read_input(string filename)
{
    std::cout << " ... reading from " << filename << std::endl;

    double pt1, theta1, phi1, pt2, theta2, phi2;
    std::fstream fin (filename.c_str());
    std::string line;
    int nlines = 0;

    std::vector<pair<p3_polar_tr, p3_polar_tr> > read_values;
    
    while (std::getline(fin, line))
    {
        std::istringstream iss(line);
        double buf; // all input data are integers
        iss >> pt1 >> theta1 >> phi1 >> pt2 >> theta2 >> phi2;
        p3_polar_tr p1;
        p3_polar_tr p2;
        
        p1.pt    = pt1;
        p1.theta = theta1;
        p1.phi   = phi1;
        
        p2.pt    = pt2;
        p2.theta = theta2;
        p2.phi   = phi2;

        read_values.push_back(make_pair(p1, p2));

        ++nlines;
    }
    
    std::cout << " ... read " << nlines << " lines" << std::endl;
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

    // double p1 = part1.pt/sin(part1.theta);
    // double p2 = part2.pt/sin(part2.theta);

    // double px1 = p1*sin(part1.theta)*cos(part1.phi);
    // double py1 = p1*sin(part1.theta)*sin(part1.phi);
    // double pz1 = p1*cos(part1.theta);

    // double px2 = p2*sin(part2.theta)*cos(part2.phi);
    // double py2 = p2*sin(part2.theta)*sin(part2.phi);
    // double pz2 = p2*cos(part2.theta);

    // double m2 = 2.*(p1*p2 - px1*px2 - py1*py2 - pz1*pz2);
    // return m2;
}


//////////////////////////////////////////////////////////////////////////////////////

template <class type_in, int out_w>
ap_fixed<out_w, 1> cos_lut (type_in input)
{
    float f_result = cos(input.to_float()); // the quantization of the input comes from the evaluation of ap.to_float()
    ap_fixed<out_w, 1> result = f_result;
    return result;
}

//////////////////////////////////////////////////////////////////////////////////////

template <class type_t, int cos_lut_w >
double minv2_digi (type_t part1, type_t part2)
{
    // double ptprod  = 2. * part1.pt.to_float() * part2.pt.to_float();
    
    // FIXME:: misses a factor of 2

    // the pt prod part
    const int nprod_w     = 2 * type_t::t_pt_w;
    const int nprod_dec_w = 2 * type_t::t_pt_dec_w;
    ap_ufixed < nprod_w , nprod_w-nprod_dec_w > ptprod = part1.pt * part2.pt;

    // the delta phi part
    // const int dphi_w = type_t::t_phi_w;
    typename type_t::phi_type  dphi        = part2.phi   - part1.phi;
    ap_fixed <cos_lut_w, 1> cos_dphi = cos_lut<type_t::phi_type, cos_lut_w> (dphi) ;

    // ap_fixed<type_t::t_theta_w> dtheta = part2.theta - part1.theta;
    // ap_fixed<type_t::t_theta_w> stheta = part2.theta + part1.theta; // check overflows?

    // ap_fixed <cos_lut_w> cos_stheta = cos_lut<type_t, cos_lut_w> ;

    // // double a1      = -1.*cos(part1.phi.to_float() - part2.phi.to_float());
    // // double a2      = (cos(part1.theta.to_float()-part2.theta.to_float()) + cos(part1.theta.to_float()+part2.theta.to_float()) - 2.) / (cos(part1.theta.to_float()+part2.theta.to_float()) - cos(part1.theta.to_float()-part2.theta.to_float()));
    // return ptprod * (a1 + a2);

    return 1.0;
}

//////////////////////////////////////////////////////////////////////////////////////


int main()
{
    // string fin_name = "BsMuMu.txt";
    // std::ofstream fout ("minv_Bsmumu.txt");

    string fin_name = "JPsiMuMu.txt";
    std::ofstream fout ("minv_JPsiMuMu.txt");

    auto trackpairs = read_input(fin_name);

    for (auto& trp : trackpairs)
    {
        // double m = sqrt(minv2(trp.first, trp.second));
        
        fw_type part1, part2;

        part1.pt    = trp.first.pt;
        part1.theta = trp.first.theta;
        part1.phi   = trp.first.phi;

        part2.pt    = trp.second.pt;
        part2.theta = trp.second.theta;
        part2.phi   = trp.second.phi;
        
        double m = sqrt(minv2_digi<fw_type> (part1, part2));
        

        fout << m << endl;
    }

    // ap_fixed<10, 3> v1 =  3.14;
    // ap_fixed<10, 3> v2 =  -3.14;
    // ap_fixed<10, 4> vd =  v2 - v1; // GOOD LOGIC : the delta accounts for the extra needed bit
    // cout << v1.to_float() << endl;
    // cout << v2.to_float() << endl;
    // cout << vd.to_float() << endl;

}