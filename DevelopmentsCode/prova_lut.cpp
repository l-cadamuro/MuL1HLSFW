#include "ap_int.h"
#include "ap_fixed.h"
#include <iostream>
#include "cos_dphi_LUT.h"
#include "cosh_deta_LUT.h"

using namespace std;

// c++ -std=c++11 -lm -o prova_lut prova_lut.cpp -I /home/zynq/software/Xilinx/Vivado/2018.3/include/

int main()
{
    cout << cos_dphi_lut [0] << " " << cos_dphi_lut[1] << " " << cos_dphi_lut[10] << endl;
    cout << cosh_deta_lut [0] << " " << cosh_deta_lut[1] << " " << cosh_deta_lut[10] << endl;
}
