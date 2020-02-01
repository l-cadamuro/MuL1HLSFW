#include <cstdio>
#include <iostream>
#include <random>
#include "src/division.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 20

// int roundDouble (double x){
//     int r = (int) (x + 0.5);
//     return r;
// }

int main()
{

    // srand (123456);
    
    for (uint itest = 0; itest < NTEST; ++itest)
    {

        std::cout << " ........ THIS IS THE TEST NUMBER .......... " << itest << std::endl;
        ap_fixed<12, 6> num = 30.456;
        ap_fixed<12, 6> den = 20.102;
        
        ap_fixed<24, 12> res = division(num, den);
        std::cout << num.to_double() << " / " << den.to_double() << " = " << res.to_double() << std::endl; 


    }
}
