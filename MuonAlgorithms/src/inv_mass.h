#ifndef INV_MASS_H
#define INV_MASS_H

#include "inv_mass_df.h"

hw_minv2over2_t inv_mass(hw_p3_t part1, hw_p3_t part2);
hw_minv2over2_t inv_mass_3body (hw_p3_t part1, hw_p3_t part2, hw_p3_t part3);

#endif // INV_MASS_H