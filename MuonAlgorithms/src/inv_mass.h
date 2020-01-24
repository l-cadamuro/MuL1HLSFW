#ifndef INV_MASS_H
#define INV_MASS_H

// the invariant of two massless particles can be expressed as:
// minv^2 = ptprod * angprod
// with 
// ptprod  = 2. * v1.Pt()*v2.Pt()
// angprod = 1./(sin(th1)*sin(th2)) - cos(phi1-phi2) - 1/(tan(th1)*tan(th2))
// where th1, phi1, th2, phi2 are the angles of the two spatial momenta, in the ROOT TLorentzVector coordinate system, i.e.
// x = p  sin(theta) cos(phi)
// y = p  sin(theta) sin(phi)
// z = p  cos(theta)
//
// therefore the implementation is:
// compute ptprod with a DSP
// compute the 3 partial sums with 3 LUT (for the cos one, input dphi)
// sum the 3 partial sums
// multiplicate all

// templated on input and output return types

template <class ret_T, class pT_T, class theta_T, class phi_T>
ret_T inv_mass_sq ()

#endif // INV_MASS_H