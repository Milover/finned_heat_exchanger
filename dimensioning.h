/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	Functions used for dimensioning of the finned heat exchanger.

\*---------------------------------------------------------------------------*/

#ifndef DIMENSIONING_H
#define DIMENSIONING_H

#include <cmath>

#include "general.h"

// * * * * * * * * * * * * * Spacing characteristics * * * * * * * * * * * * //

// compute equivalent square fin height (c_ekv)
double computeC_ekv
(
	const double d_0,			// outer fin-tube diameter
	const double delta_f,		// fin thickness
	const double psi_f,			// fin characteristic
	const double s_f			// fin spacing
)
{
	double x
	{
		0.5 * M_PI * d_0 * s_f * (psi_f - 1.0 + delta_f / s_f) +
		0.25 * M_PI * std::pow(d_0, 2)
	};
	double c_ekv
	{
		std::sqrt(x)
	};

	return c_ekv;
}


// compute equivalent fin height (l_ekv)
double computeL_ekv
(
	const double c_ekv,			// equivalent square fin height
	const double d_0			// outer fin-tube diameter
)
{
	double l_ekv {0.5 * (c_ekv - d_0)};

	return l_ekv;
}


// compute longitudinal tube spacing (s_2)
double computeS_2
(
	const double c_ekv,			// equivalent square fin height
	const double s_12,			// ratio s_1 / s_2
	const double tau,			// ratio c_f / s_2
	const double z_1			// number of tubes in transverse direction
)
{
	double x
	{
		tau * (tau + (z_1 - 0.5) * s_12)
	};
	double s_2
	{
		std::sqrt(z_1 * std::pow(c_ekv, 2) / x)
	};

	// round to 3 decimals
	s_2 = std::round(s_2 * 1000.0) / 1000.0;

	return s_2;
}


// compute transverse tube spacing (s_1)
double computeS_1
(
	const double s_12,			// ratio s_1 / s_2
	const double s_2			// longitudinal tube spacing
)
{
	double s_1 {s_12 * s_2};

	// round to 3 decimals
	s_1 = std::round(s_1 * 1000.0) / 1000.0;

	return s_1;
}


// compute diagonal tube spacing (s_3)
double computeS_3
(
	const double s_1,			// transverse tube spacing
	const double s_2			// longitudinal tube spacing
)
{
	double s_3 {std::hypot(0.5 * s_1, s_2)};

	return s_3;
}


// compute fin parameter (psi_f)
double computePsi_f
(
	const double c_ekv,			// equivalent square fin height
	const double d_0,			// outer fin-tube diameter
	const double delta_f,		// fin thickness
	const double s_f			// fin spacing
)
{
	double psi_f
	{
		2.0 * (std::pow(c_ekv, 2) - 0.25 * M_PI * std::pow(d_0, 2)) /
		(M_PI * d_0 * s_f) +
		1.0 - delta_f / s_f
	};

	return psi_f;
}


// compute optimal s_12
double computeS_12
(
	const double psi_f,			// fin characteristic
	const double Re_2			// Reynolds number (outer)
)
{
	double Phi
	{
		0.5 * std::log
		(
			(0.189 * std::log(Re_2) - 1.0) / (1.0 - 0.029 * std::log(Re_2))
		)
	};
	double s_12
	{
		1.26 / psi_f + 2.0 + Phi
	};

	return s_12;
}


// compute actual height (H_r)
double computeH_r
(
	const double s_1,			// transverse tube spacing
	const double s_2,			// longitudinal tube spacing
	const double tau,			// ratio c_f / s_2
	const double z_1			// number of tubes in transverse direction
)
{
	double H_r
	{
		tau * s_2 + (z_1 - 0.5) * s_1
	};

	return H_r;
}

// compute total heigth (H)
double computeHeight
(
	const double H_r			// actual height
)
{
	// round to 2 decimal places
	double H {std::ceil(H_r * 100.0) / 100.0};

	return H;
}


// compute total width (W)
double computeWidth
(
	const double H				// height
)
{
	double W {std::round(H * 10.0) / 10.0};

	return W;
}


// * * * * * * * * * * * * * Geometric characteristics * * * * * * * * * * * //

// compute fin area per 1m of tube (A_f)
double computeA_f
(
	const double c_ekv,			// equivalent square fin height
	const double d_0,			// outer fin-tube diameter
	const double s_f			// fin spacing
)
{
	double A_f {2.0 * (std::pow(c_ekv, 2) - 0.785 * std::pow(d_0, 2)) / s_f};

	return A_f;
}


// compute fin-tube area per 1m of tube (A_t)
double computeA_t
(
	const double d_0,			// outer fin-tube diameter
	const double delta_f,		// fin thickness
	const double s_f			// fin spacing
)
{
	double A_t {M_PI * d_0 * (1.0 - delta_f / s_f)};

	return A_t;
}


// compute total fin area per 1m of tube (A)
double computeA
(
	const double A_f,			// fin area per 1m of tube
	const double A_t			// fin-tube area per 1m of tube
)
{
	double A {A_f + A_t};

	return A;
}


// compute total fin-tube area per 1m of tube (A_ft)
double computeA_ft
(
	const double d_0			// outer fin-tube diameter
)
{
	double A_ft {M_PI * d_0};

	return A_ft;
}


// compute inner area of tube per 1m of tube (A_in)
double computeA_in
(
	const double d_in			// inner fin-tube diameter
)
{
	double A_in {M_PI * d_in};

	return A_in;
}


// compute the ratio A_f / A
double computeAA_f
(
	const double A,				// total fin area per 1m of tube
	const double A_f			// fin area per 1m of tube
)
{
	double AA_f {A_f / A};

	return AA_f;
}


// compute the ratio A_t / A
double computeAA_t
(
	const double A,				// total fin area per 1m of tube
	const double A_t			// fin-tube area per 1m of tube
)
{
	double AA_t {A_t / A};

	return AA_t;
}


// compute the ratio A / A_in
double computeAA_in
(
	const double A,				// total fin area per 1m of tube
	const double A_in			// inner area of tube per 1m of tube
)
{
	double AA_in {A / A_in};

	return AA_in;
}

// * * * * * * * * * * * * * * * * Air velocity  * * * * * * * * * * * * * * //

// compute conventional diameter (d_cl)
double computeD_cl
(
	const double d_0,			// outer fin-tube diameter
	const double delta_f,		// fin thickness
	const double l_ekv,			// equivalent fin height
	const double s_f			// fin spacing
)
{
	double d_cl {d_0 + (2.0 * l_ekv * delta_f) / s_f};

	return d_cl;
}


// compute phi_cl
double computePhi_cl
(
	const double d_cl,			// conventional diameter
	const double s_1,			// transverse tube spacing
	const double s_3			// diagonal tube spacing
)
{
	double phi_cl {(s_1 - d_cl) / (s_3 - d_cl)};

	return phi_cl;
}


// compute minimal air passage area (F)
double computeF
(
	const double d_cl,			// conventional diameter
	const double H,				// height
	const double L_c,			// length of tube projection
	const double phi_cl,
	const double W,				// width
	const double z_1			// number of tubes in transverse direction
)
{
	double F
	{
		W * H - z_1 * L_c * d_cl
	};
	if (phi_cl > 2.0)
		F *= 2.0 / phi_cl;

	return F;
}


// compute max air velocity (v_2)
double computeV_2
(
	const double F,				// minimal air passage area
	const double m_2,			// air mass flow rate
	const double rho_2			// air density
)
{
	double v_2 {m_2 / (F * rho_2)};

	return v_2;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
