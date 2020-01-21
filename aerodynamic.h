/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	Functions used for aerodynamic analysis of the finned heat exchanger.

\*---------------------------------------------------------------------------*/

#ifndef AERODYNAMIC_H
#define AERODYNAMIC_H

#include <cmath>

// * * * * * * * * * * * * * Aerodynamic resistance  * * * * * * * * * * * * //

// compute the ratio A_tot / F
double computeFA_tot
(
	const double d_0,			// outer fin-tube diameter
	const double delta_f,		// fin thickness
	const double l_f,			// equivalent fin height
	const double s_1,			// transverse tube spacing
	const double s_f			// fin spacing
)
{
	double FA_tot
	{
		M_PI * (d_0 * s_f + 2.0 * l_f * delta_f + 2.0 * l_f * (l_f + d_0)) /
		(s_1 * s_f - (d_0 * s_f + 2.0 * l_f * delta_f))
	};

	return FA_tot;
}


// compute the equivalent diameter of staggered bundle (d_eq)
double computeD_eq
(
	const double d_0,			// outer fin-tube diameter
	const double delta_f,		// fin thickness
	const double l_f,			// equivalent fin height
	const double phi_cl,
	const double s_1,			// transverse tube spacing
	const double s_f			// fin spacing
)
{
	double d_eq
	{
		2.0 * (s_f * (s_1 - d_0) - 2.0 * l_f * delta_f) / (2.0 * l_f + s_f)
	};
	if (phi_cl > 2.0)
		d_eq /= 0.5 * phi_cl;


	return d_eq;
}


// compute exponent (n_air)
double computeN_air
(
	const double FA_tot,		// the ratio A_tot / F
	const double s_12			// ratio s_1 / s_2
)
{
	double n_air
	{
		0.17 * std::pow(FA_tot, 0.25) * std::pow(s_12, 0.57) * std::exp(-0.36 * s_12)
	};

	return n_air;
}


// compute C_f_air
double computeC_f_air
(
	const double FA_tot,		// the ratio A_tot / F
	const double s_12			// ratio s_1 / s_2
)
{
	double C_f_air
	{
		2.8 * std::pow(FA_tot, 0.53) * std::pow(s_12, 1.3) * std::exp(-0.9 * s_12)
	};

	return C_f_air;
}


// compute C_Z_air
double computeC_z_air
(
	const double z_2			// number of tubes in longitudinal direction
)
{
	double C_z_air {1.0};

	if (z_2 < 6.0)
		C_z_air = std::exp(0.1 * (6.0 / z_2 - 1.0));

	return C_z_air;
}


// compute resistance coefficient (zeta_0)
double computeZeta_0
(
	const double C_f_air,
	const double C_z_air,
	const double n_air,
	const double Re_air			// Reynolds number (air)
)
{
	double zeta_0
	{
		C_f_air * C_z_air * std::pow(Re_air, -n_air)
	};

	return zeta_0;
}


// compute operating condition correction factor (C_op)
double computeC_op()
{
	double C_op {1.1};

	return C_op;
}


// compute pressure drop (deltaP)
double computeDeltaP
(
	const double C_op,			// operating condition correction factor
	const double rho_2,			// density (outer)
	const double v_2,			// velocity (outer)
	const double z_2,			// number of tubes in longitudinal direction
	const double zeta_0			// resistance coefficient
)
{
	double deltaP
	{
		C_op * zeta_0 * z_2 * 0.5 * rho_2 * std::pow(v_2, 2)
	};

	return deltaP;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
