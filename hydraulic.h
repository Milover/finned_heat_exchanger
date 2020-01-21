/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	Functions used for hydraulic analysis of the finned heat exchanger.

\*---------------------------------------------------------------------------*/

#ifndef HYDRAULIC_H
#define HYDRAULIC_H

#include <cmath>
#include <limits>
#include <map>

#include "general.h"

// * * * * * * * * * * * * * * Friction pressure loss  * * * * * * * * * * * //

// compute total coil length (L_0)
double computeL_0
(
	const double delta_H,		// housing thickness
	const double delta_W,		// weld gap/thickness
	const double L_a,			// tube end allowance
	const double L_c,			// length of tube projection
	const double L_swc,			// supply/intake collector segment length
	const double l_prot_t,		// protrusion length (collector)
	const double s_1,			// transverse tube spacing
	const double z_1			// number of tubes in transverse direction
)
{
	double L_0			// total coil length
	{
		z_1 * (L_c + 2.0 * (delta_W + delta_H + L_a))
	  + 0.5 * s_1 * M_PI * (z_1 - 1.0)
	  + 2.0 * (L_swc + l_prot_t)
	};

	return L_0;
}


// compute friction loss coefficient via Colebrook-White (lambda)
double computeLambda
(
	const double d_in,			// inner tube diameter
	double k_rough,				// absolute roughness
	double Re_1					// Reynolds number (inner)
)
{
	int maxIter {1000};

	double f_old {0.015};		// friction coefficient (old)
	double f_new {0.0};			// friction coefficient (new)
	double x;

	int iter {0};
	while
	(
		std::abs(f_old - f_new) > std::numeric_limits<double>::epsilon() &&
		iter <= maxIter
	)
	{
		x = k_rough / (3.707 * d_in) + 2.5232 / (Re_1 * std::sqrt(f_old));

		f_new = std::pow(-0.8686 * std::log(x), -2);

		f_old = f_new;
		iter++;
	}

	return f_new;
}


// compute friction pressure loss (deltaP_f)
double computeDeltaP_f
(
	const double d_in,			// inner tube diameter
	const double L_0,			// total coil length
	const double lambda,		// friction coefficient
	const double rho_1,			// density (water)
	const double v_1			// inlet velocity(water)
)
{
	double deltaP_f
	{
		lambda * L_0 / d_in * 0.5 * rho_1 * std::pow(v_1, 2)
	};

	return deltaP_f;
}


// * * * * * * * * * * * * * * * Local pressure loss * * * * * * * * * * * * //

// compute coil inlet loss coefficient (zeta_in)
double computeZeta_in
(
	const double d_in,			// inner tube diameter
	const double d_col_in		// inner collector tube diameter
)
{
	double zeta_in;			// supply collector with side inlet & radial outlet

	if
	(
		d_in / d_col_in <= 0.1
	)
		zeta_in = 0.5;
	else
		zeta_in = 0.7;

	return zeta_in;
}


// compute coil outlet loss coefficient (zeta_out)
double computeZeta_out()
{
	double zeta_out {1.1};	// intake collector with radial inlet & side outlet

	return zeta_out;
}


// compute tube bend loss coefficient (zeta_R)
double computeZeta_R		// p.218 - ISBN: 987-953-7370-02-2
(
	const double alpha,			// bend angle
	const double d_in,			// inner tube diameter
	const double R				// bend (mid-line) radius
)
{
	const std::map<double, double> lossCoeffs		// for 90° bend
	{
	//	R/d_in	zeta
		{1.0,	0.27},
		{1.5,	0.20},
		{2.0,	0.15},
		{3.0,	0.13},
		{4.0,	0.10},
		{5.0,	0.10},
		{6.0,	0.10},
		{10.0,	0.11}
	};
	const std::map<double, double> angleCoeffs		// for 0...180° bend
	{
	//	alpha	n
		{0,		0.00},
		{30,	0.4},
		{60,	0.7},
		{90,	1.0},
		{120,	1.3},
		{150,	1.5},
		{180,	1.7}
	};

	double Rd_in
	{
		R / d_in
	};
	double zeta
	{
		interpolate(lossCoeffs, Rd_in)
	};
	double n
	{
		interpolate(angleCoeffs, alpha)
	};

	double zeta_R
	{
		zeta * n
	};

	return zeta_R;
}


// compute local pressure loss (deltaP_loc)
double computeDeltaP_loc
(
	const double rho_1,			// density (water)
	const double v_1,			// inlet velocity(water)
	const double z_1,			// number of tubes in transverse direction
	const double zeta_in,		// coil inlet loss coefficient
	const double zeta_out,		// coil outlet loss coefficient
	const double zeta_R_col,	// tube-collector bend loss coefficient
	const double zeta_R_t		// tube end bend loss coefficient
)
{
	double deltaP_loc
	{
		0.5 * rho_1 * std::pow(v_1, 2) *
		(zeta_in + zeta_out + 2.0 * zeta_R_col + (z_1 - 1.0) * zeta_R_t)
	};

	return deltaP_loc;
}


// * * * * * * * * * * * * * * Collector pressure loss * * * * * * * * * * * //

// compute collector cross-section area (A_col)
double computeA_col
(
	const double d_col_in		// inner collector tube diameter
)
{
	double A_col
	{
		0.25 * M_PI * std::pow(d_col_in, 2)
	};

	return A_col;
}


// compute collector max velocity (v_col)
double computeV_col
(
	const double A_col,			// collector cross-section area
	const double m_1,			// water mass flow rate
	const double rho_1			// water density
)
{
	double v_col
	{
		m_1 / (rho_1 * A_col)
	};

	return v_col;
}


// compute supply collector coefficient (B_sup)
double computeB_sup()
{
	double B_sup {0.8};	// end inlet with full cross section

	return B_sup;
}


// compute intake collector coefficient (B_int)
double computeB_int()
{
	double B_int {2.0};	// end outlet

	return B_int;
}


// compute collector pressure loss (deltaP_col_x)
double computeDeltaP_col_x
(
	const double B,				// collector coefficient
	const double rho_1,			// collector max velocity
	const double v_col			// water density
)
{
	double deltaP_col
	{
		B * 0.5 * rho_1 * std::pow(v_col, 2)
	};

	return deltaP_col;
}


// compute total collector pressure loss (deltaP_col)
double computeDeltaP_col
(
	const double deltaP_col_sup,	// supply collector pressure loss
	const double deltaP_col_int		// supply collector pressure loss
)
{
	double deltaP_col		// Π-scheme
	{
		2.0 / 3.0 * (deltaP_col_int - deltaP_col_sup)
	};

	return deltaP_col;
}


// * * * * * * * * * * * * * * * Total pressure loss * * * * * * * * * * * * //

// compute total pressure loss (deltaP)
double computeDeltaP
(
	const double deltaP_f,			// friction pressure loss
	const double deltaP_loc,		// local pressure loss
	const double deltaP_col		// total collector pressure loss
)
{
	double deltaP
	{
		deltaP_f + deltaP_loc + deltaP_col
	};

	return deltaP;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
