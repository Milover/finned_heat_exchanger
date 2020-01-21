/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	Functions used for structural analysis of the finned heat exchanger.
	Based on EN 13445-3:2009 unless noted otherwise.

\*---------------------------------------------------------------------------*/

#ifndef STRUCTURAL_H
#define STRUCTURAL_H

#include <cmath>
#include <limits>
#include <stdexcept>

// * * * * * * * * * * * Tube stress and pressure  * * * * * * * * * * * * * //

// compute errosion thickness (C_err)
double computeC_err		// 6.6 - ISBN: 978-0-12-804397-4
(
	const double D,				// outer diameter
	const double tol_min		// max. minus tolerance
)
{
	double C_11 {tol_min};
	double C_12 {0.0};		// seamless tube
	double C_21;
	if			// D <= 0.032
	(
		D < 0.032 ||
		std::abs(D - 0.032) < std::numeric_limits<double>::epsilon()
	)
		C_21 = 0.0;
	else 
		C_21 = 5e-4;

	double C_22 {4e-4};		// min. required

	double C_err
	{
		C_11 + C_12 + C_21 + C_22
	};

	return C_err;
}


// compute analysis thickness (e_a)
double computeE_a
(
	const double C_err,			// errosion thickness
	const double delta			// thickness
)
{
	double e_a
	{
		delta - C_err
	};

	return e_a;
}


// compute min thickness (e_min)
double computeE_min
(
	const double delta,			// thickness
	const double tol_min		// max. minus tolerance
)
{
	double e_min
	{
		delta - tol_min
	};

	return e_min;
}


// compute test pressure (p_test)
double computeP_test	// 10.2.3.3.1-1 & 10.2.3.3.1-2 - EN 13445-5:2009
(
	const double f_a,			// max allowed stress of tube (test T)
	const double f_d,			// max allowed stress of tube (design T)
	const double p_1,			// design pressure
	const double p_max			// max vessel pressure
)
{
	double p_t_1
	{
		1.25 * p_1 * f_a / f_d
	};
	double p_t_2
	{
		1.43 * p_max
	};

	double p_test
	{
		std::max(p_t_1, p_t_2)
	};

	return p_test;
}


// compute tube required thickness (E_t)
double computeE_t
(
	const double D_i,			// inner diameter
	const double p_1,			// nominal gauge pressure (water)
	const double f_d			// max allowed stress of tube
)
{
	double e_t		// 7.4-1
	{
		(p_1 * D_i) / (2.0 * f_d - p_1)
	};

	return e_t;
}


// compute max tube pressure (p_max)
double computeP_max
(
	const double C_err,			// errosion thickness
	const double D,				// outer diameter
	const double D_i,			// inner diameter
	const double delta,			// wall thickness
	const double sigma_t_all	// max allowed stress
)
{
	double p_max		// 7.4-3
	{
		2.0 * sigma_t_all * (delta - C_err) / (0.5 * (D + D_i))
	};

	return p_max;
}


// * * * * * * * * * * * * * * Collector tube  * * * * * * * * * * * * * * * //

// compute collector inside curvature radius (r_is_col)
double computeR_is_col
(
	const double d_col,			// collector outer diameter
	const double e_a_col		// collector analysis thickness
)
{
	double r_is_col		// 9.5-3
	{
		0.5 * d_col - e_a_col
	};

	return r_is_col;
}


// collector isolated hole check (L_b_req)
double computeL_b_col_req
(
	const double d,				// coil tube outer diameter
	const double e_a_col,		// collector analysis thickness
	const double r_is_col		// collector inside curvature radius
)
{
	double a
	{
		0.5 * d
	};
	double l_so
	{
		std::sqrt((2.0 * r_is_col + e_a_col) * e_a_col)
	};

	double L_b_col_req	// 9.5-1
	{
		2.0 * (a + l_so)
	};

	return L_b_col_req;
}


// compute min discontinuity distance (w_min_col)
double computeW_min_col
(
	const double e_a_col,		// collector analysis thickness
	const double r_is_col		// collector inside curvature radius
)
{
	double x
	{
		(2.0 * r_is_col + e_a_col) * e_a_col
	};

	double w_min_col	// 9.7-1
	{
		std::max(0.2 * std::sqrt(x), 3.0 * e_a_col)
	};

	return w_min_col;
}


// compute max collector pressure (p_max_col)
double computeP_max_col
(
	const double d,				// tube outer diameter
	const double d_in,			// tube inner diameter
	const double d_col,			// collector outer diameter
	const double delta_t,		// tube wall thickness
	const double delta_w_col_t,	// coil tube - collector weld thickness (a)
	const double e_a_col,		// collector analysis thickness
	const double l_prot_t,		// coil tube protrusion length
	const double r_is_col,		// collector inside curvature radius
	const double sigma_col_all,	// max allowed stress of collector
	const double sigma_t_all	// max allowed stress of tube
)
{
	double l_so			// 9.5-92
	{
		std::sqrt(((d_col - 2.0 * e_a_col) + e_a_col) * e_a_col)
	};
	double l_bo			// 9.5-76
	{
		std::sqrt((d - delta_t) * delta_t)
	};
	double l_bi_eff		// 9.5-77
	{
		std::min(l_prot_t, l_bo)
	};

	double A_ps			// 9.5-94
	{
		r_is_col * (l_so + 0.5 * d)
	};
	double A_pb			// 9.5-84
	{
		0.5 * d_in * (l_bo + e_a_col)
	};
	double A_fs			// 9.5-81
	{
		l_so * e_a_col
	};
	double A_fb			// 9.5-80
	{
		delta_t * (l_bo + l_bi_eff + e_a_col)
	};
	double A_fw
	{
		std::pow(delta_w_col_t, 2)
	};

	double p_max_col	// 9.5-10
	{
		((A_fs + A_fw) * sigma_col_all + A_fb * sigma_t_all) /
		(A_ps + A_pb + 0.5 * (A_fs + A_fw + A_fb))
	};

	return p_max_col;
}


// * * * * * * * * * * * * * * Collector dish  * * * * * * * * * * * * * * * //

// compute beta
double computeBeta	// 7.5.3.5
(
	const double d_dish_in,		// dish inner diameter
	const double e,				// thickness
	const double R,				// dish radius
	const double r				// dish knuckle radius
)
{
	double beta;

	double X
	{
		r / d_dish_in
	};
	double Y
	{
		std::min(e / R, 0.04)
	};
	double Z
	{
		std::log10(1.0 / Y)
	};
	double N
	{
		1.006 - 1.0 / (6.2 + std::pow(90.0 * Y, 4))
	};

	double beta_06
	{
		N *
		(
		   - 0.3635 * std::pow(Z, 3)
		   + 2.2124 * std::pow(Z, 2)
		   - 3.2937 * Z
		   + 1.8873
		)
	};
	double beta_1
	{
		N *
		(
		   - 0.1833 * std::pow(Z, 3)
		   + 1.0383 * std::pow(Z, 2)
		   - 1.2943 * Z
		   + 0.837
		)
	};
	double beta_2
	{
		std::max(0.95 * (0.56 - 1.94 * Y - 82.5 * std::pow(Y, 2)), 0.5)
	};

	if
	(
		std::abs(X - 0.06) < std::numeric_limits<double>::epsilon()
	)
		beta = beta_06;
	else if
	(
		X > 0.06 && X < 0.1
	)
		beta = 25.0 * ((0.1 - X) * beta_06 + (X - 0.06) * beta_1);
	else if
	(
		std::abs(X - 0.1) < std::numeric_limits<double>::epsilon()
	)
		beta = beta_1;
	else if
	(
		X > 0.1 && X < 0.2
	)
		beta = 10.0 * ((0.2 - X) * beta_1 + (X - 0.1) * beta_2);
	else if
	(
		std::abs(X - 0.2) < std::numeric_limits<double>::epsilon()
	)
		beta = beta_2;

	return beta;
}


// compute collector dish required thickness (e_dish)
double computeE_dish
(
	const double d_dish_in,		// dish inner diameter
	const double p_1,			// nominal gauge pressure (water)
	const double R,				// dish radius
	const double R_p02t_dish,	// yield strength at temperature t
	const double r,				// dish knuckle radius
	const double sigma_dish_all	// max allowed stress of dish
)
{
	double e_s			// 7.5-1
	{
		(p_1 * R) / (2.0 * sigma_dish_all - 0.5 * p_1)
	};

	double e_y {e_s};		// initial
	double e_y_temp {0.0};	// dummy
	while
	(
		std::abs(e_y - e_y_temp) > std::numeric_limits<double>::epsilon()
	)
	{
		double beta
		{
			computeBeta(d_dish_in, e_y, R, r)
		};

		e_y_temp = e_y;
		e_y = beta * p_1 * (0.75 * R + 0.2 * d_dish_in) / sigma_dish_all;
	}

	double sigma_b		// 7.5-4
	{
		R_p02t_dish / 1.5
	};
	double x
	{
		p_1 / (111.0 * sigma_b) * std::pow(d_dish_in / r, 0.825)
	};
	double e_b			// 7.5-3
	{
		(0.75 * R + 0.2 * d_dish_in) * std::pow(x, 2.0 / 3.0)
	};

	double e_dish
	{
		std::min(std::min(e_s, e_y), e_b)
	};

	return e_dish;
}


// compute dish max pressure (p_max_dish)
double computeP_max_dish
(
	const double C_err_dish,	// errosion thickness
	const double d_dish_in,		// dish inner diameter
	const double delta_dish,	// dish wall thickness
	const double R,				// dish radius
	const double R_p02t_dish,	// yield strength at temperature t
	const double r,				// dish knuckle radius
	const double sigma_dish_all	// max allowed stress
)
{
	double e_a
	{
		delta_dish - C_err_dish
	};

	double p_s			// 7.5-6
	{
		2.0 * sigma_dish_all * e_a / (R + 0.5 * e_a)
	};

	double beta
	{
		computeBeta(d_dish_in, e_a, R, r)
	};
	double p_y			// 7.5-7
	{
		sigma_dish_all * e_a / (beta * (0.75 * R + 0.2 * d_dish_in))
	};

	double sigma_b		// 7.5-4
	{
		R_p02t_dish / 1.5
	};
	double x
	{
		e_a / (0.75 * R + 0.2 * d_dish_in)
	};
	double p_b			// 7.5-8
	{
		111.0 * sigma_b * std::pow(x, 1.5) * std::pow(r / d_dish_in, 0.825)
	};

	double p_max_dish
	{
		std::min(std::min(p_s, p_y), p_b)
	};

	return p_max_dish;
}


// * * * * * * * * * * * * * Coil-fin tolerance fit  * * * * * * * * * * * * //

// compute tolerance fit internal characteristic (Q_i)
double computeQ_i
(
	const double D_F,			// tolerance fit diameter
	const double D_i			// tolerance fit internal diameter
)
{
	double Q_i
	{
		D_i / D_F
	};

	return Q_i;
}


// compute tolerance fit internal characteristic (Q_o)
double computeQ_o
(
	const double D_F,			// tolerance fit diameter
	const double D_o			// tolerance fit internal diameter
)
{
	double Q_o
	{
		D_F / D_o
	};

	return Q_o;
}


// compute tolerance fit characteristic (K)
double computeK
(
	const double E_i,			// internal material Young modulus
	const double E_o,			// outer material Young modulus
	const double Q_i,			// tolerance fit internal characteristic
	const double Q_o,			// tolerance fit outer characteristic
	const double nu_i,			// internal material Poission coefficient
	const double nu_o			// outer material Poission coefficient
)
{
	double C_i
	{
		(1.0 + std::pow(Q_i, 2)) / (1.0 - std::pow(Q_i, 2))
	};
	double C_o
	{
		(1.0 + std::pow(Q_o, 2)) / (1.0 - std::pow(Q_o, 2))
	};

	double K
	{
		E_o / E_i * (C_i - nu_i) + C_o + nu_o
	};

	return K;
}


// compute allowed reduced overlap (Z_all)
double computeZ_all
(
	const double E_o,			// outer material Young modulus
	const double K,				// tolerance fit characteristic
	const double Q_i,			// tolerance fit internal characteristic
	const double Q_o,			// tolerance fit outer characteristic
	const double R_p02_i,		// internal material yield strength
	const double R_p02_o,		// outer material yield strength
	const double S				// tolerance fit safety factor
)
{
	double Z_i
	{
		(1.0 - std::pow(Q_i, 2)) / (std::sqrt(3.0) * S * E_o) * R_p02_i
	};
	double Z_o
	{
		K * (1.0 - std::pow(Q_o, 2)) / (std::sqrt(3.0) * S * E_o) * R_p02_o
	};

	double Z_all
	{
		std::min(Z_i, Z_o)
	};

	return Z_all;
}


// compute min. lower tolerance of hole (EI_min)
double computeEI_min
(
	const double D_F,			// tolerance fit diameter
	const double es,			// upper tolerance of shaft
	const double Rz_i,			// average surface roughness of internal material
	const double Rz_o,			// average surface roughness of outer material
	const double Z_all			// allowed reduced overlap
)
{
	double EI_max
	{
		es - Z_all * D_F  - 0.8 * (Rz_i + Rz_o)
	};

	return EI_max;
}


// compute max. upper tolerance of hole (ES_max)
double computeES_max
(
	const double beta_i,		// thermal expansion coefficient (internal)
	const double beta_o,		// thermal expansion coefficient (outer)
	const double D_F,			// tolerance fit diameter
	const double E_o,			// outer material Young modulus
	const double ei,			// lower tolerance of shaft
	const double K,				// tolerance fit characteristic
	const double p_contact,		// min. contact pressure
	const double Rz_i,			// average surface roughness of internal material
	const double Rz_o,			// average surface roughness of outer material
	const double S,				// tolerance fit safety factor
	const double t_max			// max. temperature
)
{
	double B_i
	{
		1.0 + beta_i * (S * t_max - 20.0)
	};
	double B_o
	{
		1.0 + beta_o * (S * t_max - 20.0)
	};

	double R
	{
		0.8 * (Rz_i + Rz_o)
	};
	double C
	{
		K * p_contact / E_o + 0.5 * (B_o - B_i)
	};

	double ES_min
	{
		(ei * (1.0 + B_i) - D_F * C - R) / (1.0 + B_o)
	};

	return ES_min;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
