/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	Functions used for thermal analysis of the finned heat exchanger.

\*---------------------------------------------------------------------------*/

#ifndef THERMAL_H
#define THERMAL_H

#include <cmath>
#include <map>

// * * * * * * * * * * * Outer heat transfer coefficient * * * * * * * * * * //

// compute chi
double computeChi
(
	const double psi_f,			// fin characteristic
	const double s_12			// ratio s_1 / s_2
)
{
	double chi {s_12 - 1.26 / psi_f - 2.0};

	return chi;
}


// compute the exponent n
double computeN
(
	const double chi,
	const double psi_f			// fin characteristic
)
{
	double n {0.7 + 0.08 * std::tanh(chi) + 0.005 * psi_f};

	return n;
}


// compute C_q
double computeC_q
(
	const double chi,
	const double psi_f			// fin characteristic
)
{
	double C_q
	{
		(1.36 - std::tanh(chi)) * (1.1 / (psi_f + 8.0) - 0.014)
	};

	return C_q;
}


// compute correction for the number of tubes in the longitudinal direction (C_z)
double computeC_z
(
	const double s_12,			// ratio s_1 / s_2
	const double z_2			// number of tubes in longitudinal direction
)
{
	double C_z;

	if (z_2 < 8.0)
	{
		if (s_12 < 2.0)
			C_z = 3.15 * std::pow(z_2, 0.05) - 2.5;
		else
			C_z = 3.5 * std::pow(z_2, 0.03) - 2.72;
	}
	else
		C_z = 1.0;


	return C_z;
}


// compute heat transfer coefficient (alpha_2)
double computeAlpha_2
(
	const double C_q,
	const double C_z,			// tube row correction coefficient
	const double d_0,			// outer fin-tube diameter
	const double lambda_2,		// heat conductance (outer)
	const double n,
	const double Pr_2,			// Prandtl number (outer)
	const double Re_2			// Reynolds number (outer)
)
{
	double alpha_2
	{
		1.13 * C_q * C_z * lambda_2 * std::pow(Re_2, n) * std::pow(Pr_2, 0.33) / d_0
	};

	return alpha_2;
}


// compute equivalent diameter (d_eq)
double computeD_eq
(
	const double c_ekv			// equivalent square fin height
)
{
	double d_eq {1.13 * c_ekv};

	return d_eq;
}


// compute equivalent fin height (l_eq)
double computeL_eq
(
	const double d_eq,			// equivalent diameter
	const double d_0			// outer fin-tube diameter
)
{
	double l_eq {0.5 * (d_eq - d_0)};

	return l_eq;
}


// compute conventional fin height (l_cl)
double computeL_cl
(
	const double d_eq,			// equivalent diameter
	const double d_0,			// outer fin-tube diameter
	const double l_eq			// equivalent fin height
)
{
	double l_cl
	{
		l_eq * (1.0 + (0.191 + 0.054 * d_eq / d_0) * std::log(d_eq / d_0))
	};

	return l_cl;
}


// compute fin parameter (beta)
double computeBeta
(
	const double alpha_2,		// heat transfer coefficient
	const double delta_f,		// fin thickness
	const double lambda_f		// fin heat conductance
)
{
	double beta
	{
		std::sqrt((2.0 * alpha_2) / (delta_f * lambda_f))
	};

	return beta;
}


// compute theoretical fin efficiency (E)
double computeE
(
	const double beta,			// fin parameter
	const double l_cl			// conventional fin height
)
{
	double E {std::tanh(beta * l_cl) / (beta * l_cl)};

	return E;
}


// compute efficiency correction factor (C_E)
double computeC_E
(
	const double beta,			// fin parameter
	const double d_eq,			// equivalent diameter
	const double d_0,			// outer fin-tube diameter
	const double l_eq			// equivalent fin height
)
{
	double C_E
	{
		1.0 - 0.016 * (d_eq / d_0 - 1.0) * (1.0 + std::tanh(2.0 * beta * l_eq - 1.0))
	};

	return C_E;
}


// compute fin shape correction factor (C_f)
double computeC_f()
{
	double C_f {1.0};

	return C_f;
}


// compute reduced heat transfer coefficient (alpha_2_r)
double computeAlpha_2_r
(
	const double AA_f,			// ratio A_f / A
	const double AA_t,			// ratio A_t / A
	const double alpha_2,		// heat transfer coefficient
	const double C_E,			// efficiency correction factor
	const double C_f,			// fin shape correction factor
	const double E				// theoretical fin efficiency
)
{
	double alpha_2_r {(AA_f * E * C_E * C_f + AA_t) * alpha_2};

	return alpha_2_r;
}


// * * * * * * * * * * Internal heat transfer coefficient  * * * * * * * * * //

// compute inlet water velocity
inline double computeV_1
(
	const double d_in,			// tube inner diameter
	const double m_1,			// water mass flow rate
	const double rho_1,			// water density
	const double z_2			// number of tubes in longitudinal direction
)
{
	double v_1 
	{
		m_1 / (rho_1 * z_2 * 0.25 * std::pow(d_in, 2) * M_PI)
	};

	return v_1;
}


// compute average internal tube wall temperature (t_in)
double computeT_in
(
	const double A_tot,			// total heat transfer surface area
	const double AA_in,			// ratio A / A_in
	const double alpha_1,		// internal heat transfer coefficient
	const double Q,				// power
	const double t_1_avg		// average temperature (water)
)
{
	double t_in
	{
		t_1_avg - Q / (A_tot / AA_in * alpha_1)
	};

	return t_in;
}


// compute viscosity correction factor (C_mu)
double computeC_mu
(
	const double mu_1,		// water viscosity
	const double t_in		// average temperature (water)
)
{
	const std::map<double, double> coeffs
	{
	//	t		mu
		{25.0,	8.8999e-4},
		{30.0,	7.9731e-4},
		{35.0,	7.1932e-4},
		{40.0,	6.5301e-4},
		{45.0,	5.9612e-4},
		{50.0,	5.4692e-4},
		{60.0,	4.6649e-4},
		{70.0,	4.04e-4},
		{80.0,	3.5446e-4}
	};

	double mu_in
	{
		interpolate(coeffs, t_in)
	};
	double C_mu
	{
		std::pow(mu_1 / mu_in, 0.11)
	};

	return C_mu;
}


// compute internal heat transfer coefficient (alpha_1)
double computeAlpha_1
(
	const double C_mu,			// viscosity correction factor
	const double d_in,			// inner fin-tube diameter
	const double lambda_1,		// heat conductance (internal)
	const double Pr_1,			// Prandtl number (internal)
	const double Re_1			// Reynolds number (internal)
)
{
	double epsilon
	{
		1.0 + 900.0 / Re_1
	};
	double zeta
	{
		std::pow((1.82 * std::log10(Re_1) - 1.64), -2)
	};

	double alpha_1
	{
		lambda_1 / d_in * (0.125 * zeta * Re_1 * Pr_1 * C_mu) /
		(epsilon + 4.5 * std::pow(zeta, 0.5) * (std::pow(Pr_1, 0.666) - 1.0))
	};

	return alpha_1;
}


// * * * * * * * * * * * Tube wall thermal resistance  * * * * * * * * * * * //

// compute tube wall thermal resistance (R_T)
double computeR_T
(
	const double delta_f,		// fin thickness
	const double delta_t,		// tube wall thickness
	const double lambda_f,		// fin conductance
	const double lambda_t,		// tube conductance
	const double R_c			// contact resistance
)
{
	double R_T {delta_f / lambda_f + R_c + delta_t / lambda_t};

	return R_T;
}


// * * * * * * * * * * Overall heat transfer coefficient * * * * * * * * * * //

// compute thermal efficiency
double computePsi()
{
	double Psi {0.95};

	return Psi;
}

// compute overall heat transfer coefficient (k)
double computeK
(
	const double AA_in,			// ratio A / A_in
	const double alpha_1,		// internal heat transfer coefficient
	const double alpha_2_r,		// reduced heat transfer coefficient
	const double Psi,			// thermal efficiency
	const double R_T			// tube wall thermal resistance
)
{
	double k {Psi / (AA_in / alpha_1 + AA_in * R_T + 1.0 / alpha_2_r)};

	return k;
}


// * * * * * * * * * * * * * * * * * Results * * * * * * * * * * * * * * * * //

// compute total heat transfer surface area (A_tot)
double computeA_tot
(
	const double deltaT,		// average temperature difference
	const double k,				// overall heat transfer coefficient
	const double Q				// power
)
{
	double A_tot {Q / (k * deltaT)};

	return A_tot;
}


// compute total heated tube length (L_t_tot)
double computeL_t_tot
(
	const double A,				// total fin area per 1m of tube
	const double A_tot			// total heat transfer surface area
)
{
	double L_t_tot {A_tot / A};

	return L_t_tot;
}


// compute total number of tubes (z)
double computeZ
(
	const double L_c,			// length of tube projection
	const double L_t_tot		// total tube length
)
{
	double z {L_t_tot / L_c};

	return z;
}


// compute number of tubes in longitudinal direction (z_2)
double computeZ_2
(
	const double z,				// total number of tubes
	const double z_1			// number of tubes in transverse direction
)
{
	double z_2 {std::ceil(z / z_1)};

	if
	(
		static_cast<int>(z_2) % 2 != 0
	)
		z_2 += 1.0;

	return z_2;
}


// * * * * * * * * * * * * * * * * Real values * * * * * * * * * * * * * * * //

// compute real total heated tube length (L_t_tot_real)
double computeL_t_tot_real
(
	const double L_c,			// length of tube projection
	const double z_1,			// number of tubes in transverse direction
	const double z_2			// number of tubes in longitudinal direction
)
{
	double L_t_tot_real
	{
		L_c * z_1 * z_2
	};

	return L_t_tot_real;
}


// compute real total heat transfer surface area (A_tot_real)
double computeA_tot_real
(
	const double A,				// total fin area per 1m of tube
	const double L_t_tot_real	// real total tube length
)
{
	double A_tot_real
	{
		A * L_t_tot_real
	};

	return A_tot_real;
}


// compute real number of tubes (z_real)
double computeZ_real
(
	const double z_1,			// number of tubes in transverse direction
	const double z_2			// number of tubes in longitudinal direction
)
{
	double z_real
	{
		z_1 * z_2
	};

	return z_real;
}
	



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
