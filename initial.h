/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	Initial and chosen parameters of the heat exchanger as well as relevant
	physical properties.

\*---------------------------------------------------------------------------*/

#ifndef INITIAL_H
#define INITIAL_H

#include <cmath>

// * * * * * * * * * * * * * * Assigned parameters * * * * * * * * * * * * * //

double Q			{55000};	// power

double t_1_in		{75};		// inlet water temperature
double t_1_out		{65};		// outlet water temperature
double t_1_avg					// average water temperature
{
	0.5 * (t_1_in + t_1_out)
};
double p_1			{0.4};		// nominal gauge pressure (water) [MPa]

double t_2_in		{-15};		// inlet air temperature
double t_2_out		{20};		// outlet air temperature
double t_2_avg					// average air temperature
{
	0.5 * (t_2_in + t_2_out)
};
double p_2			{0.0};		// nominal gauge pressure (air) [MPa]

double deltaT		{64.05};	// average heat transfer temperature difference


// * * * * * * * * * * * * * * Tube characteristics  * * * * * * * * * * * * //

// EN 10305-1
// SPH 235
double d			{0.0213};	// outer tube dimeter
double delta_t		{0.002};	// tube wall thickness
double d_in						// inner tube diameter
{
	d - 2.0 * delta_t
};
double k_rough		{5e-5};		// absolute roughness

double tol_d_t		{8e-5};		// outside diameter ±tolerance
double tol_t					// thickness ±tolerance
{
	std::min(0.1 * delta_t, 0.2)
};


// * * * * * * * * * * * * * * Fin characteristics * * * * * * * * * * * * * //

// aluminum
double s_f			{0.007};	// fin spacing
double delta_f		{0.0008};	// fin thickness
double d_0						// fin-tube diameter
{
	d + 2.0 * delta_f
};

double tol_f		{7e-5};		// fin thickness tolerance


// * * * * * * * * * * * * * * Coil characteristics  * * * * * * * * * * * * //

double delta_H		{0.003};	// housing thickness
double delta_W		{0.002};	// weld gap/thickness

double L_a			{0.01};		// tube end allowance
double L_swc		{0.05};		// supply/intake collector tube segment
double l_prot_t		{0.0};		// protrusion length (collector)

double alpha_col	{30.05};	// tube-collector bend angle
double alpha_t		{180};		// tube end bend angle


// * * * * * * * * * * * * * Collector characteristics * * * * * * * * * * * //

// EN 10216-4
// SPH 235
double d_col		{0.0761};	// outer collector tube diameter (DN 50)
double delta_col	{0.0029};	// collector tube wall thickness
double d_col_in					// inner collector tube diameter
{
	d_col - 2.0 * delta_col
};
double tol_col					// thickness ±tolerance
{
	std::min(0.1 * delta_col, 0.2)
};

double delta_w_col_t	{0.002};// coil tube - collector weld thickness (a)


// * * * * * * * * * * * * * * Dish characteristics  * * * * * * * * * * * * //

// DIN 28011
// SPH 235
double d_dish		{0.0761};	// dish outer diameter
double delta_dish	{0.003};	// dish thickness
double d_dish_in				// dish inner diameter
{
	d_dish - 2.0 * delta_dish
};
double R_dish		{d_dish};	// dish sphere radius
double r_dish					// dish knuckle radius
{
	0.1 * d_dish
};
double tol_dish		{3e-4};		// thickness -tolerance


// * * * * * * * * * * * * * * Flange characteristics  * * * * * * * * * * * //

// EN 1092-1/01/DN 50/PN 63/SPH 235


// * * * * * * * * * * * * Structural characteristics  * * * * * * * * * * * //

// SPH235
double R_m			{360};		// min tensile strength at 20°C [MPa]
double R_p02		{235};		// yield strength at 20°C [MPa]
double R_p02t		{177.94};	// yield strength at 75°C [MPa]

double f_a						// allowable stress at test temperature
{
	std::min(R_p02 / 1.5, R_m / 2.4)
};
double f_d						// allowable stress at design temperature
{
	std::min(R_p02t / 1.5, R_m / 2.4)
};
double f_test					// allowable test stress
{
	R_p02 / 1.05
};

double R_p02_f		{100};		// yield strength at 20°C (aluminum) [MPa]

double p_contact	{10};		// nominal contact pressure at working temp. [MPa]


// * * * * * * * * * * * * * * Physical properties * * * * * * * * * * * * * //

// SPH235
double lambda_t		{50};		// tube conductance
double E_t			{2.1e5};	// tube Young elasticity modulus [MPa]
double nu_t			{0.3};		// tube Poisson coefficient
double beta_t		{1.1e-5};	// tube thermal expansion coefficient
double Rz_t			{1.6e-6};	// tube average surface roughness

// aluminum
double lambda_f		{229};		// fin conductance
double E_f			{7e4};		// fin Young elasticity modulus [MPa]
double nu_f			{0.32};		// fin Poisson coefficient
double beta_f		{2.3e-5};	// fin thermal expansion coefficient
double Rz_f			{0.8e-6};	// fin average surface roughness

double R_c			{2e-4};		// thermal contact resistance

// water @ 70°C, 5bar
double rho_1		{978.0};	// water density
double mu_1			{4.04e-4};	// water viscosity
double lambda_1		{0.66334};	// water conductance
double Pr_1			{2.5502};	// water Prandtl number

double m_1			{1.314};	// water mass flow rate

// air @ 2.5°C, 1bar
double rho_2		{1.264};	// air density
double mu_2			{1.736e-5};	// air viscosity
double lambda_2		{0.0234};	// air conductance
double Pr_2			{0.7181};	// air Prandtl number

double m_2			{1.563};	// air mass flow rate


// * * * * * * * * * * * * * * * Tuning parameters * * * * * * * * * * * * * //

double z_1			{4};		// number of tubes in transverse direction
double tau			{0.6};		// the ratio c_f / s_2


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
