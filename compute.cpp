/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	Design procedure of the finned heat exchanger.

\*---------------------------------------------------------------------------*/

#include <map>

#include "initial.h"
#include "general.h"
#include "dimensioning.h"
#include "thermal.h"
#include "aerodynamic.h"
#include "hydraulic.h"
#include "structural.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{

// * * * * * * * * * * * * * * * Dimensioning  * * * * * * * * * * * * * * * //

	std::map<std::string, double> geometric;

	if (!hasKey(geometric, "psi_f"))
		geometric["psi_f"] = 20;		// initial
	if (!hasKey(geometric, "s_12_opt"))
		geometric["s_12_opt"] = 2.0;	// initial
	if (!hasKey(geometric, "s_1"))
		geometric["s_1"] = 1.0;			// dummy
	if (!hasKey(geometric, "s_2"))
		geometric["s_2"] = 1.0;			// dummy
	if (!hasKey(geometric, "z_2"))
		geometric["z_2"] = 6.0;			// initial
	if (!hasKey(geometric, "v_2"))
		geometric["v_2"] = 100.0;		// dummy

	if (!hasKey(geometric, "c_ekv"))	// initial
		geometric["c_ekv"] = computeC_ekv
		(
			d_0,
			delta_f,
			geometric["psi_f"],
			s_f
		);

	double s_1_old {0.0};
	double s_2_old {0.0};
	while
	(
		geometric["v_2"] > 5.0 ||
		(hasKey(geometric, "s_1") && geometric["s_1"] < 5.0 * d)
	)
	{
		while
		(
			std::abs(s_1_old - geometric["s_1"]) > 1e-6 ||
			std::abs(s_2_old - geometric["s_2"]) > 1e-6
		)
		{
			geometric["l_ekv"] = computeL_ekv
			(
				geometric["c_ekv"],
				d_0
			);

			s_2_old = geometric["s_2"];
			geometric["s_2"] = computeS_2
			(
				geometric["c_ekv"],
				geometric["s_12_opt"],
				tau,
				z_1
			);
			s_1_old = geometric["s_1"];
			geometric["s_1"] = computeS_1
			(
				geometric["s_12_opt"],
				geometric["s_2"]
			);

			geometric["s_3"] = computeS_3
			(
				geometric["s_1"],
				geometric["s_2"]
			);
			geometric["psi_f"] = computePsi_f
			(
				geometric["c_ekv"],
				d_0,
				delta_f,
				s_f
			);

			geometric["H_r"] = computeH_r
			(
				geometric["s_1"],
				geometric["s_2"],
				tau,
				z_1
			);
			geometric["H"] = computeHeight
			(
				geometric["H_r"]
			);
			geometric["W"] = computeWidth
			(
				geometric["H"]
			);

			geometric["A_f"] = computeA_f
			(
				geometric["c_ekv"],
				d_0,
				s_f
			);
			geometric["A_t"] = computeA_t
			(
				d_0,
				delta_f,
				s_f
			);
			geometric["A"] = computeA
			(
				geometric["A_f"],
				geometric["A_t"]
			);
			geometric["A_ft"] = computeA_ft(d_0);
			geometric["A_in"] = computeA_in(d_in);

			geometric["AA_f"] = computeAA_f
			(
				geometric["A"],
				geometric["A_f"]
			);
			geometric["AA_t"] = computeAA_t
			(
				geometric["A"],
				geometric["A_t"]
			);
			geometric["AA_in"] = computeAA_in
			(
				geometric["A"],
				geometric["A_in"]
			);

			geometric["d_cl"] = computeD_cl
			(
				d_0,
				delta_f,
				geometric["l_ekv"],
				s_f
			);
			geometric["phi_cl"] = computePhi_cl
			(
				geometric["d_cl"],
				geometric["s_1"],
				geometric["s_3"]
			);
			geometric["F"] = computeF
			(
				geometric["d_cl"],
				geometric["H"],
				geometric["W"],			// L_c
				geometric["phi_cl"],
				geometric["W"],
				z_1
			);
			geometric["v_2"] = computeV_2
			(
				geometric["F"],
				m_2,
				rho_2
			);

			geometric["s_12_opt"] = computeS_12
			(
				geometric["psi_f"],
				computeRe
				(
					d_0,
					mu_2,
					rho_2,
					geometric["v_2"]
				)
			);
		}

		// reset
		s_1_old = 0.0;
		s_2_old = 0.0;

		// increase by 1%
		geometric["c_ekv"] *= 1.01;
	}

	geometric["s_12"] = geometric["s_1"] / geometric["s_2"];


// * * * * * * * * * * * * * * * * Thermal * * * * * * * * * * * * * * * * * //

	std::map<std::string, double> thermal;

	if (!hasKey(geometric, "A_tot_real"))
		geometric["A_tot_real"] = 29.5;	// initial
	if (!hasKey(thermal, "alpha_1"))
		thermal["alpha_1"] = 3870.0;	// initial

	double alpha_1_old {0.0};
	double z_2_old {0.0};
	while
	(
		std::abs(geometric["z_2"] - z_2_old) > 1.0 ||
		std::abs(thermal["alpha_1"] - alpha_1_old) > 1.0
	)
	{
		// outer heat transfer coefficient

		thermal["chi"] = computeChi
		(
			geometric["psi_f"],
			geometric["s_12"]
		);
		thermal["n"] = computeN
		(
			thermal["chi"],
			geometric["psi_f"]
		);
		thermal["C_q"] = computeC_q
		(
			thermal["chi"],
			geometric["psi_f"]
		);
		thermal["C_z"] = computeC_z
		(
			geometric["s_12"],
			geometric["z_2"]
		);
		thermal["alpha_2"] = computeAlpha_2
		(
			thermal["C_q"],
			thermal["C_z"],
			d_0,
			lambda_2,
			thermal["n"],
			Pr_2,
			computeRe
			(
				d_0,
				mu_2,
				rho_2,
				geometric["v_2"]
			)
		);

		geometric["d_eq"] = computeD_eq
		(
			geometric["c_ekv"]
		);
		geometric["l_eq"] = computeL_eq
		(
			geometric["d_eq"],
			d_0
		);
		geometric["l_cl"] = computeL_cl
		(
			geometric["d_eq"],
			d_0,
			geometric["l_eq"]
		);

		thermal["beta"] = computeBeta
		(
			thermal["alpha_2"],
			delta_f,
			lambda_f
		);
		thermal["E"] = computeE
		(
			thermal["beta"],
			geometric["l_cl"]
		);
		thermal["C_E"] = computeC_E
		(
			thermal["beta"],
			geometric["d_eq"],
			d_0,
			geometric["l_eq"]
		);
		thermal["C_f"] = computeC_f();

		thermal["alpha_2_r"] = computeAlpha_2_r
		(
			geometric["AA_f"],
			geometric["AA_t"],
			thermal["alpha_2"],
			thermal["C_E"],
			thermal["C_f"],
			thermal["E"]
		);

		// internal heat transfer coefficient

		geometric["v_1"] = computeV_1
		(
			d_in,
			m_1,
			rho_1,
			geometric["z_2"]
		);
		thermal["t_in"] = computeT_in			// average internal tube wall temperature
		(
			geometric["A_tot_real"],
			geometric["AA_in"],
			thermal["alpha_1"],
			Q,
			t_1_avg
		);
		thermal["C_mu"] = computeC_mu			// viscosity correction factor
		(
			mu_1,
			thermal["t_in"]
		);

		alpha_1_old = thermal["alpha_1"];
		thermal["alpha_1"] = computeAlpha_1		// internal heat transfer coefficient
		(
			thermal["C_mu"],
			d_in,
			lambda_1,
			Pr_1,
			computeRe
			(
				d_in,
				mu_1,
				rho_1,
				geometric["v_1"]
			)
		);

		// tube wall thermal resistance

		thermal["R_T"] = computeR_T
		(
			delta_f,
			delta_t,
			lambda_f,
			lambda_t,
			R_c
		);

		// overall heat transfer coefficient

		thermal["Psi"] = computePsi();			// thermal efficiency

		thermal["k"] = computeK					// overall heat transfer coefficient
		(
			geometric["AA_in"],
			thermal["alpha_1"],
			thermal["alpha_2_r"],
			thermal["Psi"],
			thermal["R_T"]
		);

		// results

		geometric["A_tot"] = computeA_tot		// total heat transfer surface area
		(
			deltaT,
			thermal["k"],
			Q
		);
		geometric["L_t_tot"] = computeL_t_tot	// total heated tube length
		(
			geometric["A"],
			geometric["A_tot"]
		);
		geometric["z"] = computeZ				// number of tubes
		(
			geometric["W"],			// L_c
			geometric["L_t_tot"]
		);

		z_2_old = geometric["z_2"];
		geometric["z_2"] = computeZ_2			// number of tubes in longitudinal direction
		(
			geometric["z"],
			z_1
		);

		// real values

		geometric["L_t_tot_real"] = computeL_t_tot_real	// real total heated tube length
		(
			geometric["W"],			// L_c
			z_1,
			geometric["z_2"]
		);
		geometric["A_tot_real"] = computeA_tot_real	// real total heat transfer surface area
		(
			geometric["A"],
			geometric["L_t_tot_real"]
		);
		geometric["z_real"] = computeZ_real		// real number of tubes
		(
			z_1,
			geometric["z_2"]
		);
	}


// * * * * * * * * * * * * * * * Aerodynamic * * * * * * * * * * * * * * * * //

	std::map<std::string, double> aerodynamic;

	aerodynamic["FA_tot"] = computeFA_tot
	(
		d_0,
		delta_f,
		geometric["l_eq"],
		geometric["s_1"],
		s_f
	);
	aerodynamic["d_eq"] = computeD_eq
	(
		d_0,
		delta_f,
		geometric["l_eq"],
		geometric["phi_cl"],
		geometric["s_1"],
		s_f
	);

	aerodynamic["n_air"] = computeN_air
	(
		aerodynamic["FA_tot"],
		geometric["s_12"]
	);
	aerodynamic["C_f_air"] = computeC_f_air
	(
		aerodynamic["FA_tot"],
		geometric["s_12"]
	);
	aerodynamic["C_z_air"] = computeC_z_air
	(
		geometric["z_2"]
	);

	aerodynamic["zeta_0"] = computeZeta_0
	(
		aerodynamic["C_f_air"],
		aerodynamic["C_z_air"],
		aerodynamic["n_air"],
		computeRe
		(
			aerodynamic["d_eq"],
			mu_2,
			rho_2,
			geometric["v_2"]
		)
	);

	aerodynamic["C_op"] = computeC_op();

	aerodynamic["deltaP"] = computeDeltaP
	(
		aerodynamic["C_op"],
		rho_2,
		geometric["v_2"],
		geometric["z_2"],
		aerodynamic["zeta_0"]
	);


// * * * * * * * * * * * * * * * Hydraulic * * * * * * * * * * * * * * * * * //

	std::map<std::string, double> hydraulic;

	// friction pressure loss

	hydraulic["L_0"] = computeL_0				// total coil length
	(
		delta_H,
		delta_W,
		L_a,
		geometric["W"],			// L_c
		L_swc,
		l_prot_t,
		geometric["s_1"],
		z_1
	);
	hydraulic["lambda"]	= computeLambda			// friction coefficient
	(
		d_in,
		k_rough,
		computeRe
		(
			d_in,
			mu_1,
			rho_1,
			geometric["v_1"]
		)
	);

	hydraulic["deltaP_f"] = computeDeltaP_f		// friction pressure loss
	(
		d_in,
		hydraulic["L_0"],
		hydraulic["lambda"],
		rho_1,
		geometric["v_1"]
	);

	// local pressure loss

	hydraulic["zeta_in"] = computeZeta_in		// coil inlet loss coefficient
	(
		d_in,
		d_col_in
	);
	hydraulic["zeta_out"] = computeZeta_out();	// coil outlet loss coefficient
	hydraulic["zeta_R_col"] = computeZeta_R		// tube-collector bend loss coefficient
	(
		alpha_col,
		d_in,
		0.5 * geometric["s_1"]	// R
	);
	hydraulic["zeta_R_t"] = computeZeta_R		// tube end bend loss coefficient
	(
		alpha_t,
		d_in,
		0.5 * geometric["s_1"]	// R
	);

	hydraulic["deltaP_loc"] = computeDeltaP_loc	// local pressure loss
	(
		rho_1,
		geometric["v_1"],
		z_1,
		hydraulic["zeta_in"],
		hydraulic["zeta_out"],
		hydraulic["zeta_R_col"],
		hydraulic["zeta_R_t"]
	);

	// collector pressure loss

	hydraulic["A_col"] = computeA_col			// collector cross-section area
	(
		d_col_in
	);
	hydraulic["v_col"] = computeV_col			// collector max velocity
	(
		hydraulic["A_col"],
		m_1,
		rho_1
	);
	hydraulic["B_sup"] = computeB_sup();		// supply collector coefficient
	hydraulic["B_int"] = computeB_int();		// intake collector coefficient
	hydraulic["deltaP_col_sup"] = computeDeltaP_col_x	// supply collector pressure loss
	(
		hydraulic["B_sup"],
		rho_1,
		hydraulic["v_col"]
	);
	hydraulic["deltaP_col_int"] = computeDeltaP_col_x	// intake collector pressure loss
	(
		hydraulic["B_int"],
		rho_1,
		hydraulic["v_col"]
	);

	hydraulic["deltaP_col"] = computeDeltaP_col	// total collector pressure loss
	(
		hydraulic["deltaP_col_sup"],
		hydraulic["deltaP_col_int"]
	);

	// total pressure loss

	hydraulic["deltaP"] = computeDeltaP			// total pressure loss
	(
		hydraulic["deltaP_f"],
		hydraulic["deltaP_loc"],
		hydraulic["deltaP_col"]
	);


// * * * * * * * * * * * * * * * Structural  * * * * * * * * * * * * * * * * //

	std::map<std::string, double> structural;

	// coil tube

	structural["C_err_t"] = computeC_err	// coil tube errosion thickness
	(
		d,
		tol_t
	);
	structural["e_min_t"] = computeE_min	// coil tube analysis thickness
	(
		delta_t,
		tol_t
	);
	structural["e_a_t"] = computeE_a		// coil tube analysis thickness
	(
		structural["C_err_t"],
		delta_t
	);

	structural["e_t"] = computeE_t			// coil tube req. thickness
	(
		d_in,
		p_1,
		f_d
	);
	structural["p_max_t"] = computeP_max	// coil tube max pressure
	(
		structural["C_err_t"],
		d,
		d_in,
		delta_t,
		f_d
	);

	// collector tube

	structural["C_err_col"] = computeC_err	// collector errosion thickness
	(
		d_col,
		tol_col
	);
	structural["e_min_col"] = computeE_min	// collector analysis thickness
	(
		delta_col,
		tol_col
	);
	structural["e_a_col"] = computeE_a		// collector analysis thickness
	(
		structural["C_err_col"],
		delta_col
	);

	structural["e_col"] = computeE_t		// collector req. thickness
	(
		d_col_in,
		p_1,
		f_d
	);

	structural["r_is_col"] = computeR_is_col	// collector inside curvature radius
	(
		d_col,
		structural["e_a_col"]
	);

	structural["L_b_col_req"] = computeL_b_col_req	// collector isolated hole length
	(
		d,
		structural["e_a_col"],
		structural["r_is_col"]
	);
	structural["w_min_col"] = computeW_min_col	// collector min discontinuity distance
	(
		structural["e_a_col"],
		structural["r_is_col"]
	);

	// continuing isolated hole analysis because L_b_col_req < s_3
	// because the collector and the coil tube are SPH235

	structural["p_max_col"] = computeP_max_col	// collector max pressure
	(
		d,
		d_in,
		d_col,
		delta_t,
		delta_w_col_t,
		structural["e_a_col"],
		l_prot_t,
		structural["r_is_col"],
		f_d,
		f_d
	);

	// collector dish

	structural["C_err_dish"] = computeC_err	// dish corrosion thickness
	(
		d_dish,
		tol_dish
	);
	structural["e_min_dish"] = computeE_min	// coil tube analysis thickness
	(
		delta_dish,
		tol_dish
	);
	structural["e_a_dish"] = computeE_a		// coil tube analysis thickness
	(
		structural["C_err_dish"],
		delta_dish
	);

	structural["e_dish"] = computeE_dish	// dish required thickness
	(
		d_dish_in,
		p_1,
		R_dish,
		R_p02t,
		r_dish,
		f_d
	);
	structural["p_max_dish"] = computeP_max_dish	// dish max pressure
	(
		structural["C_err_dish"],
		d_dish_in,
		delta_dish,
		R_dish,
		R_p02t,
		r_dish,
		f_d
	);

	// test values

	structural["p_test"] = computeP_test	// test pressure
	(
		f_a,
		f_d,
		p_1,
		std::min
		(
			std::min(structural["p_max_t"], structural["p_max_col"]),
			structural["p_max_dish"]
		)
	);

	structural["e_t_test"] = computeE_t		// coil tube req. thickness
	(
		d_in,
		structural["p_test"],
		f_test
	);
	structural["e_dish_test"] = computeE_dish	// dish required thickness
	(
		d_dish_in,
		structural["p_test"],
		R_dish,
		R_p02,
		r_dish,
		f_test
	);
	structural["e_col_test"] = computeE_t	// coil tube req. thickness
	(
		d_col_in,
		structural["p_test"],
		f_test
	);

	// No flange calculation necessary, using: EN 1092-1/01/DN 50/PN 63/SPH 235

	// Contact achieved via plastic deformation of fin
/*
	// coil-fin tolerance fit
	structural["Q_i"] = computeQ_i
	(
		d,
		d_in
	);
	structural["Q_o"] = computeQ_o
	(
		d,
		d_0
	);
	structural["K"] = computeK
	(
		E_t,
		E_f,
		structural["Q_i"],
		structural["Q_o"],
		nu_t,
		nu_f
	);
	structural["Z_all"] = computeZ_all
	(
		E_f,
		structural["K"],
		structural["Q_i"],
		structural["Q_o"],
		R_p02,
		R_p02_f,
		1.1			// safety factor
	);

	structural["EI_min"] = computeEI_min
	(
		d,
		tol_d_t,
		Rz_t,
		Rz_f,
		structural["Z_all"]
	);
	structural["ES_max"] = computeES_max
	(
		beta_t,
		beta_f,
		d,
		E_f,
		-tol_d_t,
		structural["K"],
		p_contact,
		Rz_t,
		Rz_f,
		1.1,				// temperature safety factor
		structural["Z_all"]
	);
*/


	// Output

	std::cout << "\nDimensioning\n\n";
	printMap(geometric);

	std::cout << "\nThermal analysis\n\n";
	printMap(thermal);

	std::cout << "\nAerodynamic analysis\n\n";
	printMap(aerodynamic);

	std::cout << "\nHydraulic analysis\n\n";
	printMap(hydraulic);

	std::cout << "\nStructural analysis\n\n";
	printMap(structural);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
