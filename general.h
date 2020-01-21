/*---------------------------------------------------------------------------*\

	finned_heat_exchanger - Copyright (C) 2019 P. Milovic

-------------------------------------------------------------------------------
License
	See the LICENSE file for license information.

Description
	General functions used during the design of the finned heat exchanger.

\*---------------------------------------------------------------------------*/

#ifndef GENERAL_H
#define GENERAL_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <stdexcept>

// * * * * * * * * * * * * * * * * * General * * * * * * * * * * * * * * * * //

inline void printMap
(
	const std::map<std::string, double>& map
)
{
	auto totalWidth {30};

	for (const auto& [key, value] : map)
	{
		auto outputWidth {totalWidth - key.size() - 1};

		std::cout.precision(7);
		std::cout << key << ":"
				  << std::setw(outputWidth)
				  << std::fixed
				  << value << '\n';
	}
}

// check if a map has a key
inline bool hasKey
(
	const std::map<std::string, double>& map,
	const std::string& key
)
{
	auto search {map.find(key)};

	if (search == map.end())
		return false;

	return true;
}


// interpolate a value from a map
inline double interpolate
(
	const std::map<double, double>& map,
	const double x
)
{
	// range check
	if
	(
		x < map.begin()->first ||
		x > map.rbegin()->first
	)
		throw std::range_error("interpolate(): 'x' out of range");

	double y;		// interpolated value
	double X_old;	// lower bound independent variable
	double Y_old;	// lower bound value

	for (const auto& [X, Y] : map)
	{
		if (x > X)
		{
			X_old = X;
			Y_old = Y;
		}
		else
			y = Y_old
			  + (Y - Y_old)
			  * (x - X_old)
			  / (X - X_old);
	}

	return y;
}


// compute the Reynolds number
inline double computeRe
(
	const double l,				// characteristic length
	const double mu,			// dynamic viscosity
	const double rho,			// density
	const double v				// velocity
)
{
	double Re {rho * v * l / mu};

	return Re;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
