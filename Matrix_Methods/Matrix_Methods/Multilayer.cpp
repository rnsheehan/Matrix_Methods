#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions of the methods in the fresnel class

fresnel::fresnel()
{
	// Default constructor
	params_defined = false; 

	n1 = n2 = nrat = nrat_sqr = 0.0;

	theta_in = theta_t = theta_critical = theta_brewster = zero;
}

fresnel::fresnel(double n_left, double n_right)
{
	// Primary constructor
	set_params(n_left, n_right);
}

fresnel::~fresnel()
{
	T.clear(); params_defined = false; 
}

void fresnel::set_params(double n_left, double n_right)
{
	// assign values to the interface RI

	try {
		bool c1 = n_left > 0.0 ? true : false; 
		bool c2 = n_right > 0.0 ? true : false;
		bool c10 = c1 && c2; 

		if (c10) {
			n1 = n_left;
			n2 = n_right;		 
			
			nrat = n2 / n1; // nrat > 1 => external reflection, nrat < 1 => internal reflection

			nrat_sqr = template_funcs::DSQR(nrat);			
			
			theta_critical = nrat < 1.0 ? asin(nrat) : PI_2;
			
			theta_brewster = atan(nrat);

			params_defined = true; 
		}
		else {
			std::string reason;
			reason = "Error: void fresnel::set_params(double n_left, double n_right)\n";
			if (!c1 || !c2) reason += "RI values not correctly defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fresnel::set_angles(std::complex<double> angle)
{
	// compute the transmission angle based on the assigned RI values

	try {
		bool c3 = abs(angle) >= 0.0 ? true : false;

		if (params_defined && c3) {

			if (abs(angle) == 0.0) {
				theta_in = angle;
				
				theta_t = zero; 

				cos_theta_1 = cos_theta_2 = 1.0; 
			}
			else {
				theta_in = angle;

				theta_t = asin((n1 / n2)*sin(theta_in)); // transmission angle 

				cos_theta_1 = cos(theta_in); // cosine of input angle 

				cos_theta_2 = cos(theta_t); // cosine of transmission angle
			}			
		}
		else {
			std::string reason;
			reason = "Error: void fresnel::set_angles(std::complex<double> angle)\n";
			if (!c3) reason += "Input angle not correctly defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::complex<double> fresnel::reflection(bool polarisation, std::complex<double> angle)
{
	// compute the reflection coefficient for a dielectric interface

	try {
		if (params_defined) {

			set_angles(angle); 

			std::complex<double> t1, t2, numer, denom; 
			
			if (polarisation) { // TE reflection coefficient
				t1 = n1 * cos_theta_1; 
				t2 = n2 * cos_theta_2;
			}
			else { // TM reflection coefficient
				t1 = n2 * cos_theta_1;
				t2 = n1 * cos_theta_2;
			}

			numer = t1 - t2;

			denom = t1 + t2;

			if (abs(denom) > 0.0) {
				return numer / denom; 
			}
			else {
				return zero; 
			}
		}
		else {
			std::string reason;
			reason = "Error: std::complex<double> fresnel::reflection(bool polarisation, std::complex<double> angle)\n";
			reason += "interface parameters not defined\n"; 
			throw std::invalid_argument(reason);
		}			
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::complex<double> fresnel::transmission(bool polarisation, std::complex<double> angle)
{
	// compute the reflection coefficient for a dielectric interface

	try {
		if (params_defined) {

			set_angles(angle);

			std::complex<double> t1, numer, denom;

			t1 = n1 * cos_theta_1; 

			numer = 2.0 * t1; 

			if (polarisation) { // TE reflection coefficient
				denom = t1 + n2 * cos_theta_2;
			}
			else { // TM reflection coefficient
				denom = n2 * cos_theta_1 + n1 * cos_theta_2;
			}

			if (abs(denom) > 0.0) {
				return numer / denom;
			}
			else {
				return zero;
			}
		}
		else {
			std::string reason;
			reason = "Error: std::complex<double> fresnel::transmission(bool polarisation, std::complex<double> angle)\n";
			reason += "interface parameters not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fresnel::compute_T(bool polarisation, double n_left, double n_right, std::complex<double> angle)
{
	try {
		set_params(n_left, n_right); 

		std::complex<double> r12, t12, z1, z2; 

		r12 = reflection(polarisation, angle); 

		t12 = transmission(polarisation, angle);

		z1 = 1.0 / t12; z2 = r12 / t12; 

		int size = 2;
		T = vecut::zero_cmat(size, size);

		T[0][0] = z1; T[0][1] = z2; 

		T[1][0] = z2; T[1][1] = z1; 

	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<std::complex<double>>> fresnel::transition_matrix()
{
	// return the transfer matrix of a given layer

	if (params_defined) {
		return T;
	}
	else {
		int size = 2;
		return vecut::zero_cmat(size, size);
	}
}

// definitions of the methods used in the propagation class

propagation::propagation()
{
	// Default constructor
	params_defined = false; 

	lambda = n = d = 0.0; 

	phase = theta = zero;
}

propagation::propagation(double wavelength, double index, double thickness, std::complex<double> angle)
{
	set_params(wavelength, index, thickness, angle);
}

propagation::~propagation()
{
	P.clear(); params_defined = false; 
}

void propagation::set_params(double wavelength, double index, double thickness, std::complex<double> angle)
{
	try {
		bool c1 = wavelength > 0.0 ? true : false; 
		bool c2 = index > 0.0 ? true : false;
		bool c3 = thickness > 0.0 ? true : false;
		bool c4 = abs(angle) >= 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4; 

		if (c10) {
			lambda = wavelength; 
			n = index; 
			d = thickness; 
			theta = angle; 
			phase = eye * (Two_PI / wavelength) * n * d * cos(theta); 

			params_defined = true; 
		}
		else {
			std::string reason;
			reason = "Error: void propagation::set_params(double wavelength, double index, double thickness, std::complex<double> angle)\n";
			if (!c10) reason += "Input values not correctly defined\n";

			throw std::invalid_argument(reason);
		}

	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void propagation::compute_P(double wavelength, double index, double thickness, std::complex<double> angle)
{
	try {
		set_params(wavelength, index, thickness, angle);

		int size = 2; 

		P = vecut::zero_cmat(size, size); 

		P[0][0] = exp(phase); P[1][1] = exp(-1.0*phase);
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<std::complex<double>>> propagation::propagation_matrix()
{
	// return the transfer matrix of a given layer

	if (params_defined) {
		return P;
	}
	else {
		int size = 2;
		return vecut::zero_cmat(size, size);
	}
}

