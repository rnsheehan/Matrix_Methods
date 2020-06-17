#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the compute function in the spectrum namespace
// R. Sheehan 23 - 7 - 2019

void spectrum::compute(double &n_clad, double &n_sub, std::vector<std::vector<std::complex<double>>> &M, std::complex<double> &r, std::complex<double> &t)
{
	// Given a transfer matrix M(\lambda) computed for some structure at some specific wavelength
	// Compute the reflectance and transmittance coefficients when the structure described by M is placed on substrate of
	// RI n_sub and covered by material of RI n_clad
	// Computed reflectance and transmittance values will be returned through pass by reference to r and t
	// R. Sheehan 23 - 7 - 2019

	try {
		bool c1 = n_clad > 0 ? true : false; // equivalent to n_{0} in formulae
		bool c2 = n_sub > 0 ? true : false; // equivalent to n_{t} in formulae
		bool c3 = (!M.empty()) ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::complex<double> r_numer, t_numer, denom, z1, z2;

			z1 = (M[0][0] + M[0][1] * n_sub)*n_clad;

			z2 = (M[1][0] + M[1][1] * n_sub);

			denom = z1 + z2; // common denominator

			if (abs(denom) > 0.0) { // avoid division by zero
				r_numer = z1 - z2; // reflectance numerator

				t_numer = 2.0*n_clad; // transmittance numerator

				r = r_numer / denom; // reflectance

				t = t_numer / denom; // transmittance
			}
			else {
				r = t = zero;
			}
		}
		else {
			std::string reason = "Error: void spectrum::compute(double &n_clad, double &n_sub, std::vector<std::vector<std::complex<double>>> &M, std::complex<double> &r, std::complex<double> &t)\n";
			if (!c1) reason += "n_clad: " + template_funcs::toString(n_clad, 2) + " is not valid\n";
			if (!c1) reason += "n_sub: " + template_funcs::toString(n_sub, 2) + " is not valid\n";
			if (!c3) reason += "Transfer Matrix M has not been computed\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definitions of the fowles_layer class
fowles_layer::fowles_layer()
{
	// Default constructor
	defined = false;
	theta_out = R = T = 0.0; 
	r = t = zero; 
}

fowles_layer::fowles_layer(bool pol, double theta_in, double thickness, double wavelength, double n0, double n1, double n2)
{
	// constructor
	set_params(pol, theta_in, thickness, wavelength, n0, n1, n2);
}

fowles_layer::fowles_layer(const fowles_layer &lobj)
{
	// copy constructor
	*this = lobj;
}

fowles_layer::~fowles_layer()
{
	// Deconstructor
	defined = false;
	theta_out = R = T = 0.0;
	r = t = zero;
	M.clear();
}

void fowles_layer::set_params(bool pol, double theta_in, double thickness, double wavelength, double n0, double n1, double n2, bool loud)
{
	// compute the fowles_layer parameters
	// pol is the polarisation of the incoming light
	// theta_in is the wavevector input angle
	// thickness of the fowles_layer in units of nm 
	// wavelength of light in the fowles_layer in units of nm
	// n0 wavelength dependent refractive index of the cladding
	// n1 wavelength dependent refractive index of the fowles_layer
	// n2 wavelength dependent refractive index of the substrate
	// R. Sheehan 15 - 7 - 2019
	// Updated R. Sheehan 17 - 6 - 2020

	try {
		bool c1 = thickness > 0.0 ? true : false;
		bool c2 = wavelength > 0.0 ? true : false;
		bool c3 = n0 > 0.0 ? true : false;
		bool c4 = n1 > 0.0 ? true : false;
		bool c5 = n2 > 0.0 ? true : false;
		bool c6 = theta_in >= 0.0 && theta_in < PI_2 ? true : false; 
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6;

		if (c10) {
			// Declare require parameters
			int size = 2;
			double sin_theta_in, cos_theta_in, phi, cos_phi, cos_theta_out, beta, c1, c2, p0, p1, p2;
			std::complex<double> z1, z2, z3, z4, z5, A(0.0, 0.0), B(0.0, 0.0), C(0.0, 0.0), D(0.0, 0.0);

			// Compute angles if necessary
			if (fabs(theta_in) < EPS) {
				sin_theta_in = phi = theta_out = 0.0; 
				cos_theta_in = cos_phi = cos_theta_out = 1.0; 
			}
			else {
				sin_theta_in = sin(theta_in); // sine of the incidence angle
				cos_theta_in = cos(theta_in); // sine of the incidence angle
				cos_theta_out = cos(theta_out); // store cos(theta_out) as it will be used multiple times
				phi = asin((n0 / n1) * sin_theta_in); // wavevector angle inside the layer
				theta_out = asin((n0 / n2) * sin_theta_in); // it will be useful to have this as an output parameter later	
				cos_phi = cos(phi); // store cos(phi) as it will be used multiple times
				cos_theta_out = cos(theta_out); // store cos(theta_out) as it will be used multiple times
			}

			// compute parameters of elements of transfer matrix
			beta = (Two_PI * (n1 / wavelength)) * thickness; 
			p0 = pol ? n0 : 1.0 / n0;
			p1 = pol ? n1 : 1.0 / n1;
			p2 = pol ? n2 : 1.0 / n2;

			// scale by cos(angles) if necessary
			if (cos_phi != 1.0) {
				beta *= cos_phi; // beta = k_{0} n l cos(phi), k_{0} = 2 pi / lambda
				p0 *= cos_theta_in; 
				p1 *= cos_phi; 
				p2 *= cos_theta_out; 
			}
			
			// compute elements of transfer matrix
			c1 = cos(beta); c2 = sin(beta);

			z1 = -1.0*eye / p1; z2 = -1.0*eye*p1;

			A.real(c1); B = z1 * c2;
			C = z2 * c2; D.real(c1);

			M = vecut::zero_cmat(size, size);

			M[0][0] = A; M[0][1] = B;
			M[1][0] = C; M[1][1] = D;

			// compute reflectivity and transmissivity
			z3 = (A + B * p2) * p0; 
			z4 = (C + D * p2); 
			z5 = z3 + z4; 
			
			r = (z3 - z4) / z5; // reflection coefficient
			t = 2.0 * p0 / z5; // transmission coefficient

			R = template_funcs::DSQR(abs(r)); // Reflectivity
			T = (p2/p0)* template_funcs::DSQR(abs(t)); // Transmissivity

			defined = true;

			if (loud) layer_stats(); 			
		}
		else {
			std::string reason = "Error: void fowles_layer::set_params(double thickness, double wavelength, double ref_index)\n";
			if (!c1) reason += "thickness: " + template_funcs::toString(thickness, 2) + " is not valid\n";
			if (!c2) reason += "wavelength: " + template_funcs::toString(wavelength, 2) + " is not valid\n";
			if (!c3) reason += "ref_index: " + template_funcs::toString(n1, 2) + " is not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fowles_layer::layer_stats()
{
	// print some of the computed parameters to the console
	// R. Sheehan 17 - 6 - 2020

	if (defined) {
		std::cout << "Layer Statistics\n\n";
		std::cout << "Output angle: " << theta_out * RAD_TO_DEG << "\n\n"; 
		std::cout << "r: " << r << "\n"; 
		std::cout << "t: " << t << "\n\n"; 
		std::cout << "R: " << R << "\n"; 
		std::cout << "T: " << T << "\n"; 
		std::cout << "Conservation of Energy R + T: " << R+T << "\n\n"; 
		std::cout << "Transfer Matrix\n";
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++)
				std::cout << M[i][j] << " ";
			std::cout << "\n"; 
		}
		std::cout << "|M|: " << M[0][0] * M[1][1] - M[0][1] * M[1][0] << "\n\n"; 
	}

}

std::vector<std::vector<std::complex<double>>> fowles_layer::transfer_matrix()
{
	// return the transfer matrix of a given fowles_layer

	if (defined) {
		return M;
	}
	else {
		int size = 2;
		return vecut::zero_cmat(size, size);
	}
}