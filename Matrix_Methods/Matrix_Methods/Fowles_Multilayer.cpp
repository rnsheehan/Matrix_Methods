#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the compute function in the spectrum namespace
// R. Sheehan 23 - 7 - 2019

void spectrum::compute_r_t(double& p_in, double& p_out, std::vector<std::vector<std::complex<double>>>& M, std::complex<double>& r, std::complex<double>& t, double& R, double& T)
{
	// Given a transfer matrix M(\lambda) computed for some structure at some specific wavelength
	// Compute the reflectance and transmittance coefficients when the structure described by M is placed on substrate of
	// RI n_sub and covered by material of RI n_clad, this info is contained in the values p_in, p_out
	// which also contains the direction and polarisation information, see Born and Wolf for details
	// Computed reflectance and transmittance values will be returned through pass by reference to r (R) and t (T)
	// R. Sheehan 23 - 7 - 2019
	// Updated R. Sheehan 18 - 6 - 2020

	try {
		bool c1 = p_in > 0 ? true : false; 
		bool c2 = p_out > 0 ? true : false; 
		bool c3 = (!M.empty()) ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::complex<double> denom, z1, z2;

			z1 = (M[0][0] + M[0][1] * p_out )*p_in;

			z2 = (M[1][0] + M[1][1] * p_out);

			denom = z1 + z2; // common denominator

			if (abs(denom) > 0.0) { // avoid division by zero
				r = (z1 - z2) / denom; // reflectance

				t = (2.0 * p_in) / denom; // transmittance

				R = template_funcs::DSQR(abs(r)); // Reflectivity

				T = (p_out / p_in) * template_funcs::DSQR(abs(t)); // Transmissivity
			}
			else {
				r = t = zero;

				R = T = 0.0; 
			}
		}
		else {
			std::string reason = "Error: void spectrum::compute(double& p_in, double& p_out, std::vector<std::vector<std::complex<double>>>& M, std::complex<double>& r, std::complex<double>& t, double& R, double& T)\n";
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
	theta_inc = theta_layer = theta_out = R = T = p1 = p2 = p0 = 0.0;
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
	theta_inc = theta_layer = theta_out = R = T = p0 = p1 = p2 = 0.0;
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
			double beta, cb, sb;
			std::complex<double> z3, z4, z5;

			// compute parameters of elements of transfer matrix
			theta_inc = theta_in; 

			beta = (Two_PI * (n1 / wavelength)) * thickness; 

			p0 = pol ? n0 : 1.0 / n0;
			p1 = pol ? n1 : 1.0 / n1;
			p2 = pol ? n2 : 1.0 / n2;

			// scale by cos(angles) if necessary
			if ( fabs(theta_in) > EPS ) {
				theta_layer = asin( (n0 / n1) * sin(theta_inc) ); // wavevector angle inside the layer
				theta_out = asin((n0 / n2) * sin(theta_inc)); // it will be useful to have this as an output parameter later	
				beta *= cos(theta_layer); // beta = k_{0} n l cos(theta_layer), k_{0} = 2 pi / lambda
				p0 *= cos(theta_inc);
				p1 *= cos(theta_layer);
				p2 *= cos(theta_out);
			}
			else {
				theta_layer = theta_out = 0.0;
			}
			
			// compute elements of transfer matrix
			M = vecut::zero_cmat(size, size);

			cb = cos(beta); sb = sin(beta);			

			M[0][0].real(cb);					M[0][1] = (-1.0 * eye / p1) * sb;
			M[1][0] = (-1.0 * eye * p1) * sb;	M[1][1].real(cb);

			// compute reflectivity and transmissivity
			z3 = (M[0][0] + M[0][1] * p2) * p0; 
			z4 = (M[1][0] + M[1][1] * p2); 
			z5 = z3 + z4; 
			
			if (abs(z5) > 0.0) {
				r = (z3 - z4) / z5; // reflection coefficient
				t = 2.0 * p0 / z5; // transmission coefficient

				R = template_funcs::DSQR(abs(r)); // Reflectivity
				T = (p2 / p0) * template_funcs::DSQR(abs(t)); // Transmissivity

				defined = true;
			}
			else {
				r = t = zero; 
				R = T = 0.0; 
				defined = false; 
			}		

			if (loud) layer_stats(); 			
		}
		else {
			std::string reason = "Error: void fowles_layer::set_params(bool pol, double theta_in, double thickness, double wavelength, double n0, double n1, double n2, bool loud)\n";
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
		std::cout << "Incidence angle: " << theta_inc * RAD_TO_DEG << "\n"; 
		std::cout << "Layer angle: " << theta_layer * RAD_TO_DEG << "\n"; 
		std::cout << "Output angle: " << theta_out * RAD_TO_DEG << "\n\n"; 
		std::cout << "r: " << r << "\n"; 
		std::cout << "t: " << t << "\n\n"; 
		std::cout << "R: " << R << "\n"; 
		std::cout << "T: " << T << "\n"; 
		std::cout << "Conservation of Energy R + T: " << R + T << "\n\n"; 
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