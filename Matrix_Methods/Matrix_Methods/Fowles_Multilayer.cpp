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
	defined = false;
	size = 0;
}

fowles_layer::fowles_layer(double thickness, double wavelength, double ref_index)
{
	set_params(thickness, wavelength, ref_index);
}

fowles_layer::fowles_layer(const fowles_layer &lobj)
{
	*this = lobj;
}

fowles_layer::~fowles_layer()
{
	defined = false;
	M.clear();
}

void fowles_layer::set_params(double thickness, double wavelength, double ref_index)
{
	// assign the fowles_layer parameters
	// thickness of the fowles_layer in units of nm
	// wavelength of light in the fowles_layer
	// wavelength dependent refractive index of the fowles_layer
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = thickness > 0.0 ? true : false;
		bool c2 = wavelength > 0.0 ? true : false;
		bool c3 = ref_index > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			size = 2;

			double k = Two_PI * (ref_index / wavelength);
			double kl = k * thickness;
			double c1, c2;
			std::complex<double> z1, z2, A(0.0, 0.0), B(0.0, 0.0), C(0.0, 0.0), D(0.0, 0.0);

			c1 = cos(kl); c2 = sin(kl);

			z1 = -1.0*eye / ref_index; z2 = -1.0*eye*ref_index;

			A.real(c1); B = z1 * c2;
			C = z2 * c2; D.real(c1);

			M = vecut::zero_cmat(size, size);

			M[0][0] = A; M[0][1] = B;
			M[1][0] = C; M[1][1] = D;

			defined = true;
		}
		else {
			std::string reason = "Error: void fowles_layer::set_params(double thickness, double wavelength, double ref_index)\n";
			if (!c1) reason += "thickness: " + template_funcs::toString(thickness, 2) + " is not valid\n";
			if (!c2) reason += "wavelength: " + template_funcs::toString(wavelength, 2) + " is not valid\n";
			if (!c3) reason += "ref_index: " + template_funcs::toString(ref_index, 2) + " is not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<std::complex<double>>> fowles_layer::transfer_matrix()
{
	// return the transfer matrix of a given fowles_layer

	if (defined) {
		return M;
	}
	else {
		size = 2;
		return vecut::zero_cmat(size, size);
	}
}