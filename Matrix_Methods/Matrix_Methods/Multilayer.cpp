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

// definitions for the multilayer class

multilayer::multilayer()
{
	layer_mat = nullptr; substrate = nullptr; cladding = nullptr;
}

multilayer::multilayer(sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate)
{
	set_params(swp_obj, the_layer, the_cladding, the_substrate);
}

multilayer::~multilayer()
{
	M.clear(); r.clear(); t.clear();
}

void multilayer::set_params(sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate)
{
	// assign the parameters that make up the layer structure

	try {

		bool c1 = swp_obj.defined();
		bool c4 = the_layer != nullptr ? true : false;
		bool c5 = the_cladding != nullptr ? true : false;
		bool c6 = the_substrate != nullptr ? true : false;
		bool c10 = c1 && c4 && c5 && c6;

		if (c10) {
			wavelength.set_vals(swp_obj); // assign the values for the wavelength parameter space

										  // assign the material objects
			layer_mat = the_layer;
			substrate = the_substrate;
			cladding = the_cladding;

			int size = 2;
			M = vecut::zero_cmat(size, size);
		}
		else {
			std::string reason;
			reason = "Error: void multilayer::set_params(sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate)\n";
			if (!c1) reason += "swp_obj is not correct\n";
			if (!c4) reason += "the_layer has not been correctly assigned";
			if (!c5) reason += "the_cladding has not been correctly assigned";
			if (!c6) reason += "the_substrate has not been correctly assigned";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void multilayer::compute_r_t_Fowles(int n_layers, double layer_thickness, bool loud)
{
	// Compute the transfer matrix for a system with n_layers each having the same thickness
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = n_layers > 0 ? true : false;
		bool c2 = layer_thickness > 0 ? true : false;
		bool c10 = c1 && c2 && wavelength.defined();

		if (c10) {
			// sweep over all wavelengths
			double lambda, RI_layer, RI_clad, RI_sub;
			std::complex<double> reflectance, transmittance;
			std::string filename;
			filename = "Fowles_multilayer_Data_" + template_funcs::toString(n_layers) + "_Layers_" + template_funcs::toString(layer_thickness, 2) + "_Layer_Thickness.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results " + template_funcs::toString(n_layers) + " layers,  layer thickness: " + template_funcs::toString(layer_thickness, 2) + " um\n";
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {

				lambda = wavelength.get_val(i); // assign wavelength value

				layer_mat->set_wavelength(lambda); // assign wavelength to layer material object

				cladding->set_wavelength(lambda); // assign wavelength to cladding material object

				substrate->set_wavelength(lambda); // assign wavelength to substrate material object

				RI_layer = layer_mat->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				// since this is an anti-reflection film the same matrix is used for all layers for fixed wavelength
				layer this_layer(layer_thickness, lambda, RI_layer); // compute transfer matrix for the layer at this wavelength

				// get transfer matrix for all layers
				// for the anti-reflection filter multiply the single layer TM by itself n_layers times
				for (int j = 0; j < n_layers; j++) {

					if (j == 0) {
						M = this_layer.transfer_matrix();
					}
					else {
						std::vector<std::vector<std::complex<double>>> Mnow = this_layer.transfer_matrix();
						M = vecut::cmat_cmat_product(M, Mnow);
					}
				}

				// compute reflectance and transmittance of the structure at this wavelength

				RI_clad = cladding->refractive_index(); // equivalent to n_{0} in formulae

				RI_sub = substrate->refractive_index(); // equivalent to n_{t} in formulae

				spectrum::compute(RI_clad, RI_sub, M, reflectance, transmittance); // compute the wavelength dependent r, t values from the wavelength dependent M(\lambda)

				r.push_back(reflectance); // store reflectance

				t.push_back(transmittance); // store transmittance

				if (loud) std::cout << lambda << " , " << RI_clad << " , " << RI_layer << " , " << RI_sub << " , " << template_funcs::DSQR(abs(reflectance)) << " , " << template_funcs::DSQR(abs(transmittance)) << "\n";

				write << lambda << " , " << RI_clad << " , " << RI_layer << " , " << RI_sub << " , " << template_funcs::DSQR(abs(reflectance)) << " , " << template_funcs::DSQR(abs(transmittance)) << "\n";
			}

			write.close();
		}
		else {
			std::string reason;
			reason = "Error: void multilayer::compute_r_t_Fowles(int n_layers, double layer_thickness, bool loud)\n";
			if (!c1) reason += "n_layers: " + template_funcs::toString(n_layers) + " is not valid\n";
			if (!c2) reason += "layer_thickness: " + template_funcs::toString(layer_thickness, 2) + " is not valid\n";
			if (!wavelength.defined()) reason += "wavelength sweep object not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void multilayer::compute_r_t(int n_layers, double layer_thickness, bool loud)
{
	// Compute the transfer matrix for a system with n_layers each having the same thickness
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = n_layers > 0 ? true : false;
		bool c2 = layer_thickness > 0 ? true : false;
		bool c10 = c1 && c2 && wavelength.defined();

		if (c10) {
			// sweep over all wavelengths
			double lambda, RI_layer, RI_clad, RI_sub;
			std::complex<double> r_coeff, t_coeff;
			std::string filename;
			filename = "multilayer_Data_" + template_funcs::toString(n_layers) + "_Layers_" + template_funcs::toString(layer_thickness, 2) + "_Layer_Thickness.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results " + template_funcs::toString(n_layers) + " layers,  layer thickness: " + template_funcs::toString(layer_thickness, 2) + " um\n";
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {

				lambda = wavelength.get_val(i); // assign wavelength value

				layer_mat->set_wavelength(lambda); // assign wavelength to layer material object

				cladding->set_wavelength(lambda); // assign wavelength to cladding material object

				substrate->set_wavelength(lambda); // assign wavelength to substrate material object

				// compute wavelength dependent refractive index at this wavelength

				RI_clad = cladding->refractive_index(); 

				RI_layer = layer_mat->refractive_index(); 

				RI_sub = substrate->refractive_index(); 

				std::vector<std::vector<std::complex<double>>> Mnow; 

				// compute the transition matrix at the cladding - layer interface
				fresnel T12; 

				T12.compute_T(TE, RI_clad, RI_layer, zero); 

				M = T12.transition_matrix(); // M = T_{12}

				// compute the propagation matrix for the layer
				propagation P2; 

				P2.compute_P(lambda, RI_layer, n_layers * layer_thickness, zero); 

				Mnow = P2.propagation_matrix(); 

				M = vecut::cmat_cmat_product(M, Mnow); // M = T_{12} * P_{2}

				// compute the transition matrix at the layer-substrate interface
				fresnel T23;

				T23.compute_T(TE, RI_layer, RI_sub, zero);

				Mnow = T23.transition_matrix(); 

				M = vecut::cmat_cmat_product(M, Mnow); // M = T_{12} * P_{2} * T_{23}

				r_coeff = M[0][1] / M[0][0];

				t_coeff = 1.0 / M[0][0]; 

				if (loud) std::cout << lambda << " , " << RI_clad << " , " << RI_layer << " , " << RI_sub << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";

				write << lambda << " , " << RI_clad << " , " << RI_layer << " , " << RI_sub << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";
			}

			write.close();
		}
		else {
			std::string reason;
			reason = "Error: void multilayer::compute_r_t(int n_layers, double layer_thickness, bool loud)\n";
			if (!c1) reason += "n_layers: " + template_funcs::toString(n_layers) + " is not valid\n";
			if (!c2) reason += "layer_thickness: " + template_funcs::toString(layer_thickness, 2) + " is not valid\n";
			if (!wavelength.defined()) reason += "wavelength sweep object not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// definitions for the HL_stack class

HL_stack::HL_stack()
{
	layer_high = nullptr; layer_low = nullptr; substrate = nullptr; cladding = nullptr;
}

HL_stack::HL_stack(sweep &swp_obj, material *h_layer, material *l_layer, material *the_cladding, material *the_substrate)
{
	set_params(swp_obj, h_layer, l_layer, the_cladding, the_substrate);
}

HL_stack::~HL_stack()
{
	M.clear(); r.clear(); t.clear();
}

void HL_stack::set_params(sweep &swp_obj, material *h_layer, material *l_layer, material *the_cladding, material *the_substrate)
{
	// assign the parameters that make up the layer structure

	try {

		bool c1 = swp_obj.defined();
		bool c2 = h_layer != nullptr ? true : false;
		bool c4 = l_layer != nullptr ? true : false;
		bool c5 = the_cladding != nullptr ? true : false;
		bool c6 = the_substrate != nullptr ? true : false;
		bool c10 = c1 && c2 && c4 && c5 && c6;

		if (c10) {
			wavelength.set_vals(swp_obj); // assign the values for the wavelength parameter space

										  // assign the material objects
			layer_high = h_layer;
			layer_low = l_layer;
			substrate = the_substrate;
			cladding = the_cladding;

			int size = 2;
			M = vecut::zero_cmat(size, size);
		}
		else {
			std::string reason;
			reason = "Error: void multilayer::set_params(sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate)\n";
			if (!c1) reason += "swp_obj is not correct\n";
			if (!c4) reason += "the_layer has not been correctly assigned";
			if (!c5) reason += "the_cladding has not been correctly assigned";
			if (!c6) reason += "the_substrate has not been correctly assigned";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void HL_stack::compute_r_t_Fowles(int n_layers, double layer_thickness, bool loud)
{
	// Compute the transfer matrix for a system with n_layers each having the same thickness
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = n_layers > 0 ? true : false;
		bool c2 = layer_thickness > 0 ? true : false;
		bool c10 = c1 && c2 && wavelength.defined();

		if (c10) {
			// sweep over all wavelengths
			int total_layers = 2 * n_layers;
			double lambda, RI_h, RI_l, RI_clad, RI_sub;
			std::complex<double> reflectance, transmittance;
			std::string filename;
			filename = "Fowles_HL_stack_Data_" + template_funcs::toString(n_layers) + "_Layers_" + template_funcs::toString(layer_thickness, 2) + "_Layer_Thickness.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results " + template_funcs::toString(n_layers) + " layers,  layer thickness: " + template_funcs::toString(layer_thickness, 2) + " um\n";
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {

				lambda = wavelength.get_val(i); // assign wavelength value

				layer_high->set_wavelength(lambda); // assign wavelength to layer material object

				layer_low->set_wavelength(lambda); // assign wavelength to layer material object

				cladding->set_wavelength(lambda); // assign wavelength to cladding material object

				substrate->set_wavelength(lambda); // assign wavelength to substrate material object

				RI_h = layer_high->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				RI_l = layer_low->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				// since this is an anti-reflection film the same matrix is used for all layers for fixed wavelength
				layer H_layer(layer_thickness, lambda, RI_h); // compute transfer matrix for the layer at this wavelength

				layer L_layer(layer_thickness, lambda, RI_l); // compute transfer matrix for the layer at this wavelength

				std::vector<std::vector<std::complex<double>>> Mnow;

				// get transfer matrix for all layers
				// for the anti-reflection filter multiply the single layer TM by itself n_layers times
				for (int j = 0; j < total_layers; j++) {

					if (j % 2 == 0) {
						Mnow = L_layer.transfer_matrix();
					}
					else {
						Mnow = H_layer.transfer_matrix();
					}

					if (j == 0) {
						M = Mnow;
					}
					else {
						M = vecut::cmat_cmat_product(M, Mnow);
					}
				}

				// compute reflectance and transmittance of the structure at this wavelength

				RI_clad = cladding->refractive_index(); // equivalent to n_{0} in formulae

				RI_sub = substrate->refractive_index(); // equivalent to n_{t} in formulae

				spectrum::compute(RI_clad, RI_sub, M, reflectance, transmittance); // compute the wavelength dependent r, t values from the wavelength dependent M(\lambda)

				r.push_back(reflectance); // store reflectance

				t.push_back(transmittance); // store transmittance

				if (loud) std::cout << lambda << " , " << RI_clad << " , " << RI_h << " , " << RI_l << " , " << RI_sub << " , " << template_funcs::DSQR(abs(reflectance)) << " , " << template_funcs::DSQR(abs(transmittance)) << "\n";

				write << lambda << " , " << RI_clad << " , " << RI_h << " , " << RI_l << " , " << RI_sub << " , " << template_funcs::DSQR(abs(reflectance)) << " , " << template_funcs::DSQR(abs(transmittance)) << "\n";
			}

			write.close();
		}
		else {
			std::string reason;
			reason = "Error: void HL_stack::compute_r_t_Fowles(int n_layers, double layer_thickness)\n";
			if (!c1) reason += "n_layers: " + template_funcs::toString(n_layers) + " is not valid\n";
			if (!c2) reason += "layer_thickness: " + template_funcs::toString(layer_thickness, 2) + " is not valid\n";
			if (!wavelength.defined()) reason += "wavelength sweep object not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

