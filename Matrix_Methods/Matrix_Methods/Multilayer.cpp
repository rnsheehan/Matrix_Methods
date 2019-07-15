#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the classes in Multilayer.h

// Definitions of the layer class
layer::layer()
{
	defined = false;
	size = 0; 
}

layer::layer(double thickness, double wavelength, double ref_index)
{
	set_params(thickness, wavelength, ref_index);
}

layer::layer(layer &lobj)
{
	*this = lobj; 
}

layer::~layer()
{
	defined = false; 
	M.clear(); 
}

void layer::set_params(double thickness, double wavelength, double ref_index)
{
	// assign the layer parameters
	// thickness of the layer in units of nm
	// wavelength of light in the layer
	// wavelength dependent refractive index of the layer
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
			std::string reason = "Error: void layer::set_params(double thickness, double wavelength, double ref_index)\n";
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

std::vector<std::vector<std::complex<double>>> layer::transfer_matrix()
{
	// return the transfer matrix of a given layer

	if (defined) {
		return M; 
	}
	else {
		size = 2; 
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
	M.clear();
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

void multilayer::build_transfer_matrix(int n_layers, double layer_thickness, bool loud)
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
			std::complex<double> r_numer, t_numer, denom, z1, z2, reflectance, transmittance; 
			if(loud) std::cout << "Results " + template_funcs::toString(n_layers) + " layers, thickness: " + template_funcs::toString(layer_thickness, 2) + " um\n"; 
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {
				
				lambda = wavelength.get_val(i); // assign wavelength value
				
				layer_mat->set_wavelength(lambda); // assign wavelength to layer material object
				cladding->set_wavelength(lambda); // assign wavelength to cladding material object
				substrate->set_wavelength(lambda); // assign wavelength to substrate material object
				
				RI_layer = layer_mat->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				// get transfer matrix for all layers
				for (int j = 0; j < n_layers; j++) {
					layer this_layer(layer_thickness, lambda, RI_layer); // compute transfer matrix for the layer at this wavelength

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
				z1 = (M[0][0] + M[0][1] * RI_sub)*RI_clad; 
				z2 = (M[1][0] + M[1][1] * RI_sub); 
				denom = z1 + z2; // common denominator
				r_numer = z1 - z2; // reflectance numerator
				t_numer = 2.0*RI_clad; // transmittance numerator
				reflectance = r_numer / denom; // reflectance
				transmittance = t_numer / denom; // transmittance
				r.push_back(reflectance); // store reflectance
				t.push_back(transmittance); // store transmittance

				if (loud) std::cout << lambda << " , " << RI_clad << " , " << RI_layer << " , " << RI_sub << " , " << template_funcs::DSQR(abs(reflectance)) << " , " << template_funcs::DSQR(abs(transmittance)) << "\n";
			}
		}
		else {
			std::string reason;
			reason = "Error: void multilayer::build_transfer_matrix(int n_layers, double layer_thickness)\n";
			if (!c1) reason += "n_layers: " + template_funcs::toString(n_layers) + " is not valid\n";
			if (!c2) reason += "layer_thickness: " + template_funcs::toString(layer_thickness,2) + " is not valid\n";			
			if (!wavelength.defined()) reason += "wavelength sweep object not defined correctly\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}