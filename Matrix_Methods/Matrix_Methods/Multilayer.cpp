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
	// compute the transmission coefficient for a dielectric interface

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
	// compute a dielectric interface transition matrix

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
	// return a computed dielectric interface transition matrix

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
	// return the transfer matrix of a given fowles_layer

	if (params_defined) {
		return P;
	}
	else {
		int size = 2;
		return vecut::zero_cmat(size, size);
	}
}

// definitions for the layer class

layer::layer()
{
	// Default constructor

	thickness = 0.0; 
	the_mat = nullptr;
}

layer::layer(double &d, material *mat)
{
	set_params(d, mat); 
}

layer::layer(const layer &lobj)
{
	// copy constructor
	*this = lobj; 
}

void layer::set_params(double &d, material *mat)
{
	// assign parameters to the layer class
	
	try {
		bool c1 = d >= 0.0 ? true : false; 
		bool c2 = mat != nullptr ? true : false; 
		bool c10 = c1 && c2; 

		if (c10) {
			thickness = d; 
			the_mat = mat; 
		}
		else {
			std::string reason;
			reason = "Error: void layer::set_params(double d, material *mat)\n";
			if (!c1) reason += "layer thickness is not positive\n";
			if (!c2) reason += "material object has not been correctly instantiated";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double layer::get_RI(double wavelength)
{
	// compute the layer wavelength dependent refractive index
	// wavelength input in units of um

	try {
		bool c1 = wavelength > 0.0 ? true : false; 
		bool c2 = the_mat != nullptr ? true : false;
		bool c10 = c1 && c2;
		if (c10) {

			the_mat->set_wavelength(wavelength); // assign the material object refractive index
			
			return the_mat->refractive_index(); // return the computed RI value
		}
		else {
			std::string reason;
			reason = "Error: double layer::get_RI(double wavelength)\n";
			if (!c1) reason += "wavelength: " + template_funcs::toString(wavelength, 2) + " is not positive\n";
			if (!c2) reason += "material object has not been correctly instantiated";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// definitions for the multilayer class
multilayer::multilayer()
{
	// Default Constructor
}

multilayer::multilayer(sweep &swp_obj, std::vector<layer> &layer_list)
{
	// Primary constructor
	set_params(swp_obj, layer_list);
}

multilayer::~multilayer()
{
	// Deconstructor
	r.clear(); t.clear(); the_layers.clear(); 
}

void multilayer::set_params(sweep &swp_obj, std::vector<layer> &layer_list)
{
	try {
		bool c1 = swp_obj.defined();
		bool c2 = !(layer_list.empty()) ? true : false;
		bool c10 = c1 && c2; 

		if (c10) {
			wavelength = swp_obj; 

			the_layers.assign(layer_list.begin(), layer_list.end()); 

			//the_layers = layer_list; 
		}
		else {
			std::string reason;
			reason = "Error: void multilayer::set_params(sweep &swp_obj, std::list<layer> &layer_list)\n";
			if (!c1) reason += "swp_obj is not correct\n";
			if (!c2) reason += "layer_list has not been correctly assigned";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<std::complex<double>>> multilayer::transmission_matrix(double &lambda, bool polarisation, bool loud)
{
	// Compute the structure transmission matrix at a given wavelength
	// Loop over all layers in the structure and return the computed transmission matrix for that structure
	// R. Sheehan 13 - 8 - 2019

	try {
		if ( lambda > 0.0 && !the_layers.empty() ) {
			int size = 2; 

			std::vector<std::vector<std::complex<double>>> M, Msub, T, P;  // array to store transfer matrix of overall structure

			M = vecut::idn_cmat(size); 

			Msub = vecut::zero_cmat(size, size); 

			fresnel transition; 

			propagation prop; 

			// step through the layer list
			for (size_t i = 0; i < the_layers.size() - 1; i++) {
				
				// Compute the transmission matrix for a given wavelength				
				transition.compute_T(polarisation, the_layers[i].get_RI(lambda), the_layers[i + 1].get_RI(lambda), zero);

				T = transition.transition_matrix(); 

				if (i < the_layers.size() - 2) {
					// Compute the propagation matrix for a given wavelength
					prop.compute_P(lambda, the_layers[i + 1].get_RI(lambda), the_layers[i + 1].get_d(), zero);
					P = prop.propagation_matrix();
				}
				
				if (i == the_layers.size() - 2) {
					if (loud) std::cout << i << ": M * T_{" << i << " , " << i + 1 << "}\n";
					M = vecut::cmat_cmat_product(M, T);
				}
				else {
					if (loud) std::cout << i << ": M * T_{" << i << " , " << i + 1 << "} * P_{" << i + 1 << "}\n"; 
					Msub = vecut::cmat_cmat_product(T, P);
					M = vecut::cmat_cmat_product(M, Msub);
				}
			}

			if (loud)std::cout << "\n"; 

			return M; 
		}
		else {
			std::string reason;
			reason = "Error: std::vector<std::vector<std::complex<double>>> multilayer::transmission_matrix(double &lambda)\n";
			if (lambda < 0.0) reason += "input wavelength value is not correct\n";
			if (the_layers.empty()) reason += "the_layers has not been correctly assigned";
			throw std::invalid_argument(reason);
		}	
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void multilayer::compute_spectrum(bool polarisation, bool loud)
{
	// Compute the reflection and transmission spectrum for a multilayer stack
	
	try {
		if ( wavelength.defined() && !the_layers.empty() ) {
		
			double lambda; 
			std::complex<double> r_coeff, t_coeff;
			std::vector<std::vector<std::complex<double>>> M; 

			std::string filename;
			filename = "Multilayer_Spectral_Response.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results\n";

			for (int i = 0; i < wavelength.get_Nsteps(); i++) {
				lambda = wavelength.get_val(i); 

				M = transmission_matrix(lambda, polarisation); // wavelength dependent transmission matrix

				r_coeff = M[0][1] / M[0][0]; // wavelength dependent reflection coefficient

				t_coeff = 1.0 / M[0][0]; // wavelength dependent transmission coefficient

				if (loud) std::cout << lambda << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";

				write << lambda << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";
				//write << lambda << " , " << template_funcs::DSQR(abs(r_coeff)) << " , " << template_funcs::DSQR(abs(t_coeff)) << "\n";
			}	

			write.close();
		}
		else {
			std::string reason;
			reason = "Error: void multilayer::compute_spectrum(bool loud)\n";
			if (!wavelength.defined()) reason += "swp_obj is not correct\n";
			if (the_layers.empty()) reason += "the_layers has not been correctly assigned";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// definitions for the multilayer_old class

multilayer_old::multilayer_old()
{
	layer_mat = nullptr; substrate = nullptr; cladding = nullptr;
}

multilayer_old::multilayer_old(sweep &swp_obj, material *the_fowles_layer, material *the_cladding, material *the_substrate)
{
	set_params(swp_obj, the_fowles_layer, the_cladding, the_substrate);
}

multilayer_old::~multilayer_old()
{
	M.clear(); r.clear(); t.clear();
}

void multilayer_old::set_params(sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate)
{
	// assign the parameters that make up the fowles_layer structure

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
			reason = "Error: void multilayer_old::set_params(sweep &swp_obj, material *the_fowles_layer, material *the_cladding, material *the_substrate)\n";
			if (!c1) reason += "swp_obj is not correct\n";
			if (!c4) reason += "the_fowles_layer has not been correctly assigned";
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

void multilayer_old::compute_r_t_Fowles(int n_fowles_layers, double layer_thickness, bool loud)
{
	// Compute the transfer matrix for a system with n_fowles_layers each having the same thickness
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = n_fowles_layers > 0 ? true : false;
		bool c2 = layer_thickness > 0 ? true : false;
		bool c10 = c1 && c2 && wavelength.defined();

		if (c10) {
			// sweep over all wavelengths
			double lambda, RI_fowles_layer, RI_clad, RI_sub;
			std::complex<double> reflectance, transmittance;
			std::string filename;
			filename = "Fowles_multifowles_layer_Data_" + template_funcs::toString(n_fowles_layers) + "_fowles_layers_" + template_funcs::toString(layer_thickness, 2) + "_fowles_layer_Thickness.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results " + template_funcs::toString(n_fowles_layers) + " fowles_layers,  layer thickness: " + template_funcs::toString(layer_thickness, 2) + " um\n";
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {

				lambda = wavelength.get_val(i); // assign wavelength value

				layer_mat->set_wavelength(lambda); // assign wavelength to fowles_layer material object

				cladding->set_wavelength(lambda); // assign wavelength to cladding material object

				substrate->set_wavelength(lambda); // assign wavelength to substrate material object

				RI_fowles_layer = layer_mat->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				// since this is an anti-reflection film the same matrix is used for all fowles_layers for fixed wavelength
				fowles_layer this_fowles_layer(layer_thickness, lambda, RI_fowles_layer); // compute transfer matrix for the fowles_layer at this wavelength

				// get transfer matrix for all fowles_layers
				// for the anti-reflection filter multiply the single fowles_layer TM by itself n_fowles_layers times
				for (int j = 0; j < n_fowles_layers; j++) {

					if (j == 0) {
						M = this_fowles_layer.transfer_matrix();
					}
					else {
						std::vector<std::vector<std::complex<double>>> Mnow = this_fowles_layer.transfer_matrix();
						M = vecut::cmat_cmat_product(M, Mnow);
					}
				}

				// compute reflectance and transmittance of the structure at this wavelength

				RI_clad = cladding->refractive_index(); // equivalent to n_{0} in formulae

				RI_sub = substrate->refractive_index(); // equivalent to n_{t} in formulae

				spectrum::compute(RI_clad, RI_sub, M, reflectance, transmittance); // compute the wavelength dependent r, t values from the wavelength dependent M(\lambda)

				r.push_back(reflectance); // store reflectance

				t.push_back(transmittance); // store transmittance

				if (loud) std::cout << lambda << " , " << RI_clad << " , " << RI_fowles_layer << " , " << RI_sub << " , " << template_funcs::DSQR(abs(reflectance)) << " , " << template_funcs::DSQR(abs(transmittance)) << "\n";

				write << lambda << " , " << RI_clad << " , " << RI_fowles_layer << " , " << RI_sub << " , " << template_funcs::DSQR(abs(reflectance)) << " , " << template_funcs::DSQR(abs(transmittance)) << "\n";
			}

			write.close();
		}
		else {
			std::string reason;
			reason = "Error: void multifowles_layer::compute_r_t_Fowles(int n_fowles_layers, double fowles_layer_thickness, bool loud)\n";
			if (!c1) reason += "n_fowles_layers: " + template_funcs::toString(n_fowles_layers) + " is not valid\n";
			if (!c2) reason += "fowles_layer_thickness: " + template_funcs::toString(layer_thickness, 2) + " is not valid\n";
			if (!wavelength.defined()) reason += "wavelength sweep object not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void multilayer_old::compute_r_t(int n_fowles_layers, double fowles_layer_thickness, bool loud)
{
	// Compute the transfer matrix for a system with n_fowles_layers each having the same thickness
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = n_fowles_layers > 0 ? true : false;
		bool c2 = fowles_layer_thickness > 0 ? true : false;
		bool c10 = c1 && c2 && wavelength.defined();

		if (c10) {
			// sweep over all wavelengths
			double lambda, RI_fowles_layer, RI_clad, RI_sub;
			std::complex<double> r_coeff, t_coeff;
			std::string filename;
			filename = "multifowles_layer_Data_" + template_funcs::toString(n_fowles_layers) + "_fowles_layers_" + template_funcs::toString(fowles_layer_thickness, 2) + "_fowles_layer_Thickness.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results " + template_funcs::toString(n_fowles_layers) + " fowles_layers,  fowles_layer thickness: " + template_funcs::toString(fowles_layer_thickness, 2) + " um\n";
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {

				lambda = wavelength.get_val(i); // assign wavelength value

				layer_mat->set_wavelength(lambda); // assign wavelength to fowles_layer material object

				cladding->set_wavelength(lambda); // assign wavelength to cladding material object

				substrate->set_wavelength(lambda); // assign wavelength to substrate material object

				// compute wavelength dependent refractive index at this wavelength

				RI_clad = cladding->refractive_index(); 

				RI_fowles_layer = layer_mat->refractive_index(); 

				RI_sub = substrate->refractive_index(); 

				std::vector<std::vector<std::complex<double>>> Mnow; 

				// compute the transition matrix at the cladding - fowles_layer interface
				fresnel T12; 

				T12.compute_T(TE, RI_clad, RI_fowles_layer, zero); 

				M = T12.transition_matrix(); // M = T_{12}

				// compute the propagation matrix for the fowles_layer
				propagation P2; 

				P2.compute_P(lambda, RI_fowles_layer, n_fowles_layers * fowles_layer_thickness, zero); 

				Mnow = P2.propagation_matrix(); 

				M = vecut::cmat_cmat_product(M, Mnow); // M = T_{12} * P_{2}

				// compute the transition matrix at the fowles_layer-substrate interface
				fresnel T23;

				T23.compute_T(TE, RI_fowles_layer, RI_sub, T12.get_theta_t());

				Mnow = T23.transition_matrix(); 

				M = vecut::cmat_cmat_product(M, Mnow); // M = T_{12} * P_{2} * T_{23}

				r_coeff = M[0][1] / M[0][0];

				t_coeff = 1.0 / M[0][0]; 

				if (loud) std::cout << lambda << " , " << RI_clad << " , " << RI_fowles_layer << " , " << RI_sub << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";

				write << lambda << " , " << RI_clad << " , " << RI_fowles_layer << " , " << RI_sub << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";
			}

			write.close();
		}
		else {
			std::string reason;
			reason = "Error: void multifowles_layer::compute_r_t(int n_fowles_layers, double fowles_layer_thickness, bool loud)\n";
			if (!c1) reason += "n_fowles_layers: " + template_funcs::toString(n_fowles_layers) + " is not valid\n";
			if (!c2) reason += "fowles_layer_thickness: " + template_funcs::toString(fowles_layer_thickness, 2) + " is not valid\n";
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

HL_stack::HL_stack(sweep &swp_obj, material *h_fowles_layer, material *l_fowles_layer, material *the_cladding, material *the_substrate)
{
	set_params(swp_obj, h_fowles_layer, l_fowles_layer, the_cladding, the_substrate);
}

HL_stack::~HL_stack()
{
	M.clear(); r.clear(); t.clear();
}

void HL_stack::set_params(sweep &swp_obj, material *h_fowles_layer, material *l_fowles_layer, material *the_cladding, material *the_substrate)
{
	// assign the parameters that make up the fowles_layer structure

	try {

		bool c1 = swp_obj.defined();
		bool c2 = h_fowles_layer != nullptr ? true : false;
		bool c4 = l_fowles_layer != nullptr ? true : false;
		bool c5 = the_cladding != nullptr ? true : false;
		bool c6 = the_substrate != nullptr ? true : false;
		bool c10 = c1 && c2 && c4 && c5 && c6;

		if (c10) {
			wavelength.set_vals(swp_obj); // assign the values for the wavelength parameter space

										  // assign the material objects
			layer_high = h_fowles_layer;
			layer_low = l_fowles_layer;
			substrate = the_substrate;
			cladding = the_cladding;

			int size = 2;
			M = vecut::zero_cmat(size, size);
		}
		else {
			std::string reason;
			reason = "Error: void multifowles_layer::set_params(sweep &swp_obj, material *the_fowles_layer, material *the_cladding, material *the_substrate)\n";
			if (!c1) reason += "swp_obj is not correct\n";
			if (!c4) reason += "the_fowles_layer has not been correctly assigned";
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

void HL_stack::compute_r_t_Fowles(int n_fowles_layers, double fowles_layer_thickness, bool loud)
{
	// Compute the transfer matrix for a system with n_fowles_layers each having the same thickness
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = n_fowles_layers > 0 ? true : false;
		bool c2 = fowles_layer_thickness > 0 ? true : false;
		bool c10 = c1 && c2 && wavelength.defined();

		if (c10) {
			// sweep over all wavelengths
			int total_fowles_layers = 2 * n_fowles_layers;
			double lambda, RI_h, RI_l, RI_clad, RI_sub;
			std::complex<double> reflectance, transmittance;
			std::string filename;
			filename = "Fowles_HL_stack_Data_" + template_funcs::toString(n_fowles_layers) + "_fowles_layers_" + template_funcs::toString(fowles_layer_thickness, 2) + "_fowles_layer_Thickness.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results " + template_funcs::toString(n_fowles_layers) + " fowles_layers,  fowles_layer thickness: " + template_funcs::toString(fowles_layer_thickness, 2) + " um\n";
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {

				lambda = wavelength.get_val(i); // assign wavelength value

				layer_high->set_wavelength(lambda); // assign wavelength to fowles_layer material object

				layer_low->set_wavelength(lambda); // assign wavelength to fowles_layer material object

				cladding->set_wavelength(lambda); // assign wavelength to cladding material object

				substrate->set_wavelength(lambda); // assign wavelength to substrate material object

				RI_h = layer_high->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				RI_l = layer_low->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				// since this is an anti-reflection film the same matrix is used for all fowles_layers for fixed wavelength
				fowles_layer H_fowles_layer(fowles_layer_thickness, lambda, RI_h); // compute transfer matrix for the fowles_layer at this wavelength

				fowles_layer L_fowles_layer(fowles_layer_thickness, lambda, RI_l); // compute transfer matrix for the fowles_layer at this wavelength

				std::vector<std::vector<std::complex<double>>> Mnow;

				// get transfer matrix for all fowles_layers
				// for the anti-reflection filter multiply the single fowles_layer TM by itself n_fowles_layers times
				for (int j = 0; j < total_fowles_layers; j++) {

					if (j % 2 == 0) {
						Mnow = L_fowles_layer.transfer_matrix();
					}
					else {
						Mnow = H_fowles_layer.transfer_matrix();
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
			reason = "Error: void HL_stack::compute_r_t_Fowles(int n_fowles_layers, double fowles_layer_thickness)\n";
			if (!c1) reason += "n_fowles_layers: " + template_funcs::toString(n_fowles_layers) + " is not valid\n";
			if (!c2) reason += "fowles_layer_thickness: " + template_funcs::toString(fowles_layer_thickness, 2) + " is not valid\n";
			if (!wavelength.defined()) reason += "wavelength sweep object not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void HL_stack::compute_r_t(int n_fowles_layers, double fowles_layer_thickness, bool loud)
{
	// Compute the transfer matrix for a system with n_fowles_layers each having the same thickness
	// R. Sheehan 15 - 7 - 2019

	try {
		bool c1 = n_fowles_layers > 0 ? true : false;
		bool c2 = fowles_layer_thickness > 0 ? true : false;
		bool c10 = c1 && c2 && wavelength.defined();

		if (c10) {
			// sweep over all wavelengths
			int total_fowles_layers = 2 * n_fowles_layers;
			double lambda, RI_h, RI_l, RI_clad, RI_sub;
			std::complex<double> r_coeff, t_coeff;
			std::string filename;
			filename = "HL_stack_Data_" + template_funcs::toString(n_fowles_layers) + "_fowles_layers_" + template_funcs::toString(fowles_layer_thickness, 2) + "_fowles_layer_Thickness.txt";
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);
			if (loud) std::cout << "Results " + template_funcs::toString(n_fowles_layers) + " fowles_layers,  fowles_layer thickness: " + template_funcs::toString(fowles_layer_thickness, 2) + " um\n";
			for (int i = 0; i < wavelength.get_Nsteps(); i++) {

				lambda = wavelength.get_val(i); // assign wavelength value

				layer_high->set_wavelength(lambda); // assign wavelength to fowles_layer material object

				layer_low->set_wavelength(lambda); // assign wavelength to fowles_layer material object

				cladding->set_wavelength(lambda); // assign wavelength to cladding material object

				substrate->set_wavelength(lambda); // assign wavelength to substrate material object

				RI_h = layer_high->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				RI_l = layer_low->refractive_index(); // compute wavelength dependent refractive index at this wavelength

				RI_clad = cladding->refractive_index(); // equivalent to n_{0} in formulae

				RI_sub = substrate->refractive_index(); // equivalent to n_{t} in formulae

				// since this is an anti-reflection film the same matrix is used for all fowles_layers for fixed wavelength
				fresnel T12, T_hl, T_lh; 

				T12.compute_T(TE, RI_clad, RI_h, zero); 

				T_hl.compute_T(TE, RI_h, RI_l, T12.get_theta_t()); 

				T_lh.compute_T(TE, RI_l, RI_h, T_hl.get_theta_t()); 

				propagation P_h, P_l; 

				P_h.compute_P(lambda, RI_h, fowles_layer_thickness, T12.get_theta_t()); 

				P_l.compute_P(lambda, RI_l, fowles_layer_thickness, T12.get_theta_t()); 

				std::vector<std::vector<std::complex<double>>> Mnow;

				// get transfer matrix for all fowles_layers
				// for the anti-reflection filter multiply the single fowles_layer TM by itself n_fowles_layers times

				M = T12.transition_matrix(); 

				for (int j = 0; j < total_fowles_layers; j++) {
					if (j % 2 == 0) {
						Mnow = P_h.propagation_matrix(); 
						M = vecut::cmat_cmat_product(M, Mnow); 
						Mnow = T_hl.transition_matrix(); 
						M = vecut::cmat_cmat_product(M, Mnow);
					}
					else {
						Mnow = P_l.propagation_matrix(); 
						M = vecut::cmat_cmat_product(M, Mnow);
						Mnow = T_lh.transition_matrix();
						M = vecut::cmat_cmat_product(M, Mnow);
					}					
				}

				// compute reflectance and transmittance of the structure at this wavelength
				
				r_coeff = M[0][1] / M[0][0];

				t_coeff = 1.0 / M[0][0];

				if (loud) std::cout << lambda << " , " << RI_clad << " , " << RI_h << " , " << RI_l << " , " << RI_sub << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";

				write << lambda << " , " << RI_clad << " , " << RI_h << " , " << RI_l << " , " << RI_sub << " , " << abs(r_coeff) << " , " << abs(t_coeff) << "\n";
			}

			write.close();
		}
		else {
			std::string reason;
			reason = "Error: void HL_stack::compute_r_t_Fowles(int n_fowles_layers, double fowles_layer_thickness)\n";
			if (!c1) reason += "n_fowles_layers: " + template_funcs::toString(n_fowles_layers) + " is not valid\n";
			if (!c2) reason += "fowles_layer_thickness: " + template_funcs::toString(fowles_layer_thickness, 2) + " is not valid\n";
			if (!wavelength.defined()) reason += "wavelength sweep object not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

