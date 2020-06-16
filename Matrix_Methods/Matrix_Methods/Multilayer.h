#ifndef MULTILAYER_H
#define MULTILAYER_H

// Implementation of the 2*2 matrix for computing transmission / reflection spectra of multilayer structures
// Implemented here is the method described by Yariv and Yeh in "Photonics"
// Their method seems to be more widely applicable, and following from Coldren and Corzine can be extended to 
// include structures that have loss / gain elements. It also more accurately describes DBR structures
// R. Sheehan 31 - 7 - 2019

// In this implementation a "layer" comprises a transition matrix for the dielectric interface
// and a propagation matrix for propagation across the layer thickness
// A final transition matrix for the layer-substrate dielectric interface will also be required. 

// First step is to implement Fresnel equations for reflection and transmission coefficients that can handle arbitrary 
// input angles and RI ratios. Going to implement it so that it takes RI values as input

// class for computing the reflection and transmission coefficients at a dielectric interface
// propagation is assumed to take place from left to right
// R. Sheehan 31 - 7 - 2019

class fresnel {
public:
	fresnel(); 
	fresnel(double n_left, double n_right);
	~fresnel(); 

	void set_params(double n_left, double n_right);	

	inline double get_n_ratio() { return nrat;  } // return n_{2} / n_{1}
	

	inline std::complex<double> get_theta_t() { return theta_t;  } // output transmission angle
	inline std::complex<double> get_theta_critical() { return theta_critical;  } // output transmission angle
	inline std::complex<double> get_theta_brewster() { return theta_brewster;  } // output transmission angle
	inline std::complex<double> get_cos_ratio() { return cos_ratio; } // return cos(theta_{1}) / cos(theta_{1})

	std::complex<double> reflection(bool polarisation, std::complex<double> angle);
	std::complex<double> transmission(bool polarisation, std::complex<double> angle);
	
	double Reflectivity(bool polarisation, std::complex<double> angle);
	double Transmissivity(bool polarisation, std::complex<double> angle); 

	void compute_T(bool polarisation, double n_left, double n_right, std::complex<double> angle);

	std::vector<std::vector<std::complex<double>>> transition_matrix();

private:
	void set_angles(std::complex<double> angle);

private:
	bool params_defined; 

	double n1; // RI on input side, external material
	double n2; // RI on output side, internal material

	double nrat; // ratio n2/n1, nrat > 1 => external reflection, nrat < 1 => internal reflection 
	double nrat_sqr; // square of RI ratio	

	std::complex<double> theta_in; // propagation angle on input side
	std::complex<double> theta_t; // propagtaion angle on output side, angle of transmission computed using Snell's law

	// critical angle only relevant in the case of internal reflection n2 < n1, nrat < 1
	// internal reflection is real valued for angle < critical_angle, otherwise reflection is complex valued
	// nrat > 1 => external reflection, nrat < 1 => internal reflection
	std::complex<double> theta_critical; // critical angle

	// compute the brewster angle based on the input refractive index values
	// TM reflection = 0 for light incident at Brewster angle
	// Unpolarised light incident at Brewster angle is reflected in linearly polarised TE state
	// however, TE reflectivity is not great, generation of linearly polarised light in this manner is not efficient
	std::complex<double> theta_brewster; // brewster angle

	std::complex<double> cos_theta_1; // cosine of the input propagation angle 
	std::complex<double> cos_theta_2; // cosine of the transmission propagation angle 
	std::complex<double> cos_ratio; // ratio of the cosines of the angles 

	std::vector<std::vector<std::complex<double>>> T;  // array to store dielectric interface transition matrix
};

// class for computing the propagation matrix across a dielectric layer of thickness d
// R. Sheehan 31 - 7 - 2019

class propagation {

public:
	propagation(); 
	propagation(double wavelength, double index, double thickness, std::complex<double> angle); 
	~propagation(); 

	void set_params(double wavelength, double index, double thickness, std::complex<double> angle); 

	void compute_P(double wavelength, double index, double thickness, std::complex<double> angle);

	std::vector<std::vector<std::complex<double>>> propagation_matrix();

private:
	bool params_defined; 

	double lambda; // wavelength of the light in the dielectric layer
	double n; // refractive index of the dielectric layer, this should be wavelength dependent
	double d; // thickness of the dielectric layer

	std::complex<double> theta; // propagation angle in the dielectric layer
	std::complex<double> phase; // phase term in the dielectric layer

	std::vector<std::vector<std::complex<double>>> P;  // array to store dielectric layer propagation matrix
};

// implementation of layer class that can be used to compute properties of multilayer stacks
// R. Sheehan 13 - 8 -2019

class layer {
public:
	layer();
	layer(double &d, material *mat); 
	layer(const layer &lobj);

	void set_params(double &d, material *mat);

	inline double get_d() { return thickness; }

	double get_RI(double wavelength); 
private:
	double thickness; 

	material *the_mat; 
};

// implementation of a class that can be used to compute the reflection and transmission spectrum of an arbitrary multilayer stack
// the idea is to input a list of layers and then compute its r and t spectra
// the user will create the layer list and there will be one method for computing its spectra
// R. Sheehan 13 - 8 - 2019

class multilayer {
public:
	multilayer(); 
	multilayer(sweep &swp_obj, std::vector<layer> &layer_list);
	~multilayer(); 

	void set_params(sweep &swp_obj, std::vector<layer> &layer_list);

	void compute_spectrum(bool polarisation, bool loud = false);

private:
	std::vector<std::vector<std::complex<double>>> transmission_matrix(double &lambda, bool polarisation, bool loud = false);

private:
	sweep wavelength; // the sweep object defines the wavelength sweep space, all wavelength values are in units of um

	std::vector<std::complex<double>> r; // vector to hold computed reflectance values 

	std::vector<std::complex<double>> t; // vector to hold computed transmittance values

	std::vector<layer> the_layers; // list holding all info on the layers in the stack
};

// Implementation of class that is used to compute reflectance and transmittance spectrum of a multilayer stack
// in which each of the layers is composed of the same material. 

class multilayer_old {
public:
	multilayer_old();
	multilayer_old(sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate);
	~multilayer_old();

	void set_params(sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate);

	void compute_r_t_Fowles(int n_layers, double layer_thickness, bool loud = false); // make N and thickness an input parameter here. 

	void compute_r_t(int n_layers, double layer_thickness, bool loud = false);

private:
	// the sweep object defines the wavelength sweep space, all wavelength values are in units of um
	sweep wavelength;

	material *layer_mat; // object for the layer material
	material *substrate; // object for the substrate material
	material *cladding; // object for the cladding material

	std::vector<std::complex<double>> r; // vector to hold computed reflectance values 
	std::vector<std::complex<double>> t; // vector to hold computed transmittance values

	std::vector<std::vector<std::complex<double>>> M;  // array to store transfer matrix of overall structure
};

// Implementation of class that is used to compute reflectance and transmittance spectrum of a multilayer stack
// in which the layers are composed alternating materials
// R. Sheehan 23 - 7 - 2019

class HL_stack {
public:
	HL_stack();
	HL_stack(sweep &swp_obj, material *h_layer, material *l_layer, material *the_cladding, material *the_substrate);
	~HL_stack();

	void set_params(sweep &swp_obj, material *h_layer, material *l_layer, material *the_cladding, material *the_substrate);

	void compute_r_t_Fowles(int n_layers, double layer_thickness, bool loud = false); // make N and thickness an input parameter here. 

	void compute_r_t(int n_layers, double layer_thickness, bool loud = false); // make N and thickness an input parameter here. 

private:
	// the sweep object defines the wavelength sweep space, all wavelength values are in units of um
	sweep wavelength;

	material *layer_high; // object for the high RI layer material
	material *layer_low; // object for the low RI layer material
	material *substrate; // object for the substrate material
	material *cladding; // object for the cladding material

	std::vector<std::complex<double>> r; // vector to hold computed reflectance values 
	std::vector<std::complex<double>> t; // vector to hold computed transmittance values

	std::vector<std::vector<std::complex<double>>> M;  // array to store transfer matrix of overall structure
};

#endif
