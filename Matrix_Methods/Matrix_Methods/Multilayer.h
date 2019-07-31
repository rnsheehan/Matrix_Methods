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

	inline std::complex<double> get_theta_t() { return theta_t;  } // output transmission angle
	inline std::complex<double> get_theta_critical() { return theta_critical;  } // output transmission angle
	inline std::complex<double> get_theta_brewster() { return theta_brewster;  } // output transmission angle

	std::complex<double> reflection(bool polarisation, std::complex<double> angle);
	std::complex<double> transmission(bool polarisation, std::complex<double> angle);

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

#endif
