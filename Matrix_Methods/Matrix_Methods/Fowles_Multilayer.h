#ifndef FOWLES_MULTIfowles_layer_H
#define FOWLES_MULTIfowles_layer_H

// classes for implementing the matrix method for computing the reflectance / transmittance spectra of multifowles_layer thin films
// method is described in Fowles, section 4.4 and Born and Wolf, section 1.6
// R. Sheehan 15 - 7 - 2019

namespace spectrum {
	void compute(double &n_clad, double &n_sub, std::vector<std::vector<std::complex<double>>> &M, std::complex<double> &r, std::complex<double> &t);
}

// class implements a method for computing the transfer matrix of a layer of dielectric material
// assumes either TE, TM polarisation and input angle between 0 and PI/2
// R. Sheehan 17 - 6 - 2020

class fowles_layer {
public:
	fowles_layer();
	fowles_layer(bool pol, double theta_in, double thickness, double wavelength, double n0, double n1, double n2);
	fowles_layer(const fowles_layer &lobj); //copy constructor
	~fowles_layer();

	inline double get_theta_out() { return theta_out;  }
	inline double get_R() { return R;  }
	inline double get_T() { return T;  }

	void set_params(bool pol, double theta_in, double thickness, double wavelength, double n0, double n1, double n2, bool loud = false);
	void layer_stats(); 

	std::vector<std::vector<std::complex<double>>> transfer_matrix();

private:
	bool defined;
	
	double theta_out; // wavevector angle after the layer
	double R; // layer Power reflectivity
	double T; // layer Power transmissivity

	std::complex<double> r; // layer reflection coefficient
	std::complex<double> t; // layer transmission coefficient

	std::vector<std::vector<std::complex<double>>> M; // array to store transfer matrix of this fowles_layer
};

#endif
