#ifndef FOWLES_MULTIfowles_layer_H
#define FOWLES_MULTIfowles_layer_H

// classes for implementing the matrix method for computing the reflectance / transmittance spectra of multifowles_layer thin films
// method is described in Fowles, section 4.4 and Born and Wolf, section 1.6
// R. Sheehan 15 - 7 - 2019

namespace spectrum {
	void compute_r_t(double &p_in, double &p_out, std::vector<std::vector<std::complex<double>>> &M, std::complex<double> &r, std::complex<double> &t, double &R, double &T);
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
	inline double get_theta_layer() { return theta_layer;  }
	inline double get_p_in() { return p0;  }
	inline double get_p_layer() { return p1;  }
	inline double get_p_out() { return p2;  }
	inline double get_R() { return R;  }
	inline double get_T() { return T;  }
	inline std::complex<double> get_r() { return r; }
	inline std::complex<double> get_t() { return t; }

	std::vector< std::vector< std::complex<double> > > transfer_matrix();

	void set_params(bool pol, double theta_in, double thickness, double wavelength, double n0, double n1, double n2, bool loud = false);
	void layer_stats(); 

private:
	bool defined;
	
	double theta_inc; // wavevector angle at input
	double theta_layer; // wavevector angle inside the layer
	double theta_out; // wavevector angle after the layer

	double R; // layer Power reflectivity
	double T; // layer Power transmissivity
	
	double p0; // polarisation dependent wavevector scale factor at input
	double p1; // polarisation dependent wavevector scale factor inside layer
	double p2; // polarisation dependent wavevector scale factor at output

	std::complex<double> r; // layer reflection coefficient
	std::complex<double> t; // layer transmission coefficient

	std::vector<std::vector<std::complex<double>>> M; // array to store transfer matrix of this fowles_layer
};

#endif
