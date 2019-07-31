#ifndef FOWLES_MULTILAYER_H
#define FOWLES_MULTILAYER_H

// classes for implementing the matrix method for computing the reflectance / transmittance spectra of multilayer thin films
// method is described in Fowles, section 4.4
// R. Sheehan 15 - 7 - 2019

namespace spectrum {
	void compute(double &n_clad, double &n_sub, std::vector<std::vector<std::complex<double>>> &M, std::complex<double> &r, std::complex<double> &t);
}

class layer {
public:
	layer();
	layer(double thickness, double wavelength, double ref_index);
	layer(layer &lobj); //copy constructor
	~layer();

	void set_params(double thickness, double wavelength, double ref_index);

	std::vector<std::vector<std::complex<double>>> transfer_matrix();

private:
	bool defined;
	int size; // transfer matrix dimensions

	std::vector<std::vector<std::complex<double>>> M; // array to store transfer matrix of this layer
};

#endif
