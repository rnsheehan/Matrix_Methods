#ifndef FOWLES_MULTIfowles_layer_H
#define FOWLES_MULTIfowles_layer_H

// classes for implementing the matrix method for computing the reflectance / transmittance spectra of multifowles_layer thin films
// method is described in Fowles, section 4.4
// R. Sheehan 15 - 7 - 2019

namespace spectrum {
	void compute(double &n_clad, double &n_sub, std::vector<std::vector<std::complex<double>>> &M, std::complex<double> &r, std::complex<double> &t);
}

class fowles_layer {
public:
	fowles_layer();
	fowles_layer(double thickness, double wavelength, double ref_index);
	fowles_layer(const fowles_layer &lobj); //copy constructor
	~fowles_layer();

	void set_params(double thickness, double wavelength, double ref_index);

	std::vector<std::vector<std::complex<double>>> transfer_matrix();

private:
	bool defined;
	int size; // transfer matrix dimensions

	std::vector<std::vector<std::complex<double>>> M; // array to store transfer matrix of this fowles_layer
};

#endif
