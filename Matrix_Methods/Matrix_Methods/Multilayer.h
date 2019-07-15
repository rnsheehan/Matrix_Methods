#ifndef MULTILAYER_H
#define MULTILAYER_H

// class for implementing the matrix method for computing the reflectance / transmittance spectra of multilayer thin films
// method is described in Fowles, section 4.4
// R. Sheehan 15 - 7 - 2019

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

	std::vector<std::vector<std::complex<double>>> M; // transfer matrix
};



class multilayer {
public:
	multilayer(); 
	multilayer(int n_layers, sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate);
	~multilayer(); 

	void set_params(int n_layers, sweep &swp_obj, material *the_layer, material *the_cladding, material *the_substrate);

	void build_transfer_matrix(); // make N and thickness an input parameter here. 

private:
	int N; // number of layers in the structure

	double layer_thickness; // assume all layers have the same thickness for now
	
	// the sweep object defines the wavelength sweep space, all wavelength values are in units of um
	sweep wavelength;

	material *layer_mat; // object for the layer material
	material *substrate; // object for the substrate material
	material *cladding; // object for the cladding material

	//std::vector<layer> layers; // do you actually need this?
	std::vector<std::vector<std::complex<double>>> M; // transfer matrix
};

#endif
