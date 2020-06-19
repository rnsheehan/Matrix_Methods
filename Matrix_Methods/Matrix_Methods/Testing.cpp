#ifndef ATTACH_H
#include "Attach.h"
#endif

void test::matrix_test_1()
{
	// Test the matrix product calculation code
	int r1, c1, r2, c2; 

	r1 = c1 = r2 = c2 = 2; 

	// Allocate the memory for the matrices
	std::vector<std::vector<double>> M1;
	std::vector<std::vector<double>> M2;
	std::vector<std::vector<double>> M3;

	for (int i = 0; i < r1; i++) {
		M1.push_back( std::vector<double>(c1, 0.0) ); 
		M2.push_back( std::vector<double>(c2, 0.0) ); 
	}

	// Assign elements to the matrices
	M1[0][0] = 2; M1[0][1] = 6; 
	M1[1][0] = 3; M1[1][1] = 5; 

	M2[0][0] = -5; M2[0][1] = 1;
	M2[1][0] = 4; M2[1][1] = -3;

	// Result of product M1*M2
	//M3[0][0] = 14; M3[0][1] = -16; M3[1][0] = 5; M3[1][1] = -12;

	// Result of product M2*M1
	//M3[0][0] = -7; M3[0][1] = -25; M3[1][0] = -1; M3[1][1] = 9;

	M3 = vecut::mat_mat_product(M1, M2); 
	//M3 = vecut::mat_mat_product(M2, M1);

	// print the result to the screen
	/*for (size_t i = 0; i < M3.size(); i++) {
		for (size_t j = 0; j < M3[0].size(); j++)
			std::cout << M3[i][j] << " ";
		std::cout << "\n";
	}*/

	vecut::print_to_screen(M3); 
}

void test::matrix_test_2()
{
	// Test the matrix product calculation code
	int r1, c1, r2, c2;

	r1 = 5; c1 = 4; 
	r2 = 4; c2 = 7; 

	// Allocate the memory for the matrices
	std::vector<std::vector<double>> M1;
	std::vector<std::vector<double>> M2;
	std::vector<std::vector<double>> M3;

	for (int i = 0; i < r1; i++) M1.push_back(std::vector<double>(c1, 0.0));

	for (int i = 0; i < r2; i++) M2.push_back(std::vector<double>(c2, 0.0));

	// Assign elements to the matrices
	M1[0][0] = 2; M1[0][1] = 6; M1[0][2] = 1; M1[0][3] = 0;
	M1[1][0] = 3; M1[1][1] = 4; M1[1][2] = 3; M1[1][3] = 2;
	M1[2][0] = -5; M1[2][1] = -3; M1[2][2] = 8; M1[2][3] = -5;
	M1[3][0] = 1; M1[3][1] = 0; M1[3][2] = 7; M1[3][3] = -8;
	M1[4][0] = 0; M1[4][1] = 1; M1[4][2] = 0; M1[4][3] = 6;

	M2[0][0] = -5; M2[0][1] = 1; M2[0][2] = 8; M2[0][3] = -9; M2[0][4] = 6; M2[0][5] = 5; M2[0][6] = 7;
	M2[1][0] = 5; M2[1][1] = 0; M2[1][2] = 1; M2[1][3] = -11; M2[1][4] = 3; M2[1][5] = 3; M2[1][6] = 8;
	M2[2][0] = 1; M2[2][1] = 2; M2[2][2] = 3; M2[2][3] = 0; M2[2][4] = 4; M2[2][5] = 2; M2[2][6] = 9;
	M2[3][0] = 0; M2[3][1] = 7; M2[3][2] = 7; M2[3][3] = 1; M2[3][4] = 1; M2[3][5] = 0; M2[3][6] = 4;

	// Result of product M1*M2
	/*M3[0][0] = 21; M3[0][1] = 4; M3[0][2] = 25; M3[0][3] = -84; M3[0][4] = 34; M3[0][5] = 30; M3[0][6] = 71;
	M3[1][0] = 8; M3[1][1] = 23; M3[1][2] = 51; M3[1][3] = -69; M3[1][4] = 44; M3[1][5] = 33; M3[1][6] = 88;
	M3[2][0] = 18; M3[2][1] = -24; M3[3][2] = -54; M3[3][3] = 73; M3[3][4] = -12; M3[3][5] = -18; M3[3][6] = -7;
	M3[3][0] = 2; M3[3][1] = -41; M3[3][2] = -27; M3[3][3] = -17; M3[3][4] = 26; M3[3][5] = 19; M3[3][6] = 38;
	M3[3][0] = 5; M3[3][1] = 42; M3[3][2] = 43; M3[3][3] = -5; M3[3][4] = 9; M3[3][5] = 3; M3[3][6] = 32;*/

	M3 = vecut::mat_mat_product(M1, M2); 

	// print the result to the screen
	/*for (size_t i = 0; i < M3.size(); i++) {
		for (size_t j = 0; j < M3[0].size(); j++)
			std::cout << M3[i][j] << " ";
		std::cout << "\n";
	}*/

	vecut::print_to_screen(M3); 
}

void test::layer_test()
{
	// see if the layer class is working correctly
	// R. Sheehan 15 - 7 - 2019

	double angle_in = (20.0)* DEG_TO_RAD; // input angle in units of radians
	double wavelength = 1550; // wavelength in nm
	double thickness = 0.75*wavelength; // layer thickness in nm
	double cladding_index = 1.0; // RI value at some wavelength
	double layer_index = 3.45; // RI value at some wavelength
	double substrate_index = 1.45; // RI value at some wavelength

	fowles_layer l1; 
	fowles_layer l2; 

	l1.set_params(TE, angle_in, thickness, wavelength, cladding_index, layer_index, substrate_index, true); 

	// test the alternative R calculator
	double p0 = l1.get_p_in(), p2 = l1.get_p_out(), R, T;
	std::complex<double> r, t; 
	std::vector<std::vector<std::complex<double>>> M = l1.transfer_matrix(); 
	spectrum::compute_r_t(p0, p2, M, r, t, R, T); 

	std::cout << "R calc test\n";
	std::cout << "r: " << r << "\n"; 
	std::cout << "t: " << r << "\n"; 
	std::cout << "R: " << R << "\n"; 
	std::cout << "T: " << T << "\n"; 
	std::cout << "R + T: " << R + T << "\n\n"; 

	l1.set_params(TM, angle_in, thickness, wavelength, cladding_index, layer_index, substrate_index, true); 

	p0 = l1.get_p_in(), p2 = l1.get_p_out();
	M = l1.transfer_matrix();
	spectrum::compute_r_t(p0, p2, M, r, t, R, T);

	std::cout << "R calc test\n";
	std::cout << "r: " << r << "\n";
	std::cout << "t: " << r << "\n";
	std::cout << "R: " << R << "\n";
	std::cout << "T: " << T << "\n";
	std::cout << "R + T: " << R + T << "\n\n";

	// loop over the layer length and output the reflectivity
	std::string filename = "Air_Silicon_Silica_R_T.txt"; 
	//std::string filename = "Air_Silica_Silicon_R_T.txt";
	
	std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

	if (write.is_open()) {
		angle_in = (20.0) * DEG_TO_RAD; // input angle in units of radians
		
		thickness = 1000.0; 
		while (thickness < 2000.0) {
			l1.set_params(TE, angle_in, thickness, wavelength, cladding_index, layer_index, substrate_index); 
			l2.set_params(TM, angle_in, thickness, wavelength, cladding_index, layer_index, substrate_index); 
			write << thickness << " , " << l1.get_R() << " , " << l1.get_T() << " , " << l2.get_R() << " , " << l2.get_T() << "\n";
			thickness += 1.0; 
		}
		write.close(); 
	}
}

void test::AR_Coating()
{
	// Use Fowles / Born & Wolf method to compute the reflection spectrum of an AR coating
	// R. Sheehan 17 - 6 - 2020

	double angle_in = (0.0) * DEG_TO_RAD; // input angle in units of radians
	double wavelength = 550; // wavelength in nm
	
	double cladding_index = 1.0; // RI value at some wavelength
	double layer_index = 1.35; // n1 ~ sqrt(n2) for AR coating
	double substrate_index = 1.5; // RI value at some wavelength

	// layer thickness must be chosen according to n1 l = lambda / 4 for AR coating
	double thickness = wavelength / (4.0*layer_index); // layer thickness in nm

	std::cout << "AR Coating Thickness: " << thickness << "\n\n"; 

	fowles_layer l1;

	l1.set_params(TE, angle_in, thickness, wavelength, cladding_index, layer_index, substrate_index, true);

	// loop over the layer length and output the reflectivity
	std::string filename = "Air_MgF2_Glass.txt"; 
	//std::string filename = "Air_SiO2_SiN.txt"; 

	std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

	if (write.is_open()) {

		fowles_layer l1; 
		fowles_layer l2; 
		fowles_layer l3; 
		fowles_layer l4; 
		fowles_layer l5; 

		wavelength = 400.0;
		while (wavelength < 701.0) {
			thickness = 1550.0 / (4.0 * layer_index); // layer thickness in nm

			l1.set_params(TE, angle_in, thickness, wavelength, cladding_index, layer_index, substrate_index);
			l2.set_params(TE, angle_in, 2.0*thickness, wavelength, cladding_index, layer_index, substrate_index);
			l3.set_params(TE, angle_in, 3.0*thickness, wavelength, cladding_index, layer_index, substrate_index);
			l4.set_params(TE, angle_in, 4.0*thickness, wavelength, cladding_index, layer_index, substrate_index);
			l5.set_params(TE, angle_in, 5.0*thickness, wavelength, cladding_index, layer_index, substrate_index);
			
			write << wavelength << " , " << l1.get_R() << " , " << l2.get_R() << " , " << l3.get_R() << " , " << l4.get_R() << " , " << l5.get_R() << "\n";

			wavelength += 1.0;
		}

		write.close(); 
	}	
}

void test::HR_Coating()
{
	// Use Fowles Born & Wolf to compute spectrum of HR coating 
	// which comprises an even number of alternating layers of high RI and low RI material
	// each layer must have optical thickness of lambda_{0} / 4 for peak reflectivity at lambda_{0}
	// R. Sheehan 18 - 6 - 2020

	bool pol = TE; 
	
	int n_layer_pairs = 14; // minimum value here should be 2
	// min value of 2 means initial layer + 1 internal layer + substrate layer
	// value of 4 means initial layer + 3 internal layer + substrate layer
	// value of 6 means initial layer + 5 internal layer + substrate layer
	// value of 8 means initial layer + 7 internal layer + substrate layer etc..

	double angle_in = (0.0) * DEG_TO_RAD; // input angle in units of radians
	double lambda_0 = 550; // wavelength in nm

	double cladding_index = 1.0; // RI value at some wavelength
	double high_index = 2.3; 
	double low_index = 1.35; 
	double substrate_index = 2.0; // RI value at some wavelength

	// define the layer thicknesses lambda_0
	double t_high = lambda_0 / (4.0 * high_index), t_low = lambda_0 / (4.0 * low_index);
	//double t_high = lambda_0 / (4.0 * high_index), t_low = t_high;

	std::cout << "Layer thicknesses\n";
	std::cout << "High: " << t_high << "\n"; 
	std::cout << "Low: " << t_low << "\n"; 

	// sweep over all wavelengths
	int n_wl = 600, size = 2;
	double l1 = 250, l2 = 1000, p_first, p_last; 
	sweep wavelengths(n_wl, l1, l2); 

	std::string filename = "HR_Coating_" + template_funcs::toString(n_layer_pairs+1) + ".txt";
	std::ofstream write; 

	for (int i = 0; i < wavelengths.get_Nsteps(); i++) {
		
		double n0, n2, angle, R, T;

		std::complex<double> r, t; 

		std::vector<std::vector<std::complex<double>>> M = vecut::zero_cmat(size, size); 
		std::vector<std::vector<std::complex<double>>> Mnow = vecut::zero_cmat(size, size); 

		// compute the transfer matrix for the high and low RI layers
		fowles_layer high;
		fowles_layer low;

		for (int j = 0; j <= n_layer_pairs; j++) {

			if (j == 0) {
				// air-high interface
				n0 = cladding_index; 
				n2 = low_index;
				angle = angle_in; 
			}
			else if (j == n_layer_pairs) {
				//low-substrate interface
				n0 = low_index;
				n2 = substrate_index; 
				angle = low.get_theta_out();
			}
			else {
				// high-low interfaces
				n0 = low_index; 
				n2 = low_index;
				angle = low.get_theta_out(); 
			}

			high.set_params(pol, angle, t_high, wavelengths.get_val(i), n0, high_index, low_index);
			
			if (j == 0) {
				M = high.transfer_matrix(); 
				p_first = high.get_p_in();
			}
			else {
				Mnow = high.transfer_matrix();
				M = vecut::cmat_cmat_product(M, Mnow); // Mh Ml Mh
			}			

			low.set_params(pol, high.get_theta_out(), t_low, wavelengths.get_val(i), high_index, low_index, n2);
			Mnow = low.transfer_matrix();
			M = vecut::cmat_cmat_product(M, Mnow); // Mh Ml Mh Ml

			if (j == n_layer_pairs) {
				p_last = low.get_p_out();
			}
		}
		
		// M contains the product of all the transfer matrices
		// Compute the reflection coeffs etc for this wavelength
		spectrum::compute_r_t(p_first, p_last, M, r, t, R, T);

		// output the computed data
		if (i == 0) {
			write.open(filename, std::ios_base::out, std::ios_base::trunc);
		}
		else {
			write.open(filename, std::ios_base::app);
		}
		
		if (write.is_open()) {
			write << wavelengths.get_val(i) << " , " << R << " , " << T << "\n"; 
			write.close(); 
		}
	}
}

void test::AR_filter_test()
{
	// example calculation for computing reflectance of simple structure
	// for this choice of materials and layer thicknesses you are designing an anti-reflection filter
	// RI_{SiN} ~ 2, RI_{Si02} ~ 1.4 ~ \sqrt{2} for \lambda = 1.55 um

	int n_pts, n_layers;
	double start, stop, W;

	sweep WL;
	
	Air ri_air;
	SiN ri_sin;
	SiO2 ri_sio2;

	n_pts = 21; start = 1.52; stop = 1.59;
	WL.set_vals(n_pts, start, stop);

	multilayer_old calc;

	calc.set_params(WL, &ri_sin, &ri_air, &ri_sio2); 

	n_layers = 5; W = 1.55 / 4.0; 
	calc.compute_r_t(n_layers, W, true);
}

void test::high_low_test()
{
	// example calculation for computing reflectance of an alternating structure
	// R. Sheehan 23 - 7 - 2019

	int n_pts, n_layers;
	double start, stop, W;

	sweep WL;

	Air ri_air;
	SiN ri_sin;
	SiO2 ri_sio2;
	Si ri_si; 

	n_pts = 201; start = 1.52; stop = 1.59;
	WL.set_vals(n_pts, start, stop);

	HL_stack calc;

	calc.set_params(WL, &ri_sin, &ri_sio2, &ri_air, &ri_si);

	n_layers = 15; W = 1.55 / 4.0;
	calc.compute_r_t(n_layers, W);
}

void test::r_t_test()
{
	// compute the transmisson and reflection coefficients for a dielectric interface
	// R. Sheehan 31 - 7 - 2019

	bool pol = TE; 

	double n_in, n_out, angle; 

	// Need to add conservation of energy confirmation

	std::string filename = "Air_to_Glass_R_T.txt"; 
	n_out = 1.5; n_in = 1.0; angle = PI_6;

	/*std::string filename = "Glass_to_Air_R_T.txt";
	n_out = 1.0; n_in = 1.5; angle = PI_6;*/

	fresnel iface; 

	iface.set_params(n_in, n_out); 

	std::complex<double> r = iface.reflection(pol, angle); 
	std::complex<double> t = iface.transmission(pol, angle); 

	std::cout << "r: " << r << "\n";
	std::cout << "t: " << t << "\n";
	std::cout << "TE Amplitude Relation t - r: " << t - r << "\n"; // this should be unity
	std::cout << "TM Amplitude Relation (n_{2}/n_{1})t - r: " << (iface.get_n_ratio()*t) - r << "\n"; // this should be unity
	std::cout << "Power Reflection Coefficient R: " << iface.Reflectivity(pol, angle) << "\n"; 
	std::cout << "Power Transmission Coefficient T: " << iface.Transmissivity(pol, angle) << "\n";
	std::cout << "Conservation of Energy R + T: " << iface.Reflectivity(pol, angle) + iface.Transmissivity(pol, angle) << "\n"; 

	// output the air-glass coefficients
	
	std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

	if (write.is_open()) {
		
		write << "n1, " << n_in << ", n2, " << n_out << ", n1 / n2, " << n_in / n_out << "\n"; 
		write << "Critical Angle (rad), " << iface.get_theta_critical().real() << " , " << iface.get_theta_critical().imag() << "\n"; 
		write << "Brewster Angle (rad), " << iface.get_theta_brewster().real() << " , " << iface.get_theta_brewster().imag() << "\n";
		write << "angle (rad), r_{TE}, t_{TE}, r_{TM}, t_{TM}\n"; 

		int n_angle; 
		double d_angle, angle_min, angle_max;
		//double v1, v2, v3, v4;
		std::complex<double> v1, v2, v3, v4;

		n_angle = 201; 
		angle_min = 0.0; angle_max = PI_2; 
		d_angle = (angle_max - angle_min) / ( static_cast<double>(n_angle - 1) ); 

		angle = angle_min; 
		while (angle < angle_max + d_angle) {
			
			// compute and output reflection and transmission coefficients
			v1 = iface.reflection(TE, angle); 
			v2 = iface.transmission(TE, angle); 
			v3 = iface.reflection(TM, angle); 
			v4 = iface.transmission(TM, angle); 

			write << angle << " , " << v1.imag() << " , " << v2.imag() << " , " << v3.imag() << " , " << v4.imag() << "\n"; 

			// compute and output power reflection and power transmision coefficients
			/*v1 = iface.Reflectivity(TE, angle);
			v2 = iface.Transmissivity(TE, angle);
			v3 = iface.Reflectivity(TM, angle);
			v4 = iface.Transmissivity(TM, angle);

			write << angle << " , " << v1 << " , " << v2 << " , " << v3 << " , " << v4 << "\n";*/ 

			angle += d_angle; 
		}


		write.close(); 
	}
}

void test::fp_test()
{
	// Use the multilayer structure to compute the R, T spectrum of a FP cavity
	// R. Sheehan 13 - 8 - 2019

	int n_pts;
	double start, stop, W = 0.0;

	// Declarate the objects
	sweep WL;

	Air ri_air;
	SiN ri_sin;
	SiO2 ri_sio2;
	Si ri_si;

	// Fill the parameter space
	n_pts = 5; start = 1.2; stop = 1.6;
	WL.set_vals(n_pts, start, stop);

	// Create layer stack
	std::vector<layer> the_layers; 

	W = 0.0; the_layers.push_back( layer(W, &ri_air) ); 
	W = 2; the_layers.push_back( layer(W, &ri_si) );
	W = 0.0; the_layers.push_back( layer(W, &ri_sio2) );

	// Compute the r, t spectra
	multilayer compute; 

	compute.set_params(WL, the_layers); 

	compute.compute_spectrum(TE, true); 
}

void test::ar_coating()
{
	// compute the spectrum of ar stack
	// SiN layers on SiO2 layers on Si substrate

	int n_pts;
	double start, stop, W = 0.0, Wt = 0.25*1.55;

	// Declarate the objects
	sweep WL;

	Air ri_air;
	SiN ri_sin;
	SiO2 ri_sio2;
	Si ri_si;

	// Fill the parameter space
	n_pts = 51; start = 1.52; stop = 1.59;
	WL.set_vals(n_pts, start, stop);

	// Create layer stack
	std::vector<layer> the_layers;

	the_layers.push_back(layer(W, &ri_air));
	Wt *= 9;  the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(W, &ri_sin));

	// Compute the r, t spectra
	multilayer compute;

	compute.set_params(WL, the_layers);

	compute.compute_spectrum(TE, true);
}

void test::hr_coating()
{
	// compute the spectrum of hr stack
	// SiN layers on SiO2 layers on Si substrate

	int n_pts;
	double start, stop, W = 0.0, Wt = 1.55/4.0;

	// Declarate the objects
	sweep WL;

	Air ri_air;
	SiN ri_sin;
	SiO2 ri_sio2;
	Si ri_si;

	// Fill the parameter space
	n_pts = 201; start = 1.0; stop = 1.6;
	WL.set_vals(n_pts, start, stop);

	// Create layer stack
	std::vector<layer> the_layers;

	the_layers.push_back(layer(W, &ri_air));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_si));
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(W, &ri_air));

	// Compute the r, t spectra
	multilayer compute;

	compute.set_params(WL, the_layers);

	compute.compute_spectrum(TE, true);
}

void test::fp_filter()
{
	// compute the spectrum of hr stack
	// SiN layers on SiO2 layers on Si substrate

	int n_pts;
	double start, stop, W = 0.0, Wt = 0.25*1.55;

	// Declarate the objects
	sweep WL;

	Air ri_air;
	SiN ri_sin;
	SiO2 ri_sio2;
	Si ri_si;

	// Fill the parameter space
	n_pts = 101; start = 1.2; stop = 1.6;
	WL.set_vals(n_pts, start, stop);

	// Create layer stack
	std::vector<layer> the_layers;

	the_layers.push_back(layer(W, &ri_air));

	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_sin));

	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_sin));
	
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_sin));
	
	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_sin));

	the_layers.push_back(layer(Wt, &ri_sio2));
	the_layers.push_back(layer(Wt, &ri_sin));
	the_layers.push_back(layer(Wt, &ri_sin));
	the_layers.push_back(layer(Wt, &ri_sio2));
	
	the_layers.push_back(layer(Wt, &ri_sin));
	the_layers.push_back(layer(Wt, &ri_sio2));

	the_layers.push_back(layer(Wt, &ri_sin));
	the_layers.push_back(layer(Wt, &ri_sio2));

	the_layers.push_back(layer(Wt, &ri_sin));
	the_layers.push_back(layer(Wt, &ri_sio2));

	the_layers.push_back(layer(Wt, &ri_sin));
	the_layers.push_back(layer(Wt, &ri_sio2));
	
	the_layers.push_back(layer(W, &ri_si));

	// Compute the r, t spectra
	multilayer compute;

	compute.set_params(WL, the_layers);

	compute.compute_spectrum(TE, true);
}