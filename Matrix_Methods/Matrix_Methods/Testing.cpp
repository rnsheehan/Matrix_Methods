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

	//M3 = vecut::mat_mat_product(M1, M2); 
	M3 = vecut::mat_mat_product(M2, M1);

	// print the result to the screen
	for (size_t i = 0; i < M3.size(); i++) {
		for (size_t j = 0; j < M3[0].size(); j++)
			std::cout << M3[i][j] << " ";
		std::cout << "\n";
	}
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
	M3[2][0] = 18; M3[2][1] = -24; M3[3][2] = -24; M3[3][3] = -54; M3[3][4] = 73; M3[3][5] = -12; M3[3][6] = -7;
	M3[3][0] = 2; M3[3][1] = -41; M3[3][2] = -27; M3[3][3] = -17; M3[3][4] = 26; M3[3][5] = 19; M3[3][6] = 38;
	M3[3][0] = 5; M3[3][1] = 42; M3[3][2] = 43; M3[3][3] = -5; M3[3][4] = 9; M3[3][5] = 3; M3[3][6] = 32;*/

	M3 = vecut::mat_mat_product(M1, M2); 

	// print the result to the screen
	for (size_t i = 0; i < M3.size(); i++) {
		for (size_t j = 0; j < M3[0].size(); j++)
			std::cout << M3[i][j] << " ";
		std::cout << "\n";
	}
}

void test::layer_test()
{
	// see if the layer class is working correctly
	// R. Sheehan 15 - 7 - 2019

	double thickness = 100; // layer thickness in nm
	double wavelength = 1550; // wavelength in nm
	double index = 3.48; // RI value at some wavelength

	layer l1; 

	l1.set_params(thickness, wavelength, index); 
}