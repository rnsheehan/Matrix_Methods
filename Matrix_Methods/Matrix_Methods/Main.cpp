#ifndef ATTACH_H
#include "Attach.h"
#endif

int main() 
{
	//test::matrix_test_1();

	//test::r_t_test();

	bool pol = TM; 
	double angle_in = (20.0) * DEG_TO_RAD; // input angle in units of radians

	test::layer_test(pol, angle_in); 

	test::layer_test_alt(pol, angle_in);

	//test::AR_Coating(); 

	//test::HR_Coating(); 

	//test::BP_Filter(); 

	std::cout << "Press enter to close\n";
	std::cin.get();
}