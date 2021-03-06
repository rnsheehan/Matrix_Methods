#ifndef TESTING_H
#define TESTING_H

// namespace of functions to test the code

namespace test {
	void matrix_test_1(); 

	void matrix_test_2(); 

	void layer_test(bool pol, double angle_in);

	void layer_test_alt(bool pol, double angle_in); // input angle in units of radians);

	void AR_Coating();

	void HR_Coating(); 

	void BP_Filter(); 

	//void AR_filter_test(); 

	//void high_low_test(); 

	//void r_t_test(); 

	//void fp_test(); 

	//void ar_coating(); 

	//void hr_coating(); 

	//void fp_filter();
}

#endif
