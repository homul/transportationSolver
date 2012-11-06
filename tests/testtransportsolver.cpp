#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "common.h"
#include "testcommon.h"
#include "smallobjects.h"
#include "primalSolver.h"
#include "smallfunctions.h"

using namespace std;
using namespace DD;
using namespace OUT;

struct TestTransportSolver: public CppUnit::TestFixture {
CPPUNIT_TEST_SUITE( TestTransportSolver );
		CPPUNIT_TEST( test_primalSolver_min );
		CPPUNIT_TEST( test_primalSolver_max );
CPPUNIT_TEST_SUITE_END();


void PrintABbin(ofstream& fout, const SparseArray& a, const SparseArray& b,
			const Matrix& bin) {
		fout << "a:" << a << endl;
		fout << "b:" << b << endl;
		fout << "bin:" << bin << endl;
	}
	;

	template<class Solver>
	bool simpleNoPrintABbinTest(ofstream& fout, Solver& solver,
			const SparseArray& a, const SparseArray& b, const Matrix& bin) {
		fout << "solver.Solve(a,b);" << endl;
		std::pair<bool, floatType> res = solver.Solve(a, b);
		solver.PrintTestData(fout);
		fout << "result=" << res.second << endl;
		return res.first;
	}
	;

	template<class Solver>
	bool simpleTest(ofstream& fout, Solver& solver, const SparseArray& a,
			const SparseArray& b, const Matrix& bin) {
		PrintABbin(fout, a, b, bin);
		return simpleNoPrintABbinTest(fout, solver, a, b, bin);
	}
	;

	template<class Solver>
	bool simpleInitTest(ofstream& fout, Solver& solver, const SparseArray& a,
			const SparseArray& b, const Matrix& bin) {
		PrintABbin(fout, a, b, bin);
		fout << "-------------------------" << endl;
		solver.Init(a.size(), b.size(), bin);
		fout << "solver.Init(a.size(),b.size(),bin);" << endl;
		solver.PrintTestData(fout);
		fout << "-------------------------" << endl;
		return simpleNoPrintABbinTest(fout, solver, a, b, bin);
	}
	;

	void test_primalSolver_min() {

		size_t a_size = 3, b_size = 4;

		string outTst("testPrimalSolver.tst"), outChk("testPrimalSolver.chk");
		ofstream fout(outTst.c_str());

		PrimalSolver<false, Matrix, SparseArray> solver(fout);
		SparseArray a(a_size, 0.0);
		SparseArray b(b_size, 0.0);
		Matrix bin(a_size, b_size, 0.0);

		a[0] = 0.7;
		a[1] = 0.3;
		b[0] = 1.0;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		bin = ublas::trans(bin);
		fout << "Transpose: bin=ublas::trans(bin);" << endl;
		CPPUNIT_ASSERT(simpleInitTest(fout,solver,b,a,bin));

		fout << "======================================================"
				<< endl;
		CPPUNIT_ASSERT(simpleTest(fout,solver,b,a,bin));

		fout << "======================================================"
				<< endl;

		a[0] = 0.6;
		a[1] = 0.4;
		b[0] = 1.0;
		simpleTest(fout, solver, b, a, bin);

		fout << "======================================================"
				<< endl;
		bin = ublas::trans(bin);
		fout << "Transpose: bin=ublas::trans(bin);" << endl;

		b[0] = 0.5;
		b[2] = 0.5;
		bin(0, 2) = -1.0;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		b[0] = 0.45;
		b[2] = 0.55;
		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;

		b[0] = 0.4;
		b[2] = 0.6;
		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;

		b[0] = 0.3;
		b[2] = 0.7;

		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		b[0] = 0.7;
		b[2] = 0.3;
		bin(1, 0) = -1.0;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		a[0] = 0.5;
		a[1] = 0.5;
		b[0] = 0.5;
		b[2] = 0.5;
		bin(0, 2) = bin(1, 0) = 0;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		bin(0, 2) = -1;
		bin(1, 0) = -1;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));
		fout << "======================================================"
				<< endl;
		bin(0, 2) = 0;
		bin(1, 0) = 0;
		bin(0, 0) = 1;
		bin(1, 2) = 1;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout.close();

		CPPUNIT_ASSERT(TST::areEqualTXTfiles(outTst,outChk));
	}
	;

	void test_primalSolver_max() {
		size_t a_size = 3, b_size = 4;

		string outTst("testPrimalSolverMax.tst"), outChk(
				"testPrimalSolverMax.chk");
		ofstream fout(outTst.c_str());

		PrimalSolver<true, Matrix, SparseArray> solver(fout);
		SparseArray a(a_size, 0.0);
		SparseArray b(b_size, 0.0);
		Matrix bin(a_size, b_size, 0.0);

		a[0] = 0.7;
		a[1] = 0.3;
		b[0] = 1.0;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		bin = ublas::trans(bin);
		fout << "Transpose: bin=ublas::trans(bin);" << endl;
		CPPUNIT_ASSERT(simpleInitTest(fout,solver,b,a,bin));

		fout << "======================================================"
				<< endl;
		CPPUNIT_ASSERT(simpleTest(fout,solver,b,a,bin));

		fout << "======================================================"
				<< endl;

		a[0] = 0.6;
		a[1] = 0.4;
		b[0] = 1.0;
		simpleTest(fout, solver, b, a, bin);

		fout << "======================================================"
				<< endl;
		bin = ublas::trans(bin);
		fout << "Transpose: bin=ublas::trans(bin);" << endl;

		b[0] = 0.5;
		b[2] = 0.5;
		bin(0, 2) = 1.0;//!

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		b[0] = 0.45;
		b[2] = 0.55;
		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;

		b[0] = 0.4;
		b[2] = 0.6;
		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;

		b[0] = 0.3;
		b[2] = 0.7;

		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		b[0] = 0.7;
		b[2] = 0.3;
		bin(1, 0) = 1.0;//!

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		a[0] = 0.5;
		a[1] = 0.5;
		b[0] = 0.5;
		b[2] = 0.5;
		bin(0, 2) = bin(1, 0) = 0;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;
		bin(0, 2) = 1;//!
		bin(1, 0) = 1;//!

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));
		fout << "======================================================"
				<< endl;
		bin(0, 2) = 0;
		bin(1, 0) = 0;
		bin(0, 0) = -1;//!
		bin(1, 2) = -1;//!

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;

		a.resize(3);
		a.clear();
		b.resize(2);
		b.clear();
		a[0] = 1.7435310497537848e-01;
		a[1] = 5.2634502820855411e-01;
		a[2] = 2.9930186681606746e-01;
		b[0] = 4.9999828832362131e-01;
		b[1] = 5.0000171167637864e-01;
		bin.resize(3, 2);
		bin.clear();
		bin(0, 0) = -5.7719507421236260e-01;
		bin(0, 1) = -9.9572685733238553e-01;
		bin(1, 0) = -8.2402806394921058e-01;
		bin(1, 1) = -2.7476051835099258e-02;
		bin(2, 0) = -3.1663470730028798e-01;
		bin(2, 1) = -7.1173718325408974e-01;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));
		fout << "======================================================"
				<< endl;

		a.resize(5);
		a.clear();
		b.resize(5);
		b.clear();
		a[0] = 1.7435310497537848e-01;
		a[1] = 5.2634502820855411e-01;
		a[4] = 2.9930186681606746e-01;
		b[1] = 4.9999828832362131e-01;
		b[3] = 5.0000171167637864e-01;
		bin.resize(5, 5);
		bin.clear();
		bin(0, 1) = -5.7719507421236260e-01;
		bin(0, 3) = -9.9572685733238553e-01;
		bin(1, 1) = -8.2402806394921058e-01;
		bin(1, 3) = -2.7476051835099258e-02;
		bin(4, 1) = -3.1663470730028798e-01;
		bin(4, 3) = -7.1173718325408974e-01;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;

		a.resize(5);
		a.clear();
		b.resize(5);
		b.clear();
		a[1] = 0.5000043047810037;
		a[2] = 0.3145007227206573;
		a[4] = 0.185494972498339;
		b[0] = 1.8549336768616226e-01;
		b[1] = 5.0000863021898023e-01;
		b[2] = 3.1449800209485768e-01;
		bin.resize(5, 5);
		bin.clear();
		bin(1, 0) = -4.8675867402309492e-01;
		bin(1, 1) = -4.2754023332499907e-01;
		bin(1, 2) = -3.9162254305166311e-01;
		bin(2, 0) = -3.8280589756686517e-01;
		bin(2, 1) = -4.4818462964528455e-01;
		bin(2, 2) = -4.3849451487813825e-03;
		bin(4, 0) = -3.7402430101019533e-02;
		bin(4, 1) = -4.1292247754192096e-01;
		bin(4, 2) = -1.2075783830171351e-01;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		fout << "======================================================"
				<< endl;

		a.resize(5);
		a.clear();
		b.resize(5);
		b.clear();
		bin.resize(5, 5);
		bin.clear();

		a[1] = 0.5575673158533609;
		a[2] = 1.114993197639578e-14;
		a[3] = 0.0002022685480762959;
		a[4] = 0.4422304155985516;
		b[0] = 4.9999999999999906e-01;
		b[1] = 5.7567311981475813e-02;
		b[2] = 1.1149878848962176e-14;
		b[3] = 4.4243268415490600e-01;
		b[4] = 3.8636079718507516e-09;

		bin(1, 0) = -2.0928710010335180e-01;
		bin(1, 1) = -6.9043032857143805e-02;
		bin(1, 2) = -4.8170472168443013e-01;
		bin(1, 3) = -2.4935258261363608e-01;
		bin(1, 4) = -6.2896015151821083e-02;
		bin(2, 0) = -6.5427324299433889e-02;
		bin(2, 1) = -4.5734588962856021e-01;
		bin(2, 2) = -8.7675430387107386e-02;
		bin(2, 3) = -1.6126576399489576e-01;
		bin(2, 4) = -7.7213012649311227e-02;
		bin(3, 0) = -2.9230315484679453e-01;
		bin(3, 1) = -1.5189365933271762e-01;
		bin(3, 2) = -4.2349867169908184e-01;
		bin(3, 3) = -1.7548031880309820e-01;
		bin(3, 4) = -2.7961696464550539e-01;
		bin(4, 0) = -3.0708473236629957e-01;
		bin(4, 1) = -2.9169612158634517e-01;
		bin(4, 2) = -4.1554130889267721e-01;
		bin(4, 3) = -3.2917768942619566e-02;
		bin(4, 4) = -4.9921409017369807e-02;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));//degenerate sample - solved!

		fout << "======================================================"
				<< endl;

		a.resize(5);
		a.clear();
		b.resize(5);
		b.clear();
		bin.resize(5, 5);
		bin.clear();

		a[0] = 1.0945721212020032e-02;
		a[1] = 7.8500920995152549e-03;
		a[2] = 4.9999999999999994e-01;
		a[3] = 5.0636737242057930e-10;
		a[4] = 4.8120418618209737e-01;
		b[0] = 1.8795813311535212e-02;
		b[1] = 4.8120418618209743e-01;
		b[3] = 4.9999999999999994e-01;
		b[4] = 5.0636737226537033e-10;

		bin(0, 0) = -3.2094540089412844e-05;
		bin(0, 1) = -3.3488978717238166e-01;
		bin(0, 3) = -3.8749841455719358e-01;
		bin(0, 4) = -4.5889282154799105e-01;
		bin(1, 0) = -6.7508118258560136e-03;
		bin(1, 1) = -3.3595426559259850e-01;
		bin(1, 3) = -1.8055272762689402e-02;
		bin(1, 4) = -3.0566241443420872e-01;
		bin(2, 0) = -3.0921508805324096e-02;
		bin(2, 1) = -3.9383357083137777e-01;
		bin(2, 3) = -1.2068737234020949e-01;
		bin(2, 4) = -1.8499515656614451e-01;
		bin(3, 0) = -3.2840442882310805e-01;
		bin(3, 1) = -3.9982155030584965e-01;
		bin(3, 3) = -2.0073001026209911e-01;
		bin(3, 4) = -8.7864666752454204e-02;
		bin(4, 0) = -2.0557793891317114e-01;
		bin(4, 1) = -3.8355358428487259e-02;
		bin(4, 3) = -4.1237064004520452e-01;
		bin(4, 4) = -7.2320660377070620e-02;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));//degenerate sample -solved!

		fout << "======================================================"
				<< endl;

		a.resize(5);
		a.clear();
		b.resize(5);
		b.clear();
		bin.resize(5, 5);
		bin.clear();

		a[0] = 0.4999997201416438;
		a[1] = 0.4999999999999865;
		a[2] = 1.647485963204486e-09;
		a[3] = 2.782108822148093e-07;
		a[4] = 1.554179656352422e-15;
		b[0] = 6.3915970976122964e-05;
		b[2] = 4.9993608402902390e-01;
		b[3] = 1.7982985409599921e-09;
		b[4] = 4.9999999820170149e-01;
		bin(0, 0) = -8.8680443907473439e-01;
		bin(0, 2) = -3.2541941028340693e-01;
		bin(0, 3) = -7.2328323811445538e-01;
		bin(0, 4) = -2.7609481815066877e-01;
		bin(1, 0) = -2.2025188534532295e-01;
		bin(1, 2) = -2.7578421136168030e-01;
		bin(1, 3) = -1.8242431859598696e-01;
		bin(1, 4) = -1.0939130564657566e-01;
		bin(2, 0) = -6.1232874524422398e-01;
		bin(2, 2) = -1.5695524362705426e-01;
		bin(2, 3) = -1.9092470509508844e-01;
		bin(2, 4) = -7.8137388116744066e-01;
		bin(3, 0) = -1.7764167356195007e-01;
		bin(3, 2) = -1.9759737616293011e-01;
		bin(3, 3) = -5.8480153446309335e-01;
		bin(3, 4) = -2.6497888903365419e-01;
		bin(4, 0) = -8.8875147508864827e-01;
		bin(4, 2) = -5.2007304342467009e-01;
		bin(4, 3) = -2.5063470762718221e-01;
		bin(4, 4) = -5.2750310233212216e-01;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));//non-degenerate sample with very tricky change cycle


		fout << "======================================================"
				<< endl;

		a.resize(5);
		a.clear();
		b.resize(5);
		b.clear();
		bin.resize(5, 5);
		bin.clear();

		a[0] = 1.6562624859286637e-01;
		a[1] = 5.0000000000000000e-01;
		a[2] = 8.0469000374127753e-16;
		a[3] = 5.0052772984710205e-02;
		a[4] = 2.8432097842242260e-01;
		b[0] = 8.3437375140713366e-01;
		b[1] = 1.6562624859286648e-01;

		bin(0, 0) = -4.9733555503065491e-01;
		bin(0, 1) = -6.2414872256300817e-02;
		bin(0, 3) = -4.0319993202723559e-01;
		bin(1, 0) = -4.9131766869282240e-01;
		bin(1, 1) = -3.6392812610786784e-01;
		bin(1, 3) = -2.1232362241126765e-01;
		bin(2, 0) = -2.1988744950848046e-01;
		bin(2, 1) = -1.4974199475242850e-01;
		bin(2, 3) = -9.4569968103696583e-02;
		bin(3, 0) = -1.4957876929528022e-01;
		bin(3, 1) = -1.3792129030354380e-01;
		bin(3, 3) = -3.5031708835173264e-01;
		bin(4, 0) = -2.6579652785593483e-01;
		bin(4, 1) = -3.3286854640248631e-01;
		bin(4, 3) = -4.7295699267320196e-01;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));

		a[0] = 6.0489174426365568e-04;
		a[1] = 5.0000029225787046e-01;
		a[2] = 1.7628118093514564e-14;
		a[3] = 4.5648615230238321e-01;
		a[4] = 4.2908663695464956e-02;
		b[0] = 9.9226849681444906e-01;
		b[1] = 7.7312109276805995e-03;
		b[3] = 2.9225787051339848e-07;

		CPPUNIT_ASSERT(simpleTest(fout,solver,a,b,bin));

		/* Strange results - should be checked! */
		fout << "======================================================"
				<< endl;
		{
			size_t size_a = 10;
			size_t size_b = 10;

			a.resize(size_a);
			a.clear();
			b.resize(size_b);
			b.clear();
			bin.resize(size_a, size_b);
			bin.clear();

			a[0] = 0.1647713610296626;
			a[1] = 0.08491628789095651;
			a[2] = 0;
			a[3] = 0.007591662402871276;
			a[4] = 0.0692068273017841;
			a[5] = 0.1378866994980621;
			a[6] = 0.1201267121599603;
			a[7] = 0.129550931229258;
			a[8] = 0.1234157349523229;
			a[9] = 0.1625337835351224;

			b[0] = 1.6416763959202094e-01;
			b[1] = 8.4753826195601401e-02;
			b[2] = 0.0000000000000000e+0;
			b[3] = 8.5496710579625480e-03;
			b[4] = 6.8787239373389486e-02;
			b[5] = 1.3766554442553283e-01;
			b[6] = 1.2345690770621527e-01;
			b[7] = 1.2690954047312653e-01;
			b[8] = 1.2365091785770507e-01;
			b[9] = 1.6205871331844585e-01;

			for (size_t i = 0; i < size_a; ++i)
				for (size_t j = 0; j < size_b; ++j)
					bin(i, j) = 0.0;//-1.5000000000000001e+300;


			bin(0, 0) = 1.4134369362889209e-02;
			bin(1, 0) = 7.0039781535737173e+00;
			bin(1, 1) = 6.2227589603276706e-02;
			bin(3, 3) = 3.6912404267482812e-01;
			bin(4, 3) = 2.2840937434047026e+00;
			bin(4, 4) = 2.3950186750097363e-01;
			bin(5, 3) = 1.5800036457371920e+01;
			bin(5, 4) = 2.8701466661522090e+00;
			bin(5, 5) = 9.5167759676023686e-02;
			bin(6, 5) = 4.1837010499262357e+00;
			bin(6, 6) = 3.2643038066388674e-01;
			bin(7, 6) = 2.4477619821455061e+00;
			bin(7, 7) = 4.0586363818291943e-01;
			bin(8, 7) = 2.1591636305546387e+00;
			bin(8, 8) = 2.4417243903913571e-01;
			bin(9, 8) = 2.8434490099972587e+00;
			bin(9, 9) = -0.0000000000000000e+00;

			CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));
		}
		fout << "======================================================"
				<< endl;
		{
			fout << std::setprecision(16);

			size_t size_a = 7;
			size_t size_b = 7;

			a.resize(size_a);
			a.clear();
			b.resize(size_b);
			b.clear();
			bin.resize(size_a, size_b);
			bin.clear();

			a[0] =0;
			a[1] =0;
			a[2] =0.04166666666666667;
			a[3] =0;
			a[4] =0.4583333333333333;
			a[5] =0;
			a[6] =0.5;


			b[0] =0.0000000000000000e+00;
			b[1] =0.0000000000000000e+00;
			b[2] =1.0000000000000001e-01;
			b[3] =1.0000000000000001e-01;
			b[4] =2.9999999999999999e-01;
			b[5] =0.0000000000000000e+00;
			b[6] =5.0000000000000000e-01;

			bin(2,3)=-4.6000000000000000e+01;
			bin(2,4)=-1.0000000000000000e+02;
			bin(2,6)=-9.4000000000000000e+01;
			bin(4,2)=-1.0000000000000000e+02;
			bin(4,3)=-9.1000000000000000e+01;
			bin(4,6)=-4.3000000000000000e+01;
			bin(6,2)=-9.4000000000000000e+01;
			bin(6,3)=-8.4000000000000000e+01;
			bin(6,4)=-4.3000000000000000e+01;

			CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));
		}
		fout << "======================================================"<<endl;
		{
			fout << std::setprecision(16);

			size_t size_a = 4;
			size_t size_b = 4;

			a.resize(size_a);
			a.clear();
			b.resize(size_b);
			b.clear();
			bin.resize(size_a, size_b);
			bin.clear();

		a[0]=0.0;
		a[1]=0.5;
		a[2]=0.5;
		a[3]=0.0;

		b[0]=0.5;
		b[1]=0.25;
		b[2]=0.25;
		b[3]=0.0;

		bin(0,0)=-9.0807908361688199e+00;
		bin(0,1)=-3.8136225562140447e+00;
		bin(0,2)=-3.0047243314351526e+00;
		bin(0,3)=-6.2414581172361316e+00;
		bin(1,0)=-8.9799662589933575e+00;
		bin(1,1)=-2.7842199394405918e+00;
		bin(1,2)=-2.5604832071626946e+00;
		bin(1,3)=-5.0563965842856087e+00;
		bin(2,0)=3.8778798317433707e+00;
		bin(2,1)=8.9384883674040854e+00;
		bin(2,2)=9.5527275070793589e+00;
		bin(2,3)=6.3065144928668220e+00;
		bin(3,0)=-1.0308430025730482e+01;
		bin(3,1)=-4.4883860981969104e+00;
		bin(3,2)=-3.5362179051322009e+00;
		bin(3,3)=-6.3281743836254698e+00;

		CPPUNIT_ASSERT(simpleInitTest(fout,solver,a,b,bin));
		}

		fout.close();

		CPPUNIT_ASSERT(TST::areEqualTXTfiles(outTst,outChk));
	};
};
CPPUNIT_TEST_SUITE_REGISTRATION( TestTransportSolver );
