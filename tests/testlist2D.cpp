#include <cppunit/extensions/HelperMacros.h>
#include <fstream>
#include <vector>
#include "common.h"
#include "testcommon.h"
#include "list2D.h"
#include "smallobjects.h"
#include "smallfunctions.h"

using namespace std;
using namespace DD;
using namespace OUT;

struct TestList2D: public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE( TestList2D );
	CPPUNIT_TEST( test_List2D );
	CPPUNIT_TEST_SUITE_END();

	void test_List2D()
	{
		size_t xsize=3;
		size_t ysize=4;
		size_t nnz=xsize+ysize;
		List2D<floatType> lst(xsize,ysize,nnz);

		std::string foutname("testList2D.tst"), chkname("testList2D.chk");
		std::ofstream fout(foutname.c_str());

		fout << "An empty list:"<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(lst.push(0,1,10.0));
		fout << endl<< "lst.push(0,1,10.0)"<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(lst.push(0,0,9.0));
		fout << endl<< "lst.push(0,0,9.0)"<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(lst.push(0,1,11.0));
		fout << endl<<"Again lst.push(0,1,11.0)"<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(lst.push(0,2,12.0));
		fout << endl<< "lst.push(0,2,12.0)"<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(lst.push(2,3,13.0));
		fout << endl<< "lst.push(0,2,12.0)"<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(lst.push(1,2,14.0));
		fout << endl<< "lst.push(1,2,14.0)"<<endl;
		lst.PrintTestData(fout);


		CPPUNIT_ASSERT(!lst.push(1,3,15.0));
		fout << endl<< "lst.push(1,3,15.0): it should NOT change the state of the List2D."<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(lst.insert(1,3,15.0));
		fout << endl<< "lst.insert(1,3,15.0): however, insertion is still possible:"<<endl;
		lst.PrintTestData(fout);

		CPPUNIT_ASSERT(!lst.insert(0,3,15.0));
		fout << endl<< "lst.insert(0,3,15.0): it should NOT change the state of the List2D."<<endl;
		lst.PrintTestData(fout);

		List2D<floatType>::iterator it=lst.rowBegin(1);
		fout <<"it=lst.rowBegin(1), *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;
		++it;
		fout <<"++it, *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;
		it=it.changeDir();
		fout <<"it.changeDir(), *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;
		++it;
		fout <<"++it, *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;

		it=lst.colBegin(1);
		fout <<"it=lst.colBegin(1), *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;
		++it;
		fout <<"++it, *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;
		it=it.changeDir();
		fout <<"it.changeDir(), *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;
		++it;
		fout <<"++it, *it="<<*it<<", it.coordinate()="<<it.coordinate()<<endl;

		List2D<floatType>::iterator beg=lst.colBegin(2), end=lst.colEnd(2);
		CPPUNIT_ASSERT(beg!=end);
		lst.erase(beg);
		fout << endl<< "lst.erase(lst.colBegin(2))"<<endl;
		lst.PrintTestData(fout);

		beg=lst.rowBegin(1), end=lst.rowEnd(1);
		CPPUNIT_ASSERT(beg!=end);
		lst.erase(beg);
		fout << endl<< "lst.erase(lst.rowBegin(1))"<<endl;
		lst.PrintTestData(fout);

		ublas::vector<floatType> bin(xsize*ysize,1.0);
		for (size_t x=0;x<xsize;++x)
			for (size_t y=0;y<ysize;++y)
				bin(xsize*y+x)=xsize*y+x;
		fout << "bin: "<<endl<<bin<<endl;
		fout << "<lst,bin>="<<lst.inner_product1D(bin)<<endl;

		lst.rowErase(2);
		fout << "lst.rowErase(2):"<<endl;
		lst.PrintTestData(fout);

		lst.colErase(0);
		fout <<"lst.colErase(0):"<<endl;
		lst.PrintTestData(fout);

		fout.close();

		CPPUNIT_ASSERT(TST::areEqualTXTfiles(foutname,chkname));
	}
};
CPPUNIT_TEST_SUITE_REGISTRATION( TestList2D );

