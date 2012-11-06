#include "smallfunctions.h"
#include <fstream>

namespace TST{
	bool areEqualTXTfiles(const string& name1,const string& name2)
	{
		std::ifstream f1,f2;
		openfile_throw(name1,f1);
	        openfile_throw(name2,f2);

		string s1,s2;
		while (!f1.eof() && !f2.eof())
		{
			std::getline(f1,s1);
			std::getline(f2,s2);
			if (s1!=s2) 
				return false;
		}
		if (!(f1.eof() && f2.eof()))
			return false;

		return true;
	};
}
