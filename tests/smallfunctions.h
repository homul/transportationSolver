#ifndef TESTSMALLFUNCTIONS_H
#define TESTSMALLFUNCTIONS_H
#include <string>
#include <cstdlib>// for rand()

using namespace std;

namespace TST{
using std::string;
bool areEqualTXTfiles(const string& name1,const string& name2);

template <class FileStream>
void openfile_throw(const string& file_name,FileStream& file,const string& prefix="")
{
 file.open(file_name.c_str());
 if (file.fail())
  throw prefix+string("Can not open file: ")+file_name;
};

double RandomDouble(double max);

}
#endif
