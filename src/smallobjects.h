#ifndef SMALLOBJECTS_H_
#define SMALLOBJECTS_H_
#include <functional>
#include <numeric>
#include "common.h"
#include <iostream>
#include <algorithm>
#include <vector>

namespace DD{

 /*
  * copies indices of non-zero elements of a sparse array
  */
 template<class SrcIterator,class TrgIterator>
 TrgIterator copyNZindices(SrcIterator first,SrcIterator last,TrgIterator result)
 {
	 for (;first!=last;++first)
	    *result++ = first.index();
	  return result;
 }

 /****************opt_element******************************/

 template<bool maximize,class Iterator>
 class opt_element
 {
 	public:
 	//opt_element(){};
 	Iterator operator()(Iterator begin,Iterator end);
 };

 template<class Iterator>
 class opt_element <true, Iterator>
 {
 	public:
 	//opt_element(){};
 	Iterator operator()(Iterator begin,Iterator end){return std::max_element(begin,end);};
 };

 template<class Iterator>
 class opt_element <false, Iterator>
 {
 	public:
 	//opt_element(){};
 	Iterator operator()(Iterator begin,Iterator end){return std::min_element(begin,end);};
 };

 /****************opt******************************/

 template<bool maximize,class T>
 class opt
 {
 	public:
 	//opt(){};
 	const T& operator()(const T& t1,const T& t2);
 };

 template<class T>
 class opt<true, T>
 {
 	public:
 	//opt(){};
 	const T& operator()(const T& t1,const T& t2){return std::max(t1,t2);};
 };

 template<class T>
 class opt<false, T>
 {
 	public:
 	//opt(){};
 	const T& operator()(const T& t1,const T& t2){return std::min(t1,t2);};
 };

 /********************better**********************/
 template<bool maximize,class T>
 class better
 {
 	public:
 	bool operator()(const T& t1,const T& t2)const;
 };

 template<class T>
 class better<true, T>
 {
 	public:
 	bool operator()(const T& t1,const T& t2)const{return t1>t2;};
 };

 template<class T>
 class better<false, T>
 {
 	public:
 	bool operator()(const T& t1,const T& t2)const{return t1<t2;};
 };

 /********************notWorse**********************/
 template<bool maximize,class T>
 class notWorse
 {
 	public:
 	bool operator()(const T& t1,const T& t2)const;
 };

 template<class T>
 class notWorse<true, T>
 {
 	public:
 	bool operator()(const T& t1,const T& t2)const{return t1>=t2;};
 };

 template<class T>
 class notWorse<false, T>
 {
 	public:
 	bool operator()(const T& t1,const T& t2)const{return t1<=t2;};
 };

 template<class Iterator>
 floatType _Normalize(Iterator begin,Iterator end)
 {
 	if (begin==end) return 0.0;

	floatType acc=std::accumulate(begin,end,(floatType)0.0);

 	if (acc!=0.0)
 	  std::transform(begin,end,begin,std::bind1st(std::multiplies<floatType>(),1.0/acc));
 	else
 	 *begin=1.0;

 	return acc;
 };

 template<class PointerList>
 void DeallocatePointerContainer(PointerList* plist)
 {
 	for (size_t i=0;i<plist->size();++i)
 		if ((*plist)[i]!=0)
 			delete (*plist)[i];
 };


//removes almost zero entries from xarr and writes result to out.
template<class UnarySparse>
floatType  FilterSparseArray(const UnarySparse& xarr,UnarySparse* pout, floatType precision=floatTypePrecision)
 {
	UnarySparse& out=*pout;
 	typename UnarySparse::const_iterator beg=xarr.begin(),end=xarr.end();

 	floatType eps=precision;

 	typename UnarySparse::const_iterator it=std::find_if(beg,end,std::bind2nd(std::less_equal<floatType>(),eps));
 	if (it==end)
 	{
 		out=xarr;
 	}
 	else
 	{
 		out.resize(xarr.size(),false); out.clear();
 		beg=xarr.begin(), end=xarr.end();
 		for (;beg!=end;++beg)
 			if (*beg>eps)
 				out[beg.index()]=*beg;
 	}

 	return _Normalize(out.begin(),out.end());
 }

};//DD

namespace OUT{
	using std::ostream;
	using std::endl;
	using std::vector;

	template<typename ArrayType>
	ostream& operator << (ostream& logger,const vector<ArrayType>& arr)
	{
		for (size_t i=0;i<arr.size();++i)
			logger << arr[i]<<"; ";
		logger <<endl;
		return logger;
	};

	template<typename Type1,typename Type2>
	ostream& operator << (ostream& logger, const std::pair<Type1,Type2>& p)
	{
		logger <<"("<<p.first<<","<<p.second<<")";
		return logger;
	};

	template<class Iterator>
	void saveContainer(ostream& fout, Iterator begin, Iterator end)
	{
		for ( ; begin!=end; ++begin)
			fout <<std::scientific<<*begin<<"; ";
		fout << std::endl;
	}

};//OUT



#endif /* SMALLOBJECTS_H_ */
