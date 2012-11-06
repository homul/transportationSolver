#ifndef PRIMALSOLVER_H_
#define PRIMALSOLVER_H_
#include <iostream>
#include <numeric>
#include <utility>
#include <queue>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "common.h"
#include "smallobjects.h"
#include "list2D.h"

namespace DD
{
using std::ostream;
using std::cout;
using std::endl;
using std::pair;
using namespace OUT;

/*
* Use maximize = true if use want to solve maximization problem and false for a minimization one
*
* DenseMatrix represents a dense matrix type and has to
*  - contain elements of the type floatType, defined in common.h
*  and provide
*  - copy constructor
*  - floatType operator (size_t index_a, size_t index_b) to access its elements
*
*  UnarySparse represents a sparse array type and has to
*   - contain elements of the type floatType, defined in common.h
*   and provide operations
*   - size_t nnz() const returning number of non-zero elements in the vector;
*   - const_iterator begin() const
*   - iterator begin()
*   - const_iterator end() const
*   - iterator end() const, addressing non-zero elements ONLY such, that iterator.index()
*   provides an index of the elemt it addresses.
*   size_t size() const returning a maximum number of elements = the largest element index-1
*
*   All these properties are fulfilled for
*    boost::numeric::ublas;
*    typedef boost::numeric::ublas::compressed_vector<DD::floatType> SparseArray;
*    typedef boost::numeric::ublas::matrix<DD::floatType> Matrix;
*   see tests/testcommon.h
**
**/
template<bool maximize, class DenseMatrix,class UnarySparse>
class PrimalSolver
{
public:
	typedef enum{X, Y} Direction;
	typedef pair<size_t,Direction> CoordDir;
	typedef std::queue<CoordDir> Queue;
	typedef List2D<floatType> FeasiblePoint;
	typedef std::vector<floatType> UnaryDense;
	typedef std::vector<size_t> IndexArray;
	typedef std::list<typename FeasiblePoint::const_iterator> CycleList;

	PrimalSolver(ostream& fout=std::cout,floatType relativePrecision=floatTypeEps):
		_fout(fout),_pbin(0),_xsize(0),_ysize(0),_basicSolution(0,0,0),_relativePrecision(relativePrecision)
	{};

	PrimalSolver(const size_t& xsize,const size_t& ysize,const DenseMatrix& bin,ostream& fout=std::cout,floatType relativePrecision=floatTypeEps)
	: _fout(fout),_pbin(&bin),_xsize(xsize),_ysize(ysize),_basicSolution(xsize,ysize,_nnz(xsize,ysize))
	{
		Init(xsize,ysize,bin,relativePrecision);
	};

	void Init(size_t xsize,size_t ysize,const DenseMatrix& bin,floatType relativePrecision=floatTypeEps);
	std::pair<bool,floatType> Solve(const UnarySparse& xarr,const UnarySparse& yarr);//!< non-necessary near-zero elements will be removed when needed
	floatType GetPrimalValue()const{return _basicSolution.inner_product2D(*_pbin, &_primalValueNumericalPrecision);};//!< returns value of the current basic solution
	floatType GetSolution(DenseMatrix* pbin)const;
	void PrintTestData(std::ostream& fout)const;
	void PrintProblemDescription(const UnarySparse& xarr,const UnarySparse& yarr);

protected:
	bool isRestartPossible(const UnarySparse& xarr,const UnarySparse& yarr);
	bool CorrectBasicSolution(UnarySparse xarr,UnarySparse yarr);
	bool isSolutionRecalculated()const{return _recalculated;}
private:
	better<maximize,floatType> _1stbetter;
	notWorse<maximize,floatType> _1stnotWorse;

	bool _InitBasicSolution(const UnarySparse& xarr,const UnarySparse& yarr);
	bool _isOptimal(pair<size_t,size_t>* pmove);
	bool _CheckDualConstraints(const UnaryDense& xdual,const UnaryDense& ydual,pair<size_t,size_t>* pmove )const;
	CoordDir _findSingleNeighborNode(const FeasiblePoint&)const;
	void _BuildDuals(UnaryDense* pxdual,UnaryDense* pydual);
	void _FindCycle(FeasiblePoint* pfp,const pair<size_t,size_t>& move);
	void _ChangeSolution(const FeasiblePoint& fp,const pair<size_t,size_t>& move);
	bool _MovePotentials(const pair<size_t,size_t>& move);

	void _move2Neighbor(const FeasiblePoint& fp,FeasiblePoint::const_iterator &it)const;//helper function - determines, iterator direction and moves it to another element of a list fp with length 2

	static size_t _nnz(size_t xsize,size_t ysize){return xsize+ysize;}
	static floatType _FilterBound(const UnarySparse& xarr,UnarySparse& out,floatType precision);
	void _checkCounter(size_t* pcounter,const char* errmess);

	mutable floatType _primalValueNumericalPrecision;
	bool _recalculated;
	ostream& _fout;
	const DenseMatrix* _pbin;
	size_t _xsize,_ysize;
	floatType _relativePrecision;//relative precision of thresholding input marginal values
	FeasiblePoint _basicSolution;
	IndexArray _activeXbound;
	IndexArray _activeYbound;
};

template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::Init(size_t xsize,size_t ysize,const DenseMatrix& bin,floatType relativePrecision)
{
	_pbin=&bin;
	_xsize=xsize;
	_ysize=ysize;
	_basicSolution.resize(xsize,ysize,_nnz(xsize,ysize));
	_activeXbound.clear();
	_activeYbound.clear();
	_primalValueNumericalPrecision=0;
	_relativePrecision=relativePrecision;
};

template<bool maximize,class DenseMatrix,class UnarySparse>
bool PrimalSolver<maximize,DenseMatrix,UnarySparse>::isRestartPossible(const UnarySparse& xarr,const UnarySparse& yarr)
{
	if (xarr.nnz()!=_activeXbound.size())
	 return false;

	if (yarr.nnz()!=_activeYbound.size())
		return false;

	if ( (_activeXbound.size()==_xsize) && (_activeYbound.size()==_ysize))
		return true;

	size_t i=0;
	for (typename UnarySparse::const_iterator begin=xarr.begin();begin!=xarr.end();++begin)
	  if (begin.index()!=_activeXbound[i++])
		 return false;

	i=0;
	for (typename UnarySparse::const_iterator begin=yarr.begin();begin!=yarr.end();++begin)
		if (begin.index()!=_activeYbound[i++])
			 return false;

	return true;
};

template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_checkCounter (size_t* pcounter, const char* errmess)
{
	if ((*pcounter)++ < 1000)
		return;

	throw std::runtime_error(errmess);
};


template<bool maximize,class DenseMatrix,class UnarySparse>
bool PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_InitBasicSolution(const UnarySparse& xarr,const UnarySparse& yarr)
{
	UnarySparse row=xarr;
	UnarySparse col=yarr;
	//north-west corner basic solution
	_basicSolution.clear();
	typename UnarySparse::iterator rbeg=row.begin(), rend=row.end();
	typename UnarySparse::iterator cbeg=col.begin(), cend=col.end();

	size_t counter=0;
	while ((rbeg!=rend)&&(cbeg!=cend))
	{
	if (*cbeg>=*rbeg)
	{
		_basicSolution.push(rbeg.index(), cbeg.index(),*rbeg);
		(*cbeg)-=(*rbeg);
		if (rbeg!=rend)
		 ++rbeg;
		else
		 ++cbeg;
	}
	else
	{
		_basicSolution.push(rbeg.index(),cbeg.index(),*cbeg);
		(*rbeg)-=(*cbeg);
		if (cbeg!=cend)
		 ++cbeg;
		else
		 ++rbeg;
	}

	_checkCounter(&counter,"_InitBasicSolution-inifinite loop!\n");
	}

	size_t basicNum=xarr.nnz()+yarr.nnz()-1;
	if (counter==basicNum)
		return true;
	else
		return false;
};

/*
 * returns coordinate + direction of the point which is alone in its row/column.
 * e.g. (1,X) means that the column with X-coordinate equal to 1, contains a single element.
 */
template<bool maximize,class DenseMatrix,class UnarySparse>
typename PrimalSolver<maximize,DenseMatrix,UnarySparse>::CoordDir
PrimalSolver<maximize,DenseMatrix,UnarySparse>::_findSingleNeighborNode(const FeasiblePoint& fp)const
{
	IndexArray::const_iterator beg,end=_activeXbound.end();
	for (beg=_activeXbound.begin();beg!=end;++beg)
		if (fp.colSize(*beg)==1)
			return std::make_pair(*beg,X);

	end=_activeYbound.end();
	for (beg=_activeYbound.begin();beg!=end;++beg)
		if (fp.rowSize(*beg)==1)
			return std::make_pair(*beg,Y);

	return std::make_pair(MAXSIZE_T,X);
};

template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_BuildDuals(UnaryDense* pxdual,UnaryDense* pydual)
{
	UnaryDense& xdual=*pxdual;
	UnaryDense& ydual=*pydual;

	xdual.assign(_xsize,0.0);
	ydual.assign(_ysize,0.0);

	FeasiblePoint fpcopy(_basicSolution);
	CoordDir currNode=_findSingleNeighborNode(fpcopy);

	if (currNode.first==MAXSIZE_T)
		throw std::runtime_error("_BuildDuals: can not build duals: no single neighbor node available!");

	if (currNode.second==X)
	{
		currNode.second=Y;
		currNode.first=fpcopy.colBegin(currNode.first).coordinate();
	}else
	{
		currNode.second=X;
		currNode.first=fpcopy.rowBegin(currNode.first).coordinate();
	}

	Queue qu;
	qu.push(currNode);

	size_t counter=0;
	do
	{
		if (qu.front().second==Y)
		{
		 size_t y=qu.front().first;

		 typename FeasiblePoint::iterator beg=fpcopy.rowBegin(y),
										  end=fpcopy.rowEnd(y);
		 for (;beg!=end;++beg)
		 {
			 size_t x=beg.coordinate();
			 xdual[x]=(*_pbin)(x,y)-ydual[y];
			 qu.push(std::make_pair(x,X));
		 }
		 fpcopy.rowErase(y);

		}else
		{
			size_t x=qu.front().first;

			 typename FeasiblePoint::iterator beg=fpcopy.colBegin(x),
											  end=fpcopy.colEnd(x);
			 for (;beg!=end;++beg)
			 {
				 size_t y=beg.coordinate();
				 ydual[y]=(*_pbin)(x,y)-xdual[x];
				 qu.push(std::make_pair(y,Y));
			 }

			 fpcopy.colErase(x);
		}

		qu.pop();

	_checkCounter(&counter, "_BuildDuals-infinite loop!\n");
	}while (!qu.empty());

};


template<bool maximize,class DenseMatrix,class UnarySparse>
 bool PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_CheckDualConstraints(const UnaryDense& xdual,const UnaryDense& ydual,pair<size_t,size_t>* pmove )const
{
	floatType eps=(maximize ? 1.0 : -1.0)*floatTypeEps;
	floatType delta=0, precision=0;

	IndexArray::const_iterator ybeg=_activeYbound.begin(), yend=_activeYbound.end();
	for (;ybeg<yend;++ybeg)
	{
		IndexArray::const_iterator xbeg=_activeXbound.begin(), xend=_activeXbound.end();
		for (;xbeg<xend;++xbeg)
		{
			size_t x=*xbeg, y=*ybeg;
			delta=(*_pbin)(x,y)-xdual[x]-ydual[y];

			precision=(fabs((*_pbin)(x,y))+fabs(xdual[x])+fabs(ydual[y]))*eps;
			if (_1stbetter(delta,precision))
			 {
				pmove->first=x; pmove->second=y;
				return false;
			 }
		}
	}

		return true;
};


template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_FindCycle(FeasiblePoint* pfp,const pair<size_t,size_t>& move)
{
	FeasiblePoint& fp=*pfp;

	fp.insert(move.first,move.second,0);

	CoordDir cd=_findSingleNeighborNode(fp);

	size_t counter=0;//_initCounter();
	while (cd.first<MAXSIZE_T)
	{
		if (cd.second==X)
			fp.colErase(cd.first);
		else
			fp.rowErase(cd.first);

		cd=_findSingleNeighborNode(fp);

		_checkCounter(&counter,"_FindCycle-infinite loop!\n");
	}
};

//helper function - determines, iterator direction and moves it to another element of a list fp with length 2
template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_move2Neighbor(const FeasiblePoint& fp,FeasiblePoint::const_iterator &it)const
{
 FeasiblePoint::const_iterator beg=fp.rowBegin(0);
 if (it.isRowIterator())
 {
	 beg=fp.rowBegin(it.y());
 }
 else
 {
	 beg=fp.colBegin(it.x());
 }

if (beg==it)
	++it;
else
	--it;
};

template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_ChangeSolution(const FeasiblePoint& fp,const pair<size_t,size_t>& move)
{
	size_t y=0;
	for (;y<fp.ysize();++y)
	{
	  assert( (fp.rowSize(y)==2) || (fp.rowSize(y)==0) ) ;
	  if (fp.rowSize(y)!=0)
		break;
	}

	CycleList plusList, minusList;
	CycleList* pplus=&plusList, *pPlusList=0;
	CycleList* pminus=&minusList, *pMinusList=0;

	//going along the cycle to assign +/- to vertices correctly
	FeasiblePoint::const_iterator it=fp.rowBegin(y);
	pair<size_t,size_t> c0(it.x(),it.y());
	do{
		pplus->push_back(it);
		if ( (it.x()==move.first) && (it.y()==move.second) )
		{
			pMinusList=pminus; //really minus list is in *pMinusList
			pPlusList =pplus;
		}
		it=it.changeDir();
		_move2Neighbor(fp,it);
		swap(pplus,pminus);
	}while (! ( (it.x()==c0.first) && (it.y()==c0.second) ) );

	assert (pMinusList!=0);
	assert (pPlusList!=0);

	//selecting the smallest number to add...
		 floatType min=floatTypeMax;
		 CycleList::const_iterator iterMinVal, beg=pMinusList->begin(), end=pMinusList->end();
		 for (;beg!=end;++beg)
		 {
			 if ( (  **beg < min) ||
				  ( (**beg == min) &&
					    ( (beg->y() < iterMinVal->y()) ||
					    ( ((beg->y() == iterMinVal->y())
					    		 && (beg->x() < iterMinVal->x())) ) ) 	) )
			 {
				 min=**beg;
				 iterMinVal=beg;
			 }
		 }

		  //changing
		 _basicSolution.insert(move.first,move.second,0);

		  beg=pMinusList->begin(), end=pMinusList->end();
		  for (;beg!=end;++beg)
		 	 _basicSolution.buffer((*beg).index())-=min;

		  beg=pPlusList->begin(), end=pPlusList->end();
		  for (;beg!=end;++beg)
		 	 _basicSolution.buffer((*beg).index())+=min;

		  _basicSolution.erase((*iterMinVal).index());
}

template<bool maximize,class DenseMatrix,class UnarySparse>
bool PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_MovePotentials(const pair<size_t,size_t>& move)
{
	floatType ObjVal=GetPrimalValue();
	floatType primalValueNumericalPrecisionOld=_primalValueNumericalPrecision;
	FeasiblePoint fp=_basicSolution;

	_FindCycle(&fp,move);

	_ChangeSolution(fp,move);


        floatType newObjValue=GetPrimalValue();
	if ( (_1stbetter(ObjVal,newObjValue)) && (fabs(ObjVal-newObjValue)>(_primalValueNumericalPrecision+primalValueNumericalPrecisionOld) ) )
		{
			std::cerr<<_fout<<std::setprecision (std::numeric_limits<floatType>::digits10+1) << std::endl<<"ObjVal="<<ObjVal
					 <<", newObjValue="<<newObjValue
					 <<", fabs(ObjVal-newObjValue)="<<fabs(ObjVal-newObjValue)<<", _primalValueNumericalPrecision="<<_primalValueNumericalPrecision
					 << ", primalValueNumericalPrecisionOld="<< primalValueNumericalPrecisionOld <<endl;
			_fout << "Basic solution before move:" <<std::endl;
			fp.PrintTestData(_fout);
			_fout << "Move:" << move<<std::endl;

			return false;
		}

	return true;
};


template<bool maximize,class DenseMatrix,class UnarySparse>
bool PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_isOptimal(pair<size_t,size_t>* pmove)
{
	//checks current basic solution for optimality
	//1. build duals
	UnaryDense xduals,yduals;
	_BuildDuals(&xduals,&yduals);
	//2. check whether they satisfy dual constraints
	return _CheckDualConstraints(xduals,yduals,pmove);
};

template<bool maximize,class DenseMatrix,class UnarySparse>
floatType PrimalSolver<maximize,DenseMatrix,UnarySparse>::
_FilterBound(const UnarySparse& xarr,UnarySparse& out,floatType precision)
{
 	if (xarr.nnz()==0)
 		throw std::runtime_error("PrimalSolver:_FilterBound(): empty input array <xarr>!");

 	return FilterSparseArray(xarr,&out,precision);
}

template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::
PrintProblemDescription(const UnarySparse& xarr,const UnarySparse& yarr)
{
		 size_t maxprecision=std::numeric_limits<floatType>::digits10;
				_fout<< std::setprecision (maxprecision+1) << "xarr=" << xarr<< " nnz: ";
				saveContainer( _fout, xarr.begin(), xarr.end() );
				_fout<< std::setprecision (maxprecision+1) << "yarr=" << yarr <<" nnz: ";
				saveContainer( _fout, yarr.begin(), yarr.end() );
				for (typename UnarySparse::const_iterator xbeg=xarr.begin();xbeg!=xarr.end();++xbeg)
					for (typename UnarySparse::const_iterator ybeg=yarr.begin();ybeg!=yarr.end();++ybeg)
					_fout << std::setprecision (maxprecision+1)<<"; bin("<<xbeg.index()<<","<<ybeg.index()<<")="<<(*_pbin)(xbeg.index(),ybeg.index());

		_fout <<std::endl<< "Current basic solution:"<<endl;
		_basicSolution.PrintTestData(_fout);
}

template<bool maximize,class DenseMatrix,class UnarySparse>
std::pair<bool,floatType> PrimalSolver<maximize,DenseMatrix,UnarySparse>::
 Solve(const UnarySparse& xarr_ini,const UnarySparse& yarr_ini)
{
	_recalculated=false;
	UnarySparse xarr,yarr;
	_FilterBound(xarr_ini,xarr,_relativePrecision*_xsize*_ysize);
	_FilterBound(yarr_ini,yarr,_relativePrecision*_xsize*_ysize);

	if ((xarr.nnz()==0) || (yarr.nnz()==0))
		throw std::runtime_error("PrimalSolver::Solve:Empty input arrays. Is the _relativePrecision properly assigned?");

	if (isRestartPossible(xarr,yarr))
	{
		/*
		 * tries to create a non-negative solution with the same entries as in _basicSolution
		 * If it works, then this solution is automatically optimal
		 */
		if (CorrectBasicSolution(xarr,yarr))
		{
			_recalculated=true;
			return std::make_pair(true, GetPrimalValue());
		}
	}else
		{
		 _activeXbound.resize(xarr.nnz());
		 _activeYbound.resize(yarr.nnz());
		 copyNZindices(xarr.begin(),xarr.end(),_activeXbound.begin());
		 copyNZindices(yarr.begin(),yarr.end(),_activeYbound.begin());
		};

	//1. Create basic solution _basicSolution
	if (!_InitBasicSolution(xarr,yarr))
		return std::make_pair(false, GetPrimalValue());

	//2. Check optimality
	pair<size_t,size_t> move;
	bool objectiveImprovementFlag=true;
	size_t counter=0;//_initCounter();

	while ((objectiveImprovementFlag)&&(!_isOptimal(&move)))
	{
		objectiveImprovementFlag=_MovePotentials(move); //changes basic solution
		_checkCounter(&counter,"Solve-infinite loop!\n");
		if (!objectiveImprovementFlag)
		{
			PrintProblemDescription(xarr,yarr);
			throw std::runtime_error("PrimalSolver::Solve: Objective has become worse. Interrupting!");
		}
	}
	return std::make_pair(true, GetPrimalValue());
};
template<bool maximize,class DenseMatrix,class UnarySparse>
bool PrimalSolver<maximize,DenseMatrix,UnarySparse>::
CorrectBasicSolution(UnarySparse xarr,UnarySparse yarr)
{
	FeasiblePoint fp=_basicSolution;

	CoordDir cd=_findSingleNeighborNode(fp);

	while (cd.first<MAXSIZE_T)
	{
		if (cd.second==X)
		{
			FeasiblePoint::iterator setIt=fp.colBegin(cd.first);
			floatType val=_basicSolution.buffer(setIt.index())=xarr[cd.first];
			yarr[setIt.y()]-=val;
			if (yarr[setIt.y()]<0)
				return false; //!< can not create a basic solution
			fp.colErase(cd.first);
		}
		else
		{
			FeasiblePoint::iterator setIt=fp.rowBegin(cd.first);
			floatType val=_basicSolution.buffer(setIt.index())=yarr[cd.first];
			xarr[setIt.x()]-=val;
			if (xarr[setIt.x()]<0)
				return false; //!< can not create a basic solution
			fp.rowErase(cd.first);
		}

		cd=_findSingleNeighborNode(fp);
	}
	return true;
};


template<bool maximize,class DenseMatrix,class UnarySparse>
floatType PrimalSolver<maximize,DenseMatrix,UnarySparse>::GetSolution(DenseMatrix* pbin)const
{
	_basicSolution.get2DTable(pbin);
	return GetPrimalValue();
};

template<bool maximize,class DenseMatrix,class UnarySparse>
void PrimalSolver<maximize,DenseMatrix,UnarySparse>::PrintTestData(std::ostream& fout)const
{
	fout << "_relativePrecision="<<_relativePrecision<<endl;
	fout << "_xsize="<<_xsize<<", _ysize="<<_ysize<<endl;
	fout <<"_basicSolution:"<<endl;
	_basicSolution.PrintTestData(fout);
	fout <<endl<< "_activeXbound: "<<_activeXbound;
	fout <<endl<< "_activeYbound: "<<_activeYbound;
};

};

#endif /* PRIMALSOLVER_H_ */
