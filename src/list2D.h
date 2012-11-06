#ifndef LIST2D_H_
#define LIST2D_H_
#include "common.h"
#include <list>
#include <utility>
#include "smallobjects.h"
namespace DD{

using std::list;
using std::ostream;
using std::endl;
using namespace OUT;

template<class T>
class List2D
{
public:

	struct bufferElement;

	struct listElement
	{
		listElement(size_t coordinate,bufferElement* pbufElement):
			_coordinate(coordinate), _pbufElement(pbufElement)
		{};

		size_t _coordinate;
		bufferElement* _pbufElement;
	};

	typedef std::list<listElement> List1D;

	template<class Parent,class typeT>
		class iterator_template : public Parent
		{
			public:
			iterator_template(Parent it,const bufferElement* pbuffer0):Parent(it),_pbuffer0(pbuffer0){};
			typeT& operator * ()const{return this->Parent::operator *()._pbufElement->_value;}

			size_t index()const
			{
				return (this->Parent::operator *()._pbufElement)-_pbuffer0;
			}

			size_t coordinate()const{return this->Parent::operator *()._coordinate;}
			size_t x()const{return (*this->Parent::operator *()._pbufElement->_rowIterator)._coordinate;}
			size_t y()const{return (*this->Parent::operator *()._pbufElement->_colIterator)._coordinate;}

			bool isRowIterator()const{return &(*(this->Parent::operator *()._pbufElement->_rowIterator)) == &(*Parent(*this));}

			iterator_template changeDir()const{
			    if (isRowIterator())
					return iterator_template(this->Parent::operator *()._pbufElement->_colIterator,_pbuffer0);
				else
					return iterator_template(this->Parent::operator *()._pbufElement->_rowIterator,_pbuffer0);
			}

			iterator_template operator ++ (int){iterator_template it=*this; ++(*this); return it;}
			iterator_template& operator ++ (){Parent::operator ++(); return *this;}
			iterator_template& operator -- (){Parent::operator --(); return *this;}
			private:
			const bufferElement* _pbuffer0;
		};

		typedef iterator_template<typename List1D::iterator,T> iterator;
		typedef iterator_template<typename List1D::const_iterator,const T> const_iterator;

	typedef std::vector<List1D> List1DSeq;

	struct bufferElement
	{
		bufferElement(const T& val,typename List1D::iterator rowIterator,typename List1D::iterator colIterator)
		    :_value(val) {
		    if (_value != NaN()) {
			_rowIterator = rowIterator;
			_colIterator = colIterator;
		    }
		};

		bufferElement(const bufferElement &other)
		: _value(other._value)
		{
		 if (_value != NaN()) {
			 _rowIterator = other._rowIterator;
			 _colIterator = other._colIterator;
		 }
		}

		bufferElement & operator=(const bufferElement &other) {
			_value = other._value;
			 if (_value != NaN()) {
				 _rowIterator = other._rowIterator;
				 _colIterator = other._colIterator;
			 }
			 return *this;
		}

		T _value;
		typename List1D::iterator _rowIterator;
		typename List1D::iterator _colIterator;
	};

	typedef std::vector<bufferElement> Buffer;

	List2D(size_t xsize, size_t ysize, size_t nnz);
	List2D(const List2D&);
	List2D& operator = (const List2D&);

	void clear();
	/*
	 * after resizing all data is lost!
	 */
	void resize(size_t xsize, size_t ysize, size_t nnz);

	/* tries to insert to the end of lists.
	 * compares coordinates of the last element for this purpose
	 * it it does not work, returns <false>
	 */
	bool push(size_t x, size_t y, const T& val);

	/*
	 * tries to insert to the position of the last
	 * call of the function erase(). If it was not called yet -
	 * the last position in the allocated memory.
	 * If the insertion position is occupied,
	 * returns <false>
	 */
	bool insert(size_t x, size_t y, const T& val);
	void erase(iterator it);
	/*
	 * index - index in the _buffer array
	 */
	void erase(size_t index){erase(iterator(_buffer[index]._rowIterator,&_buffer[0]));}

	void rowErase(size_t y);
	void colErase(size_t x);

	size_t rowSize(size_t y)const{return _rowLists[y].size();};
	size_t xsize()const{return _colLists.size();}
	size_t colSize(size_t x)const{return _colLists[x].size();};
	size_t ysize()const{return _rowLists.size();}
	size_t nnz()const{return _buffer.size();}

	iterator rowBegin(size_t y){return iterator(_rowLists[y].begin(),&_buffer[0]);}
	const_iterator rowBegin(size_t y)const{return const_iterator(_rowLists[y].begin(),&_buffer[0]);}

	iterator rowEnd(size_t y){return iterator(_rowLists[y].end(),&_buffer[0]);}
	const_iterator rowEnd(size_t y)const{return const_iterator(_rowLists[y].end(),&_buffer[0]);}

	iterator colBegin(size_t x){return iterator(_colLists[x].begin(),&_buffer[0]);}
	const_iterator colBegin(size_t x)const{return const_iterator(_colLists[x].begin(),&_buffer[0]);}

	iterator colEnd(size_t x){return iterator(_colLists[x].end(),&_buffer[0]);}
	const_iterator colEnd(size_t x)const{return const_iterator(_colLists[x].end(),&_buffer[0]);}

	//iterator switchDirection(iterator it)const;
	template<class BinaryTable1D>
	T inner_product1D(const BinaryTable1D& bin)const;

	//pprecision - if non-zero contains an upper bound for the numerical precision of the returned  value
	template<class BinaryTable2D>
	T inner_product2D(const BinaryTable2D& bin, T* pprecision=0)const;

	template<class BinaryTable2D>
	void get2DTable(BinaryTable2D* pbin)const;

	T& buffer(size_t index){return _buffer[index]._value;}
	const T& buffer(size_t index)const{return _buffer[index]._value;}

	std::pair<bool,T> getValue(size_t x,size_t y)const;//!< not very efficient function. Implemented mainly for test purposes.
	void PrintTestData(std::ostream& fout)const;
private:
	bool _insert(size_t x, size_t y, const T& val, size_t position);
	void _copy(const List2D<T>& lst);
	static T NaN(){return std::numeric_limits<T>::max();}

	//size_t _nnz;
	size_t _insertPosition;
	size_t _pushPosition;
	List1DSeq _rowLists;
	List1DSeq _colLists;
	Buffer _buffer;
};

template<class T>
List2D<T>::List2D(size_t xsize, size_t ysize, size_t nnz):
_insertPosition(nnz-1),
_pushPosition(0),
_rowLists(ysize),
_colLists(xsize),
_buffer(nnz,bufferElement(NaN(),typename List1D::iterator(),typename List1D::iterator()))
{};

template<class T>
List2D<T>::List2D(const List2D& lst)
{
	_copy(lst);
}

template<class T>
void List2D<T>::resize(size_t xsize, size_t ysize, size_t nnz)
{
	_rowLists.assign(ysize,List1D());
	_colLists.assign(xsize,List1D());
	_buffer.assign(nnz,bufferElement(NaN(),typename List1D::iterator(),typename List1D::iterator()));
	_insertPosition=nnz-1;
	_pushPosition=0;
};


template<class T>
void List2D<T>::_copy(const List2D<T>& lst)
{
	_buffer=lst._buffer;

	_rowLists=lst._rowLists;
	typename List1DSeq::iterator itbeg=_rowLists.begin(), itend=_rowLists.end();
	for (;itbeg!=itend;++itbeg)
	{
		typename List1D::iterator beg=(*itbeg).begin(),end=(*itbeg).end();
		for (;beg!=end;++beg)
		{
			size_t offset=(*beg)._pbufElement-&(lst._buffer[0]);
			(*beg)._pbufElement= &_buffer[offset];
			(*beg)._pbufElement->_rowIterator=beg;
		}
	}

	_colLists=lst._colLists;
	itbeg=_colLists.begin(), itend=_colLists.end();
	for (;itbeg!=itend;++itbeg)
	{
		typename List1D::iterator beg=(*itbeg).begin(),end=(*itbeg).end();
		for (;beg!=end;++beg)
		{
			size_t offset=(*beg)._pbufElement-&(lst._buffer[0]);
			(*beg)._pbufElement= &_buffer[offset];
			(*beg)._pbufElement->_colIterator=beg;
		}
	}

	//_nnz=lst._nnz;
	_insertPosition=lst._insertPosition;
	_pushPosition=lst._pushPosition;
};

template<class T>
List2D<T>& List2D<T>::operator = (const List2D<T>& lst)
{
	if (this==&lst)
		return *this;

	_copy(lst);

	return *this;
}

template<class T>
bool List2D<T>::insert(size_t x,size_t y,const T& val)
{
	if (_insert(x,y,val,_insertPosition))
	{
	  _insertPosition=_buffer.size();
	  return true;
	}

	return false;
};

template<class T>
bool List2D<T>::push(size_t x, size_t y, const T& val)
{
	if (_insert(x,y,val,_pushPosition))
	{
	  ++_pushPosition;
	  //the very last position in _buffer can not be occupied de to push(), only due to insert()
	  if (_pushPosition == (_buffer.size()-1))
		  ++_pushPosition;
	  return true;
	}

	return false;
}

template<class E>
class coordLess
{
public:
 coordLess(size_t x):_x(x){}
 bool operator () (const E& e) const{return e._coordinate < _x;}
private:
 size_t _x;
};

template<class E>
class coordMore
{
public:
 coordMore(size_t x):_x(x){}
 bool operator () (const E& e) const{return e._coordinate > _x;}
private:
 size_t _x;
};

template<class T>
bool List2D<T>::_insert(size_t x, size_t y, const T& val, size_t position)
{
	assert(x<_colLists.size());
	assert(y< _rowLists.size());

	if (position >= _buffer.size())
		return false;

	bufferElement& buf=_buffer[position];
	buf._value=val;

	List1D& rowList=_rowLists[y];
	List1D& colList=_colLists[x];

	typename List1D::iterator insertPosition=std::find_if(rowList.begin(),rowList.end(),coordMore<listElement>(x));
	buf._rowIterator=rowList.insert(insertPosition,listElement(x,&buf));
	insertPosition=std::find_if(colList.begin(),colList.end(),coordMore<listElement>(y));
	buf._colIterator=colList.insert(insertPosition,listElement(y,&buf));

	return true;
};

template<class T>
void List2D<T>::erase(iterator it)
{
 _insertPosition=it.index();
 size_t x=it.x(), y=it.y();
 _rowLists[y].erase(_buffer[_insertPosition]._rowIterator);
 _colLists[x].erase(_buffer[_insertPosition]._colIterator);
 //XXX
 _buffer[_insertPosition]._value=NaN();
};

template<class T>
void List2D<T>::rowErase(size_t y)
{
	while (!_rowLists[y].empty())
		erase(iterator(_rowLists[y].begin(),&_buffer[0]));
};

template<class T>
void List2D<T>::colErase(size_t x)
{
	while (!_colLists[x].empty())
		erase(iterator(_colLists[x].begin(),&_buffer[0]));
};

template<class T>
void List2D<T>::clear()
{
	for (size_t x=0;x<_rowLists.size();++x)
		rowErase(x);

	for (size_t y=0;y<_colLists.size();++y)
		colErase(y);

	_pushPosition=0;
	_insertPosition=_buffer.size()-1;
};

template<class T>
template<class BinaryTable1D>
T List2D<T>::inner_product1D(const BinaryTable1D& bin)const
{
	T sum=0;
	for (size_t i=0; i<_colLists.size();++i)
	{
		typename List1D::const_iterator beg=_colLists[i].begin(), end=_colLists[i].end();
		for (;beg!=end;++beg)
			sum+=(*beg)._pbufElement->_value * bin(xsize()*((*beg)._coordinate)+i);
	};
	return sum;
};

template<class T>
template<class BinaryTable2D>
T List2D<T>::inner_product2D(const BinaryTable2D& bin, T* pprecision)const //DEBUG
{
	T precision_;
	T* pprecision_;
	if (pprecision!=0)
		pprecision_=pprecision;
	else
		pprecision_=&precision_;

	*pprecision_=0;

	T sum=0;
	for (size_t i=0; i<xsize();++i)
	{
		const_iterator beg=colBegin(i), end=colEnd(i);
		for (;beg!=end;++beg)
		{
			sum+=(*beg) * bin(beg.x(),beg.y());
			*pprecision_+=floatTypeEps*fabs(sum);
		}
	};
	return sum;
};

template<class T>
std::pair<bool,T> List2D<T>::getValue(size_t x,size_t y)const
{
 typename List1D::const_iterator beg=_colLists[x].begin(), end=_colLists[x].end();
 for (;beg!=end;++beg)
	if ((*beg)._coordinate==y)
		return std::make_pair(true,(*beg)._pbufElement->_value);

  return std::make_pair(false,(T)0);
};

template<class T>
void List2D<T>::PrintTestData(std::ostream& fout)const
{
	fout << "_nnz=" <<_buffer.size()<<endl;
	fout << "_insertPosition=" << _insertPosition<<endl;
	fout << "_pushPosition=" << _pushPosition<<endl;
	fout << "xsize="<<_colLists.size()<<endl;
	fout << "ysize="<<_rowLists.size()<<endl;

	std::vector<T> printBuffer(_buffer.size(),NaN());

	fout << "row Lists: "<<endl;
	for (size_t i=0; i< _rowLists.size();++i)
	{
		fout << "y="<<i<<": ";
		typename List1D::const_iterator beg=_rowLists[i].begin(), end=_rowLists[i].end();
		for (;beg!=end;++beg)
		{
			fout <<"("<<(*beg)._coordinate<<","<<(*beg)._pbufElement->_value<<")";
			printBuffer[(*beg)._pbufElement-&_buffer[0]]=(*beg)._pbufElement->_value;
		}
		fout <<endl;
	}

	fout << "column Lists: "<<endl;
	for (size_t i=0; i< _colLists.size();++i)
	{
		fout << "x="<<i<<": ";
		typename List1D::const_iterator beg=_colLists[i].begin(), end=_colLists[i].end();
		for (;beg!=end;++beg)
		{
			fout <<"("<<(*beg)._coordinate<<","<<(*beg)._pbufElement->_value<<")";
		}
		fout <<endl;
	}

	fout << "buffer: ";
	for (size_t i=0;i<printBuffer.size();++i)
		if (printBuffer[i]!=NaN())
		 fout << "("<<_buffer[i]._value<<","<<(*_buffer[i]._rowIterator)._coordinate <<","<< (*_buffer[i]._colIterator)._coordinate<<")";
		else
			fout << "(nan,nan,nan)";
	fout << endl;

};

template<class T>
template<class BinaryTable2D>
void List2D<T>::get2DTable(BinaryTable2D* pbin)const
{
	for (size_t x=0;x<xsize();++x)
		for (size_t y=0;y<ysize();++y)
			(*pbin)(x,y)=0;

	for (size_t i=0; i<xsize();++i)
	{
		const_iterator beg=colBegin(i), end=colEnd(i);
		for (;beg!=end;++beg)
			(*pbin)(beg.x(),beg.y())=(*beg);
	};
};

};

#endif /* LIST2D_H_ */
