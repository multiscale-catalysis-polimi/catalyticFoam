/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "math/OpenSMOKEVector.h"
#include "math/OpenSMOKEMatrix.h"

#if OPENSMOKE_USE_MKL == 1
	#include "mkl.h"
#elif OPENSMOKE_USE_OPENBLAS == 1
	#include "lapacke.h"
#endif

namespace OpenSMOKE
{

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(void)
	{
		SetCounters();
		Initialize(0);
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(OpenSMOKEVector<T> const& rval)
	{
		SetCounters();
		Initialize(rval.dimensions_);
		Copy(rval);
	/*
		if(rval.shadow_ == false)
		{
			Initialize(rval.dimensions_);
			if(dimensions_ != 0)
				memcpy(vector_, rval.vector_, (dimensions_+1)*sizeof(T));
		}
		else
		{
			Initialize(0);
			OpenSMOKEVector<T> aux; 
			aux = rval;
			Swap(this, &aux);
			count_--;
			whoAmI_ = aux.whoAmI_;
		}
		*/
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(const int n)
	{
		SetCounters();
		Initialize(n);
		if(dimensions_ != 0)
			memset(vector_+1, 0, dimensions_*sizeof(T));
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(const int n, const T v1, ...)
	{
		SetCounters();
		Initialize(n);
		T *w = vector_ + 1;
		va_list pointerList;
		va_start(pointerList, v1);
		*w = v1;
		for(int i=1+1;i<n+1;i++)
			*++w = va_arg(pointerList, T);
		va_end(pointerList);
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(const int n, const T* values)
	{
		SetCounters();
		Initialize(n);
		if(dimensions_ != 0)
		{
			T* w = vector_;
			w+=1;
			memcpy(w, values, n*sizeof(T));
		}
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(const int n, OpenSMOKEVector<T> const& rhs)
	{
		SetCounters();
		Initialize(n);
	
		if(n > rhs.dimensions_)
			ErrorMessage("OpenSMOKEVector(const int n, OpenSMOKEVector<T> const& rhs) requires n<=rhs.dimensions");

		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_;
			w+=1;
			v+=1;
			memcpy(w, v, n*sizeof(T));
		}
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(const int n, const int i, OpenSMOKEVector<T> const& rhs)
	{
		SetCounters();
		Initialize(n);

		if(n > rhs.dimensions_-i+1)
			ErrorMessage("OpenSMOKEVector(const int n, const int i, OpenSMOKEVector<T> const& rhs) requires n<=rhs.dimensions-i+1");
		
		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_;
			w+=1;
			v += i;
			memcpy(w, v, n*sizeof(T));
		}
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
	{
		SetCounters();
		Load(fileName, fileFormat);
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(std::istream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		SetCounters();
		Load(fInput, fileFormat);
	}

	template<typename T>
	OpenSMOKEVector<T>::OpenSMOKEVector(const std::vector<T> values)
	{
		SetCounters();
		Initialize(values.size());
		if(dimensions_ != 0)
		{
			T* w = this->vector_;
			w+=1;
			for (typename std::vector<T>::const_iterator it = values.begin(); it!=values.end(); ++it)
			{
				*w = *it;
				w++;
			}
		}
	}

	template<typename T>
	OpenSMOKEVector<T>::~OpenSMOKEVector(void)
	{
		if( matrixAsVector_ == true || subVectorAsVector_ == true )
		{
			countInScope_--;
			vector_ = 0;
			dimensions_ = 0;
			matrixAsVector_ = false;
			subVectorAsVector_ = false;
			return;
		}

		if(dimensions_ != 0)
			delete[] vector_;

		vector_ = 0;
		dimensions_ = 0;
		countInScope_--;
	}

	template<typename T>
	inline T* OpenSMOKEVector<T>::GetHandle()
	{
		return vector_ + 1;
	}
	
	template<typename T>
	inline const T* OpenSMOKEVector<T>::GetHandle() const	
	{
		return vector_ + 1;
	}

	template<typename T>
	void OpenSMOKEVector<T>::Copy(OpenSMOKEVector<T> const& orig)
	{
		if (dimensions_!=orig.dimensions_)
			ChangeDimensions(orig.Size(), this, false);
			//ErrorMessage("OpenSMOKEVector<T>::Copy(OpenSMOKEVector<T> const& orig) Dimension check failure");
						
		if(dimensions_ != 0)
			memcpy(vector_+1, orig.vector_+1, dimensions_*sizeof(T));
	}

	template<typename T>
	void OpenSMOKEVector<T>::Initialize(const int size)
	{
		matrixAsVector_		= false;
		subVectorAsVector_	= false;
		shadow_				= false;
		
		if (size < 0)
		{
			ErrorMessage("Size must be a positive integer");
		}
		else if(size == 0)
		{
			dimensions_ = 0;
			//vector_		= 0;
		}
		else
		{
			dimensions_ = size;
			vector_ = new T[size+1]; 
			if(!vector_)
				ErrorMessage("Size must be a positive integer");
		}
	}

	template<typename T>
	inline const T& OpenSMOKEVector<T>::At(const int i) const
	{
		return vector_[i];
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::GetValue(const int i) const
	{
		if( (i<1) || (i>dimensions_) )
			ErrorMessage("Vector index outside the ranges");
		return vector_[i];
	}
	
	template<typename T>
	inline void OpenSMOKEVector<T>::SetValue(const int i, const T value)
	{
		if( (i<1) || (i>dimensions_) )
			ErrorMessage("Vector index outside the ranges");
		vector_[i] = value;
	}

	template<typename T>
	inline T& OpenSMOKEVector<T>::operator() (int i)
	{
		if( (i<1) || (i>dimensions_) )
			ErrorMessage("Vector index outside the ranges");
		return vector_[i];
	}

	template<typename T>
	inline T& OpenSMOKEVector<T>::operator[] (const int i)
	{
		return vector_[i];
	}

	template<typename T>
	inline const T& OpenSMOKEVector<T>::operator[] (const int i) const
	{
		return vector_[i];
	}

	template<typename T>
	void OpenSMOKEVector<T>::operator() (const int n, OpenSMOKEVector<T> const& rhs)
	{
		Initialize(n);

		if(n>rhs.dimensions_)
			ErrorMessage("operator(const int n, OpenSMOKEVector<T> const& rhs) requires n<=rhs.dimensions");
		
		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_;
			w+=1;
			v+=1;
			memcpy(w, v, n*sizeof(T));
		}
	}

	template<typename T>
	void OpenSMOKEVector<T>::operator() (std::vector<T> const& rhs)
    {
		Initialize(rhs.size());
		std::cout << rhs.size() << std::endl;
		if(dimensions_ != 0)
		{
			T* w = this->vector_;
			w+=1;
			for (typename std::vector<T>::const_iterator it = rhs.begin(); it!=rhs.end(); ++it)
			{
				std::cout << *w << " " << *it << std::endl;
				*w= *it;
				std::cout << *w << " " << *it << std::endl;
				w++;
			}
		}
	}

	template<typename T>
	void OpenSMOKEVector<T>::operator() (const int n, const int start, OpenSMOKEVector<T> const& rhs)
	{
		Initialize(n);

		int max_n = rhs.dimensions_ - start + 1;

		if( n>max_n )
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T> const& rhs) requires n>max_n");
		else if ( max_n<1 )
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T> const& rhs) requires max_n>0");
		else if(start<1)
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T> const& rhs) requires start>=index");
		else if (start>rhs.dimensions_)
			ErrorMessage("operator(const int n, const int start, OpenSMOKEVector<T> const& rhs) requires start<=rhs.dimensions");

		if(dimensions_ != 0)
		{
			T* w = vector_;
			T* v = rhs.vector_ + start;
			w+=1;
			memcpy(w, v, n*sizeof(T));
		}
	}
        
        template<typename T>
        void OpenSMOKEVector<T>::operator = (T value)
	{
                if(value == 0 && dimensions_ > 0)
                        memset(vector_+1,0,dimensions_*sizeof(T));
                else
		{
                        for(int i=1;i<dimensions_+1;i++)
                                vector_[i] = value;
		}
	}

	template<typename T>
	OpenSMOKEVector<T>& OpenSMOKEVector<T>::operator =(OpenSMOKEVector<T> const& orig)
	{
			if (&orig!=this)
				Copy(orig);

			return *this;
	}

	template<typename T>
	OpenSMOKEVector<T>& OpenSMOKEVector<T>::operator += (OpenSMOKEVector<T> const& rval)
	{
		if(dimensions_ != rval.Size())
			ErrorMessage("operator+=(const OpenSMOKEVector<T>& rval) dimension check failure");
		if(whoAmI_ == rval.WhoAmI())
			Add(this);
		else
			Sum(dimensions_, vector_+1, rval.Vector()+1);
		return *this;
	}

	template<typename T>
	OpenSMOKEVector<T>& OpenSMOKEVector<T>::operator-=(OpenSMOKEVector<T> const& rval)
	{
		if(dimensions_ != rval.Size())
			ErrorMessage("operator-=(const OpenSMOKEVector<T>& rval) Dimension check failure");
		
		if(whoAmI_ == rval.WhoAmI())
			Sub(this);
		else
			Difference(dimensions_, vector_+1, rval.Vector() + 1);
		return *this;
	}

	
	template<typename T>
	inline OpenSMOKEVector<T>& OpenSMOKEVector<T>::operator +=(T const& rval)
	{
		Sum(dimensions_ , rval, vector_+1); 
		return *this;
	}

	template<typename T>
	inline OpenSMOKEVector<T>& OpenSMOKEVector<T>::operator -=(T const& rval)
	{
		Sum(dimensions_ , -rval, vector_+1); 
		return *this;
	}

	template<>
	inline OpenSMOKEVector<unsigned int>& OpenSMOKEVector<unsigned int>::operator -=(unsigned int const& rval)
	{
		ErrorMessage("operator-=(const OpenSMOKEVector<T>& rval) not available for unsigned int");
		return *this;
	}
	
	template<typename T>
	inline OpenSMOKEVector<T>& OpenSMOKEVector<T>::operator *=(T const& rval)
	{
		Prod(dimensions_ , rval, vector_+1); 
		return *this;
	}

	template<typename T>
	inline OpenSMOKEVector<T>& OpenSMOKEVector<T>::operator /=(T const& rval)
	{
		if (rval==0)
			ErrorMessage("operator/=(T const& rval) Division by zero");

		Prod(dimensions_ , static_cast<T>(1)/rval, vector_+1); 
		return *this;
	}

	template<>
	inline OpenSMOKEVector<int>& OpenSMOKEVector<int>::operator /=(int const& rval)
	{
		ErrorMessage("operator/=(int const& rval) not available");
		return *this;
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::SumElements() const
	{
		T sum = 0;
		T* x = vector_+1;
		for(int i=1;i<dimensions_+1;i++)
			sum += *x++;
		return sum;
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::SumAbsElements() const
	{
		T sum = 0;
		T* x = vector_+1;
		for(int i=1;i<dimensions_+1;i++)
			sum += abs(*x++);
		return sum;
	}

	template<>
	inline unsigned int OpenSMOKEVector<unsigned int>::SumAbsElements() const
	{
		unsigned int sum = 0;
		unsigned int* x = vector_ + 1;
		for (int i = 1; i<dimensions_ + 1; i++)
			sum += *x++;
		return sum;
	}

	template<typename T>
	void OpenSMOKEVector<T>::DeleteLastElements(const int n) 
	{
		if (n<=0)
			return;
		
		if(dimensions_ <= n)
		{
			OpenSMOKEVector<T> v;
			Swap(&v,this);
			return;
		}

		OpenSMOKEVector<T> v(dimensions_ - n);
		memmove(v.vector_+1, vector_+1, (dimensions_-n)*sizeof(T));
		Swap(&v,this);		
	}

	template<typename T>
	inline void OpenSMOKEVector<T>::PrintOnVideo() const
	{
		std::cout << typeid(OpenSMOKEVector<T>).name() << std::endl;
		std::cout << "Counters:      " << whoAmI_ << "/" << countInScope_ << "/" << count_ << std::endl;
		std::cout << "Size:          " << dimensions_ << std::endl;
		for(int i=1;i<dimensions_+1;i++)
			std::cout << i << " " << vector_[i] << std::endl;
		std::cout << std::endl;
	}

	template<typename T>
	void OpenSMOKEVector<T>::Load(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			std::ifstream fInput(fileName.c_str(), std::ios::in);
			if (fInput.fail())
				ErrorMessage("The " + fileName + " file does not exist");

			Load(fInput, fileFormat);
			
			fInput.close();
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			std::ifstream fInput(fileName.c_str(), std::ios::in | std::ios::binary);
			if (fInput.fail())
				ErrorMessage("The " + fileName + " binary file does not exist");

			Load(fInput, fileFormat);

			fInput.close();
		}
	}
		
	template<typename T>
	void OpenSMOKEVector<T>::Load(std::istream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			int n;
			fInput >> n;
                        Initialize(n);
			for(int i=1;i<dimensions_+1;i++)
				fInput >> vector_[i];
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			int n;

			// Reading vector size
			if(!fInput.read(reinterpret_cast<char *>(&n), sizeof(int)))
				ErrorMessage("I was unable to read from binary file (1)");
			
			// Initializing
			Initialize(n);

			// Reading vector elements
			for(int i=1;i<dimensions_+1;i++)
			if(!fInput.read(reinterpret_cast<char *>(&vector_[i]), sizeof(T)))
				ErrorMessage("I was unable to read from binary file (2)");
			else
				std::cout << n << " " << vector_[i] << std::endl;
		}
	}

		
	template<typename T>
	void OpenSMOKEVector<T>::Save(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			std::ofstream fOutput(fileName.c_str(), std::ios::out);
			if (fOutput.fail())
				ErrorMessage("The " + fileName + " file cannot be open");

			Save(fOutput, fileFormat);
			
			fOutput.close();
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			std::ofstream fOutput(fileName.c_str(), std::ios::out | std::ios::binary);
			if (fOutput.fail())
				ErrorMessage("The " + fileName + " binary file cannot be open");

			Save(fOutput, fileFormat);

			fOutput.close();
		}
	}
		
	template<typename T>
	void OpenSMOKEVector<T>::Save(std::ostream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fOutput << dimensions_ << std::endl;
			for(int i=1;i<dimensions_+1;i++)
				fOutput << std::setprecision(16) << std::scientific << vector_[i] << std::endl;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			// Writing vector size
			if(!fOutput.write( reinterpret_cast<char *>(&dimensions_), sizeof(int)))
				ErrorMessage("I was unable to write on binary file");

			// Writing vector elements
			for(int i=1;i<dimensions_+1;i++)
			if(!fOutput.write( reinterpret_cast<char *>(&vector_[i]), sizeof(T)))
				ErrorMessage("I was unable to write on binary file");
		}
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::Max(int *imax) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmax = OpenSMOKE::Max(dimensions_, vector_ + 1, imax);
		if(imax != 0)
			(*imax)+=1;
		return xmax;
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::MaxAbs(int *imax) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmax = OpenSMOKE::MaxAbs(dimensions_, vector_ + 1, imax);
		if(imax != 0) 
			(*imax)+=1;
		return xmax;
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::Min(int *imin) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmin = OpenSMOKE::Min(dimensions_, vector_ + 1, imin);
		if(imin != 0) 
			(*imin)+=1;
		return xmin;
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::MinAbs(int *imin) const
	{
		if(dimensions_ == 0)
			return 0;
		T xmin = OpenSMOKE::MinAbs(dimensions_, vector_ + 1, imin);
		if(imin != 0) 
			(*imin)+=1;
		return xmin;
	}

	template<typename T>
	void OpenSMOKEVector<T>::MinMax(int* iMin, T* min, int* iMax, T* max) const
	{
		if(dimensions_ == 0)
		{
			*iMin = 0;
			*iMax = 0;
			*min = 0;
			*max = 0;
			return;
		}

		*iMin = *iMax = 1 ;
		*min = *max = vector_[1];
	
		for(int i=1+1;i<dimensions_+1;i++)
		{
			if(vector_[i] < *min)
			{
				*min = vector_[i];
				*iMin = i;
			}
			else if(vector_[i] > *max)
			{
				*max = vector_[i];
				*iMax = i;
			}
		}
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::Norm1() const
	{
		T norm = static_cast<T>(0);
		T* x = vector_+1;
		for(int i=1;i<dimensions_+1;i++)
			norm += abs(*x++);
		return norm;
	}

	template<>
	inline unsigned int OpenSMOKEVector<unsigned int>::Norm1() const
	{
		unsigned int norm = static_cast<unsigned int>(0);
		unsigned int* x = vector_ + 1;
		for (int i = 1; i<dimensions_ + 1; i++)
			norm += *x++;
		return norm;
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::Norm2() const
	{
		T* x = vector_+1;
		return SqrtSumSqr(dimensions_, x);
	}

	template<>
	inline int OpenSMOKEVector<int>::Norm2() const
	{
		ErrorMessage("Norm2() function not available for int vectors");
		return 0;
	}

	template<>
	inline unsigned int OpenSMOKEVector<unsigned int>::Norm2() const
	{
		ErrorMessage("Norm2() function not available for unsigned int vectors");
		return 0;
	}

	template<typename T>
	inline T OpenSMOKEVector<T>::NormInf() const
	{
		return MaxAbs();
	}

	template<typename T>
	void OpenSMOKEVector<T>::Append(const T f)
	{
		OpenSMOKEVector<T> v(dimensions_+1);
		if(dimensions_ != 0)
			memcpy(v.vector_, vector_,(dimensions_+1)*sizeof(T));
		v[dimensions_+1] = f;
		Swap(&v,this);
	}

	template<typename T>
	void OpenSMOKEVector<T>::Append(OpenSMOKEVector<T> const& w)
	{
		int sz = w.dimensions_;
		OpenSMOKEVector<T> v(dimensions_ + sz);
		if(dimensions_ != 0)
			memcpy(v.vector_, vector_, (dimensions_+1)*sizeof(T));
	
		memcpy(v.vector_ + dimensions_ + 1, w.vector_ + 1, sz*sizeof(T));
		Swap(&v,this);
	}

	template<typename T>
	void OpenSMOKEVector<T>::Insert(const int i, const T f)
	{
		if(i < 1) 
			ErrorMessage("Insert(const int i, const T f) requires i>=index");
		
		if(i < dimensions_+1)
		{
			OpenSMOKEVector<T> v(dimensions_+1);
			memmove(v.vector_, vector_, i*sizeof(T));
			v[i] = f;
			memmove(v.vector_+i+1, vector_+i, (dimensions_-i+1)*sizeof(T));
			Swap(&v,this);		
		}
		else
		{
			OpenSMOKEVector<T> v(i);
			memmove(v.vector_ + 1, vector_ + 1, dimensions_ * sizeof(T));
			v[i] = f;
			Swap(&v,this);		
		}
	}

	template<typename T>
	int OpenSMOKEVector<T>::InsertElementInSortedVector(const T f)
	{
		const int j = LocateInSortedVector(f);
		Insert(j+1, f);
		return j+1;
	}

	template<typename T>
	int OpenSMOKEVector<T>::InsertElementInFirstNSortedElements(const int n, const T f)
	{
		if(n >= this->dimensions_)
		{
			OpenSMOKEVector<T> v(this->dimensions_+10);
			if(this->dimensions_ != 0)
				memcpy(v.vector_, this->vector_, (this->dimensions_+1)*sizeof(T));
			Swap(&v,this);
		}
	
		const int j = LocateInFirstNSortedElements(n,f);
		if(j != n)
		{
			memmove(this->vector_+j+2, this->vector_+j+1, (this->dimensions_-j-1)*sizeof(T));
		}
		(*this)[j+1] = f;
		return j+1;
	}

	template<typename T>
	void OpenSMOKEVector<T>::Insert(const int i, OpenSMOKEVector<T> const& w)
	{
		int sz = w.dimensions_;
		
		OpenSMOKEVector<T> v(dimensions_ + sz);
		if(i < 1) 
			ErrorMessage("Insert(const int i, OpenSMOKEVector<T> const& w) requires i>=index");

		if(i < dimensions_+1)
		{
			memmove(v.vector_, vector_,i*sizeof(T));
			memmove(v.vector_ + i, w.vector_ + 1, sz*sizeof(T));
			memmove(v.vector_ + sz + i, vector_ + i, (dimensions_ - i + 1)*sizeof(T));
			Swap(&v,this);		
		}
		else
			Append(w);
	}

	template<typename T>
	int OpenSMOKEVector<T>::LocateInSortedVector(const T f)
	{
		int n = this->dimensions_;

		if(n < 1)
			return 0;
		if(n == 1)
		{
			if(this->vector_[1] < f)	return 1;
			else                                    return 1-1;
		}
		
		int upper,lower,jm;
		lower = 0;
		upper = n + 1;
		if(this->vector_[n] > this->vector_[1])
		{
			while(upper - lower > 1)
			{
				jm = (upper + lower) >> 1;
				if(f > this->vector_[jm-1+1])
					lower = jm;
				else
					upper = jm;
			}
			return lower-1+1;
		}
		else
		{
			while(upper - lower > 1)
			{
				jm = (upper + lower) >> 1;
				if(f < this->vector_[jm-1+1])
					lower = jm;
				else
					upper = jm;
			}
			return lower-1+1;
		}
	}

	template<typename T>
	int OpenSMOKEVector<T>::LocateInFirstNSortedElements(const int position, const T f)
	{
		int n = position;
		
		if(n > this->dimensions_)
			n = this->dimensions_;
	
		if(n < 1)
			return 0;
	
		if(n == 1)
		{
			if(this->vector_[1-1+1] < f)	return 1;
			else									return 0;
		}
	
		int upper,lower;
		lower = 0;
		upper = n+1;
	
		if(this->vector_[n-1+1] > this->vector_[1-1+1])
		{
			while(upper - lower > 1)
			{
				const int jm = (upper + lower) >> 1;
				if(f > this->vector_[jm-1+1])
					lower = jm;
				else
					upper = jm;
			}
			return lower;
		}
		else
		{
			while(upper - lower > 1)
			{
				const int jm = (upper + lower) >> 1;
				if(f < this->vector_[jm-1+1])
					lower = jm;
				else
					upper = jm;
			}
			return lower;
		}
	}

	// Friend Functions
	template<typename T>
	void Swap(OpenSMOKEVector<T>* lval, OpenSMOKEVector<T>* rval)
	{
		Swap(&lval->vector_,     &rval->vector_);
		Swap(&lval->dimensions_, &rval->dimensions_);
	}

	template<typename T>
	void ChangeDimensions(const int size, OpenSMOKEVector<T>* result, bool reset)
	{
		if(size != result->dimensions_)
		{
			if(result->matrixAsVector_ == true)
				ErrorMessage("ChangeDimensions", "You can not change the dimensions when you are using a matrix as a vector");
	
			if(result->subVectorAsVector_ == true)
				ErrorMessage("ChangeDimensions", "You can not change the dimensions when you are using a subvector as a vector");
		
                        if (result->dimensions_ != 0)
                                delete[] result->vector_;

			result->Initialize(size);
		}
	
		// Resetting
		if(reset == true && result->dimensions_ != 0)
			memset(result->vector_+1, 0, (result->dimensions_)*sizeof(T) );
	}

	template<typename T>
	void Add(const OpenSMOKEVector<T>& lval, const OpenSMOKEVector<T>& rval, OpenSMOKEVector<T>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Sum(const OpenSMOKEVector<T>& lval, const OpenSMOKEVector<T>& rval, OpenSMOKEVector<T>* result)",  
							"Vector dimension check failure");
	
		if(result->whoAmI_ == lval.whoAmI_)
			(*result) += rval;
		else if(result->whoAmI_ == rval.whoAmI_)
			Add(lval,result);
		else
		{
			ChangeDimensions(lval.dimensions_, result, false);
			Sum(result->dimensions_, lval.vector_ + 1, rval.vector_ + 1, result->vector_ + 1);
		}
	}

	template<>
	void Add(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Sum(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result)",  
							"Vector dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Sum(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const float* ptlval = lval.GetHandle();
            const float* ptrval = rval.GetHandle();
				  float* ptresult = result->GetHandle();

            vsAdd( n, ptlval, ptrval, ptresult );
		}
		#else
		{
			if(result->whoAmI_ == lval.whoAmI_)
				(*result) += rval;
			else if(result->whoAmI_ == rval.whoAmI_)
				Add(lval,result);
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				Sum(result->dimensions_, lval.vector_ + 1, rval.vector_ + 1, result->vector_ + 1);
			}
		}
		#endif
	}

	template<>
	void Add(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Sum(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result)",  
							"Vector dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Sum(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const double* ptlval = lval.GetHandle();
            const double* ptrval = rval.GetHandle();
			      double* ptresult = result->GetHandle();

            vdAdd( n, ptlval, ptrval, ptresult );
		}
		#else
		{
			if(result->whoAmI_ == lval.whoAmI_)
				(*result) += rval;
			else if(result->whoAmI_ == rval.whoAmI_)
				Add(lval,result);
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				Sum(result->dimensions_, lval.vector_ + 1, rval.vector_ + 1, result->vector_ + 1);
			}
		}
		#endif
	}

	template<typename T>
	void Add(const OpenSMOKEVector<T>& lval, const T rval, OpenSMOKEVector<T>* result)
	{	
		if(result->whoAmI_ == lval.whoAmI_)
			(*result) += rval;
		else
		{
			ChangeDimensions(lval.dimensions_, result, false);
			Sum(result->dimensions_, lval.vector_ + 1, rval, result->vector_ + 1);
		}
	}


	template<typename T>
	void Add(OpenSMOKEVector<T>* lvalAndResult, const OpenSMOKEVector<T>& rval)
	{
		(*lvalAndResult) += rval;
	}

	template<typename T>
	void Add(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rvalAndResult)
	{
		(*rvalAndResult) += lval;
	}

	template<typename T>
	void Add(OpenSMOKEVector<T>* lvalRvalAndResult)
	{
		Sum(lvalRvalAndResult->dimensions_, lvalRvalAndResult->vector_+1);
	}

	template<typename T>
	void Sub(const OpenSMOKEVector<T>& lval, const OpenSMOKEVector<T>& rval, OpenSMOKEVector<T>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Difference(const OpenSMOKEVector<T>& lval, const OpenSMOKEVector<T>& rval, OpenSMOKEVector<T>* result)",
							"Dimension check failure");
		
		if(result->whoAmI_ == lval.whoAmI_)
			(*result) -= rval;
		else if(result->whoAmI_ == rval.whoAmI_)
			Sub(lval,result);
		else
		{
			ChangeDimensions(lval.dimensions_, result, false);
			Difference(result->dimensions_, lval.vector_+1, rval.vector_ + 1, result->vector_ + 1); 
		}
	}

	template<>
	void Sub(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Difference(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Difference(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const float* ptlval = lval.GetHandle();
            const float* ptrval = rval.GetHandle();
		          float* ptresult = result->GetHandle();

            vsSub( n, ptlval, ptrval, ptresult );
		}
		#else
		{
			if(result->whoAmI_ == lval.whoAmI_)
				(*result) -= rval;
			else if(result->whoAmI_ == rval.whoAmI_)
				Sub(lval,result);
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				Difference(result->dimensions_, lval.vector_+1, rval.vector_ + 1, result->vector_ + 1); 
			}
		}
		#endif
	}

	template<>
	void Sub(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Difference(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Difference(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const double* ptlval = lval.GetHandle();
            const double* ptrval = rval.GetHandle();
			      double* ptresult = result->GetHandle();

            vdSub( n, ptlval, ptrval, ptresult );
		}
		#else
		{
			if(result->whoAmI_ == lval.whoAmI_)
				(*result) -= rval;
			else if(result->whoAmI_ == rval.whoAmI_)
				Sub(lval,result);
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				Difference(result->dimensions_, lval.vector_+1, rval.vector_ + 1, result->vector_ + 1); 
			}
		}
		#endif
	}

	template<typename T>
	void Sub(OpenSMOKEVector<T>* lvalAndResult, const OpenSMOKEVector<T>& rval)
	{
		(*lvalAndResult) -= rval;
	}

	template<typename T>
	void Sub(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rvalAndResult)
	{
		if(lval.dimensions_ != rvalAndResult->dimensions_)
			ErrorMessage("Difference(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rvalAndResult)", "Dimension check failure");
		DifferenceBis(lval.dimensions_, lval.vector_+1, rvalAndResult->vector_+1);
	}

	template<typename T>
	void Sub(OpenSMOKEVector<T>* lvalRvalAndResult)
	{
		T *w = lvalRvalAndResult->vector_+1;
		for(int i=1;i<=lvalRvalAndResult->dimensions_;i++)
			*w++ = 0;
	}

	template<typename T>
	void Product(const T lval, const OpenSMOKEVector<T>& rval, OpenSMOKEVector<T>* result)
	{
    	if (rval.WhoAmI() == result->WhoAmI())
        	Prod(result->dimensions_, lval, result->vector_ + 1);
    	else
    	{
        	ChangeDimensions(rval.dimensions_, result, false);
        	Prod(rval.dimensions_, lval, rval.vector_ + 1, result->vector_ + 1);
    	}
	}

	template<>
	void Product(const float lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	if(rval.dimensions_ != result->dimensions_)
            	ChangeDimensions(rval.dimensions_, result, false);

        	result->Copy(rval); 
        	int one = 1;
        	int n = rval.dimensions_;
			float* ptresult = result->GetHandle();
        	sscal(&n, &lval, ptresult, &one);
    	}
		#else
    	if (rval.WhoAmI() == result->WhoAmI())
        	Prod(result->dimensions_, lval, result->vector_ + 1);
    	else
    	{
        	ChangeDimensions(rval.dimensions_, result, false);
        	Prod(rval.dimensions_, lval, rval.vector_ + 1, result->vector_ + 1);
    	}
		#endif
	}

	template<>
	void Product(const double lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	if(rval.dimensions_ != result->dimensions_)
            	ChangeDimensions(rval.dimensions_, result, false);

        	result->Copy(rval); 
        	int one = 1;
        	int n = rval.dimensions_;
			double* ptresult = result->GetHandle();
        	dscal(&n, &lval, ptresult, &one);
    	}
		#else
    	if (rval.WhoAmI() == result->WhoAmI())
        	Prod(result->dimensions_, lval, result->vector_ + 1);
    	else
    	{
        	ChangeDimensions(rval.dimensions_, result, false);
        	Prod(rval.dimensions_, lval, rval.vector_ + 1, result->vector_ + 1);
    	}
		#endif
	}

	template<typename T>
	void Product(const T lval, OpenSMOKEVector<T>* rvalAndResult)
	{
    	Prod(rvalAndResult->dimensions_, lval, rvalAndResult->vector_+1); 
	}

	template<>
	void Product(const float lval, OpenSMOKEVector<float>* rvalAndResult)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	int one = 1;
        	int n = rvalAndResult->dimensions_;
        	float* ptr_rvalAndResult = rvalAndResult->GetHandle();

        	sscal(&n, &lval, ptr_rvalAndResult, &one);
    	}
		#else
    	Prod(rvalAndResult->dimensions_, lval, rvalAndResult->vector_+1); 
		#endif
	}

	template<>
	void Product(const double lval, OpenSMOKEVector<double>* rvalAndResult)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	int one = 1;
        	int n = rvalAndResult->dimensions_;
        	double* ptr_rvalAndResult = rvalAndResult->GetHandle();

        	dscal(&n, &lval, ptr_rvalAndResult, &one);
    	}
		#else
    	Prod(rvalAndResult->dimensions_, lval, rvalAndResult->vector_+1); 
		#endif
	}

	template<typename T>
	void DotProduct(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("DotProduct(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, T* result)", "Dimension check failure");
		
		*result = Dot(lval.dimensions_, lval.vector_+1, rval.vector_+1); 
	}

	template<typename T>
	void UDotProduct(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("UDotProduct(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, T* result)", "Dimension check failure");
		
		*result = UDot(lval.dimensions_, lval.vector_+1, rval.vector_+1); 
	}

	template<typename T>
	T Dot(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval)
	{
		T result = 0;
		DotProduct(lval, rval, &result);
		return result;
	}

	template<typename T>
	T UDot(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval)
	{
		T result = static_cast<T>(0);
		UDotProduct(lval, rval, &result);
		return result;
	}

	template<typename T>
	void DotProduct(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, const int start, const int elements, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("Dot(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, const int start, const int elements, T* result)", "Dimension check failure");

		*result = Dot(elements, lval.vector_+start, rval.vector_+start);
	}

	template<typename T>
	void UDotProduct(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, const int start, const int elements, T* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage("UDot(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, const int start, const int elements, T* result)", "Dimension check failure");

		*result = UDot(elements, lval.vector_+start, rval.vector_+start);
	}

	template<typename T>
	void ElementByElementProduct(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, OpenSMOKEVector<T>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"ElementByElementProduct(OpenSMOKEVector<T> const& l, OpenSMOKEVector<T> const& r, OpenSMOKEVector<T>* s)",
							"Dimension check failure");

		if(rval.whoAmI_ == result->whoAmI_ || lval.whoAmI_ == result->whoAmI_)
		{
			OpenSMOKEVector<T> aux;
			ElementByElementProduct(lval, rval, &aux);	
			Swap(&aux,result);						
		}
		else
		{
			ChangeDimensions(lval.dimensions_, result, false);
			T* w = lval.vector_ + 1;
			T* v = rval.vector_ + 1;
			T* x = result->vector_ + 1;
			for(int i=1;i<=lval.dimensions_;i++)
				(*x++) = (*w++) * (*v++);
		}        
	}

	template<>
	void ElementByElementProduct(OpenSMOKEVector<float> const& lval, OpenSMOKEVector<float> const& rval, OpenSMOKEVector<float>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"ElementByElementProduct(OpenSMOKEVector<float> const& l, OpenSMOKEVector<float> const& r, OpenSMOKEVector<float>* s)",
							"Dimension check failure");
		
        #if OPENSMOKE_USE_MKL == 1
        {
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"ElementByElementProduct(OpenSMOKEVector<float> const& l, OpenSMOKEVector<float> const& r, OpenSMOKEVector<float>* s)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const float* ptlval = lval.GetHandle();
			const float* ptrval = rval.GetHandle();
			      float* ptresult = result->GetHandle();

			vsMul( n, ptlval, ptrval, ptresult );      
        }
        #else
        {
			if(rval.whoAmI_ == result->whoAmI_ || lval.whoAmI_ == result->whoAmI_)
			{
				OpenSMOKEVector<float> aux;
				ElementByElementProduct(lval, rval, &aux);	
				Swap(&aux,result);						
			}
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				float* w = lval.vector_ + 1;
				float* v = rval.vector_ + 1;
				float* x = result->vector_ + 1;
				for(int i=1;i<=lval.dimensions_;i++)
					(*x++) = (*w++) * (*v++);
			}        
        }
        #endif
	}

	template<>
	void ElementByElementProduct(OpenSMOKEVector<double> const& lval, OpenSMOKEVector<double> const& rval, OpenSMOKEVector<double>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"ElementByElementProduct(OpenSMOKEVector<double> const& l, OpenSMOKEVector<double> const& r, OpenSMOKEVector<double>* s)",
							"Dimension check failure");
		
        #if OPENSMOKE_USE_MKL == 1
        {
			if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"ElementByElementProduct(OpenSMOKEVector<double> const& l, OpenSMOKEVector<double> const& r, OpenSMOKEVector<double>* s)",
								"Vector dimension check failure");

			int n = lval.dimensions_;
			const double* ptlval = lval.GetHandle();
			const double* ptrval = rval.GetHandle();
			      double* ptresult = result->GetHandle();

			vdMul( n, ptlval, ptrval, ptresult );      
        }
        #else
        {
			if(rval.whoAmI_ == result->whoAmI_ || lval.whoAmI_ == result->whoAmI_)
			{
				OpenSMOKEVector<double> aux;
				ElementByElementProduct(lval, rval, &aux);	
				Swap(&aux,result);						
			}
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				double* w = lval.vector_ + 1;
				double* v = rval.vector_ + 1;
				double* x = result->vector_ + 1;
				for(int i=1;i<=lval.dimensions_;i++)
					(*x++) = (*w++) * (*v++);
			}        
        }
        #endif
	}

	template<typename T>
    void Division(const OpenSMOKEVector<T>& lval, const T rval, OpenSMOKEVector<T>* result)
	{
    	if (lval.WhoAmI() == result->WhoAmI())
        	Div(result->dimensions_, result->vector_ + 1, rval);
    	else
    	{
        	ChangeDimensions(lval.dimensions_, result, false);
        	Div(lval.dimensions_, lval.vector_ + 1, rval, result->vector_ + 1);
    	}
	}

	template<>
    void Division(const OpenSMOKEVector<float>& lval, const float rval, OpenSMOKEVector<float>* result)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	if(lval.Dimensions() != result->Dimensions())
            	ChangeDimensions(lval.Dimensions(), result, false);
        
        	result->Copy(lval);
        	float rval_inv = static_cast<float>(1.) / rval;
        	int one = 1;
        	int n = lval.Dimensions();
        	float* ptresult = result->GetHandle();
        	sscal(&n, &rval_inv, ptresult, &one);        
    	}
		#else
    	if (lval.WhoAmI() == result->WhoAmI())
        	Div(result->Dimensions(), result->Vector() + 1, rval);
    	else
    	{
        	ChangeDimensions(lval.Dimensions(), result, false);
        	Div(lval.Dimensions(), lval.Vector() + 1, rval, result->Vector() + 1);
    	}
		#endif
	}

	template<>
    void Division(const OpenSMOKEVector<double>& lval, const double rval, OpenSMOKEVector<double>* result)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	if(lval.Dimensions() != result->Dimensions())
            	ChangeDimensions(lval.Dimensions(), result, false);
        
        	result->Copy(lval);
        	double rval_inv = 1. / double(rval);
        	int one = 1;
        	int n = lval.Dimensions();
        	double* ptresult = result->GetHandle();
        	dscal(&n, &rval_inv, ptresult, &one);        
    	}
		#else
    	if (lval.WhoAmI() == result->WhoAmI())
        	Div(result->Dimensions(), result->Vector() + 1, rval);
    	else
    	{
        	ChangeDimensions(lval.Dimensions(), result, false);
        	Div(lval.Dimensions(), lval.Vector() + 1, rval, result->Vector() + 1);
    	}
		#endif
	}

	template<typename T>
	void Division(OpenSMOKEVector<T>* lvalAndResult, const T rval)
	{
		Div(lvalAndResult->dimensions_, lvalAndResult->vector_ + 1, rval);
	}

	template<>
	void Division(OpenSMOKEVector<float>* lvalAndResult, const float rval)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	float rval_inv = static_cast<float>(1.) / rval;
        	int one = 1;
        	int n = lvalAndResult->dimensions_;
			float* ptr_lvalAndResult = lvalAndResult->GetHandle();
        	sscal(&n, &rval_inv, ptr_lvalAndResult, &one);        
    	}
		#else
		Div(lvalAndResult->dimensions_, lvalAndResult->vector_ + 1, rval);
		#endif
	}

	template<>
	void Division(OpenSMOKEVector<double>* lvalAndResult, const double rval)
	{
		#if OPENSMOKE_USE_MKL == 1
    	{
        	double rval_inv = 1. / double(rval);
        	int one = 1;
        	int n = lvalAndResult->dimensions_;
			double* ptr_lvalAndResult = lvalAndResult->GetHandle();
        	dscal(&n, &rval_inv, ptr_lvalAndResult, &one);        
    	}
		#else
		Div(lvalAndResult->dimensions_, lvalAndResult->vector_ + 1, rval);
		#endif
	}

	template<typename T>
	void ElementByElementDivision(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, OpenSMOKEVector<T>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"ElementByElementDivision(OpenSMOKEVector<T> const& l, OpenSMOKEVector<T> const& r, OpenSMOKEVector<T>* s)",
							"Dimension check failure");
		
        #if OPENSMOKE_USE_MKL == 1
        {
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
			const T* ptrval = rval.GetHandle();
				  T* ptresult = result->GetHandle();

			vdDiv( n, ptlval, ptrval, ptresult );      
        }
        #else
        {
			if(rval.whoAmI_ == result->whoAmI_ || lval.whoAmI_ == result->whoAmI_)
			{
				OpenSMOKEVector<T> aux;
				ElementByElementDivision(lval, rval, &aux);	
				Swap(&aux,result);						
			}
			else
			{
				ChangeDimensions(lval.dimensions_, result, false);
				T* w = lval.vector_ + 1;
				T* v = rval.vector_ + 1;
				T* x = result->vector_ + 1;
				for(int i=1;i<=lval.dimensions_;i++)
					(*x++) = (*w++) / (*v++);
			}     
        }
        #endif
	}

	template<typename T>
	void OpenSMOKEVector<T>::CopyTo(T* rhs) const
	{		
		if(dimensions_ != 0)
			memcpy(rhs, vector_+1, dimensions_*sizeof(T));
	}

	template<typename T>
	void OpenSMOKEVector<T>::CopyFrom(const T* rhs)
	{		
		if(dimensions_ != 0)
			memcpy(vector_+1, rhs, dimensions_*sizeof(T));
	}
        
	template<typename T>
	void Exp(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Exp(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)", "Dimension check failure");
                
		for(int i=1;i<rval->dimensions_+1;i++)
			rval->vector_[i] = std::exp(lval.vector_[i]);          
	}

	template<>
	void Exp(OpenSMOKEVector<float> const& lval, OpenSMOKEVector<float>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Exp(OpenSMOKEVector<float> const& lval, OpenSMOKEVector<float>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const float* ptlval = lval.GetHandle();
                float* ptrval = rval->GetHandle();
                vsExp( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = std::exp(lval.vector_[i]);          
        }
        #endif
	}

	template<>
	void Exp(OpenSMOKEVector<double> const& lval, OpenSMOKEVector<double>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Exp(OpenSMOKEVector<double> const& lval, OpenSMOKEVector<double>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const double* ptlval = lval.GetHandle();
                double* ptrval = rval->GetHandle();
                vdExp( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = std::exp(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T>
	void Ln(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Ln(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdLn( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = log(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T>
	void Log10(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Log10(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdLog10( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = log10(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T>
	void Sin(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Sin(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdSin( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = std::sin(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T>
	void Cos(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Cos(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdCos( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = std::cos(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T>
	void Sqrt(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage("Sqrt(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T>* rval)", "Dimension check failure");
                
        #if OPENSMOKE_USE_MKL == 1
        {
                int n = lval.dimensions_;
                const T* ptlval = lval.GetHandle();
                T* ptrval = rval->GetHandle();
                vdSqrt( n, ptlval, ptrval );       
        }
        #else
        {
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = sqrt(lval.vector_[i]);          
        }
        #endif
	}

	template<typename T>
	void Sqr(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage(	"Sqr(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rval)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
                  T* ptrval = rval->GetHandle();

            vdSqr( n, ptlval, ptrval);
		}
		#else
		{
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = lval.vector_[i]*lval.vector_[i]; 
		}
		#endif
	}

	template<typename T>
	void Inversion(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rval)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage(	"Inversion(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rval)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
                  T* ptrval = rval->GetHandle();

            vdInv( n, ptlval, ptrval);
		}
		#else
		{
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = 1./lval.vector_[i]; 
		}
		#endif
	}

	template<typename T>
	void Pow(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rval, const double power)
	{
		if(lval.dimensions_ != rval->dimensions_)
			ErrorMessage(	"Pow(const OpenSMOKEVector<T>& lval, OpenSMOKEVector<T>* rval, const double power)",
							"Dimension check failure");
		
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
                  T* ptrval = rval->GetHandle();

            vdPowx( n, ptlval, power, ptrval);
		}
		#else
		{
                for(int i=1;i<rval->dimensions_+1;i++)
                    rval->vector_[i] = std::pow(lval.vector_[i], power);
		}
		#endif
	}

	template<typename T>
	void Pow(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, OpenSMOKEVector<T>* result)
	{
		if(lval.dimensions_ != rval.dimensions_)
			ErrorMessage(	"Pow(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, OpenSMOKEVector<T>* result)",
							"Dimension check failure");
		
		if(lval.dimensions_ != result->dimensions_)
				ErrorMessage(	"Pow(OpenSMOKEVector<T> const& lval, OpenSMOKEVector<T> const& rval, OpenSMOKEVector<T>* result)",
								"Vector dimension check failure");

        #if OPENSMOKE_USE_MKL == 1
        {
			int n = lval.dimensions_;
			const T* ptlval = lval.GetHandle();
			const T* ptrval = rval.GetHandle();
				  T* ptresult = result->GetHandle();

			vdPow( n, ptlval, ptrval, ptresult );      
        }
        #else
        {
                for(int i=1;i<rval.dimensions_;i++)
                    result->vector_[i] = std::pow(lval.vector_[i],rval.vector_[i]);     
        }
        #endif
	}

	template<typename T>
	bool operator == (const OpenSMOKEVector<T> &lval, const OpenSMOKEVector<T> &rval)
	{
		if(lval.whoAmI_ == rval.whoAmI_) 
			return true;
		
		bool flag = true;
		if(lval.dimensions_ != rval.dimensions_) 
			flag = false;
		else
		{
			if(memcmp(lval.vector_ + 1,rval.vector_ + 1, rval.dimensions_*sizeof(T)) == 0)
				flag = true;
			else flag = false;
		}
		return flag;
	}

	template<typename T>
	void Sort(OpenSMOKEVector<T> *result)
	{
		Sort(result->dimensions_, result->vector_+1); 
	}

	template<typename T>
	void Sort(OpenSMOKEVector<T> *result, OpenSMOKEVector<int>  *iS)
	{
		ChangeDimensions(result->dimensions_, iS, true);
		for(int i=1;i<=result->dimensions_;i++)
			(*iS)[i] = i;
		Sort(result->dimensions_, result->vector_ + 1,iS->vector_ + 1);
	}

	template<typename T>
	void Reorder(OpenSMOKEVector<T> *v, OpenSMOKEVector<int>  &iS)
	{
		int size =v->Size();

		if(size != iS.Size())
			ErrorMessage("void Reorder(OpenSMOKEVector<T> *v, OpenSMOKEVector<int>  &iS)", "iS has a wrong size");

		T *w = new T[size+1];
		if(!w)
			ErrorMessage("void Reorder(OpenSMOKEVector<T> *v, OpenSMOKEVector<int>  &iS)", "no enough memory");

		for(int i=1;i<=size;i++)
			w[i] = (*v)[iS[i]];
		
		Swap(&w, &v->vector_);

		delete[] w;
	}
        
        template<typename T>
	void Reverse(OpenSMOKEVector<T> *result)
	{
            const int n = result->Size();
            const int m = n/2;	
            if(n < 2)
                return;
	
            T *v = result->vector_ + 1;
            T *w = result->vector_ + n + 1 -1;
	
            for(unsigned int i=1;i<=m;i++,v++,w--)
                Swap(v,w);
	}

	// Pre-instantiations
	
		template class OpenSMOKEVector<unsigned int>;
		template class OpenSMOKEVector<int>;
		template class OpenSMOKEVector<float>;
		template class OpenSMOKEVector<double>;

		template void Add(const OpenSMOKEVector<unsigned int>& lval, const OpenSMOKEVector<unsigned int>& rval, OpenSMOKEVector<unsigned int>* result);
		template void Add(const OpenSMOKEVector<int>& lval, const OpenSMOKEVector<int>& rval, OpenSMOKEVector<int>* result);
		template void Add(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result);
		template void Add(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result);

		template void Add(const OpenSMOKEVector<unsigned int>& lval, const unsigned int rval, OpenSMOKEVector<unsigned int>* result);
		template void Add(const OpenSMOKEVector<int>& lval, const int rval, OpenSMOKEVector<int>* result);
		template void Add(const OpenSMOKEVector<float>& lval, const float rval, OpenSMOKEVector<float>* result);
		template void Add(const OpenSMOKEVector<double>& lval, const double rval, OpenSMOKEVector<double>* result);


		template void Sub(const OpenSMOKEVector<unsigned int>& lval, const OpenSMOKEVector<unsigned int>& rval, OpenSMOKEVector<unsigned int>* result);
		template void Sub(const OpenSMOKEVector<int>& lval, const OpenSMOKEVector<int>& rval, OpenSMOKEVector<int>* result);
		template void Sub(const OpenSMOKEVector<float>& lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result);
		template void Sub(const OpenSMOKEVector<double>& lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result);

		template void Sub(OpenSMOKEVector<unsigned int>* lvalAndResult, const OpenSMOKEVector<unsigned int>& rval);
		template void Sub(OpenSMOKEVector<int>* lvalAndResult, const OpenSMOKEVector<int>& rval);
		template void Sub(OpenSMOKEVector<float>* lvalAndResult, const OpenSMOKEVector<float>& rval);
		template void Sub(OpenSMOKEVector<double>* lvalAndResult, const OpenSMOKEVector<double>& rval);


		template void Product(const unsigned int lval, const OpenSMOKEVector<unsigned int>& rval, OpenSMOKEVector<unsigned int>* result);
		template void Product(const int lval, const OpenSMOKEVector<int>& rval, OpenSMOKEVector<int>* result);
		template void Product(const float lval, const OpenSMOKEVector<float>& rval, OpenSMOKEVector<float>* result);
		template void Product(const double lval, const OpenSMOKEVector<double>& rval, OpenSMOKEVector<double>* result);

		template void Product(const unsigned int lval, OpenSMOKEVector<unsigned int>* rvalAndResult);
		template void Product(const int lval, OpenSMOKEVector<int>* rvalAndResult);
		template void Product(const float lval, OpenSMOKEVector<float>* rvalAndResult);
		template void Product(const double lval, OpenSMOKEVector<double>* rvalAndResult);

		template void ElementByElementProduct(OpenSMOKEVector<unsigned int> const& lval, OpenSMOKEVector<unsigned int> const& rval, OpenSMOKEVector<unsigned int>* result);
		template void ElementByElementProduct(OpenSMOKEVector<int> const& lval, OpenSMOKEVector<int> const& rval, OpenSMOKEVector<int>* result);
		template void ElementByElementProduct(OpenSMOKEVector<float> const& lval, OpenSMOKEVector<float> const& rval, OpenSMOKEVector<float>* result);
		template void ElementByElementProduct(OpenSMOKEVector<double> const& lval, OpenSMOKEVector<double> const& rval, OpenSMOKEVector<double>* result);

		template void Exp(OpenSMOKEVector<float> const& lval, OpenSMOKEVector<float>* rval);
		template void Exp(OpenSMOKEVector<double> const& lval, OpenSMOKEVector<double>* rval);

		template float UDot(OpenSMOKEVector<float> const& lval, OpenSMOKEVector<float> const& rval);
		template double UDot(OpenSMOKEVector<double> const& lval, OpenSMOKEVector<double> const& rval);

		template unsigned int Dot(OpenSMOKEVector<unsigned int> const& lval, OpenSMOKEVector<unsigned int> const& rval);
		template int Dot(OpenSMOKEVector<int> const& lval, OpenSMOKEVector<int> const& rval);
		template float Dot(OpenSMOKEVector<float> const& lval, OpenSMOKEVector<float> const& rval);
		template double Dot(OpenSMOKEVector<double> const& lval, OpenSMOKEVector<double> const& rval);
	
		template void Reorder(OpenSMOKEVector<int> *v, OpenSMOKEVector<int>  &iS);
		template void Reorder(OpenSMOKEVector<double> *v, OpenSMOKEVector<int>  &iS);
		
		template void Sort(OpenSMOKEVector<int> *result, OpenSMOKEVector<int>  *iS);
}
