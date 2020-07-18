/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include <typeinfo>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iomanip>
#include "OpenSMOKEMatrix.h"

#if OPENSMOKE_USE_MKL == 1
	#include "mkl.h"
	#include "mkl_lapacke.h"
#elif OPENSMOKE_USE_OPENBLAS == 1
	#include "cblas.h"
	#include "lapacke.h"
#endif

namespace OpenSMOKE
{
	template<typename T>
	T* OpenSMOKEMatrix<T>::GetHandle(void)
	{
		return matrix_[0] + 1;
	}

	template<typename T>
	inline const T* OpenSMOKEMatrix<T>::GetHandle() const
	{
		return matrix_[0] + 1;
	}



	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(void)
	{
		Initialize(0,0);
	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns)
	{
		Initialize(rows,columns);
	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(OpenSMOKEMatrix<T> const &rval)
	{
		Initialize(rval.Rows(), rval.Columns());
                if(this->numRows_ != 0)
			memcpy(this->matrix_[0]+1, rval.Matrix()[0]+1, numRows_*numColumns_ * sizeof(T));

	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const T a11, ...)
	{
		Initialize(rows, columns);
		T* w = matrix_[0] + 1;
		va_list pList;
		va_start(pList, a11);
			
		*w = a11;

		for(int i=2;i<size_;i++)
			*++w = va_arg(pList, T);

		va_end(pList);
	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const T* values)
	{
		Initialize(rows, columns);
		T* w = matrix_[0] + 1;
		if(numRows_ != 0)
			memcpy(w, values, (size_ - 1) * sizeof(T));
	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(OpenSMOKEVector<T> const& rval)
	{
		Initialize(rval.dimensions_, 1);
		if(numRows_ != 0)
			memcpy(matrix_[0]+1, rval.vector_+1, numRows_*numColumns_*sizeof(T));
	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const OpenSMOKEMatrix<T> &rval)
	{
		if(rows < 1 || rows > rval.numRows_)
			ErrorMessage("OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const OpenSMOKEMatrix<T> &rval) Error in requested dimensions (rows)");

		if(columns < 1 || columns > rval.numColumns_)
			ErrorMessage("OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const OpenSMOKEMatrix<T> &rval) Error in requested dimensions (columns)");

		Initialize(rows, columns);
		for(int row = 1;row < numRows_+1;row++)
			memcpy(matrix_[row] + 1, rval.matrix_[row] + 1, numColumns_ * sizeof(T));
	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const int irow, const int jcol, const OpenSMOKEMatrix<T> &rval)
	{
		if(rows < 1 || irow < 1 || irow > rval.numRows_ || rows > (rval.numRows_ - irow + 1))
			ErrorMessage("OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const int irow, const int jcol, const OpenSMOKEMatrix<T> &rval) Error in requested dimensions (rows)");
		
		if(columns < 1 || jcol < 1 || jcol > rval.numColumns_ || columns>(rval.numColumns_ - jcol + 1))
			ErrorMessage("OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const int rows, const int columns, const int irow, const int jcol, const OpenSMOKEMatrix<T> &rval) Error in requested dimensions (columns)");

		Initialize(rows,columns);
		for(int row = 1;row < numRows_+1;row++)
			memcpy(matrix_[row] + 1, rval.matrix_[row + irow - 1] + (jcol-1+1), numColumns_*sizeof(T));
	}

	template<typename T>
	OpenSMOKEMatrix<T>::OpenSMOKEMatrix(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
	{
		Load(fileName, fileFormat);
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::Load(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
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
	void OpenSMOKEMatrix<T>::Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			int rows, columns;
			fInput >> rows;
			fInput >> columns;

			if(rows < 1 || columns < 1)
				ErrorMessage("Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat) Wrong matrix dimensions");

			Initialize(rows, columns);
                        
                        T *w = matrix_[0] + 1;
			for(int i=1;i<size_;i++)
                        {
                            fInput >> *w;
                            ++w;
                        }
                         
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			int rows, columns;

			// Reading matrix size
			if(!fInput.read(reinterpret_cast < char * > (&rows), sizeof(rows)))
				ErrorMessage("I was unable to read the number of rows from binary file");
			if(!fInput.read(reinterpret_cast < char * > (&columns), sizeof(columns)))
				ErrorMessage("I was unable to read the number of columns from binary file");
			
			if(rows < 1 || columns < 1)
				ErrorMessage("Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat) Wrong matrix dimensions");

			// Initializing
			Initialize(rows,columns);

			
	//		for(int i=1;i<numRows_+1;i++)
	//			for(int j=1;j<numColumns_+1;j++)
	//				if(!fInput.read( (char*)&matrix_[i][j], sizeof(T) ))
	//					ErrorMessage("I was unable to read from binary file");

			// Reading elements
			if(!fInput.read(reinterpret_cast < char * > (matrix_[1]),sizeof(T)*(size_)))
				ErrorMessage("I was unable to read from binary file");
		}
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::Save(const std::string fileName, const OpenSMOKE_File_Format fileFormat)
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
	void OpenSMOKEMatrix<T>::Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fOutput << numRows_    << std::endl;
			fOutput << numColumns_ << std::endl;
                        
                        T *w = matrix_[0] + 1;
			for(int i=1;i<size_;i++)
                        {
			        fOutput << std::setprecision(16) << std::scientific << *w << std::endl;
                                ++w;
                        }
                        
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			// Writing matrix size
			if(!fOutput.write( reinterpret_cast < char * > (&numRows_), sizeof(numRows_)))
				ErrorMessage("I was unable to write the number of rows on binary file");
			if(!fOutput.write( reinterpret_cast < char * > (&numColumns_), sizeof(numColumns_)))
				ErrorMessage("I was unable to write the number of rows on binary file");

			// Writing vector elements
			for(int i=1;i<numRows_+1;i++)
				for(int j=1;j<numColumns_+1;j++)
					if(!fOutput.write( reinterpret_cast < char * >(&matrix_[i][j]), sizeof(T) ))
						ErrorMessage("I was unable to write the matrix values on binary file");
		}
	}

	template<typename T>
	inline void OpenSMOKEMatrix<T>::PrintOnVideo() const
	{
		std::cout << typeid(OpenSMOKEMatrix<T>).name() << std::endl;
		std::cout << "Counters:      " << whoAmI_ << "/" << countInScope_ << "/" << count_ << std::endl;
		std::cout << "Rows:          " << numRows_ << std::endl;
		std::cout << "Columns:       " << numColumns_ << std::endl;

		for(int i=1;i<numRows_+1;i++)
		{
			for(int j=1;j<numColumns_+1;j++)
				std::cout << std::setw(16) << matrix_[i][j];
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<typename T>
	void Swap(OpenSMOKEMatrix<T> *lval, OpenSMOKEMatrix<T> *rval)
	{
		T **temp = lval->matrix_;
		lval->matrix_ = rval->matrix_;
		rval->matrix_ = temp;
		Swap(&lval->numColumns_, &rval->numColumns_);
		Swap(&lval->numRows_, &rval->numRows_);
		Swap(&lval->size_, &rval->size_);
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::Initialize(const int rows, const int columns)
	{
		count_++;
		countInScope_++;
		whoAmI_ = count_;
		if(rows < 1 || columns < 1)
		{
			numRows_ = numColumns_ = size_ = 0;
			matrix_  = 0;
			return;
		}

		numRows_ = rows;
		numColumns_ = columns;
		size_ = numRows_*numColumns_ + 1;

		matrix_ = new T *[numRows_+1];
	
		if(!matrix_)
			ErrorMessage("void OpenSMOKEMatrix<T>::Initialize(const int rows, const int columns) Memory Allocation Failure");

		matrix_[0] = new T[size_];

		if(!matrix_[0])
			ErrorMessage("void OpenSMOKEMatrix<T>::Initialize(const int rows, const int columns) Memory Allocation Failure");

		if (1 == 1)
		{
			matrix_[1] = matrix_[0];
			for(int i=2;i<=numRows_;i++)
				matrix_[i] = matrix_[i-1]+numColumns_;
		}
		else
		{

			for(int i=1; i<numRows_;i++)
				matrix_[i] = matrix_[i-1]+numColumns_;
		}

		ChangeDimensions ( std::min(numRows_, numColumns_), &ipiv, false);
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::operator = (const T c)
	{
		if(c == 0.)
			memset(matrix_[0],0,size_*sizeof(T));
		else
		{
			T *w = matrix_[0] + 1;
			for(int i = 1;i <= size_-1;i++)
				*w++ = c;
		}
	}
        
        // Assignment operator
        template<typename T>
	OpenSMOKEMatrix<T>& OpenSMOKEMatrix<T>::operator =(OpenSMOKEMatrix<T> const& orig)
        {            
            CopyPreparation(orig.numRows_, orig.numColumns_);

            if(this->numRows_ != 0)
			memcpy(this->matrix_[0]+1, orig.Matrix()[0]+1, numRows_*numColumns_ * sizeof (T));

			return *this;
        }    

	template<typename T>
	T* OpenSMOKEMatrix<T>::operator [] (const int i)
	{
		return matrix_[i];
	}

	template<typename T>
	const T* OpenSMOKEMatrix<T>::operator [] (const int i) const
	{
		return matrix_[i];
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::Deinitialize(void)
	{
            if(matrix_ == 0) return;
            
             delete[]  matrix_[0];
             delete[]  matrix_;
	}

	template<typename T>
	OpenSMOKEMatrix<T>::~OpenSMOKEMatrix(void)
	{
            Deinitialize();
            matrix_ = 0;
            size_ = 0;
            countInScope_--;
	}

	template<typename T>
	int OpenSMOKEMatrix<T>::Rows(void) const
	{
		return numRows_;
	}

	template<typename T>
	inline int OpenSMOKEMatrix<T>::Columns(void) const
	{
		return numColumns_;
	}

	template<typename T>
	inline int OpenSMOKEMatrix<T>::WhoAmI(void) const
	{
		return whoAmI_;
	}

	template<typename T>
	inline T** OpenSMOKEMatrix<T>::Matrix(void) const
	{
		return matrix_;
	}
	
	template<typename T>
	void OpenSMOKEMatrix<T>::GetRow(const int i, OpenSMOKEVector<T> *v)	
	{
		if( (i<1) || (i>numRows_) )
			ErrorMessage("void OpenSMOKEMatrix<T>::GetRow(const int i, OpenSMOKEVector<T> *v) Row out of indices");

		ChangeDimensions(numColumns_, v, false);
		memcpy(v->vector_ + 1, matrix_[i] + 1, numColumns_*sizeof(T));
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::GetColumn(const int i, OpenSMOKEVector<T> *v)
	{
			if( (i < 1) || (i >= (numColumns_+1)) )
				ErrorMessage("void OpenSMOKEMatrix<T>::GetColumn(const int i, OpenSMOKEVector<T> *v) Column out of indices");
		
			ChangeDimensions(numRows_, v, false);
			T *r = v->vector_+1;
			for(int k=1;k<numRows_+1;k++)
				*r++ = matrix_[k][i];
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::GetDiagonal(const int i, OpenSMOKEVector<T>  *v)
	{
		if(numRows_ != numColumns_)
			ErrorMessage("OpenSMOKEVector<T> OpenSMOKEMatrix<T>::GetDiagonal(const int i) const Non square matrix");

		if( i < -numRows_ + 1 || i > numColumns_ - 1 )
			ErrorMessage("OpenSMOKEVector<T> OpenSMOKEMatrix<T>::GetDiagonal(const int i) const Requested Diagonal is out of bounds");
		
		const int diag = numRows_ - abs(i);
	
		ChangeDimensions(diag, v, false);

		T *r = v->vector_+1;

		if(i >= 0)
		{
			for(int k=1; k<diag+1;k++)
				*r++ = matrix_[k][k+i];
		}
		else
		{
			for(int k=1; k<diag+1;k++)
				*r++ = matrix_[k-i][k];
		}
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::SetRow(const int j, const OpenSMOKEVector<T> &rval)
	{
		if( (j < 1 || j > numRows_) || (rval.dimensions_ != numColumns_) )
			ErrorMessage("void OpenSMOKEMatrix<T>::SetRow(const int j, const OpenSMOKEVector<T> &rval) Wrong Dimensions");
	
		memcpy(matrix_[j] + 1, rval.vector_ + 1, numColumns_ * sizeof(T));
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::SetRow(const int j, const T rval)
	{
		if( (j < 1) || (j > numRows_)  )
			ErrorMessage("void OpenSMOKEMatrix<T>::SetRow(const int j, const T rval) Wrong Dimensions");
	
		for(int k = 1; k<numColumns_+1;k++)
			matrix_[j][k] = rval;
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::SetColumn(const int j, const OpenSMOKEVector<T> &rval)
	{
		if( (j < 1 || j > numColumns_-1+1) || (rval.dimensions_ != numRows_) )
			ErrorMessage("void OpenSMOKEMatrix<T>::SetColumn(const int j, const OpenSMOKEVector<T> &rval) Wrong Dimensions");
	
		T *r = rval.vector_+1;
		for(int k = 1; k<numRows_+1;k++)
			matrix_[k][j] = *r++;
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::SetColumn(const int j, const T rval)
	{
		if( (j < 1) || (j > (numColumns_-1+1)) )
			ErrorMessage("void OpenSMOKEMatrix<T>::SetColumn(const int j, const T rval) Wrong Dimensions");
	
		for(int k = 1; k<numRows_+1;k++)
			matrix_[k][j] = rval;
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::SetMatrix(const T rval)
	{
		if(rval == 0)
			memset(matrix_[0], 0 ,size_*sizeof(T));
		else
		{
			T *w = matrix_[0] + 1;
			for(int i=1;i<=size_-1;i++)
				*w++ = rval;
		}
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::InsertRow(const int i, OpenSMOKEVector<T> &v)
	{
		if(numColumns_ == 0 && numRows_ == 0)
		{
			numColumns_ = v.Size();
			OpenSMOKEMatrix<T> A(i+1-1,numColumns_);
			memmove(A.matrix_[i] + 1,v.vector_ + 1, numColumns_*sizeof(T));
			Swap(&A,this);
			return;
		}

		if(numColumns_ != v.Size())
			ErrorMessage("void OpenSMOKEMatrix<T>::InsertRow(const int i, OpenSMOKEVector<T> &v) numColumns_ != v.Size()");
	
		const int n = numRows_ + 1;
		OpenSMOKEMatrix<T> A(n,numColumns_);

		if(i <= 1)
		{
			memmove(A.matrix_[1] + 1,v.vector_ + 1, numColumns_*sizeof(T));
			memmove(A.matrix_[1+1] + 1, matrix_[1] + 1,numColumns_*numRows_*sizeof(T));	
		}
		else if(i > numRows_-1+1)
		{
			memmove(A.matrix_[0], matrix_[0], size_*sizeof(T));		
			memmove(A.matrix_[n+1-1] + 1,v.vector_ + 1,numColumns_*sizeof(T));
		}
		else
		{
			memmove(A.matrix_[0], matrix_[0], numColumns_*(i+1-1)*sizeof(T));
			memmove(A.matrix_[i] + 1,v.vector_ + 1, numColumns_*sizeof(T));
			memmove(A.matrix_[i+1] + 1,matrix_[i] + 1, numColumns_*(numRows_ - (i+1-1) + 1)*sizeof(T));
		}
		Swap(&A,this);
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::AppendRow(OpenSMOKEVector<T> &v)
	{
		InsertRow(numRows_+1,v);
	}
	
	template<typename T>
	void OpenSMOKEMatrix<T>::RowsSum(OpenSMOKEVector<T> *sumRows)
	{
		T *a = matrix_[0] + 1;
		if (numRows_ != sumRows->dimensions_)
			ChangeDimensions(numRows_, sumRows, false);

		for (int row = 1; row < numRows_ + 1; row++)
		{
			(*sumRows)[row] = static_cast<T>(0);
			for (int j = 1; j <= numColumns_; j++)
			{
				(*sumRows)[row] += *a;
				a++;
			}
		}
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::ColumnsSum(OpenSMOKEVector<T> *sumColumns)
	{
		T *a = matrix_[0] + 1;
		if (numColumns_ != sumColumns->dimensions_)
			ChangeDimensions(numColumns_, sumColumns, false);

		for (int i = 1; i <= numColumns_; i++)
			(*sumColumns)[i] = static_cast<T>(0);

		for (int i = 1; i <= numRows_; i++)
			for (int column = 1; column < 1 + numColumns_; column++)
			{
				(*sumColumns)[column] += *a;
				a++;
			}
	}

	template<typename T>
	void OpenSMOKEMatrix<T>::CopyPreparation(const int rRows, const int rColumns)
	{
		int who = whoAmI_;
		if(numRows_ != rRows || numColumns_ != rColumns)
		{
			Deinitialize();
			Initialize(rRows, rColumns);
			count_--;
			countInScope_--;
		}
		whoAmI_ = who;
	}

	template<typename T>
	void ChangeDimensions(const int rows, const int columns, OpenSMOKEMatrix<T>* result, bool reset)
	{
		result->CopyPreparation(rows,columns);
		if(reset == true && result->size_ > 0)
			memset(result->matrix_[0], 0, result->size_ * sizeof(T));
	}


	template<typename T>
	void Add(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != B.numRows_ || A.numColumns_ != B.numColumns_)
			ErrorMessage(	"void Add(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Add(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");

		if ((C->whoAmI_ == A.whoAmI_) || (C->whoAmI_ == B.whoAmI_))
			ErrorMessage("void Add(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",
						 "The three matrices have to be different objects.");

		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
            const T* ptB = B.GetHandle();
			      T* ptC = C->GetHandle();

            vdAdd( n, ptA, ptB, ptC );
		}
        #else
		{
			Sum(B.size_ - 1, A.matrix_[1] + 1, B.matrix_[1] + 1, C->matrix_[1] + 1);
		}
        #endif                                    
	}

	template<typename T>
	void Sub(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != B.numRows_ || A.numColumns_ != B.numColumns_)
			ErrorMessage(	"void Sub(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Sub(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");

		if ((C->whoAmI_ == A.whoAmI_) || (C->whoAmI_ == B.whoAmI_))
			ErrorMessage("void Sub(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",
						 "The three matrices have to be different objects.");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
            const T* ptB = B.GetHandle();
			      T* ptC = C->GetHandle();

            vdSub( n, ptA, ptB, ptC );
		}
		#else
		{
			Difference(B.size_ - 1, A.matrix_[1] + 1, B.matrix_[1] + 1, C->matrix_[1] + 1);
		}
		#endif    
	}

	template<typename T>
	void ElementByElementProduct(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != B.numRows_ || A.numColumns_ != B.numColumns_)
			ErrorMessage(	"void ElementByElementProduct(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void ElementByElementProduct(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
                        const T* ptB = B.GetHandle();
			      T* ptC = C->GetHandle();

                        vdMul( n, ptA, ptB, ptC );
		}
		#endif
	}

	template<typename T>
	void ElementByElementDivision(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != B.numRows_ || A.numColumns_ != B.numColumns_)
			ErrorMessage(	"void ElementByElementDivision(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void ElementByElementDivision(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
                        const T* ptB = B.GetHandle();
			      T* ptC = C->GetHandle();

                        vdDiv( n, ptA, ptB, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void ElementByElementDivision(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)");
                #endif   
	}

	template<typename T>
	void Sqr(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Sqr(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdSqr( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Sqr(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif  
	}

	template<typename T>
	void Abs(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Abs(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdAbs( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Abs(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif  
	}

	template<typename T>
	void Inv(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Inv(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdInv( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Inv(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif 
	}

	template<typename T>
	void Sqrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Sqrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdSqrt( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Sqrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void InvSqrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void InvSqr(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdInvSqrt( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void InvSqrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Cbrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Cbrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdCbrt( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Cbrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void InvCbrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void InvCbrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdInvCbrt( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void InvCbrt(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Pow2o3(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Pow2o3(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdPow2o3( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Pow2o3(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Pow3o2(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Pow3o2(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdPow3o2( n, ptA, ptC);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Pow3o2(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Pow(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != B.numRows_ || A.numColumns_ != B.numColumns_)
			ErrorMessage(	"void Pow(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Pow(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
                        const T* ptB = B.GetHandle();
			      T* ptC = C->GetHandle();

                        vdPow( n, ptA, ptB, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Pow(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Pow(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Pow(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdPowx( n, ptA, b, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Pow(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Hypot(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != B.numRows_ || A.numColumns_ != B.numColumns_)
			ErrorMessage(	"void Hypot(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Hypot(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
                        const T* ptB = B.GetHandle();
			      T* ptC = C->GetHandle();

                        vdHypot( n, ptA, ptB, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Hypot(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C)");
                #endif
	}


	// Exponential and Logarithmic Functions

	template<typename T>
	void Exp(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Exp(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        vdExp( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Exp(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void ExpMinus1(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  ExpMinus1(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdExpm1( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void ExpMinus1(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Ln(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Ln(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdLn( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Ln(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Log10(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Log10(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdLog10( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Log10(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void LnPlus1(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  LnPlus1(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdLog1p( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void LnPlus1(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	// Trigonometric Functions

	template<typename T>
	void Cos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Cos(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdCos( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Cos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Sin(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Sin(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdSin( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Sin(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Tan(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Tan(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdTan( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Tan(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Acos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Acos(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdAcos( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Acos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Asin(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Asin(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdAsin( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Asin(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Atan(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Atan(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdAtan( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Atan(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void SinCos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *Sine, OpenSMOKEMatrix<T> *Cosine)
	{
		if(A.numRows_ != Sine->numRows_ || A.numColumns_ != Sine->numColumns_)
			ErrorMessage(	"void  SinCos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *Sine, OpenSMOKEMatrix<T> *Cosine)",  
							"Matrix dimension check failure");

		if(A.numRows_ != Cosine->numRows_ || A.numColumns_ != Cosine->numColumns_)
			ErrorMessage(	"void  SinCos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *Sine, OpenSMOKEMatrix<T> *Cosine)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptSine = Sine->GetHandle();
			      T* ptCosine = Cosine->GetHandle();

			vdSinCos( n, ptA, ptSine, ptCosine );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void SinCos(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *Sine, OpenSMOKEMatrix<T> *Cosine)");
                #endif
	}


	template<typename T>
	void Atan2(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> const &B, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != B.numRows_ || A.numColumns_ != B.numColumns_)
			ErrorMessage(	"void  Atan2(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T_> const &B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");

		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Atan2(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T_> const &B, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			const T* ptB = B.GetHandle();
			      T* ptC = C->GetHandle();

			vdAtan2( n, ptA, ptB, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Atan2(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> const &B, OpenSMOKEMatrix<T> *C)");
                #endif
	}


	// Hyperbolic Functions

	template<typename T>
	void Cosh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Cosh(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdCosh( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Cosh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Sinh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Sinh(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdSinh( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Sinh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Tanh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Tanh(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdTanh( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Tanh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Acosh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Acosh(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdAcosh( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Acosh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Asinh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Asinh(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdAsinh( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Asinh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Atanh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Atanh(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdAtanh( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Atanh(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}


	// Special functions

	template<typename T>
	void Erf(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Erf(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdErf( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Erf(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Erfc(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Erfc(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdErfc( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Erfc(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void ErfInv(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  ErfInv(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdErfInv( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void ErfInv(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void ErfcInv(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  ErfcInv(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdErfcInv( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void ErfcInv(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void CdfNorm(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  CdfNorm(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdCdfNorm( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void CdfNorm(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void CdfNormInv(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  CdfNormInv(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdCdfNormInv( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void CdfNormInv(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void LnGamma(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  LnGamma(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdLGamma( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void LnGamma(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Tgamma(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  TGamma(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdTGamma( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Tgamma(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	// Rounding functions

	template<typename T>
	void Floor(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Floor(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdFloor( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Floor(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Ceil(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Ceil(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdCeil( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Ceil(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Trunc(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Trunc(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdTrunc( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Trunc(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}

	template<typename T>
	void Round(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void  Round(OpenSMOKEMatrix<T> const& A, const double b, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

			vdRound( n, ptA, ptC );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Round(OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}


	
	template<typename T>
	void Product(const T alpha, OpenSMOKEMatrix<T> *C)
	{
		#if OPENSMOKE_USE_MKL == 1
		{
			int one = 1;
			int n = C->numRows_*C->numColumns_;
			T* ptC = C->GetHandle();

                        dscal( &n, &alpha, ptC, &one );
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Product(const T alpha, OpenSMOKEMatrix<T> *C)");
                #endif
	}


	template<typename T>
	void Product(const T alpha, OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)
	{
		if(A.numRows_ != C->numRows_ || A.numColumns_ != C->numColumns_)
			ErrorMessage(	"void Product(const double alpha, OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)",  
							"Matrix dimension check failure");
	
		#if OPENSMOKE_USE_MKL == 1
		{
			int one = 1;
			int n = A.numRows_*A.numColumns_;
			const T* ptA = A.GetHandle();
			      T* ptC = C->GetHandle();

                        daxpy(&n, &alpha, ptA, &one, ptC, &one);
		}
		#else
                    BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Product(const T alpha, OpenSMOKEMatrix<T> const &A, OpenSMOKEMatrix<T> *C)");
                #endif
	}


	template<typename T>
	void Product(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C,
					const bool transposea, const bool transposeb, const T alpha, const T beta)
	{
            #if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)
            {
                            CBLAS_TRANSPOSE transa;
                            CBLAS_TRANSPOSE transb;
                            const int lda = A.Columns();
                            const int ldb = B.Columns();
                            const int ldc = C->Columns();

                            int m;
                            int n;
                            int k;

                            if (transposea == false && transposeb == false)
                            {
                                    transa=CblasNoTrans;
                                    transb=CblasNoTrans;

                                    m = C->Rows();
                                    n = C->Columns();
                                    k = A.Columns();

                                    if ( m != A.Rows() || k != B.Rows() || n != B.Columns() )
                                            ErrorMessage("Product", "Product Dimension Error (NN)");
                            }
                            else if (transposea == true && transposeb == false)
                            {
                                    transa=CblasTrans;
                                    transb=CblasNoTrans;

                                    m = A.Columns();
                                    n = B.Columns();
                                    k = A.Rows();

                                    if ( k != B.Rows() || m != C->Rows() || n != C->Columns() )
                                            ErrorMessage("Product", "Product Dimension Error (TN)");
                            }
                            else if (transposea == false && transposeb == true)
                            {
                                    transa=CblasNoTrans;
                                    transb=CblasTrans;

                                    m = A.Rows();
                                    n = B.Rows();
                                    k = A.Columns();

                                    if ( m != C->Rows() || k != B.Columns() || n != C->Columns() )
                                            ErrorMessage("Product", "Product Dimension Error (NT)");
                            }
                            else if (transposea == true && transposeb == true)
                            {	
                                    transa=CblasTrans;
                                    transb=CblasTrans;

                                    m = A.Columns();
                                    n = B.Rows();
                                    k = A.Rows();

                                    if ( m != C->Rows() || k != B.Columns() || n != C->Columns() )
                                            ErrorMessage("Product", "Product Dimension Error (TT)");
                            }

                            const T* ptA = A.GetHandle();
                            const T* ptB = B.GetHandle();
                                      T* ptC = C->GetHandle();

                            cblas_dgemm(CblasRowMajor, transa, transb, m, n, k, alpha, ptA, lda, ptB, ldb, beta, ptC, ldc); 

            }
            #else
                BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Product(OpenSMOKEMatrix<T> const& A, OpenSMOKEMatrix<T> const& B, OpenSMOKEMatrix<T> *C, const bool transposea, const bool transposeb, const T alpha, const T beta)");
            #endif
	}

	template<typename T>
	void Product(OpenSMOKEMatrix<T> const& A, OpenSMOKEVector<T> const& x, OpenSMOKEVector<T> *y,
					const bool transposea, const T alpha, const T beta)
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)
		{
			CBLAS_TRANSPOSE transa = CblasNoTrans;
			const int lda = A.Columns();
			
			const int m = A.Rows();
			const int n = A.Columns();

			if (y->Size() != m || x.Size() != n)
				ErrorMessage("Matrix-Vector Product", "Product Dimension Error");

			if (transposea == true)
				transa=CblasTrans;

			const T* ptA = A.GetHandle();
			const T* ptX = x.GetHandle();
				  T* ptY = y->GetHandle();

			cblas_dgemv(CblasRowMajor, transa, m, n, alpha, ptA, lda, ptX, 1, beta, ptY, 1); 
		}
		#else
            BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void Product(OpenSMOKEMatrix<T> const& A, OpenSMOKEVector<T> const& x, OpenSMOKEVector<T> *y, const bool transposea, const T alpha, const T beta)");
        #endif
	}

	template<typename T>
	void FactorizationLU(OpenSMOKEMatrix<T> *A)
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)
		{
			
			int lda = A->Columns();

			T* ptA      = A->GetHandle();
			int* ptIpiv = A->ipiv.GetHandle();

			int status =  LAPACKE_dgetrf( CblasRowMajor, A->Rows(), A->Columns(), ptA, lda, ptIpiv );

			if (status != 0)
			{
				std::stringstream number;
				number << abs(status);
				
				if (status<0)
				{
					std::string message = "The " + number.str() + "-th had an illegal value.";
					ErrorMessage("FactorizationLU", message);
				}
				else
				{
					std::string message = "FactorizationLU: The " + number.str() + "-th element on the main diagonal is equal to 0. U is exactly singular.";
					ErrorMessage("FactorizationLU", message);
				}
			}
		}
        #else
            BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void FactorizationLU(OpenSMOKEMatrix<T> *A)");
        #endif
	}

	template<typename T>
	void SolveLU(OpenSMOKEMatrix<T>const& A, OpenSMOKEVector<T> *b, const bool transposea)
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)
		{
			char transa = 'N';
			if (transposea == true)
				transa = 'T';

			int lda = A.Columns();

			T* ptB			  = b->GetHandle();
			const T* ptA      = A.GetHandle();
			const int* ptIpiv = A.ipiv.GetHandle();

			int status = LAPACKE_dgetrs( CblasRowMajor, transa, A.Rows(), 1, ptA, lda, ptIpiv, ptB, 1 );

			if (status < 0)
			{
				std::stringstream number;
				number << abs(status);
					
				std::string message = "The " + number.str() + "-th parameter had an illegal value.";
				ErrorMessage("SolveLU", message);
			}
		}
		#else
            BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Not available operation: void SolveLU(OpenSMOKEMatrix<T>const& A, OpenSMOKEVector<T> *b, const bool transposea)");
        #endif
	}

	template class OpenSMOKEMatrix<unsigned int>;
	template class OpenSMOKEMatrix<int>;
	template class OpenSMOKEMatrix<float>;
	template class OpenSMOKEMatrix<double>;

	template void ChangeDimensions(const int rows, const int columns, OpenSMOKEMatrix<unsigned int>* result, bool reset);
	template void ChangeDimensions(const int rows, const int columns, OpenSMOKEMatrix<int>* result, bool reset);
	template void ChangeDimensions(const int rows, const int columns, OpenSMOKEMatrix<float>* result, bool reset);
	template void ChangeDimensions(const int rows, const int columns, OpenSMOKEMatrix<double>* result, bool reset);

}
