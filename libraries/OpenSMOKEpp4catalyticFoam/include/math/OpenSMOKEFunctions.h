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

#ifndef OpenSMOKE_OpenSMOKEFunctions_Hpp
#define OpenSMOKE_OpenSMOKEFunctions_Hpp

#include "OpenSMOKEStdInclude.h"
#include <boost/filesystem.hpp>
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	/**
	* Utility to print fatal error messages
	*/
	void ErrorMessage(const std::string functionName, const std::string errorMessage);

	/**
	* Utility to print fatal error messages
	*/
	int FatalErrorMessage(const std::string errorMessage);

	/**
	* Returns the MachEps (float)
	*/
	float MachEpsFloat();
	
	/**
	* Returns the MachEps (double)
	*/
	double MachEps();

	/**
	* Returns the square root of the sum of the squares of x
	* Example: z=SqrtSumSqr(n,&x) means \f$ z=todo \f$
	*/
	double SqrtSumSqr(const int n, int *x);

	/**
	* Returns the square root of the sum of the squares of x
	* Example: z=SqrtSumSqr(n,&x) means \f$ z=todo \f$
	*/
	float SqrtSumSqr(const int n, float *x);

	/**
	* Returns the square root of the sum of the squares of x
	* Example: z=SqrtSumSqr(n,&x) means \f$ z=todo \f$
	*/
	double SqrtSumSqr(const int n, double *x);

	/**
	* Returns the n-th power of x
	* Example: z=PowInt(x,n) means \f$ z=x^n \f$
	*/
	double PowInt(const double x, const int n);

	/**
	* Returns the cubic root of x
	* Example: z=CubicRoot(x) means \f$ z=x^(1/3) \f$
	*/
	double CubicRoot(const double x);
        
    /**
	* Returns the elapsed time in seconds since the process started
	*/
    double OpenSMOKEClock(void);

    /**
	* Returns the elapsed time in seconds since the OpenSMOKEDiffClock() started
	*/
	double OpenSMOKEGetCpuTime(void);

    /**
	* Returns the number of cpu clocks
	*/
	#if defined(_WIN32) || defined(_WIN64) 
		unsigned __int64 OpenSMOKEGetCpuClocks(void);
	#else
		unsigned long int OpenSMOKEGetCpuClocks(void);
	#endif

	/**
	* Returns the cpu frequency
	*/
    double OpenSMOKEGetCpuFrequency(void);

    /**
	* Returns the max cpu frequency
	*/
    double OpenSMOKEGetMaxCpuFrequency(void);

    /**
	* Returns the max cpu clock frequency
	*/
    double OpenSMOKEGetCpuClocksFrequency(void);

	std::string GetCurrentTime();
	std::string GetCurrentDate();
	void OpenSMOKE_logo(const std::string application_name, const std::string author_name);
	std::string SplitStringIntoSeveralLines(std::string source, const std::size_t width = 60, const std::string whitespace = " \t\r");
	void PrintTagOnASCIILabel(const unsigned int width, std::ostream &fOut, const std::string tag, unsigned int &counter);
	void SetXMLFile(std::ostream& xml_string);
	void CreateDirectory(const boost::filesystem::path& path_output);
	void OpenSMOKE_CheckLicense(const boost::filesystem::path executable_folder, const std::string compilation_date_os, const long maximum_license_days);
	boost::filesystem::path GetExecutableFileName(char** argv);

	void OpenInputFileXML(rapidxml::xml_document<>& doc, std::vector<char>& xml_copy, const boost::filesystem::path& file_name);

	void OpenInputFileASCII(std::ifstream &fASCII, const boost::filesystem::path input_file_ascii);



	/**
	*@brief Open the output files and checks any error
	*@param fXML the output stream used
	*@param output_file_xml XML file where the output will be written
	*/
	void OpenOutputFileXML(std::ofstream &fXML, const boost::filesystem::path output_file_xml);

	/**
	*@brief Opens the output files and checks any error
	*@param fASCII the output stream used
	*@param output_file_ascii ASCII file where the output will be written
	*/
	void OpenOutputFileASCII(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii);

	/**
	*@brief Opens an existing output files, checks any error and continue writing on it
	*@param fASCII the output stream used
	*@param output_file_ascii ASCII file where the output will be written
	*/
	void OpenOutputFileASCII_Append(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii);


	void CheckKineticsFolder(const boost::filesystem::path& path_output);

	/**
	*@brief Returns the sparsity pattern of a tridiagonal block matrix
	*@param number_blocks number of blocks
	*@param block_size size of a single block
	*@param rows row indices of non zero elements (1 index based)
	*@param cols column indices of non zero elements (1 index based)
	*/

	void SparsityPatternTridiagonal(const unsigned int number_equations, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);
	void SparsityPatternPentadiagonal(const unsigned int number_equations, const unsigned int width, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

	void SparsityPatternBlock(const unsigned int number_blocks, const unsigned int block_size,
		const std::vector<unsigned int>& rows_single, const std::vector<unsigned int>& cols_single,
		std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

	/**
	*@brief Returns the number of lines of a text file
	*@param path path of the file
	*/
	unsigned int filelines(const boost::filesystem::path path);

	/**
	*@brief Returns the number of words in a stringstream
	*@param str stringstream to analyze
	*/
	int countwords(std::stringstream& str);

	/**
	*@brief Returns the real solution of a cubic equation
	*@param a third order coefficient
	*@param b second order coefficient
	*@param c first order coefficient
	*@param d zero order coefficient
	*/
	std::vector<double> cubicrootsreal(const double a, const double b, const double c, const double d);

	/**
	*@brief Returns a normalized vector (sum of elements = 1)
	*@param x vector to normalize
	*/
	std::vector<double> normalize(std::vector<double>& x);

	/**
	*@brief Checks whether sum of elements of a vector<double> is equal to 1, without elements < 0 or > 1
	*@param x vector to analyse
	*/
	void CheckSumOfFractions(std::vector<double>& x);

	/**
	*@brief Converts temperature into [K]
	*@param temperature temperature to convert
	*@param unit original unit of measure
	*/
	double SI_temperature(const double temperature, const std::string unit);

	/**
	*@brief Converts pressure into [Pa]
	*@param pressure pressure to convert
	*@param unit original unit of measure
	*/
	double SI_pressure(const double pressure, const std::string unit);

	/**
	*@brief Converts length into [m]
	*@param length length to convert
	*@param unit original unit of measure
	*/
	double SI_length(const double length, const std::string unit);

	/**
	*@brief Converts volume into [m3]
	*@param volume volume to convert
	*@param unit original unit of measure
	*/
	double SI_volume(const double volume, const std::string unit);

	/**
	*@brief Converts velocity into [m/s]
	*@param velocity velocity to convert
	*@param unit original unit of measure
	*/
	double SI_velocity(const double velocity, const std::string unit);

	/**
	*@brief Converts time into [s]
	*@param time time to convert
	*@param unit original unit of measure
	*/
	double SI_time(const double time, const std::string unit);

	/**
	*@brief Evaluates the median value of a vector<double>
	*@param vec vector to analyse
	*/
	double median(const std::vector<double> vec);

	/**
	*@brief Evaluates the median absolute deviation of a vector<double>
	*@param vec vector to analyse
	*/
	double mad(const std::vector<double> vec);

	static const float  OPENSMOKE_MACH_EPS_FLOAT = MachEpsFloat();
	static const double OPENSMOKE_MACH_EPS_DOUBLE = MachEps();
}

#endif	// OpenSMOKE_OpenSMOKEFunctions_Hpp

