#include "stefanMaxwell.h"

OpenSMOKE::stefanMaxwell::stefanMaxwell
(
	OpenSMOKE::ThermodynamicsMap_CHEMKIN&        thermodynamicsMap,
	OpenSMOKE::TransportPropertiesMap_CHEMKIN&   transportMap
) 
:
thermodynamicsMap_(thermodynamicsMap),
transportMap_(transportMap),
dijm_(0)
{
	NC_ = thermodynamicsMap_.NumberOfSpecies();
                
	Diff_.resize(NC_,NC_),
	A_.resize(NC_-1,NC_-1);
    B_.resize(NC_-1,NC_-1);
	Dij_.resize(NC_-1,NC_-1);
    
	ChangeDimensions(NC_, &x_,     true);
	ChangeDimensions(NC_, &omega_, true);
	ChangeDimensions(NC_, &MWi_, true);
	
	inertIndex_ = 10000;

	T_ = 0.;
	P_ = 0.;
	MW_ = 0.;
}

OpenSMOKE::stefanMaxwell::stefanMaxwell
(
        OpenSMOKE::stefanMaxwell& c
)
:
thermodynamicsMap_(c.thermodynamicsMap_),
transportMap_(c.transportMap_),
dijm_(0)
{
        NC_ = thermodynamicsMap_.NumberOfSpecies();

        Diff_.resize(NC_,NC_),
        A_.resize(NC_-1,NC_-1);
        B_.resize(NC_-1,NC_-1);
        Dij_.resize(NC_-1,NC_-1);

        ChangeDimensions(NC_, &x_,     true);
        ChangeDimensions(NC_, &omega_, true);
        ChangeDimensions(NC_, &MWi_, true);

        inertIndex_ = c.inertIndex_;

        T_ = 0.;
        P_ = 0.;
        MW_ = 0.;
}


void OpenSMOKE::stefanMaxwell::setComposition(const OpenSMOKE::OpenSMOKEVectorDouble omega)
{
    	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(),MW_,omega.GetHandle());

    	for (unsigned int i=1;i<=NC_;i++)
    	{
        	omega_[i] = omega[i];
        	MWi_[i] = thermodynamicsMap_.MW(i-1);
    	}
}

void OpenSMOKE::stefanMaxwell::setComposition(Eigen::VectorXd omegaEigen)
{
	OpenSMOKEVectorDouble omega;
	ChangeDimensions(NC_, &omega, true);
	
	for (unsigned int i=1;i<=NC_;i++)
    	{
        	omega[i] = omegaEigen[i-1];
    	}
    
    	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(),MW_,omega.GetHandle());

    	for (unsigned int i=1;i<=NC_;i++)
    	{
        	omega_[i] = omega[i];
        	MWi_[i] = thermodynamicsMap_.MW(i-1);
    	}
}

void OpenSMOKE::stefanMaxwell::setTemperature(const double T)
{
    	T_= T;
}

void OpenSMOKE::stefanMaxwell::setPressure(const double P)
{
    	P_= P;
}

void OpenSMOKE::stefanMaxwell::setInertIndex(const unsigned int inertIndex)
{
    	inertIndex_= inertIndex;
}

unsigned int OpenSMOKE::stefanMaxwell::getInert()
{
    	return inertIndex_;
}

std::vector<double>  OpenSMOKE::stefanMaxwell::getdijm()
{
	return dijm_;
}

void OpenSMOKE::stefanMaxwell::solve()
{
	OpenSMOKEVectorDouble binaryDiffCoeff_(thermodynamicsMap_.NumberOfSpecies()*(thermodynamicsMap_.NumberOfSpecies()-1)*0.5);

	//- Calculate binary diffusivities
	transportMap_.MassBinaryDiffusionCoefficients(binaryDiffCoeff_.GetHandle());
	dijm_.resize(binaryDiffCoeff_.Size());

	for(unsigned int i=1;i<=binaryDiffCoeff_.Size();i++)
	{
		binaryDiffCoeff_[i] = 1./binaryDiffCoeff_[i];
		dijm_[i-1]=binaryDiffCoeff_[i];
	}

	unsigned int counter = 1;

	for(unsigned int i=0;i<NC_;i++)
	{
		Diff_(i,i) = 0.;
		
		for(unsigned int j=(i+1);j<NC_;j++)
		{
			if (j!=i)
			{
				Diff_(i,j) = binaryDiffCoeff_(counter++);
				Diff_(j,i) = Diff_(i,j);
			}
		}
	}

	//- Calculate A
	{
	unsigned int compIndexi = 0;
	unsigned int compIndexj = 0;
	
		for(unsigned int i=0;i<NC_-1;i++)
		{
		double sumDiff = 0.;

		if (i<inertIndex_)
		{
			compIndexi = i;
		}
		else
		{
			compIndexi = i+1;
		}
		
		for(unsigned int j=0;j<NC_;j++)
		{
			if(j!=compIndexi)
			sumDiff = sumDiff + x_[j+1]/Diff_(compIndexi,j)*MW_/MWi_[compIndexi+1];
		}
		
				for(unsigned int j=0;j<NC_-1;j++)
				{
			if (j<inertIndex_)
			{
				compIndexj = j;
			}
			else
			{
				compIndexj = j+1;
			}
			
			if (j==i)
			{
				A_(i,j) = -(x_[compIndexi+1]/Diff_(compIndexi,inertIndex_)*MW_/MWi_[inertIndex_+1]+sumDiff);
			}
			else
			{
				A_(i,j) = x_[compIndexi+1]*
						(1/Diff_(compIndexi,compIndexj)*MW_/MWi_[compIndexj+1]
							-1/Diff_(compIndexi,inertIndex_)*MW_/MWi_[inertIndex_+1]);
			}
				}
		}
	}

   	//- Calculate B
   	{
		unsigned int compIndexi = 0;
		unsigned int compIndexj = 0;
		
       		for(unsigned int i=0;i<NC_-1;i++)
       		{
			if (i<inertIndex_)
			{
				compIndexi = i;
			}
			else
			{
				compIndexi = i+1;
			}
			
			for(unsigned int j=0;j<NC_-1;j++)
            		{
				if (j<inertIndex_)
				{
					compIndexj = j;
				}
				else
				{
					compIndexj = j+1;
				}
				
				if (j==i)
				{
					B_(i,j) = -(
						x_[compIndexi+1]*MW_/MWi_[inertIndex_+1]
					       +(1-x_[compIndexi+1])*MW_/MWi_[compIndexi+1]);
				}
				else
				{
					B_(i,j) = x_[compIndexi+1]
						*(MW_/MWi_[compIndexj+1]-MW_/MWi_[inertIndex_+1]);
				}
            		}
        	}
    	}

    	//- Linear System Resolution
	{
		Eigen::PartialPivLU<Eigen::MatrixXd> luDec(A_);
		Dij_ = luDec.solve(B_);
	}

}



OpenSMOKE::stefanMaxwell::~stefanMaxwell()
{}
