#ifndef stefanMaxwell_H
#define stefanMaxwell_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace OpenSMOKE
{
	class stefanMaxwell
	{
		public:

			stefanMaxwell
			(
				ThermodynamicsMap_CHEMKIN&           thermodynamicsMap,
				TransportPropertiesMap_CHEMKIN&      transportMap
			);

			stefanMaxwell
			(
        			stefanMaxwell& c
			);

			
			virtual ~stefanMaxwell();

			void setComposition(const OpenSMOKEVectorDouble omega);
			
			void setComposition(Eigen::VectorXd omegaEigen);
			
			void setTemperature(const double T);

			void setPressure(const double P);
			
			void setInertIndex(const unsigned int inertIndex);
			
			unsigned int getInert();

			std::vector<double> getdijm();
			
			void solve();

			inline Eigen::MatrixXd getDij() {return Dij_;}

		private:

			ThermodynamicsMap_CHEMKIN&      thermodynamicsMap_;     //!< thermodynamic map
			TransportPropertiesMap_CHEMKIN& transportMap_;          //!< transport map

			// Binary Diffusivities from Kinetic Theory of Gases
			Eigen::MatrixXd Diff_;
			std::vector<double> dijm_;

			// Matrixes reprersenting the linear system derived from Stefan-Maxwell Equations (no Soret Effect)
			Eigen::MatrixXd A_;
			Eigen::MatrixXd B_;
			
			// Matrix of the Generalized Fick Coefficients
			Eigen::MatrixXd Dij_;

			// Molar and Massive Compisition of the cell 
			OpenSMOKEVectorDouble x_;
			OpenSMOKEVectorDouble omega_;
			
			// Molecular Weights
			OpenSMOKEVectorDouble MWi_;
			
			unsigned int NC_;
			unsigned int inertIndex_;

			double T_;
			double P_;
			double MW_;
	};
}

#include "stefanMaxwell.hpp"

#endif // DUSTYGAS_H
