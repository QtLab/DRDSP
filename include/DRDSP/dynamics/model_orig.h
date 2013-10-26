#ifndef INCLUDED_DYNAMICS_MODEL_ORIG
#define INCLUDED_DYNAMICS_MODEL_ORIG
/*
|
| Original Model
|
*/
#include "embedding.h"

namespace DRDSP {

	struct Model {
		uint32_t dimension;

		Model( uint32_t dim ) : dimension(dim) {}
		virtual VectorXd VectorField( const VectorXd& state ) = 0;
		virtual MatrixXd Partials( const VectorXd& state ) = 0;
	};
	
	struct ModelCW {
		uint32_t dimension;

		ModelCW( uint32_t dim ) : dimension(dim) {}
		virtual double VectorField( const VectorXd& state, uint32_t i ) = 0;
		virtual double Partials( const VectorXd& state, uint32_t i, uint32_t j ) = 0;
	};

	struct ModelParameterized {
		uint32_t dimension;
		uint8_t parameterDimension;

		ModelParameterized( uint32_t dim, uint32_t pDim ) : dimension(dim), parameterDimension(pDim) {}
		virtual VectorXd VectorField( const VectorXd& state, const VectorXd& parameter ) = 0;
		virtual MatrixXd Partials( const VectorXd& state, const VectorXd& parameter ) = 0;
	};
	
	struct ModelParameterizedCW {
		uint32_t dimension;
		uint8_t parameterDimension;

		ModelParameterizedCW( uint32_t dim, uint32_t pDim ) : dimension(dim), parameterDimension(pDim) {}
		virtual double VectorField( const VectorXd& state, const VectorXd& parameter, uint32_t i ) = 0;
		virtual double Partials( const VectorXd& state, const VectorXd& parameter, uint32_t i, uint32_t j ) = 0;
	};

	struct ModelEmbedded {
		Model& model;
		Embedding& embedding;

		ModelEmbedded( Model& m, Embedding& e );
		VectorXd VectorField( const VectorXd& state ) const;
		MatrixXd Partials( const VectorXd& state ) const;
	};
	
	struct ModelEmbeddedCW {
		ModelCW& model;
		EmbeddingCW& embedding;

		ModelEmbeddedCW( ModelCW& m, EmbeddingCW& e );
		double VectorField( const VectorXd& state, uint32_t i ) const;
		double Partials( const VectorXd& state, uint32_t i, uint32_t j ) const;
	};

	struct ModelParameterizedEmbedded {
		ModelParameterized& model;
		Embedding& embedding;

		ModelParameterizedEmbedded( ModelParameterized& m, Embedding& e );
		VectorXd VectorField( const VectorXd& state, const VectorXd& parameter ) const;
		MatrixXd Partials( const VectorXd& state, const VectorXd& parameter ) const;
	};
	
	struct ModelParameterizedEmbeddedCW {
		ModelParameterizedCW& model;
		EmbeddingCW& embedding;

		ModelParameterizedEmbeddedCW( ModelParameterizedCW& m, EmbeddingCW& e );
		double VectorField( const VectorXd& state, const VectorXd& parameter, uint32_t i ) const;
		double Partials( const VectorXd& state, const VectorXd& parameter, uint32_t i, uint32_t j ) const;
	};

}

#endif

