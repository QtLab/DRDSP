#ifndef INCLUDED_DYNAMICS_MODEL
#define INCLUDED_DYNAMICS_MODEL
/*
|
| Model
|
*/
#include "embedding.h"
#include <DRDSP/dynamics/dynamicalSystem.h>

namespace DRDSP {

	struct Model {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension;

		Model() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0) {}
		explicit Model( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0) {}
		explicit Model( uint32_t dim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim) {}
		Model( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim ) : wrapFunction(wrapFunc), dimension(dim) {}
		virtual VectorXd VectorField( const VectorXd& state ) = 0;
		virtual MatrixXd Partials( const VectorXd& state ) = 0;
	};
	
	struct ModelCW {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension;

		ModelCW() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0) {}
		explicit ModelCW( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0) {}
		explicit ModelCW( uint32_t dim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim) {}
		ModelCW( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim ) : wrapFunction(wrapFunc), dimension(dim) {}
		virtual double VectorField( const VectorXd& state, uint32_t i ) = 0;
		virtual double Partials( const VectorXd& state, uint32_t i, uint32_t j ) = 0;
	};

	struct ModelParameterized {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension;
		uint8_t parameterDimension;

		ModelParameterized() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0), parameterDimension(0) {}
		explicit ModelParameterized( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0), parameterDimension(0) {}
		ModelParameterized( uint32_t dim, uint32_t pDim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim), parameterDimension(pDim) {}
		ModelParameterized( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim, uint32_t pDim ) : wrapFunction(wrapFunc), dimension(dim), parameterDimension(pDim) {}
		virtual VectorXd VectorField( const VectorXd& state, const VectorXd& parameter ) = 0;
		virtual MatrixXd Partials( const VectorXd& state, const VectorXd& parameter ) = 0;
	};
	
	struct ModelParameterizedCW {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension;
		uint8_t parameterDimension;

		ModelParameterizedCW() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0), parameterDimension(0) {}
		explicit ModelParameterizedCW( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0), parameterDimension(0) {}
		ModelParameterizedCW( uint32_t dim, uint32_t pDim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim), parameterDimension(pDim) {}
		ModelParameterizedCW( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim, uint32_t pDim ) : wrapFunction(wrapFunc), dimension(dim), parameterDimension(pDim) {}
		virtual double VectorField( const VectorXd& state, const VectorXd& parameter, uint32_t i ) = 0;
		virtual double Partials( const VectorXd& state, const VectorXd& parameter, uint32_t i, uint32_t j ) = 0;
	};

	struct ModelEmbedded {
		Model& model;
		Embedding& embedding;

		ModelEmbedded( Model& m, Embedding& e );
		VectorXd VectorField( const VectorXd& state );
		MatrixXd Partials( const VectorXd& state );
	};
	
	struct ModelEmbeddedCW {
		ModelCW& model;
		EmbeddingCW& embedding;

		ModelEmbeddedCW( ModelCW& m, EmbeddingCW& e );
		double VectorField( const VectorXd& state, uint32_t i );
		double Partials( const VectorXd& state, uint32_t i, uint32_t j );
	};

	struct ModelParameterizedEmbedded {
		ModelParameterized& model;
		Embedding& embedding;

		ModelParameterizedEmbedded( ModelParameterized& m, Embedding& e );
		VectorXd VectorField( const VectorXd& state, const VectorXd& parameter );
		MatrixXd Partials( const VectorXd& state, const VectorXd& parameter );
	};
	
	struct ModelParameterizedEmbeddedCW {
		ModelParameterizedCW& model;
		EmbeddingCW& embedding;

		ModelParameterizedEmbeddedCW( ModelParameterizedCW& m, EmbeddingCW& e );
		double VectorField( const VectorXd& state, const VectorXd& parameter, uint32_t i );
		double Partials( const VectorXd& state, const VectorXd& parameter, uint32_t i, uint32_t j );
	};

}

#endif

