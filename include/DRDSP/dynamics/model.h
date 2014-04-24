#ifndef INCLUDED_DYNAMICS_MODEL
#define INCLUDED_DYNAMICS_MODEL
#include "embedding.h"
#include <DRDSP/dynamics/dynamicalSystem.h>

namespace DRDSP {

	/*!
	 * \brief Interface for a model without parameters.
	 */
	struct Model {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension; //!< Dimension of the model's state space

		Model() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0) {}
		explicit Model( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0) {}
		explicit Model( uint32_t dim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim) {}
		Model( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim ) : wrapFunction(wrapFunc), dimension(dim) {}

		
		virtual VectorXd VectorField( const VectorXd& state ) = 0; //!< The vector field
		
		virtual MatrixXd Partials( const VectorXd& state ) = 0; //!< The partial derivatives of the vector field
	};
	
	/*!
	 * \brief Interface for a model without parameters,
	 * whose derivatives are to be evaluated component-wise.
	 *
	 * This version of the model interface is for very high-dimensional systems
	 * whose matrix of partial derivatives is too large to store in memory.
	 * These systems will have their derivatives evaluated one element at a time.
	 */
	struct ModelCW {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension;

		ModelCW() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0) {}
		explicit ModelCW( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0) {}
		explicit ModelCW( uint32_t dim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim) {}
		ModelCW( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim ) : wrapFunction(wrapFunc), dimension(dim) {}
		virtual double VectorField( const VectorXd& state, uint32_t i ) = 0; //!< The ith component of the vector field
		virtual double Partials( const VectorXd& state, uint32_t i, uint32_t j ) = 0; //!< The partial derivative of the ith component of the vector field with respect to the jth element of the state
	};

	/*!
	 * \brief Interface for a model with parameters.
	 */
	struct ModelParameterized {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension; //!< Dimension of the state space
		uint8_t parameterDimension; //!< Dimension of the parameter space

		ModelParameterized() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0), parameterDimension(0) {}
		explicit ModelParameterized( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0), parameterDimension(0) {}
		ModelParameterized( uint32_t dim, uint32_t pDim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim), parameterDimension(pDim) {}
		ModelParameterized( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim, uint32_t pDim ) : wrapFunction(wrapFunc), dimension(dim), parameterDimension(pDim) {}
		virtual VectorXd VectorField( const VectorXd& state, const VectorXd& parameter ) = 0; //!< The vector field
		virtual MatrixXd Partials( const VectorXd& state, const VectorXd& parameter ) = 0; //!< The partial derivatives of the vector field (with respect to the state)
	};
	
	/*!
	 * \brief Interface for a model with parameters,
	 * whose derivatives are to be evaluated component-wise.
	 *
	 * This version of the model interface is for very high-dimensional systems
	 * whose matrix of partial derivatives is too large to store in memory.
	 * These systems will have their derivatives evaluated one element at a time.
	 */
	struct ModelParameterizedCW {
		const WrapFunction<VectorXd>* wrapFunction;
		uint32_t dimension; //!< Dimension of the state space
		uint8_t parameterDimension; //!< Dimension of the parameter space

		ModelParameterizedCW() : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(0), parameterDimension(0) {}
		explicit ModelParameterizedCW( const WrapFunction<VectorXd>* wrapFunc ) : wrapFunction(wrapFunc), dimension(0), parameterDimension(0) {}
		ModelParameterizedCW( uint32_t dim, uint32_t pDim ) : wrapFunction(&WrapFunction<VectorXd>::identity), dimension(dim), parameterDimension(pDim) {}
		ModelParameterizedCW( const WrapFunction<VectorXd>* wrapFunc, uint32_t dim, uint32_t pDim ) : wrapFunction(wrapFunc), dimension(dim), parameterDimension(pDim) {}
		virtual double VectorField( const VectorXd& state, const VectorXd& parameter, uint32_t i ) = 0; //!< The ith component of the vector field
		virtual double Partials( const VectorXd& state, const VectorXd& parameter, uint32_t i, uint32_t j ) = 0; //!< The partial derivative of the ith component of the vector field with respect to the jth element of the state
	};

	/*!
	 * \brief A model without parameters whose state space is embedded into R^n.
	 */
	struct ModelEmbedded {
		Model& model; //!< The underlying model
		Embedding& embedding; //!< An embedding into R^n

		ModelEmbedded( Model& m, Embedding& e );
		VectorXd VectorField( const VectorXd& state ); //!< The vector field on R^n
		MatrixXd Partials( const VectorXd& state ); //!< The partial derivative of the vector field on R^n with respect to the underlying state
	};

	/*!
	 * \brief A model without parameters whose state space is embedded into R^n,
	 * and whose derivates are to be evaluated component-wise.
	 */
	struct ModelEmbeddedCW {
		ModelCW& model; //!< The underlying model
		EmbeddingCW& embedding; //!< An embedding into R^n

		ModelEmbeddedCW( ModelCW& m, EmbeddingCW& e );
		double VectorField( const VectorXd& state, uint32_t i ); //!< The ith component of the vector field on R^n
		double Partials( const VectorXd& state, uint32_t i, uint32_t j ); //!< The partial derivative of the ith component of the vector field with respect to the jth element of the underlying state
	};

	/*!
	 * \brief A model with parameters whose state space is embedded into R^n.
	 */
	struct ModelParameterizedEmbedded {
		ModelParameterized& model; //!< The underlying model
		Embedding& embedding; //!< An embedding into R^n

		ModelParameterizedEmbedded( ModelParameterized& m, Embedding& e );
		VectorXd VectorField( const VectorXd& state, const VectorXd& parameter ); //!< The vector field on R^n
		MatrixXd Partials( const VectorXd& state, const VectorXd& parameter );    //!< The partial derivative of the vector field on R^n with respect to the underlying state
	};
	
	/*!
	 * \brief A model with parameters whose state space is embedded into R^n,
	 * and whose derivatives are to be evaluated component-wise.
	 */
	struct ModelParameterizedEmbeddedCW {
		ModelParameterizedCW& model; //!< The underlying model
		EmbeddingCW& embedding;      //!< An embedding into R^n

		ModelParameterizedEmbeddedCW( ModelParameterizedCW& m, EmbeddingCW& e );
		double VectorField( const VectorXd& state, const VectorXd& parameter, uint32_t i );          //!< The ith component of the vector field on R^n
		double Partials( const VectorXd& state, const VectorXd& parameter, uint32_t i, uint32_t j ); //!< The partial derivative of the ith component of the vector field with respect to the jth element of the underlying state
	};

}

#endif

