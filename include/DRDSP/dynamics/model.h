#ifndef INCLUDED_DYNAMICS_MODEL
#define INCLUDED_DYNAMICS_MODEL
#include "embedding.h"
#include <DRDSP/dynamics/dynamicalSystem.h>

namespace DRDSP {

	/*!
	 * \brief Wrapper for a Model that can be used as a SolverFunction
	 */
	template<typename Model>
	struct SolverFunctionFromModel {
		typedef typename Model::State State;
		Model model;
		explicit SolverFunctionFromModel( const Model& model ) : model(model) {}
		template<typename Time>
		State operator()( const State& x, Time t ) const {
			return model(x);
		}
	};

	/*!
	 * \brief Base for model without parameters.
	 */
	template<typename State = VectorXd>
	struct Model {
		typedef State State;
		uint32_t dimension; //!< Dimension of the model's state space

		Model() : dimension(0) {}
		explicit Model( uint32_t dim ) : dimension(dim) {}
	};

	/*!
	 * \brief Base for model with parameters.
	 */
	template<typename Model,typename Parameter = VectorXd>
	struct Family {
		typedef Model Model;
		typedef Parameter Parameter;
		uint32_t dimension,          //!< Dimension of the model's state space
		         parameterDimension; //!< Dimension of the model's parameter space

		Family() : dimension(0), parameterDimension(0) {}
		explicit Family( uint32_t dim, uint32_t pdim ) : dimension(dim), parameterDimension(pdim) {}
		/*Model operator()( const Parameter& parameter ) const {
			return Model(parameter);
		}*/
	};
	
	/*!
	 * \brief A model without parameters whose state space is embedded into R^n.
	 */
	template<typename Model,typename Embedding>
	struct ModelEmbedded {
		typedef typename Model::State State;
		Model model;         //!< The underlying model
		Embedding embedding; //!< An embedding into R^n

		ModelEmbedded() = default;
		ModelEmbedded( const Model& m, const Embedding& e ) : model(m), embedding(e) {}

		VectorXd operator()( const State& state ) const {
			return embedding.Derivative(state) * model(state);
		}

		MatrixXd Partials( const State& state ) const {

			MatrixXd embeddingPartials2, result;
			result.setZero(embedding.oDim,embedding.oDim);

			VectorXd vector = model(state);

			for(uint32_t i=0;i<embedding.eDim;++i) {
				embeddingPartials2 = embedding.Derivative2(state,i);
				result.row(i) += embeddingPartials2 * vector;
			}

			result += embedding.Derivative(state) * model.Partials(state);

			return result;
		}

	};

	/*!
	 * \brief A model without parameters whose state space is embedded into R^n,
	 * and whose derivates are to be evaluated component-wise.
	 */
	template<typename Model,typename Embedding>
	struct ModelEmbeddedCW {
		Model model;         //!< The underlying model
		Embedding embedding; //!< An embedding into R^n

		ModelEmbeddedCW() = default;
		ModelEmbeddedCW( const Model& m, const Embedding& e ) : model(m), embedding(e) {}

		double operator()( const VectorXd& state, uint32_t i ) const {
			double result = 0.0;

			for(uint32_t j=0;j<embedding.oDim;++j) {
				result += embedding.Derivative(state,i,j) * model(state,j);
			}
			return result;
		}

		double Partials( const VectorXd& state, uint32_t i, uint32_t j ) const {
			double result = 0.0;

			for(uint32_t k=0;k<embedding.oDim;++k) {
				result += embedding.Derivative2(state,i,k,j) * model(state,k);
				result += embedding.Derivative(state,i,k) * model.Partials(state,k,j);
			}
			return result;
		}
	
	};

	/*!
	 * \brief A model with parameters whose state space is embedded into R^n.
	 */
	template<typename Family,typename Embedding>
	struct FamilyEmbedded {
		typedef typename Family::Model Model;
		typedef typename Family::Parameter Parameter;
		Family family;       //!< The underlying family
		Embedding embedding; //!< An embedding into R^n

		FamilyEmbedded() = default;
		FamilyEmbedded( const Family& f, const Embedding& e ) : family(f), embedding(e) {}

		ModelEmbedded<Model,Embedding> operator()( const Parameter& parameter ) const {
			return ModelEmbedded<Model,Embedding>( family(parameter), embedding );
		}

	};

	/*!
	 * \brief A model with parameters whose state space is embedded into R^n.
	 */
	template<typename Family,typename Embedding>
	struct FamilyEmbeddedCW {
		typedef typename Family::Model Model;
		typedef typename Family::Parameter Parameter;
		Family family;       //!< The underlying family
		Embedding embedding; //!< An embedding into R^n

		FamilyEmbeddedCW() = default;
		FamilyEmbeddedCW( const Family& f, const Embedding& e ) : family(f), embedding(e) {}

		ModelEmbeddedCW<Model,Embedding> operator()( const Parameter& parameter ) const {
			return ModelEmbeddedCW<Model,Embedding>( family(parameter), embedding );
		}

	};

}

#endif
