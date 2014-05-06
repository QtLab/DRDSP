#ifndef INCLUDED_DYNAMICS_MODEL
#define INCLUDED_DYNAMICS_MODEL
#include "embedding.h"
#include "dynamicalSystem.h"

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
		State operator()( const State& x, Time ) const {
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

			MatrixXd result;
			result.setZero(embedding.eDim,embedding.oDim);

			VectorXd vector = model(state);

			for(uint32_t i=0;i<embedding.eDim;++i) {
				result.row(i) += embedding.Derivative2(state,i) * vector;
			}

			result += embedding.Derivative(state) * model.Partials(state);

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
		explicit FamilyEmbedded( const Family& f ) : family(f) {}
		explicit FamilyEmbedded( const Embedding& e ) : embedding(e) {}
		FamilyEmbedded( const Family& f, const Embedding& e ) : family(f), embedding(e) {}

		ModelEmbedded<Model,Embedding> operator()( const Parameter& parameter ) const {
			return ModelEmbedded<Model,Embedding>( family(parameter), embedding );
		}

	};

}

#endif
