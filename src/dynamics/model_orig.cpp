#include <DRDSP/dynamics/model_orig.h>

using namespace DRDSP;

ModelEmbedded::ModelEmbedded( Model& m, Embedding& e ) : model(m), embedding(e) {}

VectorXd ModelEmbedded::VectorField( const VectorXd& state ) const {
	return embedding.Derivative(state) * model.VectorField(state);
}

MatrixXd ModelEmbedded::Partials( const VectorXd& state ) const {

	MatrixXd embeddingPartials2, result;
	result.setZero(embedding.oDim,embedding.oDim);

	VectorXd vector = model.VectorField(state);

	for(uint32_t i=0;i<embedding.eDim;i++) {
		embeddingPartials2 = embedding.Derivative2(state,i);
		result.row(i) += embeddingPartials2 * vector;
	}

	result += embedding.Derivative(state) * model.Partials(state);

	return std::move(result);
}

ModelEmbeddedCW::ModelEmbeddedCW( ModelCW& m, EmbeddingCW& e ) : model(m), embedding(e) {}

double ModelEmbeddedCW::VectorField( const VectorXd& state, uint32_t i ) const {
	double result = 0.0;

	for(uint32_t j=0;j<embedding.oDim;j++) {
		result += embedding.Derivative(state,i,j) * model.VectorField(state,j);
	}
	return result;
}

double ModelEmbeddedCW::Partials( const VectorXd& state, uint32_t i, uint32_t j ) const {
	double result = 0.0;

	for(uint32_t k=0;k<embedding.oDim;k++) {
		result += embedding.Derivative2(state,i,k,j) * model.VectorField(state,k);
		result += embedding.Derivative(state,i,k) * model.Partials(state,k,j);
	}
	return result;
}

ModelParameterizedEmbedded::ModelParameterizedEmbedded( ModelParameterized& m, Embedding& e ) : model(m), embedding(e) {}

VectorXd ModelParameterizedEmbedded::VectorField( const VectorXd& state, const VectorXd& parameter ) const {
	return embedding.Derivative(state) * model.VectorField(state,parameter);
}

MatrixXd ModelParameterizedEmbedded::Partials( const VectorXd& state, const VectorXd& parameter ) const {

	MatrixXd embeddingPartials2, result;
	result.setZero(embedding.eDim,embedding.oDim);

	VectorXd vector = model.VectorField(state,parameter);

	for(uint32_t i=0;i<embedding.eDim;i++) {
		embeddingPartials2 = embedding.Derivative2(state,i);
		result.row(i) += embeddingPartials2 * vector;
	}

	result += embedding.Derivative(state) * model.Partials(state,parameter);

	return std::move(result);
}

ModelParameterizedEmbeddedCW::ModelParameterizedEmbeddedCW( ModelParameterizedCW& m, EmbeddingCW& e ) : model(m), embedding(e) {}

double ModelParameterizedEmbeddedCW::VectorField( const VectorXd& state, const VectorXd& parameter, uint32_t i ) const {
	double result = 0.0;

	for(uint32_t j=0;j<embedding.oDim;j++) {
		result += embedding.Derivative(state,i,j) * model.VectorField(state,parameter,j);
	}
	return result;
}

double ModelParameterizedEmbeddedCW::Partials( const VectorXd& state, const VectorXd& parameter, uint32_t i, uint32_t j ) const {
	double result = 0.0;

	for(uint32_t k=0;k<embedding.oDim;k++) {
		result += embedding.Derivative2(state,i,k,j) * model.VectorField(state,parameter,k);
		result += embedding.Derivative(state,i,k) * model.Partials(state,parameter,k,j);
	}
	return result;
}
