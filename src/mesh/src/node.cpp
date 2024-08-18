#include "mesh/include/node.h"
#include "iniCondition/include/GlobalData.h"
#include "utility/include/log.h"


Node::Node(const Coordinate& coor, const size_t index, const std::shared_ptr<Mesh>& ptr, Nodetype type)
	: nodecoor_(coor), nodeindex_(index), nodetype_(type), meshptr_(ptr),primitiveValue_(4, 0.0), conservedValue_(4, 0.0),
	 coordTrans_ { /* ksi_x_= */ 0.0, /* ksi_y_= */ 0.0, /* eta_x_= */ 0.0, /* eta_y_= */ 0.0, /* jacobian_= */ 0.0 }
{
}

Node& Node::operator=(const Node& other) {
	if (this != &other) {
		nodecoor_ = other.nodecoor_;
		nodeindex_ = other.nodeindex_;
		primitiveValue_ = other.primitiveValue_;
		conservedValue_ = other.conservedValue_;
		coordTrans_ = other.coordTrans_;
	}
	return *this;
}

bool Node::operator==(const Node& other) const {
	return (nodecoor_ == other.nodecoor_) && (nodeindex_ == other.nodeindex_);
}

void Node::toConservedform()
{
	const double gamma = GlobalData::GetDouble("refGamma");
	const auto& density = r();
	const auto& xVel = u();
	const auto& yVel = v();
	//const auto& zVel = w();
	const auto& pressure = p();

	conservedValue_[0] = density;
	conservedValue_[1] = density * xVel;
	conservedValue_[2] = density * yVel;
	//conservedValue_[3] = density * zVel;
	conservedValue_[3] = pressure / (gamma - 1) + 0.5 * density * (xVel * xVel + yVel * yVel /*+ zVel * zVel*/);

}

void Node::toPrimitiveform()
{
	const double gamma = GlobalData::GetDouble("refGamma");
	const auto& cons0 = rho();
	const auto& cons1 = rhou();
	const auto& cons2 = rhov();
	//const auto& cons3 = rhow();
	const auto& cons4 = E();

	double density = cons0;
	double xVel = cons1 / cons0;
	double yVel = cons2 / cons0;
	//double zVel = cons3 / cons0;
	double pressure =  (gamma - 1) * (cons4 - 0.5 * density * (xVel * xVel + yVel * yVel /*+ zVel * zVel*/));
	primitiveValue_.at(0) = density;
	primitiveValue_.at(1) = xVel;
	primitiveValue_.at(2) = yVel;
	//primitiveValue_.at(3) = zVel;
	primitiveValue_.at(3) = pressure;
	if (pressure < 0.0)
	{
		spdlog::error("the pressure is negative !!!");
		throw std::runtime_error("Error the pressure is negative.");
	}
	//checkNaN(pressure);
}

void Node::setCoordTrans(double inksi_x, double inksi_y, double ineta_x, double ineat_y, double inJacobian)
{
	coordTrans_ = { inksi_x, inksi_y, 0.0, ineta_x, ineat_y,0.0, inJacobian };
}


void Node::checkNaN(const double value)
{
	if (std::isnan(value))
	{
		spdlog::error("the result is NaN !!!");
		throw std::runtime_error("Errorthe result is NaN.");
	}
}