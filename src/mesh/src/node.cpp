#include "mesh/include/node.h"
#include "iniCondition/include/GlobalData.h"

Node::Node(const Coordinate& coor, const size_t index) : nodecoor_(coor), nodeindex_(index)
{
	initValues();
}

Node& Node::operator=(const Node& other) {
	if (this != &other) {
		nodecoor_ = other.nodecoor_;
		nodeindex_ = other.nodeindex_;
	}
	return *this;
}

bool Node::operator==(const Node& other) const {
	return (nodecoor_ == other.nodecoor_) && (nodeindex_ == other.nodeindex_);
}

void Node::toConservedform()
{
	const double gamma = GlobalData::GetDouble("refGamma");
	auto& density = primitiveValue_.at("rho");
	auto& xVel = primitiveValue_.at("u");
	auto& yVel = primitiveValue_.at("v");
	auto& zVel = primitiveValue_.at("w");
	auto& pressure = primitiveValue_.at("p");

	conservedValue_["rho"] = density;
	conservedValue_["rho_u"] = density * xVel;
	conservedValue_["rho_v"] = density * yVel;
	conservedValue_["rho_w"] = density * zVel;
	conservedValue_["E"] = pressure / (gamma - 1) + 0.5 * density * (xVel * xVel + yVel * yVel + zVel * zVel);

}

void Node::toPrimitiveform()
{
	const double gamma = GlobalData::GetDouble("refGamma");
	auto& cons0 = conservedValue_.at("rho");
	auto& cons1 = conservedValue_.at("rho_u");
	auto& cons2 = conservedValue_.at("rho_v");
	auto& cons3 = conservedValue_.at("rho_w");
	auto& cons4 = conservedValue_.at("E");

	double density = cons0;
	double xVel = cons1 / cons0;
	double yVel = cons2 / cons0;
	double zVel = cons3 / cons0;
	double pressure = pressure = (gamma - 1) * (cons4 - 0.5 * density * (xVel * xVel + yVel * yVel + zVel * zVel));
	primitiveValue_.at("rho") = density;
	primitiveValue_.at("u") = xVel;
	primitiveValue_.at("v") = yVel;
	primitiveValue_.at("w") = zVel;
	primitiveValue_.at("p") = pressure;

}

void Node::initValues()
{
	primitiveValue_ = {
	{"rho", 0.0},
	{"u", 0.0},
	{"v", 0.0},
	{"w", 0.0},
	{"p", 0.0}
	};
	conservedValue_ = {
		{"rho", 0.0},
		{"rho_u", 0.0},
		{"rho_v", 0.0},
		{"rho_w", 0.0},
		{"E", 0.0}
	};
}