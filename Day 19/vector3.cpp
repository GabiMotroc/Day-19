#include "vector3.hpp"

#include "LinearAlgebra.hpp"


vector3::vector3(const double x, const double y, const double z) : x(x), y(y), z(z) {};

vector3::vector3() : x(0), y(0), z(0)
{
}

bool vector3::operator==(const vector3 & otherPoint) const
{
	return areSame(this->x, otherPoint.x) and areSame(this->y, otherPoint.y) and areSame(this->z, otherPoint.z);
}