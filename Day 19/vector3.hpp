#pragma once
#include <ostream>
#include <vector>

class vector3
{
public:
	double x, y, z;

	vector3(const double x, const double y, const double z);
	vector3(const vector3& point) = default;
	vector3();

	bool operator==(const vector3& otherPoint) const;


	friend std::ostream& operator<<(std::ostream& os, const vector3& point)
	{
		os << point.x << "," << point.y << "," << point.z << ' ';
		return os;
	}
};