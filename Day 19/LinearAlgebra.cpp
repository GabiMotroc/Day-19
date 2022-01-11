#include "LinearAlgebra.hpp"

#include <cassert>
#include <ostream>
#include <unordered_set>
#include <vector>

constexpr double epsilon = 100000 * DBL_EPSILON;

template<>
struct std::hash<vector3>
{
	std::size_t operator()(const vector3& k) const noexcept
	{
		using std::size_t;
		using std::hash;

		// Compute individual hash values for first,
		// second and third and combine them using XOR
		// and bit shifting:

		return ((hash<double>()(k.x)
			^ (hash<double>()(k.y) << 1)) >> 1)
			^ (hash<double>()(k.z) << 1);
	}
};

bool areSame(const double a, const double b)
{
	return fabs(a - b) < epsilon;
}


double getDistance(const vector3& a, const vector3& b)
{
	auto result = pow(b.x - a.x, 2) + pow(b.y - a.y, 2) + pow(b.z - a.z, 2);
	result = sqrt(result);
	return result;
}

std::vector<std::vector<double>> transformPointTo2DArray(const vector3& point)
{
	std::vector<std::vector<double>> result =
	{
		{point.x},
		{point.y},
		{point.z}
	};
	return result;
}

inline std::vector<std::vector<double>> generateMatrix(int row, int col) {
	std::vector<std::vector<double>> ret;
	ret.reserve(row);
	for (int i = 0; i < row; i++)
		ret.emplace_back(col);
	return ret;
}

std::vector<std::vector<double>> multiplyMatrices(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b)
{
	const int n = a.size();     // a rows
	const int m = a[0].size();  // a cols
	const int p = b[0].size();  // b cols

	std::vector <std::vector<double>> c(n, std::vector<double>(p, 0));
	for (auto j = 0; j < p; ++j)
	{
		for (auto k = 0; k < m; ++k)
		{
			for (auto i = 0; i < n; ++i)
			{
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return c;
}

vector3 multiplyMatrices(const std::vector<std::vector<double>>& a, const vector3& point)
{
	const auto b = transformPointTo2DArray(point);

	const auto result = multiplyMatrices(a, b);

	const auto result2 = vector3{ result[0][0], result[1][0], result[2][0] };
	return result2;
}

std::vector<std::vector<double>> multiplyWithScalar(const std::vector<std::vector<double>>& vector, double aux)
{
	auto result = vector;

	for (auto& line : result)
	{
		for (auto& value : line)
		{
			value *= aux;
		}
	}

	return result;
}

vector3 findCenterOfGravity(const std::unordered_set<vector3>& vectors)
{
	vector3 result;

	for (const auto& point : vectors)
	{
		result.x += point.x;
		result.y += point.y;
		result.z += point.z;
	}

	result.x = result.x / static_cast<double>(vectors.size());
	result.y = result.y / static_cast<double>(vectors.size());
	result.z = result.z / static_cast<double>(vectors.size());

	return result;
}

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& a) {
	const int rows = a.size();

	if (rows == 0) return { {} };

	const int cols = a[0].size();

	std::vector<std::vector<double>> r(cols, std::vector<double>(rows));

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			r[j][i] = a[i][j];
		}
	}

	return r;
}

std::vector<std::vector<double>> getRotationMatrix2(const vector3& first, const vector3& second)
{
	auto result = generateMatrix(3, 3);

	for (int i = 0; i < 3; ++i)
	{
		result[0][i] = first.x;
		result[1][i] = first.y;
		result[2][i] = first.z;
	}

	for (int i = 0; i < 3; ++i)
	{
		if (areSame(fabs(result[i][0]), fabs(second.x)))
			result[i][0] /= second.x;
		else
			result[i][0] = 0;
	}

	for (int i = 0; i < 3; ++i)
	{
		if (areSame(fabs(result[i][1]), fabs(second.y)))
			result[i][1] /= second.y;
		else
			result[i][1] = 0;
	}

	for (int i = 0; i < 3; ++i)
	{
		if (areSame(fabs(result[i][2]), fabs(second.z)))
			result[i][2] /= second.z;
		else
			result[i][2] = 0;
	}

	result = transpose(result);

	if (const auto rotated = multiplyMatrices(result, first); !(second == rotated))
	{
		throw 100;
	}

	return result;
}

std::vector<std::vector<double>> getRotationMatrix3(const std::unordered_set<vector3>& a, const std::unordered_set<vector3>& destination)
{
	std::unordered_set<vector3> firstVectors;
	std::unordered_set<vector3> secondVectors;

	const auto centroidA = findCenterOfGravity(a);
	const auto centroidB = findCenterOfGravity(destination);

	for (const auto& point : a)
	{
		firstVectors.insert(vector3{ point.x - centroidA.x, point.y - centroidA.y, point.z - centroidA.z });
	}

	for (const auto& point : destination)
	{
		secondVectors.insert(vector3{ point.x - centroidB.x, point.y - centroidB.y, point.z - centroidB.z });
	}

	std::vector<std::vector<double>> result;

	for (const auto& firstVector : firstVectors)
	{
		for (const auto& secondVector : secondVectors)
		{
			if (areSame(fabs(firstVector.x), fabs(secondVector.x))
				and areSame(fabs(firstVector.y), fabs(secondVector.y))
				and areSame(fabs(firstVector.z), fabs(secondVector.z)))
			{
				return getRotationMatrix2(firstVector, secondVector);
			}
			if (areSame(fabs(firstVector.x), fabs(secondVector.y))
				and areSame(fabs(firstVector.y), fabs(secondVector.z))
				and areSame(fabs(firstVector.z), fabs(secondVector.x)))
			{
				return getRotationMatrix2(firstVector, secondVector);
			}
			if (areSame(fabs(firstVector.x), fabs(secondVector.z))
				and areSame(fabs(firstVector.y), fabs(secondVector.x))
				and areSame(fabs(firstVector.z), fabs(secondVector.y)))
			{
				return getRotationMatrix2(firstVector, secondVector);
			}
			if (areSame(fabs(firstVector.x), fabs(secondVector.x))
				and areSame(fabs(firstVector.y), fabs(secondVector.z))
				and areSame(fabs(firstVector.z), fabs(secondVector.y)))
			{
				return getRotationMatrix2(firstVector, secondVector);
			}
			if (areSame(fabs(firstVector.x), fabs(secondVector.z))
				and areSame(fabs(firstVector.y), fabs(secondVector.y))
				and areSame(fabs(firstVector.z), fabs(secondVector.x)))
			{
				return getRotationMatrix2(firstVector, secondVector);
			}
			if (areSame(fabs(firstVector.x), fabs(secondVector.y))
				and areSame(fabs(firstVector.y), fabs(secondVector.x))
				and areSame(fabs(firstVector.z), fabs(secondVector.z)))
			{
				return getRotationMatrix2(firstVector, secondVector);
			}
		}
	}

	return result;
}

vector3 computeTranslation(const std::unordered_set<vector3>& a, const std::unordered_set<vector3>& b, const std::vector<std::vector<double>>& rotationMatrix)
{
	const auto firstPoint = findCenterOfGravity(a);

	auto secondPoint = findCenterOfGravity(b);
	secondPoint = (multiplyMatrices(rotationMatrix, secondPoint));

	const auto shiftMatrix = vector3
	{
		(firstPoint.x - secondPoint.x),
		(firstPoint.y - secondPoint.y),
		(firstPoint.z - secondPoint.z)
	};
	return shiftMatrix;
}


double manhattanDistance(const vector3& a, const vector3& b)
{
	double result = round(fabs(a.x - b.x));
	result += round(fabs(a.y - b.y));
	result += round(fabs(a.z - b.z));
	return result;
}