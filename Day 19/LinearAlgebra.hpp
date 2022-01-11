#pragma once
#include <unordered_set>
#include <vector>
#include "vector3.hpp"

bool areSame(const double a, const double b);

vector3 transform2DArrayIntoPoint(const std::vector<std::vector<double>>& a);

double getDistance(const vector3& a, const vector3& b);

vector3 normalizeVector(const vector3& a);

vector3 getVector(std::vector<std::pair<vector3, vector3>>::const_reference pair);

vector3 getCrossProduct(const vector3& a, const vector3& b);

double getDotProduct(const vector3& a, const vector3& b);

std::vector<std::vector<double>> transformPointTo2DArray(const vector3& point);

inline std::vector<std::vector<double>> generateMatrix(int row, int col);

std::vector<std::vector<double>> multiplyMatrices(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b);

vector3 multiplyMatrices(const std::vector<std::vector<double>>& a, const vector3& point);

std::vector<std::vector<double>> addMatrices(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b);

std::vector<std::vector<double>> multiplyWithScalar(const std::vector<std::vector<double>>& vector, double aux);

/**
 * \brief Mda
 * \param first The first vector
 * \param second The second vector
 * \return Returns the rotation matrix necessary to align a with b
 */
std::vector<std::vector<double>> getRotationMatrix(const vector3& first, const vector3& second);

vector3 findCenterOfGravity(const std::vector<std::pair<vector3, vector3>>& vectors);

std::vector<std::vector<double>> getRotationMatrix2(const vector3& first, const vector3& second);

std::vector<std::vector<double>> getRotationMatrix3(const std::vector<std::pair<vector3, vector3>>& a, const std::vector<std::pair<vector3, vector3>>& b);

std::vector<std::vector<double>> getRotationMatrix3(const std::unordered_set<vector3>& a, const std::unordered_set<vector3>& destination);

vector3 computeTranslation(const std::vector<std::pair<vector3, vector3>>& a, const std::vector<std::pair<vector3, vector3>>& b, const std::vector<std::vector<double>>& rotationMatrix);

vector3 computeTranslation(const std::unordered_set<vector3>& a, const std::unordered_set<vector3>& b, const std::vector<std::vector<double>>& rotationMatrix);

void getCofactor(std::vector<std::vector<double>> mat, std::vector<std::vector<double>> temp, int p, int q, int n);

bool isDeterminantEqualTo1(const std::vector<std::vector<double>>& matrix);

double manhattanDistance(const vector3& a, const vector3& b);