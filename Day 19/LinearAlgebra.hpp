#pragma once
#include <unordered_set>
#include <vector>
#include "vector3.hpp"

bool areSame(const double a, const double b);

double getDistance(const vector3& a, const vector3& b);

std::vector<std::vector<double>> transformPointTo2DArray(const vector3& point);

inline std::vector<std::vector<double>> generateMatrix(int row, int col);

std::vector<std::vector<double>> multiplyMatrices(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b);

vector3 multiplyMatrices(const std::vector<std::vector<double>>& a, const vector3& point);

std::vector<std::vector<double>> multiplyWithScalar(const std::vector<std::vector<double>>& vector, double aux);

vector3 findCenterOfGravity(const std::unordered_set<vector3>& vectors);

std::vector<std::vector<double>> getRotationMatrix2(const vector3& first, const vector3& second);

std::vector<std::vector<double>> getRotationMatrix3(const std::unordered_set<vector3>& a, const std::unordered_set<vector3>& destination);

vector3 computeTranslation(const std::unordered_set<vector3>& a, const std::unordered_set<vector3>& b, const std::vector<std::vector<double>>& rotationMatrix);

double manhattanDistance(const vector3& a, const vector3& b);