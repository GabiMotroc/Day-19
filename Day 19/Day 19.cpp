// Day 19.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "LinearAlgebra.hpp"

#pragma region Structs

struct Scanner
{
	int id{};
	bool used = false;
	std::vector<vector3> points;
};

template<>
struct std::hash<vector3>
{
	std::size_t operator()(const vector3& k) const noexcept
	{
		using std::size_t;
		using std::hash;
		using std::string;

		// Compute individual hash values for first,
		// second and third and combine them using XOR
		// and bit shifting:

		return ((hash<double>()(k.x)
				^ (hash<double>()(k.y) << 1)) >> 1)
			^ (hash<double>()(k.z) << 1);
	}
};

struct PairHash {
	inline std::size_t operator()(const std::pair<vector3, vector3>& v) const {
		const auto h1 = std::hash<vector3>()(v.first);
		const auto h2 = std::hash<vector3>()(v.second);
		return h1 ^ h2;
	}
};

#pragma endregion

std::vector<Scanner> scanners;
std::unordered_set<vector3> solution;

std::ifstream fin("19.txt");
std::ifstream fin2("solution.txt");

vector3 splitStringInPoint(std::string s)
{
	std::vector<double> aux;
	const std::string delimiter = ",";
	
	size_t pos;

	while ((pos = s.find(delimiter)) != std::string::npos) 
	{
		std::string token = s.substr(0, pos);
		aux.push_back(stoi(token));
		s.erase(0, pos + delimiter.length());
	}
	aux.push_back(stoi(s));

	return vector3{ aux[0], aux[1], aux[2] };
}

void readInput()
{
	std::string s;

	Scanner aux;
	int counter = 0;
	while(getline(fin, s))
	{
		if (s.empty())
		{
			aux.id = counter++;
			scanners.push_back(aux);
			aux.points.clear();
			continue;
		}

		if (s[1] == '-')
			continue;

		aux.points.push_back(splitStringInPoint(s));
	}
	aux.id++;
	scanners.push_back(aux);

}

void readSolution()
{
	std::string s;
	while (getline(fin2, s))
	{
		const auto aux = splitStringInPoint(s);
		solution.insert(aux);
	}
}

std::unordered_map<double, std::unordered_set<std::pair<vector3, vector3>, PairHash>> getDistancesOfScanner(std::vector<Scanner>::const_reference scanner1)
{
	std::unordered_map<double, std::unordered_set<std::pair<vector3, vector3>, PairHash>> result;
	for (std::size_t i = 0; i < scanner1.points.size(); ++i)
	{
		for (std::size_t j = i + 1; j < scanner1.points.size(); ++j)
		{
			if(i != j)
			{
				const auto aux = getDistance(scanner1.points[i], scanner1.points[j]);
				result[aux].insert(std::make_pair(scanner1.points[i], scanner1.points[j]));
			}
		}
	}
	return result;
}

Scanner connect2Scanners(std::vector<Scanner>::const_reference scanners1, std::vector<Scanner>::const_reference scanners2, const std::vector<std::vector<double>>& rotationMatrix, const vector3& shift)
{
	std::unordered_set<vector3> result;

	for (const auto& point : scanners1.points)
	{
		result.insert(point);
	}

	for (const auto& point : scanners2.points)
	{
		auto aux = multiplyMatrices(rotationMatrix, point);
		aux.x += shift.x;
		aux.y += shift.y;
		aux.z += shift.z;

		aux.x = round(aux.x);
		aux.y = round(aux.y);
		aux.z = round(aux.z);

		if (solution.find(aux) == solution.end())
		{
			//throw std::exception("done an oopsie");
		}

		result.insert(aux);
	}

	return Scanner{scanners1.id, true, std::vector<vector3>(result.begin(), result.end())};
}

std::pair<std::vector<std::pair<vector3, vector3>>, std::vector<std::pair<vector3, vector3>>>
extractCommonPoints(
	const std::unordered_map<double, std::unordered_set<std::pair<vector3, vector3>, PairHash>>& distance1,
	std::unordered_map<double, std::unordered_set<std::pair<vector3, vector3>, PairHash>>& distance2,
	int counter)
{

	std::unordered_map<vector3, int> pointFrequency1;
	std::unordered_map<vector3, int> pointFrequency2;

	std::vector<std::pair<vector3, vector3>> vectors1;
	std::vector<std::pair<vector3, vector3>> vectors2;

	for (auto& [distance, points] : distance1)
	{
		if (distance2.find(distance) != distance2.end())
		{
			if (auto aux = distance2[distance]; true/* points.size() == 1 and aux.size() == 1*/)
			{
				const auto firstVector = *(points.begin());
				const auto secondVector = *(aux.begin());
				vectors1.push_back(firstVector);
				vectors2.push_back(secondVector);
			}
		}
	}

	if (vectors1.size() != vectors2.size())
	{

		if (vectors1.size() < vectors2.size())
		{
			vectors2.erase(
				std::remove_if(
					vectors2.begin(),
					vectors2.end(),
					[&vectors1](auto a)
					{

					}
				),
				vectors2.end()
						);
		}
	}

	if (counter == 66)
		return std::make_pair(vectors1, vectors2);

	for (const auto & [first, second] : vectors1)
	{
		pointFrequency1[first]++;
		pointFrequency1[second]++;
	}

	for (const auto & [first, second] : vectors2)
	{
		pointFrequency2[first]++;
		pointFrequency2[second]++;
	}

	vectors1.erase(
		std::remove_if(
			vectors1.begin(),
			vectors1.end(),
			[&pointFrequency1](const std::pair<vector3, vector3>& pair)
			{
				return (pointFrequency1[pair.first] < 11 or pointFrequency1[pair.second] < 11);
			}
		),
		vectors1.end()
				);

	vectors2.erase(
		std::remove_if(
			vectors2.begin(),
			vectors2.end(),
			[&pointFrequency2](const std::pair<vector3, vector3>& pair)
			{
				return (pointFrequency2[pair.first] < 11 or pointFrequency2[pair.second] < 11);
			}
		),
		vectors2.end()
				);


	return std::make_pair(vectors1, vectors2);
}

Scanner doWork(std::vector<Scanner>::const_reference scanner1, std::vector<Scanner>::const_reference scanner2)
{
	auto distance1 = getDistancesOfScanner(scanner1);
	auto distance2 = getDistancesOfScanner(scanner2);

	int counter = 0;

	for (const auto& value : distance1)
	{
		if(distance2.find(value.first) != distance2.end())
		{
			counter++;
		}
	}

	// HACK just an assumption
	if (counter >= 66)
	{
		// TODO implement the rotations and reverence point calculations
		auto [vectors1, vectors2] =
			extractCommonPoints(distance1, distance2, counter);

		std::sort(vectors1.begin(), vectors1.end(),
			[](const std::pair<vector3, vector3>& a, const std::pair<vector3, vector3>& b) -> bool
			{
				return getDistance(a.first, a.second) < getDistance(b.first, b.second);
			});

		std::sort(vectors2.begin(), vectors2.end(),
			[](const std::pair<vector3, vector3>& a, const std::pair<vector3, vector3>& b) -> bool
			{
				return getDistance(a.first, a.second) < getDistance(b.first, b.second);
			});

		const auto rotationMatrix = getRotationMatrix3(vectors2, vectors1);

		for (const auto& matrix : rotationMatrix)
		{
			for (const double& value : matrix)
			{
				std::cout << value << ' ';
			}
			std::cout << '\n';
		}
		/*
		//auto i = 0;
		//auto vectorA = getVector(vectors1[i]);
		//auto vectorB = getVector(vectors2[i]);
		//i++;

		//while ((areSame(fabs(vectorA.x), fabs(vectorA.y))
		//	or areSame(fabs(vectorA.y), fabs(vectorA.z))
		//	or areSame(fabs(vectorA.x), fabs(vectorA.z))
		//	or areSame(fabs(vectorA.x), 0)
		//	or areSame(fabs(vectorA.y), 0)
		//	or areSame(fabs(vectorA.z), 0)
		//	/*or isDeterminantEqualTo1(rotationMatrix) == false)
		//	and i < 66)
		//{
		//	vectorA = getVector(vectors1[i]);
		//	vectorB = getVector(vectors2[i]);
		//	i++;
		//}

		//auto rotationMatrix = std::vector<std::vector<double>>();

		//try
		//{
		//	rotationMatrix = getRotationMatrix2(vectorB, vectorA);
		//}
		//catch (int e)
		//{
		//	std::cout << e;
		//}
		*/
		const auto shiftMatrix = computeTranslation(vectors1, vectors2, rotationMatrix);

		auto result = connect2Scanners(scanner1, scanner2, rotationMatrix, shiftMatrix);

		return result;
	}

	return Scanner{ -1, false, std::vector<vector3>() };
}

Scanner solve1()
{
	Scanner scanner;
	bool firstPass = true;
	for (auto i = 0; i < scanners.size() - 1; ++i)
	{
		for (int j = i + 1; j < scanners.size(); ++j)
		{
			if ( firstPass == true and i != j and scanners[j].used == false)
			{
				const Scanner aux = doWork(scanners[i], scanners[j]);
				if (aux.id != -1)
				{
					scanner = aux;
					firstPass = false;
					scanners[i].used = true;
					scanners[j].used = true;
				}
			}
			else if(scanners[j].used == false and firstPass == false and i != j)
			{
				const Scanner aux = doWork(scanner, scanners[j]);
				if (aux.id != -1)
				{
					scanner = aux;
					scanners[j].used = true;
				}
			}
		}
	}

	// TODO implement the math for finding the number of beacons
	return scanner;
}

int main()
{
	readInput();
	readSolution();


	const auto scanner = solve1();
	for (const auto& _scanner : scanners)
	{
		std::cout << _scanner.id << ": " << _scanner.used << '\n';
	}

	std::cout << scanner.points.size() << '\n';

	auto debug = scanner.points;
	std::sort(debug.begin(), debug.end(),
		[](const vector3& a, const vector3& b) -> bool
		{
			return a.x * a.x + a.y * a.y + a.z * a.z < b.x* b.x + b.y * b.y + b.z * b.z;
		});

	for (const auto& point :debug)
	{
		std::cout << point << '\n';
	}
}
