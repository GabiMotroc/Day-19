// Day 19.cpp : This file contains the 'main' function. Program execution begins and ends there.
// https://adventofcode.com/2021/day/19

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <optional>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "LinearAlgebra.hpp"

struct PairHash;
using std::vector;
using std::unordered_set;
using std::pair;
using std::unordered_map;
using std::string;

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

struct PairHash {
	inline std::size_t operator()(const std::pair<vector3, vector3>& v) const {
		const auto h1 = std::hash<vector3>()(v.first);
		const auto h2 = std::hash<vector3>()(v.second);
		return h1 ^ h2;
	}
};


#pragma region Structs

struct Distance
{
	double distance;
	vector3 a, b;
};

struct Scanner
{
	int id{};
	bool used = false;
	vector<vector3> points;
	std::optional<int> reference;

	unordered_map<double, pair<vector3, vector3>> distances;
	vector<vector<double>> rotationMatrix;
	vector3 shift;

	pair<unordered_set<vector3>, unordered_set<vector3>> findCommonBeacons(Scanner const& s);

	void computeTransformation(Scanner const& s, const pair<unordered_set<vector3>, unordered_set<vector3>>&);
	void computeDistances();
};

pair<unordered_set<vector3>, unordered_set<vector3>> Scanner::findCommonBeacons(Scanner const& destination)
{
	unordered_set<vector3> pointsIn;
	unordered_set<vector3> pointsInDestination;


	for (auto& [key, value] : distances)
	{
		if(destination.distances.find(key) != destination.distances.end())
		{
			const auto& [fst, snd] = destination.distances.at(key);
			pointsInDestination.insert(fst);
			pointsInDestination.insert(snd);

			pointsIn.insert(value.first);
			pointsIn.insert(value.second);
		}
	}

	return std::make_pair(pointsInDestination, pointsIn);
}

void Scanner::computeTransformation(Scanner const& s, const pair<unordered_set<vector3>, unordered_set<vector3>>& mapping)
{
	rotationMatrix = getRotationMatrix3(mapping.second, mapping.first);
	shift = computeTranslation(mapping.first, mapping.second, rotationMatrix);

	reference = s.id;
}

void Scanner::computeDistances()
{
	for (auto i = 0ULL; i < points.size(); ++i)
	{
		for (auto j = i + 1; j < points.size(); ++j)
		{
			if (i != j)
			{
				const auto aux = getDistance(points[i], points[j]);
				distances[aux] = std::make_pair(points[i], points[j]);
			}
		}
	}
}


#pragma endregion

unordered_set<vector3> solution;

std::ifstream fin("19.txt");
std::ifstream fin2("solution.txt");

vector3 splitStringInPoint(string s)
{
	vector<double> aux;
	const string delimiter = ",";
	
	size_t pos;

	while ((pos = s.find(delimiter)) != string::npos) 
	{
		string token = s.substr(0, pos);
		aux.push_back(stoi(token));
		s.erase(0, pos + delimiter.length());
	}
	aux.push_back(stoi(s));

	return vector3{ aux[0], aux[1], aux[2] };
}

vector<Scanner> readInput()
{
	string s;

	vector<Scanner> scanners;

	Scanner aux;
	int counter = 0;
	while(getline(fin, s))
	{
		if (s.empty())
		{
			aux.id = counter++;
			aux.computeDistances();
			scanners.push_back(aux);
			aux.points.clear();
			aux.distances.clear();
			continue;
		}

		if (s[1] == '-')
			continue;

		aux.points.push_back(splitStringInPoint(s));
	}
	aux.id++;
	aux.computeDistances();
	scanners.push_back(aux);

	scanners[0].reference = 0;

	scanners[0].rotationMatrix =
	{
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	};

	scanners[0].shift = vector3{ 0,0,0 };

	return scanners;

}

void readSolution()
{
	string s;
	while (getline(fin2, s))
	{
		const auto aux = splitStringInPoint(s);
		solution.insert(aux);
	}
}

unordered_map<double, unordered_set<std::pair<vector3, vector3>, PairHash>> getDistancesOfScanner(vector<Scanner>::const_reference scanner1)
{
	unordered_map<double, unordered_set<pair<vector3, vector3>, PairHash>> result;
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

Scanner connect2Scanners(vector<Scanner>::const_reference scanners1, vector<Scanner>::const_reference scanners2, const vector<vector<double>>& rotationMatrix, const vector3& shift)
{
	unordered_set<vector3> result;

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

	return Scanner{scanners1.id, true, vector<vector3>(result.begin(), result.end())};
}

pair<vector<pair<vector3, vector3>>, vector<pair<vector3, vector3>>>
extractCommonPoints(
	const unordered_map<double, unordered_set<pair<vector3, vector3>, PairHash>>& distance1,
	unordered_map<double, unordered_set<pair<vector3, vector3>, PairHash>>& distance2,
	int counter)
{

	unordered_map<vector3, int> pointFrequency1;
	unordered_map<vector3, int> pointFrequency2;

	unordered_map<double, pair<vector3, vector3>> vectors1;
	unordered_map<double, pair<vector3, vector3>> vectors2;

	for (auto& [distance, points] : distance1)
	{
		if (distance2.find(distance) != distance2.end())
		{
			if (auto aux = distance2[distance];
				points.size() == 1 and aux.size() == 1)
			{
				const auto firstVector = *(points.begin());
				const auto secondVector = *(aux.begin());
				vectors1[distance] = firstVector;
				vectors2[distance] = secondVector;
			}

		}
	}

	for (auto it = vectors1.begin(); it != vectors1.end();)
	{
		if (vectors2.find(it->first) == vectors2.end())
			it = vectors1.erase(it);
		else
			++it;
	}

	for (auto it = vectors2.begin(); it != vectors2.end();)
	{
		if (vectors1.find(it->first) == vectors1.end())
			it = vectors2.erase(it);
		else
			++it;
	}
	/*vectors1.erase(
		std::remove_if(
			vectors1.begin(),
			vectors1.end(),
			[&vectors2](const std::pair<double, std::pair<vector3, vector3>>& pair)
			{
				return (vectors2.find(pair.first) == vectors2.end()) == false;
			}
		),
		vectors1.end()
				);*/



	if (counter == 66)
	{
		vector<pair<vector3, vector3>> result1;
		result1.reserve(vectors1.size());

		for (const auto& [_, second] : vectors1)
		{
			result1.push_back(second);
		}

		vector<pair<vector3, vector3>> result2;
		result2.reserve(vectors2.size());

		for (const auto& [_, second] : vectors2)
		{
			result2.push_back(second);
		}

		return make_pair(result1, result2);
	}

	for (const auto& [first, second] : vectors1)
	{
		pointFrequency1[second.first]++;
		pointFrequency1[second.second]++;
	}

	for (const auto& [first, second] : vectors2)
	{
		pointFrequency2[second.first]++;
		pointFrequency2[second.second]++;
	}

	//std::map<int, int> mostProbableNumberOfPoints;
	//for (const auto & [_, second] : pointFrequency1)
	//{
	//	mostProbableNumberOfPoints[second]++;
	//}
	///*for (const auto& [_, second] : pointFrequency2)
	//{
	//	mostProbableNumberOfPoints[second]++;
	//}*/

	//auto max = std::max_element(mostProbableNumberOfPoints.begin(), mostProbableNumberOfPoints.end(),
	//	[](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
	//		return p1.second < p2.second; });

	//for (auto it = vectors1.begin(); it != vectors1.end();)
	//{
	//	if (pointFrequency1[it->second.first] != max->second - 1 or pointFrequency1[it->second.second] != max->second - 1)
	//		it = vectors1.erase(it);
	//	else
	//		++it;
	//}

	//for (auto it = vectors2.begin(); it != vectors2.end();)
	//{
	//	if (pointFrequency2[it->second.first] != max->second - 1 or pointFrequency2[it->second.second] != max->second - 1)
	//		it = vectors2.erase(it);
	//	else
	//		++it;
	//}


	for (auto it = vectors1.begin(); it != vectors1.end();)
	{
		if (vectors2.find(it->first) == vectors2.end())
			it = vectors1.erase(it);
		else
			++it;
	}

	for (auto it = vectors2.begin(); it != vectors2.end();)
	{
		if (vectors1.find(it->first) == vectors1.end())
			it = vectors2.erase(it);
		else
			++it;
	}
	{
		vector<pair<vector3, vector3>> result1;
		result1.reserve(vectors1.size());
		for (const auto& [_, second] : vectors1)
		{
			result1.push_back(second);
		}

		vector<pair<vector3, vector3>> result2;
		result2.reserve(vectors2.size());
		for (const auto& [_, second] : vectors2)
		{
			result2.push_back(second);
		}

		return make_pair(result1, result2);
	}
}

Scanner doWork(vector<Scanner>::const_reference scanner1, vector<Scanner>::const_reference scanner2)
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

		if(vectors1.empty() or vectors2.empty())
			return Scanner{ -1, false, vector<vector3>() };

		stable_sort(vectors1.begin(), vectors1.end(),
			[](const pair<vector3, vector3>& a, const pair<vector3, vector3>& b) -> bool
			{
				return getDistance(a.first, a.second) < getDistance(b.first, b.second);
			});

		stable_sort(vectors2.begin(), vectors2.end(),
			[](const pair<vector3, vector3>& a, const pair<vector3, vector3>& b) -> bool
			{
				return getDistance(a.first, a.second) < getDistance(b.first, b.second);
			});

		const auto rotationMatrix = getRotationMatrix3(vectors2, vectors1);

		const auto shiftMatrix = computeTranslation(vectors1, vectors2, rotationMatrix);

		auto result = connect2Scanners(scanner1, scanner2, rotationMatrix, shiftMatrix);

		return result;
	}

	return Scanner{ -1, false, vector<vector3>() };
}

Scanner solve1()
{
	vector<Scanner> scanners = readInput();

	Scanner scanner;
	bool firstPass = true;
	for (auto i = 0L; i < scanners.size() - 1; ++i)
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

void solve2()
{
	vector<Scanner> scanners = readInput();
	vector known{ scanners[0] };

	while (known.size() < scanners.size())
	{
		const auto size = known.size();
		for (auto i = 0ULL; i < size; ++i)
		{
			auto const& scanner = known[i];
			for (auto & value : scanners)
			{
				if(value.reference.has_value())
					continue;
				if(auto mapping = value.findCommonBeacons(scanner); mapping.first.size() == 12 and mapping.second.size() == 12)
				{
					value.computeTransformation(scanner, mapping);
					known.push_back(value);
					break;
				}
			}
		}
	}

	unordered_set<vector3> uniquePoints;

	vector<pair<vector<vector<double>>, vector3>> newTransformations(scanners.size());
	newTransformations.reserve(scanners.size());

	for ( auto i = 1ULL; i < scanners.size(); i++)
	{

		auto a = scanners[i].reference;

		if (a.value() == 0)
			continue;

		std::vector<int> references;
		references.push_back(i);

		while(a.value() != 0)
		{
			references.push_back(a.value());
			a = scanners[a.value()].reference;
		}

		//std::cout << references.size();

		vector3 newShift = scanners[references[0]].shift;
		vector<vector<double>> newRotationMatrix = scanners[references[0]].rotationMatrix;

		for (int j = 1; j < references.size(); ++j)
		{
			auto rotMat = scanners[references[j]].rotationMatrix;
			newRotationMatrix = multiplyMatrices(rotMat, newRotationMatrix);

			auto b = scanners[references[j]].shift;
			newShift = multiplyMatrices(scanners[references[j]].rotationMatrix, newShift);

			newShift.x += b.x;
			newShift.y += b.y;
			newShift.z += b.z;
			newShift.x = round(newShift.x);
			newShift.y = round(newShift.y);
			newShift.z = round(newShift.z);
		}

		
		newTransformations[i].first = newRotationMatrix;
		newTransformations[i].second = newShift;
	}

	for (auto i = 1ULL; i < scanners.size(); ++i)
	{
		if (newTransformations[i].first.empty() == false) 
		{
			scanners[i].rotationMatrix = newTransformations[i].first;
			scanners[i].shift = newTransformations[i].second;
		}
	}

	for (int i = 0; i < scanners.size(); ++i)
	{
		
		for (const auto & point : scanners[i].points)
		{
			if (i != 0) 
			{
				auto a = multiplyMatrices(scanners[i].rotationMatrix, point);
				auto b = scanners[i].shift;

				a.x += b.x;
				a.y += b.y;
				a.z += b.z;

				a.x = round(a.x);
				a.y = round(a.y);
				a.z = round(a.z);

				if (solution.find(a) == solution.end())
					std::cout << "";

				uniquePoints.insert(a);
			}
			else
			{
				if (solution.find(point) == solution.end())
					std::cout << "";

				uniquePoints.insert(point);
			}
		}
	}

	std::cout << uniquePoints.size();
}

int main()
{
	//readInput();
	readSolution();

	//const auto scanner = solve1();
	solve2();
	//std::cout << scanner.points.size() << '\n';
}
