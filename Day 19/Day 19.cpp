// Day 19.cpp : This file contains the 'main' function. Program execution begins and ends there.
// https://adventofcode.com/2021/day/19

#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <optional>
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

	pair<unordered_set<vector3>, unordered_set<vector3>> findCommonBeacons(Scanner const& destination);

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
	scanners.reserve(30);

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

void solve()
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

		vector3 newShift = scanners[references[0]].shift;
		vector<vector<double>> newRotationMatrix = scanners[references[0]].rotationMatrix;

		for (auto j = 1ULL; j < references.size(); ++j)
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

	for (int i = 0ULL; i < scanners.size(); ++i)
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

	double max = 0;
	for (const auto & scanner1 : scanners)
	{
		for (const auto & scanner2 : scanners)
		{
			if (const auto aux = manhattanDistance(scanner1.shift, scanner2.shift); aux > max)
				max = aux;
		}
	}
	std::cout << uniquePoints.size() << '\n' << max << '\n';
}

int main()
{
	using namespace std::chrono;

	readSolution();
	const auto start = high_resolution_clock::now();
	solve();
	const auto stop = high_resolution_clock::now();

	const auto elapsed = duration_cast<milliseconds>(stop - start).count();
	std::cout << std::setprecision(4) << elapsed << "ms\n";

}
