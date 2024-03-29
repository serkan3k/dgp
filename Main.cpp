#define HAVE_SINGLEPRECISION_MATH
#define _CRT_SECURE_NO_WARNINGS
#define COIN_DLL
#define SOWIN_DLL
#define HAVE_INT8_T
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Eigen/Dense>
#include "glm/glm.hpp"
#include "glm/ext/scalar_constants.inl"
#include "GMM.h"

#include "Mesh.h"
#include "Painter.h"
#include <string>
#include <cstring>
#include <queue>
#include <chrono>
#include <functional>
#include <cmath>
#include <fstream>
#include <stack>
#include <set>
#include <algorithm>
#include <random>

//using Eigen::MatrixXd;
using namespace Eigen;

class Ray
{
public:
	Ray();
	Ray(glm::vec3 origin, glm::vec3 direction) : Origin(origin), Direction(direction){}
	~Ray();
	glm::vec3 Origin,
			Direction;
};

Ray::Ray(){}

Ray::~Ray(){}

class RayHitInfo
{
public:
	RayHitInfo();
	~RayHitInfo();
	Ray Ray;
	glm::vec3 Normal,
		Wo,
		Wi,
		Point;
	float T, U, V;
	bool IsHit;
	int HitMaterialID;
};

RayHitInfo::RayHitInfo(){}


RayHitInfo::~RayHitInfo(){}


class BBox
{
public:
	glm::vec3 Min, Max, Center;
	BBox();
	~BBox();
	bool Intersect(const Ray& ray) const;
};

BBox::BBox(){}
BBox::~BBox(){}
bool BBox::Intersect(const Ray& ray) const
{
	const glm::vec3 inverseDir = 1.0f / ray.Direction;
	const glm::vec3 invMin = (Min - ray.Origin) * inverseDir;
	const glm::vec3 invMax = (Max - ray.Origin) * inverseDir;
	const glm::vec3 t0 = glm::min(invMin, invMax);
	const glm::vec3 t1 = glm::max(invMin, invMax);
	float tmin = t0.x;
	float tmax = t1.x;
	tmin = std::max(tmin, std::max(t0.y, t0.z));
	tmax = std::min(tmax, std::min(t1.y, t1.z));
	return tmin <= tmax;
}

class Object
{
public:
	Object();
	~Object();
	Triangle* Triangle;
	BBox BoundingBox;
};

Object::Object() {};
Object::~Object() {};


class BVH
{
public:
	BVH();
	BVH(std::vector<Object*>& triangles, const Mesh* mesh, int currentAxis);
	~BVH();
	bool Intersect(Ray& ray, const Mesh* mesh, double& tmin, RayHitInfo& rayHitInfo);
	BBox CalculateBoundingBox();
	BVH * LeftNode;
	BVH * RightNode;
	BBox BoundingBox;
	Triangle * ShapeObject;
	bool IsLeaf;
};

BVH::BVH()
{
	LeftNode = nullptr;
	RightNode = nullptr;
	ShapeObject = nullptr;
	BoundingBox = BBox();
}
BVH::~BVH(){}

bool compareBBoxX(Object* o1, Object* o2) {
	return (o1->BoundingBox.Center.x < o2->BoundingBox.Center.x);
}

bool compareBBoxY(Object* o1, Object* o2) {
	return (o1->BoundingBox.Center.y < o2->BoundingBox.Center.y);
}

bool compareBBoxZ(Object* o1, Object* o2) {
	return (o1->BoundingBox.Center.z < o2->BoundingBox.Center.z);
}

BVH::BVH(std::vector<Object*>& triangles, const Mesh* mesh, int currentAxis)
{
	if(triangles.size() == 1)
	{
		ShapeObject = triangles[0]->Triangle;
		BoundingBox = triangles[0]->BoundingBox;
		IsLeaf = true;
		return;
	}
	IsLeaf = false;
	int axis = (currentAxis % 3);
	auto trianglesSorted = triangles;
	if(axis == 0){
		std::sort(trianglesSorted.begin(), trianglesSorted.end(), compareBBoxX);
	}
	else if(axis == 1){
		std::sort(trianglesSorted.begin(), trianglesSorted.end(), compareBBoxY);
	}
	else{
		std::sort(trianglesSorted.begin(), trianglesSorted.end(), compareBBoxZ);
	}
	std::vector<Object*> leftObjects, rightObjects;
	
	for (auto it = trianglesSorted.begin(); it != trianglesSorted.begin() + trianglesSorted.size() / 2; ++it)
	{
		leftObjects.push_back(*it);	
	}
	for (auto it = trianglesSorted.begin() + trianglesSorted.size() / 2; it != trianglesSorted.end(); ++it) {
		rightObjects.push_back(*it);
	}
	LeftNode = new BVH(leftObjects, mesh, axis + 1);
	RightNode = new BVH(rightObjects, mesh, axis + 1);
	BoundingBox.Min.x = std::min(LeftNode->BoundingBox.Min.x, RightNode->BoundingBox.Min.x);
	BoundingBox.Min.y = std::min(LeftNode->BoundingBox.Min.y, RightNode->BoundingBox.Min.y);
	BoundingBox.Min.z = std::min(LeftNode->BoundingBox.Min.z, RightNode->BoundingBox.Min.z);
	BoundingBox.Max.x = std::max(LeftNode->BoundingBox.Max.x, RightNode->BoundingBox.Max.x);
	BoundingBox.Max.y = std::max(LeftNode->BoundingBox.Max.y, RightNode->BoundingBox.Max.y);
	BoundingBox.Max.z = std::max(LeftNode->BoundingBox.Max.z, RightNode->BoundingBox.Max.z);
	BoundingBox.Center.x = (BoundingBox.Min.x + BoundingBox.Max.x) * 0.5f;
	BoundingBox.Center.y = (BoundingBox.Min.y + BoundingBox.Max.y) * 0.5f;
	BoundingBox.Center.z = (BoundingBox.Min.z + BoundingBox.Max.z) * 0.5f;
	return;
}

bool BVH::Intersect(Ray &ray, const Mesh* mesh, double &tmin, RayHitInfo& rayHitInfo)
{
	if (BoundingBox.Intersect(ray))
	{
		if (IsLeaf)
		{
			const glm::vec3 v1(mesh->verts[ShapeObject->v1i]->coords[0],
						mesh->verts[ShapeObject->v1i]->coords[1],
						mesh->verts[ShapeObject->v1i]->coords[2]);
			const glm::vec3 v2(mesh->verts[ShapeObject->v2i]->coords[0],
				mesh->verts[ShapeObject->v2i]->coords[1],
				mesh->verts[ShapeObject->v2i]->coords[2]);
			const glm::vec3 v3(mesh->verts[ShapeObject->v3i]->coords[0],
				mesh->verts[ShapeObject->v3i]->coords[1],
				mesh->verts[ShapeObject->v3i]->coords[2]);
			const glm::vec3 v1v2 = v2 - v1;
			const glm::vec3 v1v3 = v3 - v1;
			const glm::vec3 p = glm::cross(ray.Direction, v1v3);
			const double d = glm::dot(v1v2, p);
			if(abs(d) < 0){ return false; }
			const glm::vec3 t = ray.Origin - v1;
			const double u = glm::dot(t, p) * (1.0 / d);
			if(u < 0.0f || u > 1.0f){return false;}
			const glm::vec3 q = glm::cross(t, v1v2);
			const double v = glm::dot(ray.Direction, q) * (1.0f / d);
			if (v < 0.0f || u + v > 1.0f) { return false; }
			double dist = glm::dot(v1v3, q) * (1.0 / d);
			if (dist < 0) { return false; }
			if(dist <= tmin)
			{
				glm::vec3 normalAtIntersection = -glm::normalize(glm::cross(v1v2, v1v3));	//inwards facing normal
				float dotBetween = glm::dot(ray.Direction, normalAtIntersection);
				if(glm::degrees(glm::acos(dotBetween)) >= 90.0f)	// termination condition for same direction facing normals
				{
					// add termination condition for rays (from paper)	
					tmin = dist;
					rayHitInfo.T = tmin;
					return true;
				}
			}
			return false;
			//return ShapeObject->Shape->Intersect(ray, tmin, rayHitInfo);
		}
		RayHitInfo leftHitInfo, rightHitInfo;
		bool leftHit = LeftNode == nullptr ? false : LeftNode->Intersect(ray, mesh, tmin, leftHitInfo);
		bool rightHit = RightNode == nullptr ? false : RightNode->Intersect(ray, mesh, tmin, rightHitInfo);
		if (leftHit && rightHit)
		{
			rayHitInfo = leftHitInfo.T < rightHitInfo.T ? leftHitInfo : rightHitInfo;
			tmin = rayHitInfo.T;
			return true;
		}
		rayHitInfo = leftHit ? leftHitInfo : rayHitInfo;
		rayHitInfo = rightHit ? rightHitInfo : rayHitInfo;
		if (leftHit || rightHit)
		{
			tmin = rayHitInfo.T;
		}
		return (leftHit || rightHit);
	}
	return false;
}

BBox BVH::CalculateBoundingBox()
{
	return BBox();
}



int main(int, char ** argv)
{

	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);
	SoSeparator * root = new SoSeparator;
	root->ref();
	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();
	// load mesh
	char* x = (char*)malloc(strlen("man.off") + 1); 
	strcpy(x, "man.off");
	mesh->loadOff(x);

	const int numVertices = mesh->verts.size();

	vector<glm::vec3> normals(mesh->verts.size());
	vector<glm::vec3> centroids(mesh->verts.size());
	vector<glm::vec3> directions(mesh->verts.size());
	vector<glm::vec3> negatedNormals(mesh->verts.size());
	vector<vector<Ray>> rays(mesh->verts.size());

	//vector<glm::vec3> normals(mesh->tris.size());
	//vector<glm::vec3> centroids(mesh->tris.size());
	//vector<glm::vec3> directions(mesh->tris.size());
	//vector<glm::vec3> negatedNormals(mesh->tris.size());
	//vector<vector<Ray>> rays(mesh->tris.size());

	for (unsigned i = 0; i < normals.size(); ++i) {
		normals[i] = glm::vec3(0,0,0);
	}
	std::vector<BBox> boundingBoxes(mesh->tris.size());
	std::vector<Object*> objects(mesh->tris.size());
#pragma omp parallel for
	for (int i = 0; i < mesh->tris.size(); ++i) 
	{
		int v1Index = (mesh->tris[i]->v1i);
		int v2Index = (mesh->tris[i]->v2i);
		int v3Index = (mesh->tris[i]->v3i);
		glm::vec3 v1(mesh->verts[v1Index]->coords[0], mesh->verts[v1Index]->coords[1], mesh->verts[v1Index]->coords[2]);
		glm::vec3 v2(mesh->verts[v2Index]->coords[0], mesh->verts[v2Index]->coords[1], mesh->verts[v2Index]->coords[2]);
		glm::vec3 v3(mesh->verts[v3Index]->coords[0], mesh->verts[v3Index]->coords[1], mesh->verts[v3Index]->coords[2]);
		glm::vec3 v1v2 = v2 - v1;
		glm::vec3 v1v3 = v3 - v1;
		glm::vec3 p = glm::cross(v1v2, v1v3);
		normals[v1Index] += p;
		normals[v2Index] += p;
		normals[v3Index] += p;
		BBox triangleBox; // bounding box for acceleration
		triangleBox.Center = glm::vec3((v1.x + v2.x + v3.x) / 3.0,
			(v1.y + v2.y + v3.y) / 3.0,
			(v1.z + v2.z + v3.z) / 3.0f);
		triangleBox.Min = glm::vec3(std::min(std::min(v1.x, v2.x), v3.x),
			std::min(std::min(v1.y, v2.y), v3.y),
			std::min(std::min(v1.z, v2.z), v3.z));
		triangleBox.Max = glm::vec3(std::max(std::max(v1.x, v2.x), v3.x),
			std::max(std::max(v1.y, v2.y), v3.y),
			std::max(std::max(v1.z, v2.z), v3.z));
		boundingBoxes[i] = triangleBox;
		objects[i] = new Object();
		objects[i]->Triangle = mesh->tris[i];
		objects[i]->BoundingBox = boundingBoxes[i];
	}

	BVH bvhNode(objects, mesh, 0);

	for (unsigned i = 0; i < normals.size(); ++i) {
		normals[i] = glm::normalize(normals[i]);
	}

	std::mt19937 generator;
	std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
#pragma omp parallel for
	for(int i = 0; i < mesh->verts.size(); ++i)
	{
		glm::vec3 normal(0, 0, 0);
		glm::vec3 vx(mesh->verts[i]->coords[0], mesh->verts[i]->coords[1], mesh->verts[i]->coords[2]);
		centroids[i] = vx;
		normal = glm::normalize(normal);
		negatedNormals[i] = -normals[i];
		rays[i].resize(30);
		for(unsigned j = 0; j < rays[i].size(); ++j)
		{
			// rejection sampling, over the unit sphere, reject the ones that are not in cone
			bool accepted = false;
			while(!accepted)
			{
				// take samples over the unit sphere
				float z = distribution(generator) * 2.0f - 1.0f;	// z uniformly distributed btw [-1, 1]
				float t = distribution(generator) * 2.0f * glm::pi<float>();	// t uniformly distributed btw [0, 2pi)
				float r = sqrt(1.0f - z * z);
				float xx = r * cos(t);
				float y = r * sin(t);
				glm::vec3 sampledVec(xx, y, z);
				float dotProduct = glm::dot(glm::normalize(negatedNormals[i]), glm::normalize(sampledVec));
				float angleBetween = glm::acos(dotProduct);	// vector lengths are 1 since they are normals & normalized
				if(glm::degrees(angleBetween) <= 60.0f)	//accept the sample
				{
					rays[i][j].Direction = glm::normalize(sampledVec);
					accepted = true;
				}
			}
			rays[i][j].Origin = centroids[i] + rays[i][j].Direction * (float)1e-6;	// some epsilon for robustness in intersection tests
		}
		directions[i] = rays[i][0].Direction;
		centroids[i] = rays[i][0].Origin;
	}
	vector<vector<float>> rayDistances(rays.size());
	#pragma omp parallel for
	for(int i = 0; i < rays.size(); ++i)
	{
		for(unsigned j = 0; j < rays[i].size(); ++j)
		{
			//bool isHit = false;
			double tmin = DBL_MAX;
			RayHitInfo info;
			if(bvhNode.Intersect(rays[i][j], mesh, tmin, info))
			{
				//isHit = true;
				rayDistances[i].push_back(tmin);
			}
			//for(unsigned k = 0; k < mesh->tris.size(); ++k)
			//{
			//	const glm::vec3 v1(mesh->verts[mesh->tris[k]->v1i]->coords[0],
			//		mesh->verts[mesh->tris[k]->v1i]->coords[1],
			//		mesh->verts[mesh->tris[k]->v1i]->coords[2]);
			//	const glm::vec3 v2(mesh->verts[mesh->tris[k]->v2i]->coords[0],
			//		mesh->verts[mesh->tris[k]->v2i]->coords[1],
			//		mesh->verts[mesh->tris[k]->v2i]->coords[2]);
			//	const glm::vec3 v3(mesh->verts[mesh->tris[k]->v3i]->coords[0],
			//		mesh->verts[mesh->tris[k]->v3i]->coords[1],
			//		mesh->verts[mesh->tris[k]->v3i]->coords[2]);
			//	const glm::vec3 v1v2 = v2 - v1;
			//	const glm::vec3 v1v3 = v3 - v1;
			//	const glm::vec3 p = glm::cross(rays[i][j].Direction, v1v3);
			//	const double d = glm::dot(v1v2, p);
			//	if(abs(d) < 0){ continue; }
			//	const glm::vec3 t = rays[i][j].Origin - v1;
			//	const double u = glm::dot(t, p) * (1.0 / d);
			//	if(u < 0.0f || u > 1.0f){continue;}
			//	const glm::vec3 q = glm::cross(t, v1v2);
			//	const double v = glm::dot(rays[i][j].Direction, q) * (1.0f / d);
			//	if (v < 0.0f || u + v > 1.0f) { continue; }
			//	double dist = glm::dot(v1v3, q) * (1.0 / d);
			//	if (dist < 0) { continue; }
			//	if(dist <= tmin)
			//	{
			//		glm::vec3 normalAtIntersection = -glm::normalize(glm::cross(v1v2, v1v3));	//inwards facing normal
			//		float dotBetween = glm::dot(rays[i][j].Direction, normalAtIntersection);
			//		if(glm::degrees(glm::acos(dotBetween)) >= 90.0f)	// termination condition for same direction facing normals
			//		{
			//			// add termination condition for rays (from paper)	
			//			isHit = true;
			//			tmin = dist;
			//		}
			//	}		
			//}
			//if(isHit)
			//{
			//	rayDistances[i].push_back(tmin);
			//}
		}
	}
	vector<float> sdf;
	for(unsigned i = 0; i < rayDistances.size(); ++i)
	{
		if(rayDistances[i].empty())
		{
			sdf.push_back(0.0f);
			continue;
		}
		else if(rayDistances[i].size() == 1)
		{
			sdf.push_back(rayDistances[i][0]);
			continue;
		}
		std::sort(rayDistances[i].begin(), rayDistances[i].end());
		float mean = 0;
		float variance = 0;
		for(unsigned j = 0; j < rayDistances[i].size(); ++j)
		{
			mean += rayDistances[i][j] / (float)rayDistances[i].size();
		}
		for(unsigned j = 0; j < rayDistances[i].size(); ++j)
		{
			variance += (rayDistances[i][j] - mean) * (rayDistances[i][j] - mean);
		}
		variance /= (float)rayDistances[i].size();
		float stdDev = glm::sqrt(variance);
		float rightMax = rayDistances[i][rayDistances[i].size() / 2] + stdDev; // 1 std dev away from median
		float leftMax = std::max(0.0f, rayDistances[i][rayDistances[i].size() / 2] - stdDev);
		float sum = 0.0f;
		float numElementsInSum = 0.0f;
		for(unsigned j = 0; j < rayDistances[i].size(); ++j)
		{
			if(rayDistances[i][j] <= rightMax && rayDistances[i][j] >= leftMax)
			{
				numElementsInSum += 1.0f;
				sum += rayDistances[i][j];
			}
		}
		sdf.push_back(sum / numElementsInSum);
	}
	float minSdf = FLT_MAX;
	float maxSdf = FLT_MIN;
	for(unsigned i = 0; i < sdf.size(); ++i){
		if(sdf[i] <= minSdf){
			minSdf = sdf[i];
		}
		if (sdf[i] >= maxSdf) {
			maxSdf = sdf[i];
		}
	}
	vector<float> nsdf;
	const float sdfAlpha = 4.0f;
	ofstream nsdfFile;
	nsdfFile.open("nsdf.csv");
	for(unsigned i = 0; i < sdf.size(); ++i)
	{
		float normalized = glm::log(((sdf[i] - minSdf) / (maxSdf - minSdf)) * sdfAlpha + 1.0f) / glm::log(sdfAlpha + 1.0f);
		nsdf.push_back(normalized);
		nsdfFile << normalized << std::endl;
	}
	nsdfFile.close();

	int dimension = 1;
	int numData = nsdf.size();
	int numIterations = 200;
	int numGaussianComponents = 3;

	Gaussian_Mixture_Model GMM = Gaussian_Mixture_Model("full", dimension, numGaussianComponents);
	double **data = new double*[numData];
	for (int i = 0; i < numData; ++i) {
		data[i] = new double[dimension];
		data[i][0] = nsdf[i];
	}
	for (int i = 0; i < numIterations; ++i) {
		double logLikelihood;
		if (i == 0) GMM.Initialize(numData, data);
		logLikelihood = GMM.Expectaion_Maximization(numData, data);
		if ((i + 1) % 10 == 0) { 
			cout << i + 1 << " , " << logLikelihood << std::endl;
		}
	}
	cout << std::endl <<"mean" << std::endl;
	for (int i = 0; i < numGaussianComponents; ++i) {
		for (int j = 0; j < dimension; ++j) {
			cout << GMM.mean[i][j];
		}
		cout << std::endl;
	}
	ofstream gmmfile;
	gmmfile.open("gmmresult.txt");
	std::vector<int> nsdfSegments(nsdf.size());
	/*for (int j = 0; j < numGaussianComponents; ++j) {
		for (int i = 0; i < numData; ++i) {
			if (GMM.Classify(data[i]) == j) {
				nsdfSegments[i] = GMM.Classify(data[i]);
				gmmfile << data[i][0] << " , " << GMM.Classify(data[i]) << std::endl;
 			}
		}
	}*/
	for (int i = 0; i < numData; ++i) {
		for (int j = 0; j < numGaussianComponents; ++j) {	
			if (GMM.Classify(data[i]) == j) {
				nsdfSegments[i] = GMM.Classify(data[i]);
				gmmfile << data[i][0] << " , " << j << std::endl;
			}
		}
	}
	gmmfile.close();
	for (int i = 0; i < numData; ++i) {
		delete[] data[i];
	}
	delete[] data;
	
	vector<int> histogramBins(50);
	for (unsigned i = 0; i < histogramBins.size(); ++i) {
		histogramBins[i] = 0;
	}
	const double binStep = 1.0 / (double)histogramBins.size();
	for (unsigned i = 0; i < nsdf.size(); ++i) {
		auto idx = std::min(glm::floor(nsdf[i] / binStep), histogramBins.size() - 1.0);
		histogramBins[idx] += 1;
	}
	ofstream histogramFile;
	histogramFile.open("histogram.csv");
	for (unsigned i = 0; i < histogramBins.size(); ++i) {
		histogramFile << (i + 1) * binStep << "," << histogramBins[i] << std::endl;
	}
	histogramFile.close();

	/*
	cout << "------------------" << endl << "Dijkstra" << endl << "------------------" << endl;
#pragma region dijkstraQuery
	
	int dijkstraQueryFirst = -1;
	int dijkstraQuerySecond = -1;
	while (dijkstraQueryFirst < 0 || dijkstraQueryFirst >= numVertices) {
		cout << "Enter first vertex index, [0," << numVertices - 1 << "] :";
		cin >> dijkstraQueryFirst;
		if(dijkstraQueryFirst < 0 || dijkstraQueryFirst >= numVertices){
			cout << "Invalid index, try again: ";
		}
	}
	while (dijkstraQuerySecond < 0 || dijkstraQuerySecond >= numVertices) {
		cout << "Enter second vertex index, [0," << numVertices-1 << "] :";
		cin >> dijkstraQuerySecond;
		if (dijkstraQuerySecond < 0 || dijkstraQuerySecond >= numVertices) {
			cout << "Invalid index, try again: ";
		}
	}
	cout << dijkstraQueryFirst << " " << dijkstraQuerySecond << endl;
	
#pragma endregion 
*/
#pragma region array
	/*
	distances.clear();
	distances.resize(numVertices);
	for(auto &d:distances){
		d.resize(numVertices);
	}
	parents.clear();
	parents.resize(numVertices);
	t0 = chrono::high_resolution_clock::now();
	for (int i = 0; i < numVertices; ++i) {
		float * arr = nullptr;
		arr = new float[numVertices];
		for(int j = 0; j < numVertices; ++j){
			if (i == j) arr[j] = 0;
			else arr[j] = FLT_MAX;
		}
		std::vector<int> parent(numVertices, -1);
		parent[i] = i;
		std::vector<bool> visited(numVertices, false);
		while (true)
		{
			float minDist = FLT_MAX;
			int minDistIndex = -1;
			for(int k = 0; k < numVertices; k++){
				if(!visited[k] && arr[k] < minDist){
					minDist = arr[k];
					minDistIndex = k;
				}
			}
			if (minDistIndex == -1) break;
			int u = minDistIndex;
			visited[u] = true;
			for (int l = 0; l < mesh->verts[u]->vertList.size(); ++l){
				int v = mesh->verts[u]->vertList[l];
				auto v1 = mesh->verts[u]->coords;
				auto v2 = mesh->verts[v]->coords;
				float weight = sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) +
					(v1[1] - v2[1]) * (v1[1] - v2[1]) +
					(v1[2] - v2[2]) * (v1[2] - v2[2]));
				if (!visited[v] && (arr[v] > arr[u] + weight)){
					arr[v] = arr[u] + weight;
					parent[v] = u;
				}
			}
		}
		for(int m = 0; m < numVertices; ++m){
			distances[i][m] = arr[m];
		}
		parents[i] = parent;
		delete [] arr;
	}
	t1 = chrono::high_resolution_clock::now();
	duration= chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << "Array: " << duration << " seconds"  << endl;
	*/
#pragma endregion 
	
#pragma region minHeap
	/*
	std::vector<std::vector<float>> distances(numVertices);
	std::vector<std::vector<int>> parents(numVertices);
	chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
	for (int i = 0; i < numVertices; ++i) {
		std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<>> pq;
		int source = i;
		std::vector<int> parent(numVertices, -1);
		parent[source] = source;
		std::vector<float> dist(numVertices, FLT_MAX);
		dist[source] = 0;
		std::vector<bool> visited(numVertices, false);
		visited[source] = true;
		pq.push(std::make_pair(0, source));
		while (!pq.empty())
		{
			auto top = pq.top();
			int u = top.second;
			visited[u] = true;
			pq.pop();
			for (int j = 0; j < mesh->verts[u]->vertList.size(); ++j)
			{
				int v = mesh->verts[u]->vertList[j];
				auto v1 = mesh->verts[u]->coords;
				auto v2 = mesh->verts[v]->coords;
				float weight = sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) +
					(v1[1] - v2[1]) * (v1[1] - v2[1]) +
					(v1[2] - v2[2]) * (v1[2] - v2[2]));
				if (!visited[v] && (dist[v] > dist[u] + weight))
				{
					dist[v] = dist[u] + weight;
					pq.push(make_pair(dist[v], v));
					parent[v] = u;
				}
			}
		}
		distances[i] = dist;
		parents[i] = parent;
	}
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << "Dijkstra: " << duration << " seconds" << endl;
	*/
#pragma endregion 
	/*
	std::vector<int> boundaryIndices;
	const int numEdges = mesh->edges.size();
	const int numTris = mesh->tris.size();
	const auto edges = mesh->edges;
	const auto tris = mesh->tris;
	const auto verts = mesh->verts;
	std::vector<bool> isVertexBoundary(numVertices, false);
	int i = 0;
	int belongsTo = 0;
	t0 = chrono::high_resolution_clock::now();
	//#pragma omp parallel for private (i) 
	for (i = 0; i < numEdges; ++i) {
		belongsTo = 0;
		const int ev1 = edges[i]->v1i;
		const int ev2 = edges[i]->v2i;
		//#pragma omp parallel for reduction(+: belongsTo)
		for (int j = 0; j < numTris; ++j) {
			const int tv1 = tris[j]->v1i;
			const int tv2 = tris[j]->v2i;
			const int tv3 = tris[j]->v3i;
			if (ev1 == tv1) {
				if (ev2 == tv2 || ev2 == tv3) {
					belongsTo++;
					continue;
				}
			}
			else if (ev1 == tv2) {
				if (ev2 == tv1 || ev2 == tv2) {
					belongsTo++;
					continue;
				}
			}
			else if (ev1 == tv3) {
				if (ev2 == tv1 || ev2 == tv3) {
					belongsTo++;
					continue;
				}
			}
			if (ev2 == tv1) {
				if (ev1 == tv2 || ev1 == tv3) {
					belongsTo++;
				}
			}
			else if (ev2 == tv2) {
				if (ev1 == tv1 || ev1 == tv3) {
					belongsTo++;
				}
			}
			else if (ev2 == tv3) {
				if (ev1 == tv1 || ev1 == tv2) {
					belongsTo++;
				}
			}
		}
		if (belongsTo == 1) {
			if (!isVertexBoundary[ev1]) boundaryIndices.push_back(ev1);
			if (!isVertexBoundary[ev2]) boundaryIndices.push_back(ev2);
			isVertexBoundary[ev1] = true;
			isVertexBoundary[ev2] = true;
		}
	}
	t1 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << "Boundary edges: " << duration << " seconds" << endl;

	std::vector<int> firstSetIndices;
	std::set<int> firstSet;
	firstSetIndices.push_back(boundaryIndices[0]);
	firstSet.insert(boundaryIndices[0]);
	while(firstSetIndices.size() < boundaryIndices.size())
	{
		const auto top = firstSetIndices[firstSetIndices.size() - 1];
		const auto neighbours = mesh->verts[top]->vertList;
		std::vector<int> candidateVertices;
		for(int i = 0; i < boundaryIndices.size(); ++i){
			if (boundaryIndices[i] == top) continue;
			bool isVertexCandidate = false;
			for(int j = 0; j < neighbours.size(); ++j){
				if(neighbours[j] == boundaryIndices[i]){
					candidateVertices.push_back(boundaryIndices[i]);
					break;
				}
			}
		}
		float minDist = FLT_MAX;
		int nextIdx = -1;
		if (candidateVertices.size() == 0) break;
		for(int i = 0; i < candidateVertices.size(); ++i){
			if(distances[top][candidateVertices[i]] < minDist){
				bool isContained = false;
				for(int j = 0; j < firstSetIndices.size(); ++j){
					if(firstSetIndices[j] == candidateVertices[i]){
						isContained = true;
						break;
					}
				}
				if (isContained == false) {
					minDist = distances[top][candidateVertices[i]] < minDist;
					nextIdx = candidateVertices[i];
				}
			}
		}
		if (nextIdx == -1) break;
		firstSetIndices.push_back(nextIdx);
		firstSet.insert(nextIdx);
	}
	vector<int> secondSetIndices;
	for(int i = 0; i < boundaryIndices.size(); ++i){
		if(firstSet.find(boundaryIndices[i]) == firstSet.end()){
			secondSetIndices.push_back(boundaryIndices[i]);
		}
	}
	bool useSymmetricTriangles = false;
	if (firstSetIndices.size() != boundaryIndices.size()) {
		int symmetric = -1;
		while (symmetric != 0 && symmetric != 1) {
			std::cout << "The mesh contains holes. Should I use symmetric triangles? 0 -> No, 1 -> Yes" << std::endl;
			cin >> symmetric;
		}
		if (symmetric == 0) useSymmetricTriangles = false;
		else if (symmetric == 1) useSymmetricTriangles = true;
		float bboxFirstMin[3], bboxFirstMax[3], bboxSecondMin[3], bboxSecondMax[3];
		bboxFirstMin[0] = FLT_MAX; bboxFirstMin[1] = FLT_MAX; bboxFirstMin[2] = FLT_MAX;
		bboxFirstMax[0] = FLT_MIN; bboxFirstMax[1] = FLT_MIN; bboxFirstMax[2] = FLT_MIN;
		for (int i = 0; i < firstSetIndices.size(); ++i) {
			const auto coords = mesh->verts[firstSetIndices[i]]->coords;
			bboxFirstMin[0] = min(bboxFirstMin[0], coords[0]);
			bboxFirstMin[1] = min(bboxFirstMin[1], coords[1]);
			bboxFirstMin[2] = min(bboxFirstMin[2], coords[2]);
			bboxFirstMax[0] = max(bboxFirstMax[0], coords[0]);
			bboxFirstMax[1] = max(bboxFirstMax[1], coords[1]);
			bboxFirstMax[2] = max(bboxFirstMax[2], coords[2]);
		}
		bboxSecondMin[0] = FLT_MAX; bboxSecondMin[1] = FLT_MAX; bboxSecondMin[2] = FLT_MAX;
		bboxSecondMax[0] = FLT_MIN; bboxSecondMax[1] = FLT_MIN; bboxSecondMax[2] = FLT_MIN;
		for (int i = 0; i < secondSetIndices.size(); ++i) {
			const auto coords = mesh->verts[firstSetIndices[i]]->coords;
			bboxSecondMin[0] = min(bboxSecondMin[0], coords[0]);
			bboxSecondMin[1] = min(bboxSecondMin[1], coords[1]);
			bboxSecondMin[2] = min(bboxSecondMin[2], coords[2]);
			bboxSecondMax[0] = max(bboxSecondMax[0], coords[0]);
			bboxSecondMax[1] = max(bboxSecondMax[1], coords[1]);
			bboxSecondMax[2] = max(bboxSecondMax[2], coords[2]);
		}

		const auto bboxFirstXSqr = pow(bboxFirstMin[0] - bboxFirstMax[0], 2.0);
		const auto bboxFirstYSqr = pow(bboxFirstMin[1] - bboxFirstMax[1], 2.0);
		const auto bboxFirstZSqr = pow(bboxFirstMin[2] - bboxFirstMax[2], 2.0);
		const auto bboxFirstLen = sqrt(bboxFirstXSqr + bboxFirstYSqr + bboxFirstZSqr);

		const auto bboxSecondXSqr = pow(bboxSecondMin[0] - bboxSecondMax[0], 2.0);
		const auto bboxSecondYSqr = pow(bboxSecondMin[1] - bboxSecondMax[1], 2.0);
		const auto bboxSecondZSqr = pow(bboxSecondMin[2] - bboxSecondMax[2], 2.0);
		const auto bboxSecondLen = sqrt(bboxSecondXSqr + bboxSecondYSqr + bboxSecondZSqr);
		//
		int abandonSecondaryVertices = -1;
		std::cout << "Should I abandon secondary vertices? 0 -> No, 1 -> Yes" << endl;
		while(abandonSecondaryVertices != 0 && abandonSecondaryVertices != 1)
		{
			cin >> abandonSecondaryVertices;
		}
		if (bboxFirstLen > bboxSecondLen)
		{
			boundaryIndices = firstSetIndices;
			if (abandonSecondaryVertices == 1) {
				for (int j = 0; j < secondSetIndices.size(); ++j) {
					isVertexBoundary[secondSetIndices[j]] = false;
				}
			}
		}
		else
		{
			boundaryIndices = secondSetIndices;
			if (abandonSecondaryVertices == 1) {
				for (int j = 0; j < firstSetIndices.size(); ++j) {
					isVertexBoundary[firstSetIndices[j]] = false;
				}
			}
		}
		
	}
	
	MatrixXd w(numVertices, numVertices), xx(numVertices, 1), bx(numVertices, 1),
		xy(numVertices, 1), by(numVertices, 1);
	t0 = chrono::high_resolution_clock::now();

	int mode = -1;
	std::cout << "Enter mode: 0 -> Uniform, 1 -> Harmonic, 2 -> Mean Value" << endl;
	while(mode != 0 && mode != 1 && mode != 2)
	{
		cin >> mode;
	}
	if (mode == 0) {
#pragma region uniform
		for (int i = 0; i < numVertices; i++) {
			for (int j = 0; j < numVertices; j++) {
				if (isVertexBoundary[i]) {
					w(i, j) = i == j ? 1 : 0;
					continue;
				}
				if (i == j) {
					w(i, j) = (double)verts[i]->vertList.size() * -1.0;
					continue;
				}
				w(i, j) = 0;
				for (const auto &k : verts[i]->vertList) {
					if (k == j) {
						w(i, j) = 1;
						break;
					}
				}
			}
		}
#pragma endregion
	}
	else if (mode == 1) {
#pragma region harmonic
		std::vector<int> nonBoundaryIndices;
		for (int i = 0; i < numVertices; ++i) {
			if (!isVertexBoundary[i]) {
				nonBoundaryIndices.push_back(i);
			}
			for (int j = 0; j < numVertices; ++j) {
				if (isVertexBoundary[i]) {
					if (i == j) {
						w(i, j) = 1;
					}
					else {
						w(i, j) = 0;
					}
				}
				else {
					w(i, j) = 0;
				}
			}
		}
		for (int i = 0; i < nonBoundaryIndices.size(); ++i) {
			const auto vertexId = nonBoundaryIndices[i];
			const auto neighbouringEdges = mesh->verts[vertexId]->edgeList;
			const auto neighbouringTris = mesh->verts[vertexId]->triList;
			if (neighbouringTris.size() < 2) {
				std::cout << "error on vertex " << vertexId << " : non-boundary vertex belongs to less than 2 triangles" << std::endl;
				continue;
			}
			for (int j = 0; j < neighbouringEdges.size(); ++j) {
				std::vector<int> triangleIndices;
				const int ev1 = mesh->edges[neighbouringEdges[j]]->v1i;
				const int ev2 = mesh->edges[neighbouringEdges[j]]->v2i;
				for (int k = 0; k < neighbouringTris.size(); ++k) {
					const int tv1 = mesh->tris[neighbouringTris[k]]->v1i;
					const int tv2 = mesh->tris[neighbouringTris[k]]->v2i;
					const int tv3 = mesh->tris[neighbouringTris[k]]->v3i;
					if (ev1 == tv1) {
						if (ev2 == tv2 || ev2 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
							continue;
						}
					}
					else if (ev1 == tv2) {
						if (ev2 == tv1 || ev2 == tv2) {
							triangleIndices.push_back(neighbouringTris[k]);
							continue;
						}
					}
					else if (ev1 == tv3) {
						if (ev2 == tv1 || ev2 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
							continue;
						}
					}
					if (ev2 == tv1) {
						if (ev1 == tv2 || ev1 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
						}
					}
					else if (ev2 == tv2) {
						if (ev1 == tv1 || ev1 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
						}
					}
					else if (ev2 == tv3) {
						if (ev1 == tv1 || ev1 == tv2) {
							triangleIndices.push_back(neighbouringTris[k]);
						}
					}
				}
				if (triangleIndices.size() == 2) {
					auto firstTriangleIdx = mesh->tris[triangleIndices[0]];
					auto secondTriangleIdx = mesh->tris[triangleIndices[1]];
					int firstTriangleOtherVertexIdx = -1;
					int secondTriangleOtherVertexIdx = -1;
					if (firstTriangleIdx->v1i == ev1) {
						if (firstTriangleIdx->v2i == ev2) {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
						}
						else {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
						}
					}
					else if (firstTriangleIdx->v2i == ev1) {
						if (firstTriangleIdx->v1i == ev2) {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
						}
						else {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
						}
					}
					else {
						if (firstTriangleIdx->v2i == ev2) {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
						}
						else {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
						}
					}
					if (secondTriangleIdx->v1i == ev1) {
						if (secondTriangleIdx->v2i == ev2) {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v3i;
						}
						else {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v2i;
						}
					}
					else if (secondTriangleIdx->v2i == ev1) {
						if (secondTriangleIdx->v1i == ev2) {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v3i;
						}
						else {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v1i;
						}
					}
					else {
						if (secondTriangleIdx->v2i == ev2) {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v1i;
						}
						else {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v2i;
						}
					}
					auto x1 = mesh->verts[ev1]->coords[0] - mesh->verts[firstTriangleOtherVertexIdx]->coords[0];
					auto y1 = mesh->verts[ev1]->coords[1] - mesh->verts[firstTriangleOtherVertexIdx]->coords[1];
					auto z1 = mesh->verts[ev1]->coords[2] - mesh->verts[firstTriangleOtherVertexIdx]->coords[2];
					auto x2 = mesh->verts[ev2]->coords[0] - mesh->verts[firstTriangleOtherVertexIdx]->coords[0];
					auto y2 = mesh->verts[ev2]->coords[1] - mesh->verts[firstTriangleOtherVertexIdx]->coords[1];
					auto z2 = mesh->verts[ev2]->coords[2] - mesh->verts[firstTriangleOtherVertexIdx]->coords[2];
					auto dot = x1 * x2 + y1 * y2 + z1 * z2;
					auto lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
					auto lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
					auto angle = acos(dot / sqrt(lenSq1 * lenSq2));
					if (angle < 0) angle += 2 * M_PI;
					auto cotangent1 = cos(angle) / sin(angle);
					x1 = mesh->verts[ev1]->coords[0] - mesh->verts[secondTriangleOtherVertexIdx]->coords[0];
					y1 = mesh->verts[ev1]->coords[1] - mesh->verts[secondTriangleOtherVertexIdx]->coords[1];
					z1 = mesh->verts[ev1]->coords[2] - mesh->verts[secondTriangleOtherVertexIdx]->coords[2];
					x2 = mesh->verts[ev2]->coords[0] - mesh->verts[secondTriangleOtherVertexIdx]->coords[0];
					y2 = mesh->verts[ev2]->coords[1] - mesh->verts[secondTriangleOtherVertexIdx]->coords[1];
					z2 = mesh->verts[ev2]->coords[2] - mesh->verts[secondTriangleOtherVertexIdx]->coords[2];
					dot = x1 * x2 + y1 * y2 + z1 * z2;
					lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
					lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
					angle = acos(dot / sqrt(lenSq1 * lenSq2));
					if (angle < 0) angle += 2 * M_PI;
					auto cotangent2 = cos(angle) / sin(angle);
					if (vertexId == ev1) {
						w(ev1, ev2) = (cotangent1 + cotangent2) * 0.5f;
					}
					else {
						w(ev2, ev1) = (cotangent1 + cotangent2) * 0.5f;
					}
				}
				else {
					std::cout << "Vertex " << vertexId << " : non-boundary hole vertex belongs to " << triangleIndices.size() << " triangles" << std::endl;
					if (useSymmetricTriangles) {
						std::cout << "But I'll override with symmetric triangle" << std::endl;
						auto firstTriangleIdx = mesh->tris[triangleIndices[0]];
						int firstTriangleOtherVertexIdx = -1;
						if (firstTriangleIdx->v1i == ev1) {
							if (firstTriangleIdx->v2i == ev2) {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
							}
							else {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
							}
						}
						else if (firstTriangleIdx->v2i == ev1) {
							if (firstTriangleIdx->v1i == ev2) {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
							}
							else {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
							}
						}
						else {
							if (firstTriangleIdx->v2i == ev2) {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
							}
							else {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
							}
						}
						auto x1 = mesh->verts[ev1]->coords[0] - mesh->verts[firstTriangleOtherVertexIdx]->coords[0];
						auto y1 = mesh->verts[ev1]->coords[1] - mesh->verts[firstTriangleOtherVertexIdx]->coords[1];
						auto z1 = mesh->verts[ev1]->coords[2] - mesh->verts[firstTriangleOtherVertexIdx]->coords[2];
						auto x2 = mesh->verts[ev2]->coords[0] - mesh->verts[firstTriangleOtherVertexIdx]->coords[0];
						auto y2 = mesh->verts[ev2]->coords[1] - mesh->verts[firstTriangleOtherVertexIdx]->coords[1];
						auto z2 = mesh->verts[ev2]->coords[2] - mesh->verts[firstTriangleOtherVertexIdx]->coords[2];
						auto dot = x1 * x2 + y1 * y2 + z1 * z2;
						auto lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
						auto lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
						auto angle = acos(dot / sqrt(lenSq1 * lenSq2));
						if (angle < 0) {
							angle += 2 * M_PI;
						}
						auto cotangent1 = cos(angle) / sin(angle);
						auto cotangent2 = cos(angle) / sin(angle);
						if (vertexId == ev1) {
							w(ev1, ev2) = (cotangent1) * 0.5f;
						}
						else {
							w(ev2, ev1) = (cotangent1) * 0.5f;
						}
					}
				}
			}
		}
		for (int i = 0; i < numVertices; ++i)
		{
			if (isVertexBoundary[i]) continue;
			double sum = 0.0;
			for (int j = 0; j < numVertices; ++j)
			{
				sum += w(i, j);
			}
			w(i, i) = -sum;
		}
#pragma endregion
	}
	else if (mode == 2) {
#pragma region meanvalueweights
		std::vector<int> nonBoundaryIndices;
		for (int i = 0; i < numVertices; ++i) {
			if (!isVertexBoundary[i]) {
				nonBoundaryIndices.push_back(i);
			}
			for (int j = 0; j < numVertices; ++j) {
				if (isVertexBoundary[i]) {
					if (i == j) {
						w(i, j) = 1;
					}
					else {
						w(i, j) = 0;
					}
				}
				else {
					w(i, j) = 0;
				}
			}
		}
		for (int i = 0; i < nonBoundaryIndices.size(); ++i) {
			const auto vertexId = nonBoundaryIndices[i];
			const auto neighbouringEdges = mesh->verts[vertexId]->edgeList;
			const auto neighbouringTris = mesh->verts[vertexId]->triList;
			if (neighbouringTris.size() < 2) {
				std::cout << "error on vertex " << vertexId << " : non-boundary vertex belongs to less than 2 triangles" << std::endl;
				continue;
			}
			for (int j = 0; j < neighbouringEdges.size(); ++j) {
				std::vector<int> triangleIndices;
				const int ev1 = mesh->edges[neighbouringEdges[j]]->v1i;
				const int ev2 = mesh->edges[neighbouringEdges[j]]->v2i;
				for (int k = 0; k < neighbouringTris.size(); ++k) {
					const int tv1 = mesh->tris[neighbouringTris[k]]->v1i;
					const int tv2 = mesh->tris[neighbouringTris[k]]->v2i;
					const int tv3 = mesh->tris[neighbouringTris[k]]->v3i;
					if (ev1 == tv1) {
						if (ev2 == tv2 || ev2 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
							continue;
						}
					}
					else if (ev1 == tv2) {
						if (ev2 == tv1 || ev2 == tv2) {
							triangleIndices.push_back(neighbouringTris[k]);
							continue;
						}
					}
					else if (ev1 == tv3) {
						if (ev2 == tv1 || ev2 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
							continue;
						}
					}
					if (ev2 == tv1) {
						if (ev1 == tv2 || ev1 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
						}
					}
					else if (ev2 == tv2) {
						if (ev1 == tv1 || ev1 == tv3) {
							triangleIndices.push_back(neighbouringTris[k]);
						}
					}
					else if (ev2 == tv3) {
						if (ev1 == tv1 || ev1 == tv2) {
							triangleIndices.push_back(neighbouringTris[k]);
						}
					}
				}
				if (triangleIndices.size() == 2) {
					auto firstTriangleIdx = mesh->tris[triangleIndices[0]];
					auto secondTriangleIdx = mesh->tris[triangleIndices[1]];
					int firstTriangleOtherVertexIdx = -1;
					int secondTriangleOtherVertexIdx = -1;
					if (firstTriangleIdx->v1i == ev1) {
						if (firstTriangleIdx->v2i == ev2) {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
						}
						else {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
						}
					}
					else if (firstTriangleIdx->v2i == ev1) {
						if (firstTriangleIdx->v1i == ev2) {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
						}
						else {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
						}
					}
					else {
						if (firstTriangleIdx->v2i == ev2) {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
						}
						else {
							firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
						}
					}
					if (secondTriangleIdx->v1i == ev1) {
						if (secondTriangleIdx->v2i == ev2) {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v3i;
						}
						else {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v2i;
						}
					}
					else if (secondTriangleIdx->v2i == ev1) {
						if (secondTriangleIdx->v1i == ev2) {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v3i;
						}
						else {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v1i;
						}
					}
					else {
						if (secondTriangleIdx->v2i == ev2) {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v1i;
						}
						else {
							secondTriangleOtherVertexIdx = secondTriangleIdx->v2i;
						}
					}
					float x1, y1, z1, x2, y2, z2, angle, tangent1, tangent2;
					if (vertexId == ev1)
					{
						x1 = mesh->verts[secondTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev1]->coords[0];
						y1 = mesh->verts[secondTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev1]->coords[1];
						z1 = mesh->verts[secondTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev1]->coords[2];
						x2 = mesh->verts[ev2]->coords[0] - mesh->verts[ev1]->coords[0];
						y2 = mesh->verts[ev2]->coords[1] - mesh->verts[ev1]->coords[1];
						z2 = mesh->verts[ev2]->coords[2] - mesh->verts[ev1]->coords[2];
						auto dot = x1 * x2 + y1 * y2 + z1 * z2;
						auto lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
						auto lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
						angle = acos(dot / sqrt(lenSq1 * lenSq2));
						if (angle < 0) angle += 2 * M_PI;
						tangent1 = sin(angle * 0.5f) / cos(angle * 0.5f);
						x1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev1]->coords[0];
						y1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev1]->coords[1];
						z1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev1]->coords[2];
						x2 = mesh->verts[ev2]->coords[0] - mesh->verts[ev1]->coords[0];
						y2 = mesh->verts[ev2]->coords[1] - mesh->verts[ev1]->coords[1];
						z2 = mesh->verts[ev2]->coords[2] - mesh->verts[ev1]->coords[2];
						dot = x1 * x2 + y1 * y2 + z1 * z2;
						lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
						lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
						angle = acos(dot / sqrt(lenSq1 * lenSq2));
						if (angle < 0) angle += 2 * M_PI;
						tangent2 = sin(angle * 0.5f) / cos(angle * 0.5f);
					}
					else
					{
						x1 = mesh->verts[secondTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev2]->coords[0];
						y1 = mesh->verts[secondTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev2]->coords[1];
						z1 = mesh->verts[secondTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev2]->coords[2];
						x2 = mesh->verts[ev1]->coords[0] - mesh->verts[ev2]->coords[0];
						y2 = mesh->verts[ev1]->coords[1] - mesh->verts[ev2]->coords[1];
						z2 = mesh->verts[ev1]->coords[2] - mesh->verts[ev2]->coords[2];
						auto dot = x1 * x2 + y1 * y2 + z1 * z2;
						auto lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
						auto lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
						angle = acos(dot / sqrt(lenSq1 * lenSq2));
						if (angle < 0) angle += 2 * M_PI;
						tangent1 = sin(angle * 0.5f) / cos(angle * 0.5f);
						x1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev2]->coords[0];
						y1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev2]->coords[1];
						z1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev2]->coords[2];
						x2 = mesh->verts[ev1]->coords[0] - mesh->verts[ev2]->coords[0];
						y2 = mesh->verts[ev1]->coords[1] - mesh->verts[ev2]->coords[1];
						z2 = mesh->verts[ev1]->coords[2] - mesh->verts[ev2]->coords[2];
						dot = x1 * x2 + y1 * y2 + z1 * z2;
						lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
						lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
						angle = acos(dot / sqrt(lenSq1 * lenSq2));
						if (angle < 0) angle += 2 * M_PI;
						tangent2 = sin(angle * 0.5f) / cos(angle * 0.5f);
					}
					auto xsq = pow(mesh->verts[ev1]->coords[0] - mesh->verts[ev2]->coords[0], 2.0);
					auto ysq = pow(mesh->verts[ev1]->coords[1] - mesh->verts[ev2]->coords[1], 2.0);
					auto zsq = pow(mesh->verts[ev1]->coords[2] - mesh->verts[ev2]->coords[2], 2.0);
					auto len = sqrt(xsq + ysq + zsq);
					if (vertexId == ev1) {
						w(ev1, ev2) = (tangent1 + tangent2) / (2.0f * len);
					}
					else {
						w(ev2, ev1) = (tangent1 + tangent2) / (2.0f * len);
					}
				}
				else {
					std::cout << "Vertex " << vertexId << " : non-boundary vertex belongs to " << triangleIndices.size() << " triangles" << std::endl;
					if (useSymmetricTriangles)
					{
						std::cout << "But I'll override with a symmetric triangle" << std::endl;
						auto firstTriangleIdx = mesh->tris[triangleIndices[0]];
						int firstTriangleOtherVertexIdx = -1;
						if (firstTriangleIdx->v1i == ev1) {
							if (firstTriangleIdx->v2i == ev2) {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
							}
							else {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
							}
						}
						else if (firstTriangleIdx->v2i == ev1) {
							if (firstTriangleIdx->v1i == ev2) {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v3i;
							}
							else {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
							}
						}
						else {
							if (firstTriangleIdx->v2i == ev2) {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v1i;
							}
							else {
								firstTriangleOtherVertexIdx = firstTriangleIdx->v2i;
							}
						}
						float x1, y1, z1, x2, y2, z2, angle, tangent1, tangent2;
						if (vertexId == ev1)
						{
							x1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev1]->coords[0];
							y1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev1]->coords[1];
							z1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev1]->coords[2];
							x2 = mesh->verts[ev2]->coords[0] - mesh->verts[ev1]->coords[0];
							y2 = mesh->verts[ev2]->coords[1] - mesh->verts[ev1]->coords[1];
							z2 = mesh->verts[ev2]->coords[2] - mesh->verts[ev1]->coords[2];
							auto dot = x1 * x2 + y1 * y2 + z1 * z2;
							auto lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
							auto lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
							angle = acos(dot / sqrt(lenSq1 * lenSq2));
							if (angle < 0) angle += 2 * M_PI;
							tangent1 = sin(angle * 0.5f) / cos(angle * 0.5f);
							x1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev1]->coords[0];
							y1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev1]->coords[1];
							z1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev1]->coords[2];
							x2 = mesh->verts[ev2]->coords[0] - mesh->verts[ev1]->coords[0];
							y2 = mesh->verts[ev2]->coords[1] - mesh->verts[ev1]->coords[1];
							z2 = mesh->verts[ev2]->coords[2] - mesh->verts[ev1]->coords[2];
							dot = x1 * x2 + y1 * y2 + z1 * z2;
							lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
							lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
							angle = acos(dot / sqrt(lenSq1 * lenSq2));
							if (angle < 0) angle += 2 * M_PI;
							tangent2 = sin(angle * 0.5f) / cos(angle * 0.5f);
						}
						else
						{
							x1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev2]->coords[0];
							y1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev2]->coords[1];
							z1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev2]->coords[2];
							x2 = mesh->verts[ev1]->coords[0] - mesh->verts[ev2]->coords[0];
							y2 = mesh->verts[ev1]->coords[1] - mesh->verts[ev2]->coords[1];
							z2 = mesh->verts[ev1]->coords[2] - mesh->verts[ev2]->coords[2];
							auto dot = x1 * x2 + y1 * y2 + z1 * z2;
							auto lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
							auto lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
							angle = acos(dot / sqrt(lenSq1 * lenSq2));
							if (angle < 0) angle += 2 * M_PI;
							tangent1 = sin(angle * 0.5f) / cos(angle * 0.5f);
							x1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[0] - mesh->verts[ev2]->coords[0];
							y1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[1] - mesh->verts[ev2]->coords[1];
							z1 = mesh->verts[firstTriangleOtherVertexIdx]->coords[2] - mesh->verts[ev2]->coords[2];
							x2 = mesh->verts[ev1]->coords[0] - mesh->verts[ev2]->coords[0];
							y2 = mesh->verts[ev1]->coords[1] - mesh->verts[ev2]->coords[1];
							z2 = mesh->verts[ev1]->coords[2] - mesh->verts[ev2]->coords[2];
							dot = x1 * x2 + y1 * y2 + z1 * z2;
							lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
							lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
							angle = acos(dot / sqrt(lenSq1 * lenSq2));
							if (angle < 0) angle += 2 * M_PI;
							tangent2 = sin(angle * 0.5f) / cos(angle * 0.5f);
						}
						auto xsq = pow(mesh->verts[ev1]->coords[0] - mesh->verts[ev2]->coords[0], 2.0);
						auto ysq = pow(mesh->verts[ev1]->coords[1] - mesh->verts[ev2]->coords[1], 2.0);
						auto zsq = pow(mesh->verts[ev1]->coords[2] - mesh->verts[ev2]->coords[2], 2.0);
						auto len = sqrt(xsq + ysq + zsq);
						if (vertexId == ev1) {
							w(ev1, ev2) = (tangent1 + tangent2) / (2.0f * len);
						}
						else {
							w(ev2, ev1) = (tangent1 + tangent2) / (2.0f * len);
						}
					}
				}
			}
		}
		for (int i = 0; i < numVertices; ++i)
		{
			if (isVertexBoundary[i]) continue;
			double sum = 0.0;
			for (int j = 0; j < numVertices; ++j)
			{
				sum += w(i, j);
			}
			w(i, i) = -sum;
		}
#pragma endregion
	}
	t1 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << "Creating W matrix: " << duration << " seconds" << endl;
	for (int i = 0; i < numVertices; ++i) {
		xx(i, 0) = 0; xy(i, 0) = 0; bx(i, 0) = 0; by(i, 0) = 0; // init
	}

	t0 = chrono::high_resolution_clock::now();
	std::vector<std::pair<float, float>> diskPoints;
	auto stepSize = M_PI * 2.0 / (double)(boundaryIndices.size());
	double currentPointAngle = 0;
	//while (currentPointAngle <= M_PI * 2.0)
	while (diskPoints.size() != boundaryIndices.size())
	{
		diskPoints.push_back(std::make_pair(std::cos(currentPointAngle), std::sin(currentPointAngle)));
		currentPointAngle += stepSize;
	}
	//int currentDiskPoint = diskPoints.size() - 1;
	int currentDiskPoint = 0;
	int selectedIndex = boundaryIndices[0];
	std::vector<bool> isVertexBoundaryBackup = isVertexBoundary;
	bx(selectedIndex, 0) = diskPoints[0].first;
	by(selectedIndex, 0) = diskPoints[0].second;
	isVertexBoundary[selectedIndex] = false;
	currentDiskPoint++;
	int minBoundaryIdx = -1;
	while(currentDiskPoint < boundaryIndices.size() )
	{
		float minDist = FLT_MAX;
		minBoundaryIdx = -1;
		auto neighbours = mesh->verts[selectedIndex]->vertList;
		for(int i = 0; i < neighbours.size(); ++i)
		{
			auto neighbourIdx = neighbours[i];
			if(isVertexBoundary[neighbourIdx])
			{
				if(distances[selectedIndex][neighbourIdx] < minDist)
				{
					minDist = distances[selectedIndex][neighbourIdx];
					minBoundaryIdx = neighbourIdx;
					isVertexBoundary[neighbourIdx] = false;
				}
				
			}
		}
		if (minBoundaryIdx != -1) {
			bx(minBoundaryIdx, 0) = diskPoints[currentDiskPoint].first;
			by(minBoundaryIdx, 0) = diskPoints[currentDiskPoint].second;
			selectedIndex = minBoundaryIdx;
		}
		else
		{
			bx(selectedIndex, 0) = diskPoints[currentDiskPoint].first;
			by(selectedIndex, 0) = diskPoints[currentDiskPoint].second;
		}
		currentDiskPoint++;
	}
	currentDiskPoint = 0;
	t1 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << "Mapping on disk: " << duration << " seconds" << endl;

	// checking for degenerate cases for debugging, namely bx and by values
	for(int i = 0; i < numVertices; ++i)
	{
		if(isVertexBoundaryBackup[i] && (bx(i, 0) == 0.0 && by(i, 0) == 0.0))
		{
			std::cout << "WARNING BOUNDARY VERTEX " << i << " HAS 0 VAL" << std::endl;
		}
		else if(isVertexBoundaryBackup[i] == false && (bx(i, 0) != 0.0 || by(i, 0) != 0))
		{
			std::cout << "WARNING NON-BOUNDARY VERTEX " << i << " HAS NONZERO VAL" << std::endl;
		}
	}

	//std::ofstream file("matrices.txt");

	//file << "bx" << std::endl;
	//file << bx << std::endl;
	//file << "by" << std::endl;
	//file << by << std::endl;
	//file << "w" << std::endl;
	//file << w << std::endl;;

	currentDiskPoint = 0;
	auto winverse = w.inverse();
	t0 = chrono::high_resolution_clock::now();
	//	xx = w.bdcSvd(ComputeThinU | ComputeThinV).solve(bx);
		//xx = w.colPivHouseholderQr().solve(bx);
	xx = winverse * bx;
		//xx = w.ldlt().solve(bx);
	t1 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << "Solving xx: " << duration << " seconds" << endl;
	t0 = chrono::high_resolution_clock::now();
	//xy = w.colPivHouseholderQr().solve(by);
	xy = winverse * by;
	//xy = w.bdcSvd(ComputeThinU | ComputeThinV).solve(by);
	//xy = w.ldlt().solve(by);
	t1 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << "Solving xy: " << duration << " seconds" << endl;
	//file << "x" << std::endl;
	//file << xx << std::endl;
	//file << "y" << std::endl;
	//file << xy << std::endl;
	//file.close();
	/*
	FILE * pFile;
	pFile = fopen("dijkstra_out.txt", "w");

	for(int i = 0; i < numVertices; ++i)
	{
		for(int j = 0; j < numVertices; j++)
		{
			fprintf(pFile, "%g ",distances[i][j]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	*/
	//int p1 = 0, p2 = 200;
	/*
	int p1 = -1, p2 = -1;
	while (p1 < 0 || p1 >= numVertices) {
		cout << "Enter first vertex index for query, [0," << numVertices - 1 << "] :";
		cin >> p1;
		if (p1 < 0 || p1 >= numVertices) {
			cout << "Invalid index, try again: ";
		}
	}
	while (p2 < 0 || p2 >= numVertices) {
		cout << "Enter second vertex index for query, [0," << numVertices - 1 << "] :";
		cin >> p2;
		if (p2 < 0 || p2 >= numVertices) {
			cout << "Invalid index, try again: ";
		}
	}
	
	std::vector<int> shortestPathVertices;
	shortestPathVertices.push_back(p2);	// add dest node
	int p2Parent = parents[p1][p2];
	while(p2Parent != p1){
		shortestPathVertices.push_back(p2Parent);
		p2Parent = parents[p1][p2Parent];
	}
	shortestPathVertices.push_back(p1);	// add source node
	cout << "Shortest path vertices for queried indices " << p1 << ", " << p2 << endl;
	for(const auto &spv:shortestPathVertices){
		cout << spv << " ";
	}
	cout << endl;
	
#pragma region fps
	cout << "------------------" << endl << "FPS" << endl << "------------------" << endl;
	t0 = chrono::high_resolution_clock::now();
	std::vector<int> fpsVertices;
	int randomIndex = rand() % numVertices;
	fpsVertices.push_back(randomIndex);
	float maxDist = FLT_MIN;
	int maxDistIndex = -1;
	for (int i = 0; i < numVertices; ++i) {
		if (distances[randomIndex][i] > maxDist) {
			maxDist = distances[randomIndex][i];
			maxDistIndex = i;
		}
	}
	fpsVertices.push_back(maxDistIndex);
	cout << "Random seed and farthest vertex to it are: ";
	for (int i = 0; i < fpsVertices.size(); ++i) {
		cout << fpsVertices[i] << " ";
	}
	cout << endl;
	const int numSamples = 100;
	while (fpsVertices.size() < numSamples) {
		std::vector<pair<int, int>> associations;	// which vertex is associated with whom	(i, j)
		for (int i = 0; i < numVertices; ++i) {
			float minDist = FLT_MAX;
			int minIndex = -1;
			for (int j = 0; j < fpsVertices.size(); ++j) {
				if (i == fpsVertices[j]) {
					minIndex = -1;
					minDist = FLT_MAX;
					break;
				}
				if (distances[i][fpsVertices[j]] < minDist) {
					minDist = distances[i][fpsVertices[j]];
					minIndex = fpsVertices[j];
				}
			}
			if (minIndex != -1) {
				associations.push_back(make_pair(i, minIndex));
			}		// pair.first is in the elements to be sampled, pair.second are already sampled
		}
		// get the argmax out of all
		float maxGeoDist = FLT_MIN;
		int maxIndex = -1;
		for (int i = 0; i < associations.size(); ++i) {
			if (distances[associations[i].first][associations[i].second] > maxGeoDist) {
				maxIndex = associations[i].first;
				maxGeoDist = distances[associations[i].first][associations[i].second];
			}
		}
		fpsVertices.push_back(maxIndex);
	}
	mesh->samples = fpsVertices;
	t1 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << endl << "Sampling points: " << duration << " seconds" << endl;
#pragma endregion

#pragma region geodesic isocurves
	cout << "------------------" << endl << "Geodesic Isocurves" << endl << "------------------" << endl;
	int seedIndex = -1;
	while (seedIndex < 0 || seedIndex >= numVertices) {
		cout << "Enter the seed index, [0," << numVertices - 1 << "] :";
		cin >> seedIndex;
		if (seedIndex < 0 || seedIndex >= numVertices) {
			cout << "Invalid seed index, try again: ";
		}
	}
	int k = -1;
	while (k < 0 || k >= numVertices) {
		cout << "Enter the number of bins, [0," << numVertices - 1 << "] :";
		cin >> k;
		if (k < 0 || k >= numVertices) {
			cout << "Invalid bin count, try again: ";
		}
	}
	t0 = chrono::high_resolution_clock::now();
	float maxVertexDist = FLT_MIN;
	for(int i = 0; i < numVertices; ++i){
		if(distances[seedIndex][i] > maxVertexDist){
			maxVertexDist = distances[seedIndex][i];
		}
	}
	float d = maxVertexDist / (float)k;
	int numTris = mesh->tris.size();
	auto triangles = mesh->tris;
	std::vector<float> histogramBins(k);
	std::vector<std::vector<pair<std::vector<float>, std::vector<float>>>> isoCurveLines(k);
	for(int i = 1; i <= k; ++i){
		float radius = i * d;
		float isoCurveLength = 0.0f;
		for(int j = 0; j < numTris; ++j){
			auto distV1 = distances[seedIndex][triangles[j]->v1i];
			auto distV2 = distances[seedIndex][triangles[j]->v2i];
			auto distV3 = distances[seedIndex][triangles[j]->v3i];
			if ((distV1 < radius && distV2 < radius && distV3 < radius)
				|| (distV1 > radius && distV2 > radius && distV3 > radius)) continue;
			std::vector<int> lt, gt;
			distV1 <= radius ? lt.push_back(triangles[j]->v1i) : gt.push_back(triangles[j]->v1i);
			distV2 <= radius ? lt.push_back(triangles[j]->v2i) : gt.push_back(triangles[j]->v2i);
			distV3 <= radius ? lt.push_back(triangles[j]->v3i) : gt.push_back(triangles[j]->v3i);
			if (lt.empty() || gt.empty()) continue;
			float a1, a2, g0, g1, g2;
			std::vector<float> p1(3), p2(3);
			float * v0, * v1, * v2;
			if(lt.size() < gt.size()){
				g0 = distances[seedIndex][lt[0]];
				g1 = distances[seedIndex][gt[0]];
				g2 = distances[seedIndex][gt[1]];
				a1 = fabs(radius - g0) / fabs(g1 - g0);
				a2 = fabs(radius - g0) / fabs(g2 - g0);
				v0 = mesh->verts[lt[0]]->coords;
				v1 = mesh->verts[gt[0]]->coords;
				v2 = mesh->verts[gt[1]]->coords;
				p1[0] = (1.0f - a1) * v0[0] + a1 * v1[0];
				p1[1] = (1.0f - a1) * v0[1] + a1 * v1[1];
				p1[2] = (1.0f - a1) * v0[2] + a1 * v1[2];
				p2[0] = (1.0f - a2) * v0[0] + a2 * v2[0];
				p2[1] = (1.0f - a2) * v0[1] + a2 * v2[1];
				p2[2] = (1.0f - a2) * v0[2] + a2 * v2[2];
			}
			else if(lt.size() > gt.size()){
				g0 = distances[seedIndex][lt[0]];
				g1 = distances[seedIndex][lt[1]];
				g2 = distances[seedIndex][gt[0]];
				a1 = fabs(radius - g0) / fabs(g2 - g0);
				a2 = fabs(radius - g1) / fabs(g2 - g1);
				v0 = mesh->verts[lt[0]]->coords;
				v1 = mesh->verts[lt[1]]->coords;
				v2 = mesh->verts[gt[0]]->coords;
				p1[0] = (1.0f - a1) * v0[0] + a1 * v2[0];
				p1[1] = (1.0f - a1) * v0[1] + a1 * v2[1];
				p1[2] = (1.0f - a1) * v0[2] + a1 * v2[2];
				p2[0] = (1.0f - a2) * v1[0] + a2 * v2[0];
				p2[1] = (1.0f - a2) * v1[1] + a2 * v2[1];
				p2[2] = (1.0f - a2) * v1[2] + a2 * v2[2];
			}
			float p1p2dist = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
								(p1[1] - p2[1]) * (p1[1] - p2[1]) + 
								(p1[2] - p2[2]) * (p1[2] - p2[2]));
			isoCurveLength += p1p2dist;
			isoCurveLines[i - 1].push_back(make_pair(p1, p2));
		}
		histogramBins[i - 1] = isoCurveLength;
	}
	t1 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::duration<float>>(t1 - t0).count();
	std::cout << endl << "Geodesic isocurve: " << duration << " seconds " << endl;
#pragma endregion
	float globalMaxDist = FLT_MIN;
	for (int i = 0; i < distances.size(); ++i) {
		for (int j = 0; j < distances[i].size(); ++j) {
			if (distances[i][j] > globalMaxDist) {
				globalMaxDist = distances[i][j];
			}
		}
	}
	*/
	
	//root->addChild( painter->getShapeSep(mesh) );
	//root->addChild(painter->getSdfShapeSep(mesh, nsdf));
	root->addChild(painter->getSdfSegmentedShapeSep(mesh, nsdf, nsdfSegments));
	//root->addChild(painter->getRayCastRaysShapeSep(mesh, centroids, directions));
	int visualization = 1;
	while (visualization <= 0 || visualization > 4) {
		cout << endl << "Select the visualization: 1 -> Dijkstra, 2 -> Geodesic Isocurves, 3 -> Farthest Point Sampling, 4 -> Boundary Vertices :";
		cin >> visualization;
		if (visualization <= 0 || visualization > 3) {
			cout << "Invalid visualization query, try again: " << endl;
		}
	}
	if(visualization == 1){
		//root->addChild(painter->getShortestPathSep(mesh, shortestPathVertices));	// visualization for shortest path vertices
	}
	else if(visualization == 2){
		//root->addChild(painter->getGeodesicIsoCurveSep(mesh, isoCurveLines, histogramBins, seedIndex));
	}
	else if (visualization == 3) {
		root->addChild(painter->getSpheresSep(mesh, 0, 0, 1.0f)); // visualization for sampled points
	}
	else if(visualization == 4)
	{
		//mesh->samples = boundaryIndices;
		//root->addChild(painter->getSpheresSep(mesh, 0, 0, 1.0f));
		//root->addChild(painter->getParametrizedMeshSep(mesh, xx, xy));

	}
	
	
	//viewer->setSize(SbVec2s(800, 600));
	viewer->setSize(SbVec2s(1280, 760));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
