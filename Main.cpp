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


#include "Mesh.h"
#include "Painter.h"
#include <string>
#include <cstring>
#include <queue>
#include <chrono>
#include <functional>
#include <set>
#include <math.h>

using Eigen::MatrixXd;

int main(int, char ** argv)
{
	MatrixXd m(2, 2);
	m(0, 0) = 3;
	m(1, 0) = 2.5;
	m(0, 1) = -1;
	m(1, 1) = m(1, 0) + m(0, 1);
	std::cout << m << std::endl;

	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);
	SoSeparator * root = new SoSeparator;
	root->ref();
	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();
	// load mesh
	char* x = (char*)malloc(strlen("facem.off") + 1); 
	strcpy(x, "facem.off");
	mesh->loadOff(x);

	const int numVertices = mesh->verts.size();
	std::vector<int> boundaryIndices;
	std::set<int> boundaryVertices;
	const int numEdges = mesh->edges.size();
	const int numTris = mesh->tris.size();
	const auto edges = mesh->edges;
	const auto tris = mesh->tris;
	const auto verts = mesh->verts;
	for(int i = 0; i < numEdges; ++i){
		int belongsTo = 0;
		const int ev1 = edges[i]->v1i;
		const int ev2 = edges[i]->v2i;
		for(int j = 0; j < numTris; ++j){
			const int tv1 = tris[j]->v1i;
			const int tv2 = tris[j]->v2i;
			const int tv3 = tris[j]->v3i;
			if(ev1 == tv1){
				if(ev2 == tv2 || ev2 == tv3){
					belongsTo++;
					continue;
				}
			}
			else if(ev1 == tv2){
				if(ev2 == tv1 || ev2 == tv2){
					belongsTo++;
					continue;
				}
			}
			else if(ev1 == tv3){
				if(ev2 == tv1 || ev2 == tv3){
					belongsTo++;
					continue;
				}
			}
			if(ev2 == tv1){
				if(ev1 == tv2 || ev1 == tv3){
					belongsTo++;
					continue;
				}
			}
			else if(ev2 == tv2){
				if(ev1 == tv1 || ev1 == tv3){
					belongsTo++;
					continue;
				}
			}
			else if(ev2 == tv3){
				if(ev1 == tv1 || ev1 == tv2){
					belongsTo++;
					continue;
				}
			}
		}
		if(belongsTo == 1){
			boundaryVertices.insert(ev1);
			boundaryVertices.insert(ev2);
		}
	}
	for(auto it = boundaryVertices.begin(); it != boundaryVertices.end(); ++it){
		boundaryIndices.push_back(*it);
	}
	std::vector<std::pair<float, float>> diskPoints;
	auto stepSize = M_PI * 2.0 / (double)boundaryIndices.size();
	double currentPointAngle = 0;
	while(currentPointAngle < M_PI * 2.0)
	{
		diskPoints.push_back(std::make_pair(std::sin(currentPointAngle), std::cos(currentPointAngle)));
		currentPointAngle += stepSize;
	}
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
	/*
#pragma region minHeap
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
	std::cout << "Min heap: " << duration << " seconds" << endl;
#pragma endregion 
	
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

	//int p1 = 0, p2 = 200;
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
	root->addChild( painter->getShapeSep(mesh) );
	int visualization = 4;
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
		mesh->samples = boundaryIndices;
		root->addChild(painter->getSpheresSep(mesh, 0, 0, 1.0f));
	}
	
	
	viewer->setSize(SbVec2s(800, 600));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
