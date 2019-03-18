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


#include "Mesh.h"
#include "Painter.h"
#include <string>
#include <cstring>
#include <queue>
#include <chrono>

int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);

	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator * root = new SoSeparator;
	root->ref();

	//stuff to be drawn on screen must be added to the root
//	SoCone * cone = new SoCone;
//	root->addChild(cone);

	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();

	char* x = (char*)malloc(strlen("0.off") + 1); 
	strcpy(x, "0.off");
	mesh->loadOff(x);
	//mesh->createCube(20.0f);


	cout << "my (verts[4]) 1-ring neighborhood is: \n";
	for (int nv = 0; nv < mesh->verts[4]->vertList.size(); nv++)
		cout << mesh->verts[4]->vertList[nv] << " neighbb\n";

	cout << mesh->verts[4]->coords[0] << ", " << mesh->verts[4]->coords[1] << ", " << mesh->verts[4]->coords[2] << endl;
	cout << "my (verts[4]) 1-ring neighborhood is: \n";
	for (int ne = 0; ne < mesh->verts[4]->edgeList.size(); ne++)
		if (mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v1i == 4)
			cout << mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v2i << " nnnnnnnnn\n";
		else
			cout << mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v1i << " nnnnnnnnn\n";

//		cout << mesh->verts[4]->vertList[nv] << " neighbb\n";

	const int numVertices = mesh->verts.size();
#pragma region dijkstraQuery
	/*
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
	*/
#pragma endregion 

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
	std::cout << "Min heap: " << duration << endl;
#pragma endregion 
	
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
	std::cout << "Array: " << duration << endl;
	*/
#pragma endregion 

	FILE * pFile;
	pFile = fopen("out.txt", "w");

	for(int i = 0; i < numVertices; ++i)
	{
		for(int j = 0; j < numVertices; j++)
		{
			fprintf(pFile, "%g ",distances[i][j]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);

	int p1 = -1, p2 = -1;
	while (p1 < 0 || p1 >= numVertices) {
		cout << "Enter first vertex index, [0," << numVertices - 1 << "] :";
		cin >> p1;
		if (p1 < 0 || p1 >= numVertices) {
			cout << "Invalid index, try again: ";
		}
	}
	while (p2 < 0 || p2 >= numVertices) {
		cout << "Enter second vertex index, [0," << numVertices - 1 << "] :";
		cin >> p2;
		if (p2 < 0 || p2 >= numVertices) {
			cout << "Invalid index, try again: ";
		}
	}

	root->addChild( painter->getShapeSep(mesh) );

	SoSeparator *sphereRoot = new SoSeparator;
	SoTransform *sphereTransform = new SoTransform;
	sphereTransform->translation.setValue(17., 17., 0.);
	sphereTransform->scaleFactor.setValue(8., 8., 8.);
	sphereRoot->addChild(sphereTransform);

	auto sphereMaterial = new SoMaterial;
	sphereMaterial->diffuseColor.setValue(.8, .8, .8);
	sphereRoot->addChild(sphereMaterial);
	sphereRoot->addChild(new SoSphere);
	root->addChild(sphereRoot);

	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
