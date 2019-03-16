#define HAVE_SINGLEPRECISION_MATH
#define _CRT_SECURE_NO_WARNINGS
#define COIN_DLL
#define SOWIN_DLL
#define HAVE_INT8_T
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

#include "Mesh.h"
#include "Painter.h"
#include <string>
#include <cstring>
#include <queue>

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


	cout << "my (verts[4]) 1-ring neighborhood is: \n";
	for (int ne = 0; ne < mesh->verts[4]->edgeList.size(); ne++)
		if (mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v1i == 4)
			cout << mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v2i << " nnnnnnnnn\n";
		else
			cout << mesh->edges[   mesh->verts[4]->edgeList[ne]   ]->v1i << " nnnnnnnnn\n";

//		cout << mesh->verts[4]->vertList[nv] << " neighbb\n";

	const int numVertices = mesh->verts.size();
	std::priority_queue<float, std::vector<float>, std::greater<>> pq;
	std::vector<std::vector<int>> parent;
	parent.resize(numVertices);
	for(int i = 0; i < numVertices; ++i)
	{
		parent[i].resize(numVertices);
		for(int j = 0; j < numVertices; ++j)
		{
			parent[i][j] = -1;
		}
	}
	

	root->addChild( painter->getShapeSep(mesh) );


	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
