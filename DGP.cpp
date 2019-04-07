#define SOWIN_DLL
#define COIN_DLL
#define HAVE_INT8_T

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCone.h>

int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);

	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator * root = new SoSeparator;
	root->ref();

	//stuff to be drawn on screen must be added to the root
	SoCone * cone = new SoCone;
	root->addChild(cone);

	viewer->setSceneGraph(root);
	viewer->show();
//
	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
