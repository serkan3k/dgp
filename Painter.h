#pragma once
#define HAVE_SINGLEPRECISION_MATH
#define _CRT_SECURE_NO_WARNINGS
#define COIN_DLL
#define HAVE_INT8_T
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoSphere.h>

//#include <Inventor/nodes/SoCone.h>

#include "Mesh.h"


class Painter
{
public:
	SoSeparator* getShapeSep(Mesh* mesh);
	SoSeparator * getSpheresSep(Mesh * mesh, float deltaX, float deltaY, float scale);
};
