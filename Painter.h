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
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoText2.h>

#include <Eigen/Dense>
//#include <Inventor/nodes/SoCone.h>

#include "Mesh.h"
#include "glm/vec3.hpp"
using namespace Eigen;

class Painter
{
public:
	SoSeparator* getShapeSep(Mesh* mesh);
	SoSeparator* getSdfShapeSep(Mesh* mesh, std::vector<float> &normalizedSdf);
	SoSeparator* getSdfSegmentedShapeSep(Mesh* mesh, std::vector<float> &normalizedSdf, std::vector<int> &nsdfSegment);
	SoSeparator* getRayCastRaysShapeSep(Mesh* mesh, std::vector<glm::vec3> rayOrigins, std::vector<glm::vec3> rayDirections);
	SoSeparator * getSpheresSep(Mesh * mesh, float deltaX, float deltaY, float scale);
	SoSeparator * getShortestPathSep(Mesh * mesh, const vector<int> &shortestPathVertices);
	SoSeparator * getGeodesicIsoCurveSep(Mesh * mesh, const std::vector<std::vector<pair<std::vector<float>, std::vector<float>>>> &isoCurves, const std::vector<float> &histogramBins, const int seedVertex);
	SoSeparator * getParametrizedMeshSep(Mesh * mesh, const MatrixXd &xx, const MatrixXd &xy);
};
