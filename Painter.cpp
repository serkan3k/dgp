#define HAVE_SINGLEPRECISION_MATH
#include "Painter.h"

SoSeparator* Painter::getShapeSep(Mesh* mesh)
{
	SoSeparator* res = new SoSeparator();

	//transformation
	//not needed

	//color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(1, 1, 1); //paint all vertices with this color
	//mat->transparency = 0.5f : 0.0f; //0 makes it completely opaque, the default

	bool youWantToPaintEachVertexDifferently = false;
	if (youWantToPaintEachVertexDifferently)
		for (int i = 0; i < (int) mesh->verts.size(); i++) //i = 0 obj->color above overwritten here
			mat->diffuseColor.set1Value(i, mesh->verts[i]->color); //vert color according to its x-y-z coord (for mesh1) and to the transferred color (for mesh2)


	res->addChild(mat);

	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints); //Gouraud shading

	if (youWantToPaintEachVertexDifferently)
	{
		SoMaterialBinding* materialBinding = new SoMaterialBinding; //for 2+ diffuse color usage on the same mesh
		materialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
		res->addChild(materialBinding);
	}

	//shape
	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->verts.size(); c++)
		coords->point.set1Value(c, mesh->verts[c]->coords[0], mesh->verts[c]->coords[1], mesh->verts[c]->coords[2]);
	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	for (int c = 0; c < mesh->tris.size(); c++)
	{
		faceSet->coordIndex.set1Value(c*4, mesh->tris[c]->v1i);
		faceSet->coordIndex.set1Value(c*4 + 1, mesh->tris[c]->v2i);
		faceSet->coordIndex.set1Value(c*4 + 2, mesh->tris[c]->v3i);
		faceSet->coordIndex.set1Value(c*4 + 3, -1);

		if (youWantToPaintEachVertexDifferently)
		{
			int t = c;
			int nt = mesh->tris.size();
			faceSet->materialIndex.set1Value(0 + 4*nt, mesh->tris[c]->v1i);
			faceSet->materialIndex.set1Value(1 + 4*nt, mesh->tris[c]->v2i);
			faceSet->materialIndex.set1Value(2 + 4*nt, mesh->tris[c]->v3i);
		}
	}
	res->addChild(coords);
	res->addChild(faceSet);

	return res;
}

SoSeparator* Painter::getSpheresSep(Mesh* mesh, float deltaX, float deltaY, float scale)
{
	//returns a set of spheres to highlight each mesh.samples[i]

	SoSeparator* spheresSep = new SoSeparator();

	float radius = 0.25f;
	//float radius = 2.0f;

	for (int i = 0; i < (int)mesh->samples.size(); i++)
	{
		//1 sphere for this sample
		SoSeparator* sphere1Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		tra->translation.setValue(scale*mesh->verts[mesh->samples[i]]->coords[0] + deltaX, scale*mesh->verts[mesh->samples[i]]->coords[1] + deltaY, scale*mesh->verts[mesh->samples[i]]->coords[2]);
		sphere1Sep->addChild(tra);

		//material
		SoMaterial* ma = new SoMaterial;
		if (i == 0)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.0f, 0.7f));
		else if (i == 1)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.0f, 0.0f));
		else if (i == 2)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.7f, 0.0f));
		else if (i == 3)
			ma->diffuseColor.setValue(SbColor(0.7f, 0.0f, 0.7f));
		else if (i == 4)
			ma->diffuseColor.setValue(SbColor(0.7f, 0.7f, 0.0f));
		else
			ma->diffuseColor.setValue(SbColor(0.7f, 0.0f, 0.0f));

		sphere1Sep->addChild(ma);

		//shape
		SoSphere* sph1 = new SoSphere();
		sph1->radius = radius;
		sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above

		spheresSep->addChild(sphere1Sep);
	}

	return spheresSep;
}

SoSeparator * Painter::getShortestPathSep(Mesh * mesh, const vector<int> &shortestPathVertices)
{
	SoSeparator* thickEdgeSep = new SoSeparator;
	//material
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.set1Value(0, 1.0f, 0.0f, 0.0f);
	thickEdgeSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 3.0f;	thickEdgeSep->addChild(sty);

	//shape
	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;

	//assumes no edge in sedges is removed
	for(int i = 0; i < shortestPathVertices.size() - 1; ++i){
		SbVec3f end1 = mesh->verts[shortestPathVertices[i]]->coords;
		SbVec3f end2 = mesh->verts[shortestPathVertices[i + 1]]->coords;
		co->point.set1Value(2 * i, end1);
		co->point.set1Value(2 * i + 1, end2);
	}
	for(int i = 0; i < shortestPathVertices.size()-1; ++i){
		ils->coordIndex.set1Value(3 * i, 2 * i);
		ils->coordIndex.set1Value(3 * i + 1, 2 * i + 1);
		ils->coordIndex.set1Value(3 * i + 2, -1);
	}
	thickEdgeSep->addChild(co);	thickEdgeSep->addChild(ils);

	SoSeparator* sphere1SepStart = new SoSeparator;	
	SoTransform *sphereTransform = new SoTransform;
	sphereTransform->translation.setValue(
		mesh->verts[shortestPathVertices[0]]->coords[0], 
		mesh->verts[shortestPathVertices[0]]->coords[1], 
		mesh->verts[shortestPathVertices[0]]->coords[2]);
	sphereTransform->scaleFactor.setValue(1, 1, 1);
	sphere1SepStart->addChild(sphereTransform);
	auto sphereMaterial = new SoMaterial;
	sphereMaterial->diffuseColor.setValue(0.0f, 1.0f, 0.0f);
	sphere1SepStart->addChild(sphereMaterial);
	SoSphere* sph1 = new SoSphere();
	sph1->radius = 2.0f;
	sphere1SepStart->addChild(sph1);
	
	SoSeparator* sphere1SepEnd = new SoSeparator;
	SoTransform *sphereTransformEnd = new SoTransform;
	sphereTransformEnd->translation.setValue(
		mesh->verts[shortestPathVertices[shortestPathVertices.size()-1]]->coords[0], 
		mesh->verts[shortestPathVertices[shortestPathVertices.size() - 1]]->coords[1],
		mesh->verts[shortestPathVertices[shortestPathVertices.size() - 1]]->coords[2]);
	sphereTransformEnd->scaleFactor.setValue(1, 1, 1);
	sphere1SepEnd->addChild(sphereTransformEnd);
	auto sphereMaterialEnd = new SoMaterial;
	sphereMaterialEnd->diffuseColor.setValue(0.0f, 1.0f, 0.0f);
	sphere1SepEnd->addChild(sphereMaterialEnd);
	SoSphere* sph2 = new SoSphere();
	sph2->radius = 2.0f;
	sphere1SepEnd->addChild(sph2);

	thickEdgeSep->addChild(sphere1SepStart);
	thickEdgeSep->addChild(sphere1SepEnd);
	return thickEdgeSep;
}

SoSeparator* Painter::getGeodesicIsoCurveSep(Mesh* mesh,
	const std::vector<std::vector<pair<std::vector<float>, std::vector<float>>>>& isoCurves,
	const std::vector<float>& histogramBins, const int seedVertex)
{
	SoSeparator* res = new SoSeparator();
	for (int i = 0; i < isoCurves.size(); ++i) {
		if (isoCurves[i].size() == 0) continue;
		SoSeparator* thickEdgeSep = new SoSeparator;
		//material
		SoMaterial* ma = new SoMaterial;
		ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 0.0f);
		thickEdgeSep->addChild(ma);
		SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 3.0f;	thickEdgeSep->addChild(sty);

		//shape
		SoIndexedLineSet* ils = new SoIndexedLineSet;
		SoCoordinate3* co = new SoCoordinate3;

		//assumes no edge in sedges is removed
		for (int j = 0; j < isoCurves[i].size() - 1; ++j) {
			SbVec3f end1;// = isoCurves[i][j];//mesh->verts[shortestPathVertices[i]]->coords;
			end1.setValue(isoCurves[i][j].first[0], isoCurves[i][j].first[1], isoCurves[i][j].first[2]);
			SbVec3f end2; // = //mesh->verts[shortestPathVertices[i + 1]]->coords;
			end2.setValue(isoCurves[i][j].second[0], isoCurves[i][j].second[1], isoCurves[i][j].second[2]);
			co->point.set1Value(2 * j, end1);
			co->point.set1Value(2 * j + 1, end2);
		}
		for (int j = 0; j < isoCurves[i].size() - 1; ++j) {
			ils->coordIndex.set1Value(3 * j, 2 * j);
			ils->coordIndex.set1Value(3 * j + 1, 2 * j + 1);
			ils->coordIndex.set1Value(3 * j + 2, -1);
		}
		thickEdgeSep->addChild(co);	thickEdgeSep->addChild(ils);
		res->addChild(thickEdgeSep);
	}

	SoSeparator* sphere1SepStart = new SoSeparator;
	SoTransform *sphereTransform = new SoTransform;
	sphereTransform->translation.setValue(
		mesh->verts[seedVertex]->coords[0],
		mesh->verts[seedVertex]->coords[1],
		mesh->verts[seedVertex]->coords[2]);
	sphereTransform->scaleFactor.setValue(1, 1, 1);
	sphere1SepStart->addChild(sphereTransform);
	auto sphereMaterial = new SoMaterial;
	// COLOR THE SEED VERTEX
	float distanceTotal = 0;
	float maxDist = FLT_MIN;
	float minDist = FLT_MAX;
	for (unsigned int i = 0; i < histogramBins.size(); ++i) {
		if (histogramBins[i] == 0) continue;
		if (histogramBins[i] > maxDist) {
			maxDist = histogramBins[i];
		}
		if (histogramBins[i] < minDist) {
			minDist = histogramBins[i];
		}
		distanceTotal += histogramBins[i];
	}
	distanceTotal /= (float)histogramBins.size();
	float rColor = (distanceTotal -minDist) / (maxDist - minDist);	// red color is normalized value w.r.t. max histogram value
	if (distanceTotal < minDist) rColor = 0.0f;
	sphereMaterial->diffuseColor.setValue(rColor, 0.0f, 0.0f);
	// COLOR THE SEED VERTEX
	sphere1SepStart->addChild(sphereMaterial);
	SoSphere* sph1 = new SoSphere();
	sph1->radius = 2.0f;
	sphere1SepStart->addChild(sph1);
	res->addChild(sphere1SepStart);


	return res;
}


/* stuff below are from my old projects; should run fine and be useful in your development

if (drawThickEdges) //draw thick edges (may be useful in geodesic path drawing)
	{
		SoSeparator* thickEdgeSep = new SoSeparator;
		//material
		SoMaterial* ma = new SoMaterial;
		ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 1.0f);
		thickEdgeSep->addChild(ma);
		SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 5.0f;	thickEdgeSep->addChild(sty);

		//shape
		SoIndexedLineSet* ils = new SoIndexedLineSet;
		SoCoordinate3* co = new SoCoordinate3;

		//assumes no edge in sedges is removed
		for (unsigned int se = 0; se < mesh->sedges.size(); se++)
		{
			SbVec3f end1 = mesh->verts[ mesh->sedges[se]->v1i ]->coords + SbVec3f(deltaX, 0.0f, 0.0f),
					end2 = mesh->verts[ mesh->sedges[se]->v2i ]->coords + SbVec3f(deltaX, 0.0f, 0.0f);
			co->point.set1Value(2*se, end1);
			co->point.set1Value(2*se + 1, end2);
		}

		for (unsigned int ci = 0; ci < mesh->sedges.size(); ci++)
		{
			ils->coordIndex.set1Value(3*ci, 2*ci);	ils->coordIndex.set1Value(3*ci + 1, 2*ci + 1);
			ils->coordIndex.set1Value(3*ci + 2, -1); //end this edge with -1
		}
		thickEdgeSep->addChild(co);	thickEdgeSep->addChild(ils);
		obj->sep->addChild(thickEdgeSep);
	}
	
	
SoSeparator* Painter::get1PointSep(ScreenObject* obj, int pnt, int drawWhat, float deltaX, float deltaY, float scale)
{
	//renders only 1 pnt in blue, w/ drawWhat = 1 for spectral coords, = 2 for spatial coords, = 5 for coord written here

	Mesh* mesh = obj->getMesh();

	SoSeparator* pntSep = new SoSeparator;
	//material
	SoMaterial* mat = new SoMaterial;
	if (mesh->targetMesh)
		mat->diffuseColor.setValue(SbColor(1.0f, 0.0f, 0.0f)); //red
	else
		mat->diffuseColor.setValue(SbColor(0.0f, 1.0f, 0.0f)); //green
if (pnt == 594) mat->diffuseColor.setValue(SbColor(1.0f, 0.0f, 1.0f)); //magenta
//if (pnt == 6916) mat->diffuseColor.setValue(SbColor(0.0f, 1.0f, 1.0f));

	pntSep->addChild(mat);
	SoDrawStyle* style = new SoDrawStyle;
	style->pointSize = 17.0f;
	pntSep->addChild(style);
	
	//shape
	SoVertexProperty* vp = new SoVertexProperty;
	if (drawWhat == 2)
		vp->vertex.set1Value(0, scale*mesh->verts[pnt]->coords[0]+deltaX, scale*mesh->verts[pnt]->coords[1]+deltaY, scale*mesh->verts[pnt]->coords[2]);
	else if (drawWhat == 5)
		vp->vertex.set1Value(0, 0.016721f, -0.000984876f, 0.0f);
	else
		vp->vertex.set1Value(0, scale*mesh->verts[pnt]->spectralK[0]+deltaX, scale*mesh->verts[pnt]->spectralK[1]+deltaX, scale*mesh->verts[pnt]->spectralK[2]+deltaX);
	SoPointSet* pSet = new SoPointSet;
	pSet->numPoints = 1;
	pSet->vertexProperty = vp;
	pntSep->addChild(pSet);

//cout << pnt << " ------> " << mesh->verts[pnt]->matchIdx << endl;
	return pntSep;
}

SoSeparator* Painter::getPointsSep(Mesh* mesh, SbColor c)
{
	//renders grid points, i.e. voxel centers

	SoSeparator* pntsSep = new SoSeparator;

	//material
	SoMaterial* mat = new SoMaterial;
	mat->diffuseColor.set1Value(0, c); //SbColor(199.0f/255.0f, 166.0f/255.0f, 1.0f));
	pntsSep->addChild(mat);
	SoDrawStyle* style = new SoDrawStyle;
	style->pointSize = 7.0f;
	pntsSep->addChild(style);

	//shape
	SoVertexProperty* vp = new SoVertexProperty;	
	int nPnts = (int) mesh->verts.size(), nAdds = 0;
	for (int p = 0; p < nPnts; p++)
	{
//if (selection[p]->center[1] > 120) continue; //just draw allowed voxels see its volume/thickness better
		vp->vertex.set1Value(nAdds++, mesh->verts[p]->coords);
	}
	SoPointSet* pSet = new SoPointSet;
	pSet->numPoints = nAdds;
	pSet->vertexProperty = vp;
	pntsSep->addChild(pSet);

	return pntsSep;
}


*/
