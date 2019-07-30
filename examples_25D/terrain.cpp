#include <Fade_2D.h>

using namespace std;
using namespace GEOM_FADE25D;
void readTerrain(std::vector<Point2>& vInPoints,const char* filename);
void isoContours(std::vector<Triangle2*>& vTriangles);
void getPolygon(Fade_2D* pDt,std::vector<Segment2>& vSegments);




int terrain_main()
{
	std::cout<<"\n* Fade2.5D Demo - triangulation, constraints and contours"<<std::endl;
	std::cout<<"---------------------------------------------------------"<<std::endl<<std::endl;

	// 1) Create a new 2.5D triangulation. 
	std::cout<<"\n* Triangulate 2.5D points"<<std::endl;
	Fade_2D* pDt=new Fade_2D();

	// 2) Read terrain points from an ASCII file (xyz-coordinates)
	std::vector<Point2> vInputPoints;
	readTerrain(vInputPoints,"../examples_25D/gaeta_small.xyz");


	// NEW since version 1.46: You may want to prune the input points to 
	// represent the model more efficiently
	cout<<"\nOriginal number of points: "<<vInputPoints.size()<<endl;
	EfficientModel em(vInputPoints);
	vInputPoints.clear();
	double maxError(1.0);
	em.extract(maxError,vInputPoints);
	cout<<"Efficient number of points: "<<vInputPoints.size()<<endl<<endl;


	// 3) Insert the points
	pDt->insert(vInputPoints);
	// [X] Done - Show what we have:

	// Write a wavefront *.obj file
	pDt->writeObj("gaeta.obj");

	// Write a visualization for (technically up-to-date) web browsers
	pDt->writeWebScene("gaeta_webScene");

	// Write a file for the Geomview viewer
	pDt->showGeomview("gaeta_geomview.list");

	// ** Test constraint insertion:
	std::cout<<"\n* Testing constraint insertion"<<std::endl;

	// 1) Create a polygon
	std::vector<Segment2> vSegments;
	getPolygon(pDt,vSegments);

	// 2) Insert it as constraint, strategy can be:
	//    +) CIS_CONSTRAINED_DELAUNAY: without subdivision
	//    +) CIS_CONFORMING_DELAUNAY_SEGMENT_LEVEL: subdivision at segment level
	//    +) CIS_CONFORMING_DELAUNAY: at triangulation level

	//ConstraintInsertionStrategy cis(CIS_CONSTRAINED_DELAUNAY);
	//ConstraintInsertionStrategy cis(CIS_CONFORMING_DELAUNAY_SEGMENT_LEVEL);
	ConstraintInsertionStrategy cis(CIS_CONFORMING_DELAUNAY);

	pDt->createConstraint(vSegments,cis);
	pDt->applyConstraintsAndZones();

	// 3) Write a visualization for (technically up-to-date) web browsers
	//    The constraint polygon adapts to the heights of the terrain. Thus
	//    it is only visible in wireframe mode.
	pDt->writeWebScene("gaeta_constraint_webScene");

	// 4) Write a file for the Geomview viewer
	pDt->showGeomview("gaeta_constraint_geomview.list");

	// 5) Write a postscript file that shows a 2D version of the triangulation
	pDt->show("gaeta_2d_postscript.ps");




	// ** Test ISO-contour computation
	std::cout<<"\n* Testing ISO-contour computation"<<std::endl;

	// 1) Create a global zone
	Zone2* pZone(pDt->createZone(NULL,ZL_GLOBAL));
	if(pZone==NULL)
	{
		std::cout<<"No zone, something went wrong."<<std::endl;
		return 1;
	}

	// 2) Get its triangles
	std::vector<Triangle2*> vTriangles;
	pZone->getTriangles(vTriangles);
	std::cout<<"Number of triangles: "<<vTriangles.size()<<std::endl;

	// 3) Compute ISO contours
	isoContours(vTriangles);


	std::cout<<"Normal program end"<<std::endl;
	delete pDt;
	return 0;
}


// Read an ASCII file containing xyz coordinates
void readTerrain(std::vector<Point2>& vInPoints,const char* filename)
{
	std::ifstream f(filename);
	if(f.is_open())
	{
		std::cout<<"reading "<<filename<<std::endl;
	}
	else
	{
		std::cout<<filename<<" is unreadable"<<std::endl;
		exit(1);
	}

	// For min/max statistics
	double minx(DBL_MAX),miny(DBL_MAX),minz(DBL_MAX);
	double maxx(-DBL_MAX),maxy(-DBL_MAX),maxz(-DBL_MAX);

    double x,y,z;
	while(f>>x>>y>>z)
	{
		vInPoints.push_back(Point2(x,y,z));
		if(x<minx) minx=x;
		if(y<miny) miny=y;
		if(z<minz) minz=z;
		if(x>maxx) maxx=x;
		if(y>maxy) maxy=y;
		if(z>maxz) maxz=z;
	}

	std::cout<<vInPoints.size()<<" points"<<std::endl;
	std::cout<<"RANGE X:  "<<minx<<" to "<<maxx<<endl;
	std::cout<<"RANGE Y:  "<<miny<<" to "<<maxy<<endl;
	std::cout<<"RANGE Z:  "<<minz<<" to "<<maxz<<endl;
}


// For demonstration purposes we create a polygon that can be
// inserted in the 2.5D triangulation.
void getPolygon(Fade_2D* pDt,std::vector<Segment2>& vSegments)
{
	// Initially the z coordinate is set to 0
	std::vector<Point2> vPoints;
	vPoints.push_back(Point2(2391300,4569100,0));
	vPoints.push_back(Point2(2402180,4571644,0));
	vPoints.push_back(Point2(2400700,4562800,0));
//	vPoints.push_back(Point2(2390000,4569100,0));
//	vPoints.push_back(Point2(2404300,4569100,0));
//	vPoints.push_back(Point2(2398300,4563800,0));

	// Adapt the z coordinates to the existing surface
	pDt->applyConstraintsAndZones(); // If constraints exist already then they should be inserted before to get correct heights
	for(size_t i=0;i<vPoints.size();++i)
	{
		Point2& p(vPoints[i]);
		double height;
		bool bHasHeight(pDt->getHeight(p.x(),p.y(),height));

		if(bHasHeight)
		{
			// Replace the point with one that has a height value
			vPoints[i]=Point2(p.x(),p.y(),height);
		}
		else
		{
			// Projection to the xy-plane has not hit a triangle
			std::cout<<"Warning, this point can't be projected to a triangle, no height: "<<std::endl<<p<<std::endl;
		}
	}

	// Create segments
	vSegments.push_back(Segment2(vPoints[0],vPoints[1]));
	vSegments.push_back(Segment2(vPoints[1],vPoints[2]));
	vSegments.push_back(Segment2(vPoints[2],vPoints[0]));

	cout<<"getPolygon: "<<endl;
	for(std::vector<Segment2>::iterator it(vSegments.begin());it!=vSegments.end();++it)
	{
		cout<<*it<<endl;
	}
	cout<<endl;
}


// The isoContours(..) function below computes ISO lines at 3 different heights
void isoContours(std::vector<Triangle2*>& vTriangles)
{
	// Create an instance of the IsoContours class and 3 height values
	IsoContours isoManager(vTriangles);
	double minZ(isoManager.getMinHeight());
	double maxZ(isoManager.getMaxHeight());
	vector<double> vHeights;
	vHeights.push_back(minZ+(maxZ-minZ)*1/4);
	vHeights.push_back(minZ+(maxZ-minZ)*2/4);
	vHeights.push_back(minZ+(maxZ-minZ)*3/4);


	// Open a result file
	std::ofstream f("contours.txt");
	if(!f.is_open())
	{
		std::cout<<"File contours.txt can't be written"<<std::endl;
		return;
	}
	std::cout<<"isoContours(): Writing result file contours.txt"<<std::endl;

	// For all height values
	for(size_t i=0;i<vHeights.size();++i)
	{
		double height(vHeights[i]);

		// The result will be stored in a vector of polygons or polylines
		std::vector<std::vector<Segment2> > vvIsoContours;

		// IsoContours::getContours() intersects the triangulation with a horizontal
		// plane at height z. Certain height values lead to degenerate intersections
		// and IsoContours::getContours() returns false in this case. Then a different
		// height value must be chosen (small random perturbation is sufficient).
		while(!isoManager.getContours(height,vvIsoContours,true))
		{
			double rnd=(1e-4*rand()/(RAND_MAX+1.0));
			std::cout<<"Degenerate intersection at height "<<height<<std::endl;
			height+=rnd;
		}

		f<<"\n\n\n** Number of contours at height "<<height<<": "<<vvIsoContours.size()<<std::endl<<std::endl;

		// Write the result
		for(size_t i=0;i<vvIsoContours.size();++i)
		{
			std::vector<Segment2>& vContour(vvIsoContours[i]);
			f<<"Contour no. "<<i<<" at z="<<height<<" consists of "<<vContour.size()<<" segments"<<std::endl;
			Point2 sourcePoint(vContour[0].getSrc());
			Point2 targetPoint(vContour.back().getTrg());
			if(sourcePoint==targetPoint)
			{
				f<<"the contour is a closed polygon"<<std::endl;
			}
			else
			{
				f<<"the contour is an open polyline"<<std::endl;
			}

			for(std::vector<Segment2>::iterator it(vContour.begin());it!=vContour.end();++it)
			{
				Segment2& seg(*it);
				f<<seg<<std::endl;
			}
			f<<std::endl;
		}
	}
}





