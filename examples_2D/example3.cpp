#include <Fade_2D.h>
#include <stdlib.h>
#include <stdio.h>
using namespace GEOM_FADE2D;
using namespace std;


int example3_main()
{
	std::cout<<"\n";
	std::cout<<"example3: Constraints - Enforce constraint edges\n";
	std::cout<<"* Constrained Delaunay\n";
	std::cout<<"* Conforming Delaunay\n\n";

	// 1) Generate some input points
	std::vector<Point2> vInputPoints;
	vInputPoints.push_back(Point2(-100,-100));
	vInputPoints.push_back(Point2(+100,+100));
	vInputPoints.push_back(Point2(-50,-70));
	vInputPoints.push_back(Point2(-50,-30));
	vInputPoints.push_back(Point2(50,70));
	vInputPoints.push_back(Point2(50,30));

	// 2) Triangulate the points and show
	Fade_2D dt;
	dt.insert(vInputPoints);
	dt.show("example3_noConstraints.ps");

	// 3) Prepare a vector of one or more edges to be enforced
	std::vector<Segment2> vSegments;
	vSegments.push_back(Segment2(vInputPoints[0],vInputPoints[1]));

	// 4) Use one of two basic constraint insertion strategies:
	//    CIS_CONSTRAINED_DELAUNAY inserts the whole constraint edge.
	//    In contrast, CIS_CONFORMING_DELAUNAY subdivides the constraint
	//    edge in order to maintain the Delaunay property.

	//    a) Constrained Delaunay - the segment will not be subdivided,
	//       except it intersects an existing vertex which is different
	//       from its endpoints.
	ConstraintGraph2* pCG=dt.createConstraint(vSegments,CIS_CONSTRAINED_DELAUNAY);

	//    b) Conforming Delaunay - the segment will be subdivided such
	//       that its subsegments appear automatically in the Delaunay
	//       triangulation. Better triangle quality results but possibly
	//       much more triangles.
	//ConstraintGraph2* pCG=dt.createConstraint(vSegments,CIS_CONFORMING_DELAUNAY);

	// The very last step that establishes the conforming constraints (if any) in 
	// the triangulation is a call to applyConstraintsAndZones(). If afterwards 
	// insertions are made, then conforming constraints may disappear and 
	// applyConstraintsAndZones() must be repeated. CIS_CONSTRAINED_DELAUNAY
	// does not require applyConstraintsAndZones().
	dt.applyConstraintsAndZones();

	// Get the points of the (possibly subdivided) constraint edge
	// and visualize.
	Visualizer2 vis2("example3_withConstraints.ps");
	dt.show(&vis2);

	vector<Point2*> vPointsOfConstraintEdge;
	pCG->getPolygonVertices(vPointsOfConstraintEdge);
	for(size_t i=0;i+1<vPointsOfConstraintEdge.size();++i)
	{
		Point2* p0(vPointsOfConstraintEdge[i]);
		Point2* p1(vPointsOfConstraintEdge[i+1]);
		vis2.addObject(Segment2(*p0,*p1),Color(1,0,0,0.01));
		vis2.addObject(*p0,Color(0,0,1,0.1));
		vis2.addObject(*p1,Color(0,0,1,0.1));
	}
	vis2.writeFile();

	return 0;
}



