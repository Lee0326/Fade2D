// Copyright (C) Geom Software e.U, Bernhard Kornberger, Graz/Austria
//
// This file is part of the Fade2D library. The student license is free
// of charge and covers personal non-commercial research. Licensees
// holding a commercial license may use this file in accordance with
// the Commercial License Agreement.
//
// This software is provided AS IS with NO WARRANTY OF ANY KIND,
// INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE.
//
// Please contact the author if any conditions of this licensing are
// not clear to you.
//
// Author: Bernhard Kornberger, bkorn (at) geom.at
// http://www.geom.at

#pragma once

#include "common.h"
#if GEOM_PSEUDO3D==GEOM_TRUE

#include "Point2.h"



namespace GEOM_FADE25D {

struct EMData; // FWD


/** \brief EfficientModel
*
* Point clouds from terrains collected by scanners are often
* oversampled. The present class aims to reduce these point clouds
* in a controlled way such that the resulting surface keeps a certain
* user specified quality.\n
*
* \note EfficientModel is pre-released. But quite some further
* functionality like automatic line recognition and simplification
* exists already under the hood and will be completed soon. Thus
* the interface of the present class may change in the near future.
*/
class CLASS_DECLSPEC EfficientModel
{
public:
	explicit EfficientModel(const std::vector<Point2>& vPoints);
	~EfficientModel();

	void go();

/** \brief Extract a subset of points
*
* This method extracts a subset of the original point cloud that
* represents the model more efficiently. Thereby the original and
* the simplified model cover the same area.
*
* @param maxError is the maximum height difference between the original
* points and the simplified model.
* @param [out] vEfficientPointsOut is used to return a subset of the
* original points that represents the model more efficiently.
*/
	void extract(double maxError,std::vector<Point2>& vEfficientPointsOut);

protected:
	void part1_extractFC();
	void part2_setWeights();

	void sortVtx(std::vector<Point2*>& vVtx);
	int insertKeepError(double maxErr,std::vector<Point2*>& vA,std::vector<Point2*>& vB);
	void insertMinHull();
	void show(const std::string& name);

	EMData* pEMData;
private:
	EfficientModel(const EfficientModel&);
};

} // NAMESPACE

// up to here: if GEOM_PSEUDO3D==GEOM_TRUE
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif
