
/********************************************************************
* Copyright (c) 2014 北京超维创想信息技术有限公司
* All rights reserved.
* @file    GlbSAH.h
* @brief 该kdtree树结构 使用 surface area heuristic (SAH)来定位   
* splitting plane的位置
* Reference:http://graphics.cs.illinois.edu
* SAH (Surface Area Heuristic) Kd-tree implementation.

* Reference:
* Ingo Wald and Vlastimil Havran. On building fast kd-trees for ray tracing, and on doing that in
* O(N*log N). In Proceedings of the 2006 IEEE Symposium on Interactive Ray Tracing, 2006. (accepted
* for publication, minor revision pending).
* @version 1.0
* @author  xt
* @date    2016-09-12 11:30
*********************************************************************/


#ifndef _GLB_SAH_H_
#define _GLB_SAH_H_

#include <osg/BoundingBox>

#define Ct_DEFAULT 15.0 //
#define Ci_DEFAULT 20.0 //
#define emptyBonus_DEFAULT 0.0


namespace GlbGlobe
{
	class GlbSAH 

	{
	public:
		GlbSAH(double Ct = Ct_DEFAULT,
			double Ci = Ci_DEFAULT,
			double emptyBonus = emptyBonus_DEFAULT);

		double calculateSAH(const osg::BoundingBox &nodeExtent, char axis,
			unsigned int  nA, unsigned int nB, double position) const;

		// convenient function
		double operator()(const osg::BoundingBox &nodeExtent, char axis,
			unsigned int nA, unsigned int nB, double position) const;

		// parameters in SAH formula
		const double m_Ci;          // triangle intersection cost
		const double m_Ct;	         // node traversal cost
		const double m_emptyBonus;  // takes value between 0 and 1
		// closer to 1 => encourage empty space splitting
		// closer to 0 => discourage
	};

}


#endif 
