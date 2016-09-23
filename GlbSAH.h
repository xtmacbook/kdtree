
/********************************************************************
* Copyright (c) 2014 ������ά������Ϣ�������޹�˾
* All rights reserved.
* @file    GlbSpatialKdTree.h
* @brief ��kdtree���ṹ ʹ�� surface area heuristic (SAH)����λ   
* splitting plane��λ��
*	��ϸ�ο�:http://graphics.cs.illinois.edu

* @version 1.0
* @author  xt
* @date    2016-09-12 11:30
*********************************************************************/


#ifndef _GLB_SAH_H_
#define _GLB_SAH_H_

#include <osg/BoundingBox>

#define Ct_DEFAULT 15.0 //���ṹ����ʱ�� 
#define Ci_DEFAULT 20.0 //���ṹ��ʱ��
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
