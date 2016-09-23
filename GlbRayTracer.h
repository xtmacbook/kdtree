/********************************************************************
* Copyright (c) 2014 北京超维创想信息技术有限公司
* All rights reserved.
* @file    GlbRayTracer.h
* @ 光线追踪 参考: 
*	Stackless KD-Tree Traversal for High Performance GPU Ray Tracing
*   Volume 26, Number 3, 2007
*
*
* @version 1.0
* @author  xt
* @date    2016-09-12 11:30
*********************************************************************/

#include <osg/Vec3>

class GLbKdTree;


namespace GlbGlobe
{

	struct Ray
	{
		Ray(osg::Vec3 orig, osg::Vec3 dir):origin(orig),direction(dir)
		{ }
	
		osg::Vec3 origin;
		osg::Vec3 direction;

	};

	osg::Vec3 get_ray_trace(const Ray&ray,GLbKdTree*kdtree,unsigned int depth);

}
