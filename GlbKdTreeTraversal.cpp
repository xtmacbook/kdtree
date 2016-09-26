#include "GlbSpatialKdTree.h"
#include "GlbRayTracer.h"
#include "../comm/xMath.h"

extern const double INFINITY;
extern const double kdTreeEpsilon;
//AABB 射线进入点位置
static double get_enter_distance(const osg::BoundingBox &box, const GlbGlobe::Ray &r)
{
	double enter_distance = -1e100; 
	for (int i = 0; i < 3; i++)
	{
		double dist;
		if (r.direction[i] > kdTreeEpsilon)
			dist = (box._min[i] - r.origin[i]) / r.direction[i];
		else if (r.direction[i] < -kdTreeEpsilon)
			dist = (box._max[i] - r.origin[i]) / r.direction[i];
		else 
			continue;

		if (dist > enter_distance) //三个轴向最短距离
			enter_distance = dist;
	}
	return enter_distance;
}

//AABB 射线进出点位置
static double get_exit_distance(const osg::BoundingBox &box, const GlbGlobe::Ray &r) 
{
	double exit_distance = 1e100;
	for (int i = 0; i < 3; i++) 
	{
		double dist;
		if (r.direction[i] > kdTreeEpsilon)
			dist = (box._max[i] - r.origin[i]) / r.direction[i];
		else if (r.direction[i] < -kdTreeEpsilon)
			dist = (box._min[i] - r.origin[i]) / r.direction[i];
		else
			continue;

		if (dist < exit_distance)
			exit_distance = dist;
	}
	return exit_distance;
}


bool  GlbGlobe::GLbKdTree::GetIntersecting(const KDTNode*node, const Vec3&origin ,const Vec3&dir,Vec3&intersectionP)
{
	if(!node->is_leaf())
	{
		const KDTNode * leaf = node->backtrack_leaf(origin);
		if(leaf == NULL)
			return false;
		else
			return GetIntersecting(leaf,origin,dir,intersectionP);
	}
	else
	{
		double dist;
		Ray ray(origin,dir);

		const std::vector<int>& trisV = *(node->triangleIndices);

		for (unsigned int i = 0; i < trisV.size(); i++) 
		{
			Triangle &tri = meshTriangles[i];

			if(tri.getRayIntersection(origin,dir,intersectionP))
			{
				return true;
			}
		}

		/* no intersection found */
		double exit_distance = get_exit_distance(node->box, ray);
		ray.origin = ray.origin +  ray.direction * (exit_distance + 0.00001) ;

		return GetIntersecting(node->parent,ray.origin,ray.direction,intersectionP);
	}
}

static unsigned int Test = 0;//test

bool  GlbGlobe::GLbKdTree::GetIntersectingM(const KDTNodeM*node, const Vec3&origin ,const Vec3&dir,Vec3&intersectionP)
{
	Test++; //test
	unsigned int nowLeve = node->leve; //test

	if(!node->is_leaf())
	{
		const KDTNodeM * leaf = node->backtrack_leaf(origin);
		if(leaf == NULL)
			return false;
		else
			return GetIntersectingM(leaf,origin,dir,intersectionP);
	}
	else
	{
		double dist;
		Ray ray(origin,dir);

		const unsigned int *trIndices = node->tri_indices;

		for (unsigned int i = 0; i < node->num_tris; i++) 
		{
			Triangle &tri = meshTriangles[trIndices[i]];

			if(tri.getRayIntersection(origin,dir,intersectionP))
			{
				return true;
			}
		}

		/* no intersection found */
		double exit_distance = get_exit_distance(node->box, ray);
		ray.origin = ray.origin +  ray.direction * (exit_distance + 0.00001) ;

		return GetIntersectingM(node->parent,ray.origin,ray.direction,intersectionP);
	}
}

bool GlbGlobe::GLbKdTree::GetIntersectPoint(const Vec3&ori,const Vec3&dir,Vec3&intersectionP)
{
#ifndef KDTREE_NEIGHBORLINKS
	{
		// sah 和 media space都没有使用
		if(!treeRoot && !treeRootM ) return false;

		
		BoundingBox& bound = (sahUse) ? treeRoot->box :treeRootM->box;

		Ray ray(ori,dir);
		//线起始点不在树内部
		if(!bound.contains(ori,kdTreeEpsilon))
		{
			double enter_distance = get_enter_distance(bound, ray);
			double exit_distance  = get_exit_distance(bound, ray);

			if (enter_distance > exit_distance)
				return false;
			else
				ray.origin = ray.origin +  ray.direction * (enter_distance + kdTreeEpsilon);
		}

		if(sahUse)
		{
			return GetIntersecting(treeRoot,ray.origin,ray.direction,intersectionP);
		}
		else
		{
			return GetIntersectingM(treeRootM,ray.origin,ray.direction,intersectionP);
		}
	}
#else
	{
		return false;
	}
#endif

}

 