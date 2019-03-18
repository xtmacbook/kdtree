#include "GlbSpatialKdTree.h"


#include "../../comm/xMath.h"
extern const double INFINITYM;
extern const double kdTreeEpsilon;
static const double RayEpsilonOffset = 0.0001;

//AABB
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

		if (dist > enter_distance) //»˝∏ˆ÷·œÚ◊Ó∂Ãæ‡¿Î
			enter_distance = dist;
	}
	return enter_distance;
}

//AABB
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


bool  GlbGlobe::GLbKdTree::RayTravAlgSEQ(const KDTNodeM*node, const Ray&ray,Vec3&intersectionP)
{
	if(!node->is_leaf())
	{
		const KDTNodeM * leaf = node->backtrack_leaf(ray.origin);
		if(leaf == NULL)
			return false;
		else
			return RayTravAlgSEQ(leaf,ray,intersectionP);
	}
	else
	{
		Ray r(ray); //tmp

		const unsigned int *trIndices = node->tri_indices;

        double bestT = INFINITYM;

		for (unsigned int i = 0; i < node->num_tris; i++) 
		{
			Triangle &tri = meshTriangles[trIndices[i]];
            double t = -1.0;
            double u = 0.0;
            double v = 0.0;

            if(RayTriangleIntersect(ray.origin,ray.direction,tri.vertex[0],tri.vertex[1],tri.vertex[2],t,u,v))
            {
                if(t < bestT)
                {
                    bestT = t;
                    intersectionP =   tri.vertex[0] * (1 - u - v) +
                                      tri.vertex[1] * u+
                                      tri.vertex[2] * v;
                }
            }
		}
        if(bestT != INFINITYM) return  true;

		/* no intersection found */
		double exit_distance = get_exit_distance(node->box, ray);
		r.origin = r.origin +  r.direction * (exit_distance + RayEpsilonOffset) ;

		return RayTravAlgSEQ(node->parent,r,intersectionP);
	}
}

bool GlbGlobe::GLbKdTree::RayTracer(const Ray&ray,Vec3&intersectionP)
{
#ifndef KDTREE_NEIGHBORLINKS
	{
		// sah ∫Õ media space∂º√ª”– π”√
		if(!treeRoot && !treeRootM ) return false;

        if(sahUse)
        {
            return RayTravAlgRECA(treeRoot,ray,intersectionP);
        }
        else
        {
            const BoundingBox& bound = treeRootM->box;

            Ray r(ray);
            
            if(!bound.contains(ray.origin,kdTreeEpsilon))
            {
                double enter_distance = get_enter_distance(bound, ray);
                double exit_distance  = get_exit_distance(bound, ray);

                if (enter_distance > exit_distance)
                    return false;
                else
                    r.origin = r.origin +  r.direction * (enter_distance + RayEpsilonOffset);
            }

            RayTravAlgSEQ(treeRootM,r,intersectionP);
        }

	}
#else
	{
		return false;
	}
#endif

}

 bool GlbGlobe::GLbKdTree:: RayTravAlgRECA(const KDTNode* node,const Ray& ray,Vec3&intersectionP)
 {
    double a = 0.0; //entry point signed distances
    double b = 0.0; //exit point signed  distances

    double t ; // signed distance to the splitting plane

    a = get_enter_distance(node->box, ray);
    b = get_exit_distance(node->box, ray);

    return  true;

 }