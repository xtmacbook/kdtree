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

		if (dist > enter_distance) //
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
		// sah 和 media space都没有使用
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

 bool GlbGlobe::GLbKdTree::RayTravAlgRECA(const KDTNode* node,const Ray& ray,Vec3&intersectionP)
 {
    double A = 0.0; //entry point signed distances
    double B = 0.0; //exit point signed  distances

    double t ; // signed distance to the splitting plane

    if( !RayAABB(node->box._min, node->box._max, ray.origin, ray.direction, A, B))
    {
        return false;
    }

    //stack to avoid recursive calls
    StackElem stack[TERMINATION_CRITERIA_D_MAX];
    int stackPtr = 0; /*pointer to the stack*/

    //push the initial values onto the stack
    
    stack[stackPtr++] = StackElem(node,A,B);

    KDTNode * farChild,*nearChild;
    const KDTNode* currNode;

    while (stackPtr != 0)
    {
        /*pop values from the stack*/
        StackElem& ele = stack[--stackPtr];

        currNode = ele.node;

        while(!currNode->is_leaf())
        {
            A = ele.a;
            B = ele.b;

            const GlbGlobe::BoxEdge * splitting = currNode->splitEdge;
            double diff = currNode->right->box._min[splitting->axis] - ray.origin[splitting->axis];

            // the signed distance to splitting plane
            t = diff / ray.direction[splitting->axis];

            if(diff > 0.0)
            {
                nearChild = currNode->left;
                farChild  = currNode->right;
            }
            else
            {
                nearChild = currNode->right;
                farChild  = currNode->left;
            }

            if((t > B) || (t < 0.0))
            {
                currNode = nearChild;
            }
            else
            {
                if(t < A)
                {
                    currNode = farChild;
                }
                else
                {
                    stack[stackPtr++] = StackElem(farChild,t,B);
                    currNode = nearChild;
                    B = t;
                }
            }
        }

        //then is leaf
        std::vector<int>& tris = *(currNode->triangleIndices);

        for(unsigned int i = 0;i < tris.size();i++)
        {
            Triangle& tri = meshTriangles[tris[i]];

            if(tri.getRayIntersection(ray.origin,ray.direction,intersectionP))
            {
                return true;
            }
        }

    }

    return  false;

 }

 bool GlbGlobe::GLbKdTree::RayTravAlgRECB(const KDTNode * node,const Ray&ray,Vec3&intersectionP)
 {
    double a,b;
    double t;

    if( !RayAABB(node->box._min, node->box._max, ray.origin, ray.direction, a, b))
    {
        return false;
    }

    StackElemA stack[TERMINATION_CRITERIA_D_MAX];

    /*point to the far child node and current node*/
    KDTNode * farChild,*currNode;
    //currNode = node;


 }

