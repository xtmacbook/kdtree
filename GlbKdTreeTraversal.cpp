#include "GlbSpatialKdTree.h"


#include "../../comm/xMath.h"

//test
#include <iostream.h>
extern const double INFINITYM;
extern const double kdTreeEpsilon;
extern const double KD_TREE_EPSILON = 0.000001;
static const double RayEpsilonOffset = 0.0001;

static const GlbGlobe::GLbKdTree * localKdTreePtr;
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


static bool treeLeafTrace(const GlbGlobe::KDTNodeM*node,const GlbGlobe::Ray ray,
                          const GlbGlobe::GLbKdTree*tree,GlbGlobe::Vec3&intersectionP)
{
    const unsigned int *trIndices = node->tri_indices;

    double bestT = INFINITYM;

    const GlbGlobe::Triangle * meshTriangles = tree->getMeshTriangles();

    for (unsigned int i = 0; i < node->num_tris; i++)
        {

        const  GlbGlobe::Triangle &tri = meshTriangles[trIndices[i]];

        double t = -1.0;
        double u = 0.0;
        double v = 0.0;

        if(RayTriangleIntersect(ray.origin,ray.direction,tri.vertex[0],tri.vertex[1],
                                tri.vertex[2],t,u,v))
            {
            if(t < bestT)
                {
                //test
                std::cout << "kdtree trignle :"
                    << "(" << tri.vertex[0].x() << " " << tri.vertex[0].y() << " " << tri.vertex[0].z() << ")" << std::endl
                    << "(" << tri.vertex[1].x() << " " << tri.vertex[1].y() << " " << tri.vertex[1].z() << ")" << std::endl
                    << "(" << tri.vertex[2].x() << " " << tri.vertex[2].y() << " " << tri.vertex[2].z() << ")" << std::endl;

                bestT = t;
                intersectionP =   tri.vertex[0] * (1 - u - v) +
                tri.vertex[1] * u+
                tri.vertex[2] * v;
                }
            }
        }
    return (bestT != INFINITYM);

}

static bool treeLeafTrace(const GlbGlobe::KDTNode*node,const GlbGlobe::Ray ray,
                          const GlbGlobe::GLbKdTree*tree,GlbGlobe::Vec3&intersectionP)
{
    std::vector<int>& tris = *(node->triangleIndices);

    double bestT = INFINITYM;

    const GlbGlobe::Triangle * meshTriangles = tree->getMeshTriangles();

    for(unsigned int i = 0;i < tris.size();i++)
        {
        double t = 0.0,u = 0.0,v = 0.0;

        const GlbGlobe::Triangle& tri = meshTriangles[tris[i]];

        if(RayTriangleIntersect(ray.origin,ray.direction,tri.vertex[0],
                                tri.vertex[1],tri.vertex[2],t,u,v))
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
    return (bestT != INFINITYM);
}
//template <typedef T>
static bool stackLessNeighLink(GlbGlobe::KDTNodeM*node,const GlbGlobe::Ray& ray,
                               double&t_entry,double&t_exit)
{
    bool intersection = false;

    double t_entry_prev = - INFINITYM;

    while(t_entry < t_exit && t_entry > t_entry_prev)
    {
        t_entry_prev = t_entry;

        GlbGlobe::Vec3 p_entry = ray.origin + (ray.direction *  t_entry );

        while(!node->is_leaf())
        {
            node = node->isPointToLeftOfSplittingPlane(p_entry)?node->left:node->right;
        }
        //the leaf
        const unsigned int * trisIndices = node->tri_indices;

        const GlbGlobe::Triangle * meshTris = localKdTreePtr->getMeshTriangles();

        for(unsigned int i = 0; i < node->num_tris;i++)
        {
            const GlbGlobe::Triangle&tri = meshTris[trisIndices[i]];

            double t,u,v;
            if(RayTriangleIntersect(ray.origin, ray.direction, tri.vertex[0],
                                        tri.vertex[1], tri.vertex[2], t, u, v))
            {
                if(t < t_exit)
                {
                    intersection = true;
                    t_exit = t;
                }
            }
        }

        double tmp_t_near,tmp_t_far;
        bool intersectsNodeBox = RayAABB(node->box._min, node->box._max, ray.origin, ray.direction, tmp_t_near, tmp_t_far);
        if(intersectsNodeBox)
        {
            t_entry = tmp_t_far;
        }
        else
        {
            break;
        }

        GlbGlobe::Vec3 p_exit = ray.origin + (ray.direction * t_entry);
        node = node->getNeighboringNode(p_exit);

        if(node == NULL)
        {
            return false;
        }
    }
}

bool  GlbGlobe::GLbKdTree::RayTravAlgSEQ(const KDTNodeM*node, const Ray&ray,Vec3&intersectionP)
{
#ifndef KDTREE_NEIGHBORLINKS
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
    #else
    return NULL;
    #endif
}

bool GlbGlobe::GLbKdTree::RayTracer(const Ray&r,Vec3&intersectionP)
{
    if(!sahUse)
    {
        const BoundingBox& bound = treeRootM->box;
        Ray ray(r);

        #if 0 // ok
        if(!bound.contains(ray.origin,kdTreeEpsilon))
        {
            double enter_distance = get_enter_distance(bound, ray);
            double exit_distance  = get_exit_distance(bound, ray);

            if (enter_distance > exit_distance) return false;
            else ray.origin = ray.origin +  ray.direction * (enter_distance + RayEpsilonOffset);
        }
        return RayTravAlgSEQ(treeRootM,ray,intersectionP);
        #else
        double t_near,t_far;
        bool intersects_root_node_bounding_box = RayAABB(bound._min, bound._max,ray.origin ,ray.direction, t_near, t_far);

        if ( intersects_root_node_bounding_box )
        {
            localKdTreePtr = this;
            bool hit = stackLessNeighLink( treeRootM, ray, t_near, t_far );
            if ( hit )
            {
                intersectionP = ray.origin  + ( ray.direction * t_far );
            }
            return hit;
        }
        else
        {
            return false;
        }
        #endif
    }
    else
    {

    }
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
        double bestT = INFINITYM;

        for(unsigned int i = 0;i < tris.size();i++)
            {
            double t = 0.0,u = 0.0,v = 0.0;

            Triangle& tri = meshTriangles[tris[i]];

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
        if(bestT != INFINITYM)
            {
            return true;
            }
        }

    return  false;

}

bool GlbGlobe::GLbKdTree::RayTravAlgRECA(const KDTNodeM* node,const Ray& ray,Vec3&intersectionP)
{
    double A = 0.0; //entry point signed distances
    double B = 0.0; //exit point signed  distances

    double t ; // signed distance to the splitting plane

    if( !RayAABB(node->box._min, node->box._max, ray.origin, ray.direction, A, B))
        {
        return false;
        }

        //stack to avoid recursive calls
    StackElemM stack[TERMINATION_CRITERIA_D_MAX];
    int stackPtr = 0; /*pointer to the stack*/

        //push the initial values onto the stack

    stack[stackPtr++] = StackElemM(node,A,B);

    KDTNodeM * farChild,*nearChild;
    const KDTNodeM* currNode;

    while (stackPtr != 0)
        {
        /*pop values from the stack*/
        StackElemM& ele = stack[--stackPtr];

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
                    stack[stackPtr++] = StackElemM(farChild,t,B);
                    currNode = nearChild;
                    B = t;
                    }
                }
            }

            //then is leaf
        const unsigned int *trIndices = currNode->tri_indices;

        double bestT = INFINITYM;

        for (unsigned int i = 0; i < currNode->num_tris; i++)
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
    KDTNode * farChild;
    const KDTNode *currNode;

    currNode = node;

    int enPt = 0; /*setup initial entry point*/
    stack[enPt].t = a;

    if(a >= 0.0)/*a ray with external origin*/
        {
        stack[enPt].pb = ray.origin + ray.direction * a;
        }
    else /*a ray with internal origin*/
        {
        stack[enPt].pb = ray.origin;
        }

    /*setup initial exit point in the stack*/

    int exPt = 1;
    stack[exPt].t = b;
    stack[exPt].pb = ray.origin + ray.direction * b;
    stack[exPt].node = NULL;

    while(true)
        {
        /*loop until a leaf is found*/
        while(!currNode->is_leaf())
            {
            double splitVal = currNode->splitEdge->splitPlanePosition;
            Axes axis = currNode->splitEdge->axis;

            //nextAxis: x->y y->z z->x
            //preAxis : x->z y->x z->y

            Axes nextAxis = (Axes)(((int)axis + 1) % 3);
            Axes preAxis  = (axis == X_axis)? Z_axis:((axis == Y_axis)?Z_axis:Y_axis);

            if(stack[enPt].pb[axis] <= splitVal)
            {
                if(stack[enPt].pb[axis] <= splitVal)
                {
                    currNode = currNode->left;
                    continue;
                }
                if(stack[exPt].pb[axis] == splitVal)
                {
                    currNode = currNode->right;
                    continue;
                }
                farChild = currNode->right;
                currNode = currNode->left;
            }
            else
            {
                if(splitVal < stack[exPt].pb[axis])
                {
                    currNode = currNode->right;
                    continue;
                }
                farChild = currNode->left;
                currNode = currNode->right;
            }

            t = (splitVal - ray.origin[axis]) / ray.direction[axis];

            /*setup the exit point*/
            int tmp = exPt;
            exPt++;

            if(exPt == enPt)
            {
                exPt++;
            }

            //push values onto the stack
            stack[exPt].prev = tmp;
            stack[exPt].t = t;
            stack[exPt].node = farChild;
            stack[exPt].pb[axis] = splitVal;
            stack[exPt].pb[nextAxis] = ray.origin[nextAxis] + t * ray.direction[nextAxis];
            stack[exPt].pb[preAxis] = ray.origin[preAxis] + t * ray.direction[preAxis];

            }
            //then the leaf

            if(/*any intersection*/1)
            {
                return true;
            }

            //pop from the stack
            enPt = exPt;
            currNode = stack[exPt].node;

            exPt = stack[enPt].prev;

        }

        return false;
}



