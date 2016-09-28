#include "GlbSpatialKdTree.h"

#include "../comm/xMath.h"

static const GlbGlobe::GLbKdTree * localKdTreePtr;

//AABB entry
static double get_enter_distance(const osg::BoundingBox &box, const GlbGlobe::Ray &r)
{
	double enter_distance = -1e100;
	for (int i = 0; i < 3; i++)
	{
		double dist;
		if (r.direction[i] > KD_TREE_EPSILON)
			dist = (box._min[i] - r.origin[i]) / r.direction[i];
		else if (r.direction[i] < -KD_TREE_EPSILON)
			dist = (box._max[i] - r.origin[i]) / r.direction[i];
		else
			continue;

		if (dist > enter_distance) 
			enter_distance = dist;
	}
	return enter_distance;
}

//AABB exit
static double get_exit_distance(const osg::BoundingBox &box, const GlbGlobe::Ray &r)
{
	double exit_distance = 1e100;
	for (int i = 0; i < 3; i++)
	{
		double dist;
		if (r.direction[i] > KD_TREE_EPSILON)
			dist = (box._max[i] - r.origin[i]) / r.direction[i];
		else if (r.direction[i] < -KD_TREE_EPSILON)
			dist = (box._min[i] - r.origin[i]) / r.direction[i];
		else
			continue;

		if (dist < exit_distance)
			exit_distance = dist;
	}
	return exit_distance;
}

//定位包含该点的叶子节点,每次都从根节点开始
//template<typename T>
static const GlbGlobe::KDTNodeM* LocateLeaf(const const GlbGlobe::KDTNodeM* node,const GlbGlobe::Vec3&point)
{
	const GlbGlobe::KDTNodeM * currNode = node;

	if(!currNode->box.contains(point,KD_TREE_EPSILON)) return NULL;

	while(!currNode->is_leaf())
	{
		GlbGlobe::Axes axis = currNode->splitEdge->axis;

		//if(point[axis]  < currNode->right->box._min[axis])
		if(point[axis]  < currNode->splitEdge->splitPlanePosition)
		{
			currNode = currNode->left;
		}
		else
		{
			currNode = currNode->right;
		}

	}

	return currNode;
}

//template <typedef T>
static bool stackLessNeighLink(GlbGlobe::KDTNodeM*node,const GlbGlobe::Ray& ray,double&t_entry,double&t_exit)
{
	bool intersection = false;

	double t_entry_prev = - KDTREEDOUBLEINFINITYM;

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
		bool intersectsNodeBox = RayAABB(node->box._min, node->box._max, ray.origin, ray.direction,
			tmp_t_near, tmp_t_far);
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

	return intersection;
}

/* Sequential ray traversal algorithm
* Reference :Heuristic Ray Shooting Algorithms by Vlastimil Havran
*/

template <typename T>
static bool RayTravAlgSEQ(const T*rootNode, const GlbGlobe::Ray&ray,GlbGlobe::Vec3&intersectionP)
{
	double a,b; /*entry/exit 点 距离射线原点的距离*/

	GlbGlobe::Vec3 point;

	if(!RayAABB(rootNode->box._min,rootNode->box._max,ray.origin,ray.direction,a,b))
	{
		return false;
	}

	const T * currentNode =  rootNode;

	//原点在box内部
	if(a < 0.0)
	{
		point = ray.origin;
	}
	else
	{
		point = ray.origin + ray.direction *(a +  RAY_EPSILON_OFFSET);
	}

	currentNode = LocateLeaf(rootNode,point);

	static GlbGlobe::Vec3 prePoint = point;

	while(currentNode && currentNode->is_leaf())
	{
		if(RayAABB(currentNode->box._min,currentNode->box._max,ray.origin,
			ray.direction,a,b))
		{
			if(currentNode->treeLeafTrace(ray,intersectionP,localKdTreePtr))
				return true;
		}
		point = ray.origin + ray.direction * (b + 0.0001);

		if(prePoint == point)
		{
			return false;
		}
		prePoint = point;
		currentNode = LocateLeaf(rootNode,point);
	}

	return false;
}

//parent
static bool RayTravParent(const GlbGlobe::KDTNodeM*node, const GlbGlobe::Ray&ray,GlbGlobe::Vec3&intersectionP)
{
#ifndef KDTREE_NEIGHBORLINKS
	if(!node->is_leaf())
	{
		const GlbGlobe::KDTNodeM * leaf = node->backtrack_leaf(ray.origin);
		if(leaf == NULL)
			return false;
		else
			return RayTravParent(leaf,ray,intersectionP);
	}
	else
	{
		GlbGlobe::Ray r(ray); //tmp

		const unsigned int *trIndices = node->tri_indices;

		double bestT = KDTREEDOUBLEINFINITYM;

		const GlbGlobe::Triangle* meshTriangles = localKdTreePtr->getMeshTriangles();

		for (unsigned int i = 0; i < node->num_tris; i++)
		{
			const GlbGlobe::Triangle &tri = meshTriangles[trIndices[i]];
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
		if(bestT != KDTREEDOUBLEINFINITYM) return  true;

		/* no intersection found */
		double exit_distance = get_exit_distance(node->box, ray);
		r.origin = r.origin +  r.direction * (exit_distance +  RAY_EPSILON_OFFSET) ;

		return RayTravParent(node->parent,r,intersectionP);
	}
#else
	return NULL;
#endif
}

static bool RayTravAlgRECA(const GlbGlobe::KDTNodeM* node,const GlbGlobe::Ray& ray,GlbGlobe::Vec3&intersectionP)
{
	double A = 0.0; //进入AABB的进入点距离(有效距离 包括负值)
	double B = 0.0; //出AABB的点距离(有效距离 包括负值)

	if( !RayAABB(node->box._min, node->box._max, ray.origin, ray.direction, A, B))
	{
		return false;
	}

	double t ; // 射线原点距离splittting plane的有效距离

	//堆栈用来避免递归
	GlbGlobe::StackElem<GlbGlobe::KDTNodeM> stack[TERMINATION_CRITERIA_D_MAX];
	int stackPtr = 0; /*指向堆栈*/

	stack[stackPtr++] = GlbGlobe::StackElem<GlbGlobe::KDTNodeM>(node,A,B);

	GlbGlobe::KDTNodeM * farChild,*nearChild;
	const GlbGlobe::KDTNodeM* currNode;

	while (stackPtr != 0)
	{
		GlbGlobe::StackElem<GlbGlobe::KDTNodeM>& ele = stack[--stackPtr];

		currNode = ele.node;

		A = ele.a;
		B = ele.b;

		while(!currNode->is_leaf())
		{
			const GlbGlobe::BoxEdge * splitting = currNode->splitEdge;
			//double diff = currNode->right->box._min[splitting->axis] - ray.origin[splitting->axis];
			double diff = splitting->splitPlanePosition - ray.origin[splitting->axis];

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
					stack[stackPtr++] = GlbGlobe::StackElem<GlbGlobe::KDTNodeM>(farChild,t,B);
					currNode = nearChild;
					B = t;
				}
			}
		}

		//叶子节点
		const unsigned int *trIndices = currNode->tri_indices;

		double bestT = KDTREEDOUBLEINFINITYM;

		const GlbGlobe::Triangle* meshTriangles = localKdTreePtr->getMeshTriangles();

		for (unsigned int i = 0; i < currNode->num_tris; i++)
		{
			const GlbGlobe::Triangle &tri = meshTriangles[trIndices[i]];
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

		if(bestT != KDTREEDOUBLEINFINITYM) return  true;

	}

	return  false;

}

/*
* Reference :Heuristic Ray Shooting Algorithms by Vlastimil Havran
*/
/*只是依靠进入点、出点和splitting plane位置判断当前点位置*/
template <typename T>
static bool RayTravAlgRECB(const T * node,const GlbGlobe::Ray&ray,GlbGlobe::Vec3&intersectionP)
{
	double A = 0.0; //进入AABB的进入点距离(有效距离 包括负值)
	double B = 0.0; //出AABB的点距离(有效距离 包括负值)

	double t;

	if( !RayAABB(node->box._min, node->box._max, ray.origin, ray.direction, A, B))
	{
		return false;
	}

	GlbGlobe::StackElemA<T> stack[TERMINATION_CRITERIA_D_MAX];

	/*指向当前和远处的节点*/
	T * farChild;
	const T *currNode = node;

	int enPt = 0; /*进入点堆栈指针*/
	stack[enPt].t = A;

	if(A >= 0.0)/*射线原点在AABB外部*/
	{
		stack[enPt].pb = ray.origin + ray.direction * A;
	}
	else /*射线原点在AABB内部*/
	{
		stack[enPt].pb = ray.origin;
	}

	/*将初始化的出点放入堆栈*/
	int exPt = 1;
	stack[exPt].t = B;
	stack[exPt].pb = ray.origin + ray.direction * B;
	stack[exPt].node = NULL;

	while(currNode && currNode != NULL)
	{
		/*循环到叶子*/
		while(!currNode->is_leaf())
		{
			double splitVal = currNode->splitEdge->splitPlanePosition;
			GlbGlobe::Axes axis = currNode->splitEdge->axis;

			//nextAxis: x->y y->z z->x
			//preAxis : x->z y->x z->y
			GlbGlobe::Axes nextAxis = (GlbGlobe::Axes)(((int)axis + 1) % 3);
			GlbGlobe::Axes preAxis  = (axis == GlbGlobe::X_axis)? GlbGlobe::Z_axis:((axis == GlbGlobe::Y_axis)?
							GlbGlobe::Z_axis:GlbGlobe::Y_axis);

			if(stack[enPt].pb[axis] <= splitVal)
			{
				if(stack[exPt].pb[axis] <= splitVal)
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
			//射线原点到splitting plane距离
			t = (splitVal - ray.origin[axis]) / ray.direction[axis];

			/*设置出点*/
			int tmp = exPt;
			exPt++;

			if(exPt == enPt)
			{
				exPt++;
			}

			//将farChild放入堆栈,继续遍历currentNode(nearChild)
			stack[exPt].prev = tmp;
			stack[exPt].t = t;
			stack[exPt].node = farChild;
			stack[exPt].pb[axis] = splitVal;
			stack[exPt].pb[nextAxis] = ray.origin[nextAxis] + t * ray.direction[nextAxis];
			stack[exPt].pb[preAxis]  = ray.origin[preAxis]  + t * ray.direction[preAxis];

		}

		//叶子
		if(currNode->is_leaf())
		{
			if( currNode->treeLeafTrace(ray,intersectionP,localKdTreePtr))
				return true;
		}

		//从堆栈弹出
		enPt = exPt;
		currNode = stack[exPt].node;

		exPt = stack[enPt].prev;

	}

	return false;
}

bool GlbGlobe::GLbKdTree::RayTracer(const Ray&r,Vec3&intersectionP,int t)
{
	if(!sahUse)
	{
		const BoundingBox& bound = treeRootM->box;
		Ray ray(r);
		if(t == 0)
		{
			/*parent node traversal*/
			if(!bound.contains(ray.origin,0.00000001))
			{
				double enter_distance = get_enter_distance(bound, ray);
				double exit_distance  = get_exit_distance(bound, ray);

				if (enter_distance > exit_distance) return false;
				else ray.origin = ray.origin +  ray.direction * (enter_distance + RAY_EPSILON_OFFSET);
			}
			
			localKdTreePtr = this;
			return RayTravParent(treeRootM,ray,intersectionP);
		}
 
		if(t == 1)
		  {
			// NeightLink
			double t_near,t_far;
			if ( RayAABB(bound._min, bound._max,ray.origin ,ray.direction, t_near, t_far) )
			{
				localKdTreePtr = this;

				if ( stackLessNeighLink( treeRootM, ray, t_near, t_far ) )
				{
					intersectionP = ray.origin  + ( ray.direction * t_far );
					return true;
				}
			}
			return false;
		  }
 
		if(t == 2)
		{
			localKdTreePtr = this;
			return RayTravAlgSEQ(treeRootM,ray,intersectionP);
		}
	
		if(t == 3)
		{
			localKdTreePtr = this;
			return RayTravAlgRECA(treeRootM,ray,intersectionP);
		}
		
		if(t == 4)
		{
			localKdTreePtr = this;
			return RayTravAlgRECB(treeRootM,ray,intersectionP);
		}
	}
	else
	{

	}

	return false;
}
