
#include "GlbSpatialKdTree.h"
#include <limits>

#include <osg/TriangleFunctor>
#include "../../comm/xMath.h"

using namespace GlbGlobe;

extern const double INFINITYM = std::numeric_limits<double>::max();
extern const double kdTreeEpsilon = 1E-9;


#define MAX(a, b) ((a) < (b) ? (b) : (a))
#define MAX3(a, b, c) MAX( MAX(a ,b) ,c)
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MIN3(a, b, c) MIN( MIN(a, b) ,c)


BoxEdge::BoxEdge()
{
	splitPlanePosition = INFINITYM;
}

Ray::Ray(const Ray&r)
{ origin = r.origin; direction = r.direction; }

Triangle::Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
{
	vertex[0] = v0;
	vertex[1] = v1;
	vertex[2] = v2;

	bound._max[0] = MAX3(vertex[0][0], vertex[1][0], vertex[2][0]);
	bound._min[0] = MIN3(vertex[0][0], vertex[1][0], vertex[2][0]);

	bound._max[1] = MAX3(vertex[0][1], vertex[1][1], vertex[2][1]);
	bound._min[1] = MIN3(vertex[0][1], vertex[1][1], vertex[2][1]);

	bound._max[2] = MAX3(vertex[0][2], vertex[1][2], vertex[2][2]);
	bound._min[2] = MIN3(vertex[0][2], vertex[1][2], vertex[2][2]);

}

bool Triangle::getRayIntersection(const Vec3&origin,const Vec3&dir,Vec3&intersectionPoint)
{

	 double t, u, v;

	if( RayTriangleIntersect(origin,dir,vertex[0],vertex[1],vertex[2],t,u,v))
	{

		intersectionPoint = vertex[0] * (1 - u - v) + vertex[1] * u + vertex[2] * v;

		return true;
	}

	 return false;
}

TriangleM::TriangleM(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
{
	vertex[0] = v0;
	vertex[1] = v1;
	vertex[2] = v2;
}

KDTNode::KDTNode(void)
{
	left = NULL;
	right = NULL;

	splitEdge = new BoxEdge();

#ifdef KDTREE_NEIGHBORLINKS
	for ( int i = 0; i < 6; ++i )
	{
		ropes[i] = NULL;
	}
#endif

}

KDTNode::~KDTNode(void)
{

	if ( left ) 
	{
		delete left;
		left = NULL;
	}
	if ( right ) 
	{
		delete right;
		right = NULL;
	}

	if(splitEdge)
	{
		delete splitEdge;
		splitEdge = NULL;
	}

#ifdef KDTREE_NEIGHBORLINKS
	for ( int i = 0; i < 6; ++i )
	{
		ropes[i] = NULL;
	}
#endif
}

KDTNodeM::KDTNodeM(void)
{
	left = NULL;
	right = NULL;

	splitEdge = new BoxEdge();

	tri_indices = NULL;
	num_tris = 0;

#ifdef KDTREE_NEIGHBORLINKS
	for ( int i = 0; i < 6; ++i )
	{
		ropes[i] = NULL;
	}
#else
	parent = NULL;
#endif

}
KDTNodeM::~KDTNodeM(void)
{

	if ( left ) 
	{
		delete left;
		left = NULL;
	}
	if ( right ) 
	{
		delete right;
		right = NULL;
	}

	if(splitEdge)
	{
		delete splitEdge;
		splitEdge = NULL;
	}

#ifdef KDTREE_NEIGHBORLINKS
	for ( int i = 0; i < 6; ++i )
	{
		ropes[i] = NULL;
	}
#endif
	if(tri_indices)
	{
		delete [] tri_indices;
		tri_indices = NULL;
	}
}

const KDTNodeM* KDTNodeM::backtrack_leaf(const Vec3 &point) const
{
	if (box.contains(point,kdTreeEpsilon))
	{
		return find_leaf(point);		
	}
	else if (parent == NULL)
	{
		return NULL;
	}
	else
	{
		return parent->backtrack_leaf(point);
	}
}

const KDTNodeM * KDTNodeM::find_leaf(const Vec3 &point) const
{
	//递归寻找最终的叶节点
	if (is_leaf())
	{
		return this;
	}
	else if (point[splitEdge->axis] < splitEdge->splitPlanePosition)
	{
		return left->find_leaf(point);
	}
	else
	{
		return right->find_leaf(point);
	}
}

GLbKdTree::GLbKdTree(bool t):sahUse(t)
{
	treeRoot =  NULL;
	treeRootM = NULL;
}


GlbGlobe::GLbKdTree::~GLbKdTree(void )
{
	if(treeRoot)
	{
		delete treeRoot;
		treeRoot = NULL;
	}

	if(treeRootM)
	{
		delete treeRootM;
		treeRootM = NULL;
	}
}



std::vector<Triangle> triangles;

struct TriangleIntersector
{
	 

	TriangleIntersector()
	{
		 
	}

	inline void operator () (const osg::Vec3& v1,const osg::Vec3& v2,const osg::Vec3& v3, bool treatVertexDataAsTemporary)
	{
		triangles.push_back(Triangle(v1,v2,v3));
	}

};

unsigned int  GLbKdTree::GetMeshTriangleAndVertexs(const osg::Drawable*drawable)
{
	if(!drawable) return false;

	//unsigned int drawableNum = geometry->getdraw
	osg::TriangleFunctor<TriangleIntersector> ti;
	
	drawable->accept(ti);

	meshTriangles = new Triangle[triangles.size()];

	std::vector<Triangle>::iterator iter = triangles.begin();
	unsigned int triIndice = 0;
	for(;iter != triangles.end();iter++)
	{
		meshTriangles[triIndice++] = *iter;
	}

	triangleSize = triangles.size();

	return triangles.size();
}


void GLbKdTree::expandBoundBox(const Vec3&v,Vec3&_max,Vec3&_min)
{
	if ( v.x() < _min.x() ) 
	{
		_min.x() = v.x();
	}
	if ( v.y() < _min.y ()) 
	{
		_min.y() = v.y();
	}
	if ( v.z() < _min.z() )
	{
		_min.z() = v.z();
	}
	if ( v.x() > _max.x() )
	{
		_max.x() = v.x();
	}
	if ( v.y() > _max.y() )
	{
		_max.y() = v.y();
	}
	if ( v.z() > _max.z() )
	{
		_max.z() = v.z();
	}
}


BoundingBox GLbKdTree::computeTightFittingBoundingBox(unsigned int num_tris,unsigned int *tri_indices)
{
	// Compute bounding box for input mesh.
	Vec3 _max = Vec3( -INFINITYM, -INFINITYM, -INFINITYM );
	Vec3 _min = Vec3( INFINITYM, INFINITYM, INFINITYM );

	for ( unsigned int i = 0; i < num_tris; ++i )
	{
		unsigned triIndice = tri_indices[i];

		Triangle&tri = meshTriangles[triIndice];

		expandBoundBox(tri.vertex[0], _max,_min);
		expandBoundBox(tri.vertex[1], _max,_min);
		expandBoundBox(tri.vertex[2], _max,_min);
	}

	BoundingBox bbox;
	bbox._min = _min;
	bbox._max = _max;

	return bbox;
}

 





