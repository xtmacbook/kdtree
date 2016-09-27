

#include <limits>
#include <osg/TriangleFunctor>
#include <osg/Geode>
#include <osg/Group>
#include "GlbSpatialKdTree.h"

#include "../comm/xMath.h"

using namespace GlbGlobe;


BoxEdge::BoxEdge()
{
	splitPlanePosition = KDTREEDOUBLEINFINITYM;
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

bool Triangle::getRayIntersection(const Vec3&origin,const Vec3&dir,Vec3&intersectionPoint)const
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
#ifndef KDTREE_NEIGHBORLINKS

	if (box.contains(point,KD_TREE_EPSILON))
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
#else
	return NULL;
#endif
}

const bool KDTNodeM::isPointToLeftOfSplittingPlane(const Vec3&p)const
{
	if ( splitEdge->axis == X_axis )
	{
		return ( p.x() < splitEdge->splitPlanePosition );
	}
	else if ( splitEdge->axis == Y_axis )
	{
		return ( p.y() < splitEdge->splitPlanePosition );
	}
	else if ( splitEdge->axis == Z_axis )
	{
		return ( p.z() < splitEdge->splitPlanePosition );
	}
	else
	{
		return false;
	}
}

KDTNodeM* KDTNodeM::getNeighboringNode(const Vec3&p)const
{
#ifdef KDTREE_NEIGHBORLINKS
	// Check left face.
	if ( fabs( p.x() - box._min.x() ) < KD_TREE_EPSILON )
	{
		return ropes[FLeft];
	}
	// Check front face.
	else if ( fabs( p.z() - box._max.z() ) < KD_TREE_EPSILON )
	{
		return ropes[FFront];
	}
	// Check right face.
	else if ( fabs( p.x() - box._max.x() ) < KD_TREE_EPSILON )
	{
		return ropes[FRight];
	}
	// Check back face.
	else if ( fabs( p.z() - box._min.z() ) < KD_TREE_EPSILON )
	{
		return ropes[FBack];
	}
	// Check top face.
	else if ( fabs( p.y() - box._max.y() ) < KD_TREE_EPSILON )
	{
		return ropes[FTop];
	}
	// Check bottom face.
	else if ( fabs( p.y() - box._min.y() ) < KD_TREE_EPSILON )
	{
		return ropes[FBottom];
	}
	// p should be a point on one of the faces of this node's bounding box, but in this case, it isn't.
	else
	{
		// std::cout << "ERROR: Node neighbor could not be returned." << std::endl;
		return NULL;
	}
#else
	return NULL;
#endif
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

GLbKdTree::GLbKdTree(bool t,bool r):sahUse(t),rope(r)
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

const Triangle * GlbGlobe::GLbKdTree::getMeshTriangles(void)const
{
	return meshTriangles;
}

unsigned int GLbKdTree::GetMeshTriangleAndVertexs(const osg::Node* mesh)
{
	if(!mesh) return 0;

	if(mesh->asGroup())
	{
		const osg::Group * meshGroup = mesh->asGroup();

		for(unsigned int i = 0;i <meshGroup->getNumChildren();i++)
		{
			if(meshGroup->getChild(i)->asGeode())
			{
				const osg::Geode::DrawableList& list =  meshGroup->getChild(i)->asGeode()->getDrawableList();

				for(unsigned int j = 0;j < list.size();j++)
				{
					GetMeshTriangleAndVertexsFromDrawable(list.at(j).get());
				}
			}
		}
	}

	meshTriangles = new Triangle[triangles.size()];

	std::vector<Triangle>::iterator iter = triangles.begin();
	unsigned int triIndice = 0;
	for(;iter != triangles.end();iter++)
	{
		meshTriangles[triIndice++] = *iter;
	}

	triangleSize = triangles.size();

	triangles.clear();

	return triangleSize;
}

struct TriangleIntersector
{


	TriangleIntersector()
	{

	}

	inline void operator () (const osg::Vec3d& v1,const osg::Vec3d& v2,const osg::Vec3d& v3, bool treatVertexDataAsTemporary)
	{
		triangles.push_back(Triangle(v1,v2,v3));
	}

};

unsigned int  GLbKdTree::GetMeshTriangleAndVertexsFromDrawable(const osg::Drawable*drawable)
{
	if(!drawable) return false;

	//unsigned int drawableNum = geometry->getdraw
	osg::TriangleFunctor<TriangleIntersector> ti;

	drawable->accept(ti);

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
	Vec3 _max = Vec3( -KDTREEDOUBLEINFINITYM, -KDTREEDOUBLEINFINITYM, -KDTREEDOUBLEINFINITYM );
	Vec3 _min = Vec3( KDTREEDOUBLEINFINITYM, KDTREEDOUBLEINFINITYM, KDTREEDOUBLEINFINITYM );


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







