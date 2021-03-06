/********************************************************************
* Copyright (c) 2014 北京超维创想信息技术有限公司
* All rights reserved.
* @file    GlbSpatialKdTree.h
* @brief 该kdtree树结构 使用 surface area heuristic (SAH)来定位   
* splitting plane的位置
*	http://dcgi.felk.cvut.cz/home/havran/publications.html
*   2000 14. V. Havran: "Heuristic Ray Shooting Algorithms", Ph.D. Thesis, November 2000
*   主要用到 Geometric probability (空间几何概率)
* @version 1.0
* @author  xt
* @date    2016-09-12 11:30
*********************************************************************/
#pragma once

#ifndef _GLB_SPATIAL_KD_TREE_
#define _GLB_SPATIAL_KD_TREE_

#include <osg/Vec3d>
#include <osg/BoundingBox>
#include <osg/Geometry>
#include <vector>

#include "GlbSAH.h"

#define TERMINATION_CRITERIA_D_MAX 16
#define TERMINATION_CRITERIA_N_MAX 100

#define KDTREE_NEIGHBORLINKS
#define KDTREE_SAH_CONSTRUCT

const double   KDTREEDOUBLEINFINITYM			= std::numeric_limits<double>::max();
const unsigned int KDTREEINTINFINITY			= std::numeric_limits<unsigned int >::max();
const double   KD_TREE_EPSILON					= 0.00001;
const double   RAY_EPSILON_OFFSET				= 0.0001;

namespace GlbGlobe
{
	struct BoxEdge;
	class GLbKdTree;

	typedef osg::Vec3d Vec3;
	typedef osg::BoundingBoxd BoundingBox;

	typedef std::vector<BoxEdge > v_BoxEdge;
	typedef std::vector<v_BoxEdge> vv_BoxEdge;
	typedef std::vector<BoxEdge* > vp_BoxEdge;
	typedef std::vector<vp_BoxEdge > vvp_BoxEdge;

	enum SegmentPT
	{
		START = 0, END
	};

	//splitting plane Axes
	enum Axes
	{
		X_axis = 0, Y_axis, Z_axis
	};

	//splitting plane Faces
	enum Faces
	{
		FLeft = 0x00, FRight, FFront, FBack, FBottom, FTop
	};

	struct Triangle
	{
	public:
		Triangle(){};

		Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);

		bool getRayIntersection(const Vec3&origin,const Vec3&dir,Vec3&intersectionPoint)const;

		BoundingBox bound;

		Vec3 vertex[3];
	};

	struct TriangleM
	{
	public:
		TriangleM(){};

		TriangleM(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);

		Vec3 vertex[3];
	};

	struct BoxEdge
	{
	public:

		BoxEdge();

		BoxEdge(double f, unsigned int t, SegmentPT edgeType, Axes axis) :
		splitPlanePosition(f), triangleIndex(t), edgeType(edgeType), axis(axis) {}

		bool operator<(const BoxEdge &RHS) const
		{
			if (this->splitPlanePosition == RHS.splitPlanePosition) 
			{
				if (this->triangleIndex != RHS.triangleIndex) 
				{
					return this->triangleIndex < RHS.triangleIndex;
				} else 
				{  
					return (int)(this->edgeType) < (int)(RHS.edgeType);
				}
			} else {
				return this->splitPlanePosition < RHS.splitPlanePosition;
			}
		}

		SegmentPT edgeType;
		Axes axis;

		unsigned int triangleIndex;

		double splitPlanePosition;

		friend std::ostream &operator<<(std::ostream &o, const BoxEdge &edge);
		friend std::ostream &operator<<(std::ostream &o, const BoxEdge *edge);
	};

	struct Ray
	{
		Ray(const Ray&r);

		Ray(Vec3 orig, Vec3 dir):origin(orig),direction(dir)
		{ }

		Vec3 origin;
		Vec3 direction;

	};

	struct KDTNode
	{

		KDTNode(void);

		~KDTNode(void);

		//bool is_root()const;
		inline bool is_leaf()const {return (left == NULL);}

		bool treeLeafTrace(const Ray ray,Vec3&intersectionP,const GLbKdTree*tree);

		struct KDTNode * left;
		struct KDTNode * right;

#ifdef KDTREE_NEIGHBORLINKS
		struct KDTNode* ropes[6];
#endif // KDTREE_NEIGHBORLINKS

		BoxEdge *splitEdge; //splitting plane

		std::vector<int> *triangleIndices;

		BoundingBox box;  //box


	};

	struct KDTNodeM
	{

		KDTNodeM(void);

		~KDTNodeM(void);

		//bool is_root()const;
		inline bool is_leaf()const {return left == NULL;}

		const KDTNodeM* backtrack_leaf(const Vec3 &point)const;
		const KDTNodeM * find_leaf(const Vec3 &point)const ;
		const bool isPointToLeftOfSplittingPlane(const Vec3&point)const;

		KDTNodeM * getNeighboringNode(const Vec3&exitP)const;

		bool treeLeafTrace(const Ray ray,Vec3&intersectionP,const GLbKdTree*localKdTreePtr)const;

		struct KDTNodeM * left;
		struct KDTNodeM * right;

#ifdef KDTREE_NEIGHBORLINKS
		struct KDTNodeM* ropes[6];
#else
		struct KDTNodeM * parent;
#endif // KDTREE_NEIGHBORLINKS

		BoxEdge *splitEdge; ////splitting plane

		unsigned int num_tris;
		unsigned int *tri_indices;

		BoundingBox box;  //box 
	};

	template <typename T>
	struct StackElem
	{
		StackElem(){}
		StackElem( const T* n,double aa,double bb):
		node(n),a(aa),b(bb)
		{}
		const T * node;
		double a;
		double b;
	};

	template <typename T>
	struct StackElemA
	{
		StackElemA(){}

		const T* node;//指向子节点
		double t; //  进 /出 有效距离
		Vec3 pb;  // 进 / 出 点;
		int prev; // the pointer to the pre stack item
	};

	class GLbKdTree
	{
	public:

		GLbKdTree(bool constructType,bool r);

		~GLbKdTree(void);

		unsigned int GetMeshTriangleAndVertexs(const osg::Node* mesh);

		unsigned int  GetMeshTriangleAndVertexsFromDrawable(const osg::Drawable*drawable);

		bool ConstructKdTree(unsigned int num_tris,const BoundingBox&bounds);

		/* ray tracer */
		bool RayTracer(const Ray&ray,Vec3&intersectionP, int  t);

		const Triangle * getMeshTriangles(void)const;
	protected:

	protected:
		/*
		(median space)
		*/
		KDTNodeM* ConstructTreeMedianSpaceSplit(unsigned int num_tris,const BoundingBox& bounds);
		/*
		SAH
		*/
		KDTNode* ConstructTreeSAHSplit(unsigned int num_tris,const BoundingBox& bounds);	

		//////////////////////////////////////////////////////////////////////////
		//neighbor links

		/* 

		*/
		KDTNode* FindSingleNeighborLink(KDTNode *node, Faces face, KDTNode *rootNode);

		KDTNode* CreateNeighborLinksTree(KDTNode *node, KDTNode *subtree, Faces face);

		void BuildNeighborLinksTree(KDTNode *node, Faces face);

		void BuildRopeStructure();	

	private:

		/*
		TighFitting AABBox
		*/
		BoundingBox computeTightFittingBoundingBox( unsigned int num_tris,unsigned int *tri_indices );
		void expandBoundBox(const Vec3&v,Vec3&max,Vec3&min);

		KDTNodeM* constructTreeMedianSpaceSplit(unsigned int num_tris,unsigned int *tri_indices,
			const BoundingBox& bounds,unsigned int curr_depth,KDTNodeM*parent);

		KDTNode * buildTree_boxEdges(const BoundingBox& nodeExtent, vv_BoxEdge& boxEdgeList,
			int maxDepth);
		
	private:

		bool     sahUse; //true: sah false :median space
		bool     rope;
		KDTNode  *          treeRoot;
		KDTNodeM *          treeRootM;

		Triangle*           meshTriangles;
		unsigned int        triangleSize;

#ifdef KDTREE_SAH_CONSTRUCT
		const GlbSAH       sah;
#endif

	};
}


#endif