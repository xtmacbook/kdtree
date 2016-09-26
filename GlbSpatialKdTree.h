/********************************************************************
* Copyright (c) 2014 北京超维创想信息技术有限公司
* All rights reserved.
* @file    GlbSpatialKdTree.h
* @brief 该kdtree树结构 使用 surface area heuristic (SAH)来定位   
* splitting plane的位置
*	详细参考:http://dcgi.felk.cvut.cz/home/havran/publications.html
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

/*
	ad hoc termination criteria
	是树结构构建过程中 ，终止构建的阀值
	dMax 和 Nmax
*/

#define TERMINATION_CRITERIA_D_MAX 16  //最深树节点(16 + 2)
#define TERMINATION_CRITERIA_N_MAX 100   //叶节点

//#define KDTREE_NEIGHBORLINKS        //是否使用neighborLinks遍历树
#define KDTREE_SAH_CONSTRUCT        //构建SAH树


namespace GlbGlobe
{
	struct BoxEdge;

	typedef osg::Vec3d Vec3;
	typedef osg::BoundingBoxd BoundingBox;

	typedef std::vector<BoxEdge > v_BoxEdge;
	typedef std::vector<v_BoxEdge> vv_BoxEdge;
	typedef std::vector<BoxEdge* > vp_BoxEdge;
	typedef std::vector<vp_BoxEdge > vvp_BoxEdge;

	//线段上点类型，起点 、终点
	enum SegmentPT
	{
		START = 0, END
	};

	//splitting plane是哪个坐标轴
	enum Axes
	{
		X_axis = 0, Y_axis, Z_axis, No_axis
	};

	//节点AABB的6个面
	enum Faces
	{
		FLeft = 0x00, FRight, FFront, FBack, FBottom, FTop
	};

	struct Triangle
	{
	public:
		Triangle(){};
		// Constructor with arguments
		Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);

		bool getRayIntersection(const Vec3&origin,const Vec3&dir,Vec3&intersectionPoint);

		BoundingBox bound;

		Vec3 vertex[3];
	};
	
	//用在中点构建
	struct TriangleM
	{
	public:
		TriangleM(){};
		// Constructor with arguments
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

		unsigned int triangleIndex; //对应的三角形索引

		double splitPlanePosition; //splitting plane的位置 

		friend std::ostream &operator<<(std::ostream &o, const BoxEdge &edge);
		friend std::ostream &operator<<(std::ostream &o, const BoxEdge *edge);
	};

	struct KDTNode
	{

		KDTNode(void);

		~KDTNode(void);

		//bool is_root()const;
		inline bool is_leaf()const {return (left == NULL) && (right == NULL);}
		const KDTNode* backtrack_leaf(const Vec3 &point)const;
		const KDTNode * find_leaf(const Vec3 &point)const ;
		
		struct KDTNode * left;  //子节点 child
		struct KDTNode * right;

#ifdef KDTREE_NEIGHBORLINKS
		struct KDTNode* ropes[6]; ///* 6个平面的链接 */
#else
		struct KDTNode * parent;
#endif // KDTREE_NEIGHBORLINKS

		BoxEdge *splitEdge; ////splitting plane

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

		struct KDTNodeM * left;  //子节点 child
		struct KDTNodeM * right;

#ifdef KDTREE_NEIGHBORLINKS
		struct KDTNode* ropes[6]; ///* 6个平面的链接 */
#else
		struct KDTNodeM * parent;
#endif // KDTREE_NEIGHBORLINKS

		BoxEdge *splitEdge; ////splitting plane

		unsigned int num_tris;     //内部三角形
		unsigned int *tri_indices;  //三角形索引

		//test
		unsigned int leve;

		BoundingBox box;  //box 
	};

	class GLbKdTree
	{
	public:

		GLbKdTree(bool constructType);

		~GLbKdTree(void);

		unsigned int  GetMeshTriangleAndVertexs(const osg::Drawable*geometry);

		bool ConstructKdTree(unsigned int num_tris,const BoundingBox&bounds);

		/* 遍历   */
		bool GetIntersectPoint(const Vec3&origin,const Vec3&dir,Vec3&intersectionP);
	protected:

		/*
			空间中点划分 (median space)
		*/
		KDTNodeM* ConstructTreeMedianSpaceSplit(unsigned int num_tris,const BoundingBox& bounds);
		/*
			SAH 构造
		*/
		KDTNode* ConstructTreeSAHSplit(unsigned int num_tris,const BoundingBox& bounds);	

		
		//////////////////////////////////////////////////////////////////////////
		//neighbor links

		/* 
		node :节点
		face:该节点的一个面
		rootNode:顶级节点
		return :该节点的面对应的下一个节点
		*/
		KDTNode* FindSingleNeighborLink(KDTNode *node, Faces face, KDTNode *rootNode);

		KDTNode* CreateNeighborLinksTree(KDTNode *node, KDTNode *subtree, Faces face);

		void BuildNeighborLinksTree(KDTNode *node, Faces face);
	
		void BuildRopeStructure();	

	private:

		/*
			构建TighFitting AABBox
		*/
		BoundingBox computeTightFittingBoundingBox( unsigned int num_tris,unsigned int *tri_indices );
		void expandBoundBox(const Vec3&v,Vec3&max,Vec3&min);

		KDTNodeM* ConstructTreeMedianSpaceSplit(unsigned int num_tris,unsigned int *tri_indices,
			const BoundingBox& bounds,unsigned int curr_depth,KDTNodeM*parent);

		KDTNode * buildTree_boxEdges(const BoundingBox& nodeExtent, vv_BoxEdge& boxEdgeList,
				int maxDepth,KDTNode*parent);

		//traveral
		bool GetIntersecting(const KDTNode*node, const Vec3&origin ,const Vec3&dir,Vec3&intersectionP);
		bool GetIntersectingM(const KDTNodeM*node, const Vec3&origin ,const Vec3&dir,Vec3&intersectionP);
	private:

		bool     sahUse; //true 中点 false sah

		KDTNode * treeRoot;
		KDTNodeM *treeRootM;

		Triangle*     meshTriangles;
		unsigned int triangleSize;
	
#ifdef KDTREE_SAH_CONSTRUCT
		 const GlbSAH sah;
#endif

	};
}


#endif