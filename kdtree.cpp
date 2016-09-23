
#include <osg/Material>
#include <osg/Geode>
#include <osg/BlendFunc>
#include <osg/Depth>
#include <osg/PolygonOffset>
#include <osg/MatrixTransform>
#include <osg/Camera>
#include <osg/RenderInfo>
#include <osg/ShapeDrawable>  
#include <osg/ComputeBoundsVisitor>  
#include <osg/BoundingBox>  
#include <osg/BoundingSphere>  
#include <osg/AnimationPath> 
#include <osg/PositionAttitudeTransform>
#include <osg/VertexProgram>
#include <osg/FragmentProgram>
#include <osg/Plane>

#include <OpenThreads/Mutex>

#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osgViewer/CompositeViewer>
#include <osgViewer/ViewerEventHandlers>
#include <osgWidget/PdfReader>

#include <osgDB/ReadFile>
#include <osgDB/WriteFile>


#include <osgGA/NodeTrackerManipulator>
#include <osgGA/TrackballManipulator>
#include <osgGA/KeySwitchMatrixManipulator>
#include <osgGA/TerrainManipulator>
#include <osgGA/StateSetManipulator>

#include <osgUtil/Optimizer>
#include <osgUtil/CullVisitor>
#include <osgUtil/SmoothingVisitor>
#include <osgUtil/PolytopeIntersector>
#include <osgUtil/LineSegmentIntersector>

#include "../kdtree/GlbSpatialKdTree.h"
#include "../comm/Geometry.h"

class PickHandler : public osgGA::GUIEventHandler {
public:

	PickHandler(osg::Node * n,GlbGlobe::GLbKdTree& tree):
	  _mx(0.0f),
		  _my(0.0f),
	  _node(n),
	  kdtree(tree){}

	  ~PickHandler() {}

	  bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa)
	  {
		  osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
		  if (!view) return false;

		  switch(ea.getEventType())
		  {
		  case(osgGA::GUIEventAdapter::PUSH):
			  {
				  _mx = ea.getX();
				  _my = ea.getY();
				  break;
			  }
		  case(osgGA::GUIEventAdapter::RELEASE):
			  {
				  if (_mx==ea.getX() && _my==ea.getY())
				  {
					  pick(view, ea);
				  }
				  break;
			  }
		  default:
			  break;
		  }
		  return false;
	  }

	  void pick(osgViewer::View* view, const osgGA::GUIEventAdapter& event)
	  {

		  double _x = _mx;
		  double _y = _my;

		  _y = view->getCamera()->getViewport()->height() - _y;

		  osg::Vec3d vStart(_x,_y,0.0);
		  osg::Vec3d vEnd(_x,_y,1.0);	


		  osg::Camera* p_camera = view->getCamera();
		  osg::Matrixd VPW = p_camera->getViewMatrix() *
			  p_camera->getProjectionMatrix() *
			  p_camera->getViewport()->computeWindowMatrix();
		  osg::Matrixd inverseVPW;
		  inverseVPW.invert(VPW);
		  vStart = vStart * inverseVPW;
		  vEnd = vEnd * inverseVPW;


		  osg::Timer timer;
		  osg::Timer_t t1=0,t2=0,t3 = 0;


		  osg::Vec3 osgInterP;
		  osg::Vec3d kdTreeInterP;
#if 0
		  osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector = new osgUtil::LineSegmentIntersector(vStart, vEnd);
		  osgUtil::IntersectionVisitor intersectVisitor( intersector.get());

		  t1 = timer.tick();
		  _node->accept(intersectVisitor);

		  if(intersector->containsIntersections())
		  {
			  osgInterP = intersector->getIntersections().begin()->getWorldIntersectPoint();
		  }
		  t2 = timer.tick();

		  //osgTime = osg::Timer::instance()->delta_s(t1,t2);
		  std::cout<<"OSG Completed in: "<<osg::Timer::instance()->delta_s(t1,t2)<<std::endl;
		  std::cout << "OSG Point: " << osgInterP.x() << " " << osgInterP.y() << " " << osgInterP.z() << std::endl;
#endif
		  t2 = timer.tick();

		 if( kdtree.GetIntersectPoint(vStart,vEnd - vStart,kdTreeInterP))
		 {

		 }
		 t3 = timer.tick();
		 std::cout<<"kdtree Completed in: "<<osg::Timer::instance()->delta_s(t2,t3)<<std::endl;
		 std::cout << "kdtree Point: " << kdTreeInterP.x() << " " << kdTreeInterP.y() << " " << kdTreeInterP.z() << std::endl;
	  }

	  float _mx, _my;
	  osg::Node * _node;
	  GlbGlobe::GLbKdTree& kdtree;

};


int main()
{

	osgViewer::Viewer viewer;
	viewer.addEventHandler(new osgViewer::StatsHandler);

	osg::ref_ptr<osg::Group> root = new osg::Group;

#if 1
	std::string meshFileName = "D:/OpenSceneGraph-Data/objs/bunny.osg";
	osg::ref_ptr<osg::Node> sponaza = osgDB::readNodeFile(meshFileName);
	osg::Drawable * drawable = sponaza->asGroup()->getChild(0)->asGeode()->getDrawable(0);
#else
	osg::ref_ptr<osg::Node> sponaza = openGLBox(osg::Vec3(0.0,0.0,0.0),osg::Vec3(2,4,5));
	osg::Drawable * drawable = sponaza->asGeode()->getDrawable(0);
#endif
	GlbGlobe::GLbKdTree kdtree(false);

	unsigned int num = kdtree.GetMeshTriangleAndVertexs(drawable);

	kdtree.ConstructKdTree(num,drawable->getBound());

	root->addChild(sponaza);

	viewer.setSceneData(root);

	viewer.addEventHandler(new PickHandler(root,kdtree));

	viewer.run();

	return 0;
}

//task_scheduler_init init(task_scheduler_init::deferred);

 