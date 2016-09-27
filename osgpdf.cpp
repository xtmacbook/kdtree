/************************************************************************/
/* test h grid                                                                      */
/************************************************************************/


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
#include <osg/Point>

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
#include "comm/KeyboardHandler.h"

/************************************************************************/

osg::MatrixTransform * m0;
osg::MatrixTransform * m1;
void testkey()
{
    osg::ComputeBoundsVisitor cbVisitor;
    m1->accept(cbVisitor);
    osg::BoundingBox &bb = cbVisitor.getBoundingBox();

}

int mainx()
{

    osg::ref_ptr<osg::Group> root = new osg::Group;
    m0 = new osg::MatrixTransform;
    m1 = new osg::MatrixTransform;

    osg::Matrix mt0 = osg::Matrix::translate(osg::Vec3(2,3,4));
    osg::Matrix mt1 = osg::Matrix::scale(osg::Vec3(0.5,0.5,0.5));

    m0->setMatrix(mt0);
    m1->setMatrix(mt1);


    osg::Node * cow = osgDB::readNodeFile("/Users/glp/Documents/osgResource/OpenSceneGraph-Data-3.0.0/cow.osg");

    m0->addChild(m1);
    m1->addChild(cow);

    root->addChild(m0);

    osgViewer::Viewer viewer;
    viewer.setUpViewInWindow(10, 10, 800, 600);
    viewer.addEventHandler(new osgViewer::StatsHandler);
    viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet() ));
    keyboardEventHandler * kf = new keyboardEventHandler;
    kf->addFunction('b',testkey);
//    kf->addFunction('r',romveCOw);
    viewer.addEventHandler(kf);
    
    viewer.setSceneData(root);
    
    return viewer.run();
}


///////////////////////////////////////


#include "kdtree/kdtree/GlbSpatialKdTree.h"
#include "comm/Geometry.h"


osg::Vec3d SStart(0.0148538, -0.0633929, 0.0756582);
osg::Vec3d Send( -0.00060602 ,0.0688803 ,0.0783454);

class PickHandler : public osgGA::GUIEventHandler
{
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

         //std::cout << "start: " << vStart.x() << " "  << vStart.y() << " " << vStart.z()<< std::endl;
        // std::cout << "End: " << vEnd.x() << " "  << vEnd.y() << " " << vEnd.z()<< std::endl;




		  osg::Vec3 osgInterP;
		  osg::Vec3d kdTreeInterP;
#if 1
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
		 // std::cout<<"OSG Completed in: "<<osg::Timer::instance()->delta_s(t1,t2)<<std::endl;
		  std::cout << "OSG Point: " << osgInterP.x() << " " << osgInterP.y() << " " << osgInterP.z() << std::endl;
#endif
		  t2 = timer.tick();
         osg::Vec3d dir = (vEnd - vStart);
         dir.normalize();

        if( kdtree.RayTracer(GlbGlobe::Ray(vStart,dir),kdTreeInterP))
        {

        }
    t3 = timer.tick();
    //std::cout<<"kdtree Completed in: "<<osg::Timer::instance()->delta_s(t2,t3)<<std::endl;
    std::cout << "kdtree Point: " << kdTreeInterP.x() << " " << kdTreeInterP.y() << " " << kdTreeInterP.z() << std::endl;
    }

    float _mx, _my;
    osg::Node * _node;
    GlbGlobe::GLbKdTree& kdtree;

};



osg::Node * hubby;
GlbGlobe::GLbKdTree * kdtreeTest;
void testKdtree()
{
    osg::Vec3 osgInterP;
    osg::Vec3d kdTreeInterP;

    osg::ref_ptr<osgUtil::LineSegmentIntersector> intersector =
        new osgUtil::LineSegmentIntersector(SStart, Send);
		  osgUtil::IntersectionVisitor intersectVisitor( intersector.get());


		  hubby->accept(intersectVisitor);

		  if(intersector->containsIntersections())
            {
              osgInterP = intersector->getIntersections().begin()->getWorldIntersectPoint();
              std::cout << "OSG Point: " << osgInterP.x() << " " << osgInterP.y() << " " << osgInterP.z() << std::endl;
            }


    osg::Vec3d dir = (Send - SStart);
    dir.normalize();
    if( kdtreeTest->RayTracer(GlbGlobe::Ray(SStart,dir),kdTreeInterP))
    {
        std::cout << "kdtree Point: " << kdTreeInterP.x() << " " << kdTreeInterP.y() << " " << kdTreeInterP.z() << std::endl;
        //0.0137945 -0.0543293 0.0758423
    }



}

int main()
{

    osgViewer::Viewer viewer;
    viewer.setUpViewInWindow(10, 10, 800, 600);
    viewer.addEventHandler(new osgViewer::StatsHandler);

    osg::ref_ptr<osg::Group> root = new osg::Group;

#if 1
    std::string meshFileName = "/Users/glp/Documents/osgResource/OpenSceneGraph-Data-3.0.0/hubby.osg";
    //std::string meshFileName = "/Users/glp/Downloads/teapot/teapot.obj";
    osg::ref_ptr<osg::Node> sponaza = osgDB::readNodeFile(meshFileName);
    hubby = sponaza;
    osg::Drawable * drawable = sponaza->asGroup()->getChild(0)->asGeode()->getDrawable(0);
#else
    osg::Vec3 center (0.0,0.0,0.0);
    osg::Vec3 radio(2.0,4.0,5.0);
    osg::ref_ptr<osg::Node> sponaza = openGLBox(center,radio);
    osg::Drawable * drawable = sponaza->asGeode()->getDrawable(0);
#endif
    GlbGlobe::GLbKdTree kdtree(false,true);
    kdtreeTest =  &kdtree;

    unsigned int num = kdtree.GetMeshTriangleAndVertexs(drawable);
    
    kdtree.ConstructKdTree(num,drawable->getBound());
    
    root->addChild(sponaza);
    
    viewer.setSceneData(root);

#if 1
    viewer.addEventHandler(new PickHandler(root,kdtree));
#else
    keyboardEventHandler * kf = new keyboardEventHandler;
    kf->addFunction('a',testKdtree );
    viewer.addEventHandler(kf);
#endif
    viewer.run();
    
    return 0;
}

    //task_scheduler_init init(task_scheduler_init::deferred);


