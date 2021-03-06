/* 
 * heliangliang (heliangliang.bj@gmail.com), USTB, 2014.05.14, CopyRight Reserved
 *
 * The SPH algrithm is implemented based on these papers:
 * 1. Smoothed particle hydrodynamics (Monaghan,J.J. 1992)			[Mon92]
 * 2. Smoothed particles: a new paradigm for animating highly deformable
 *	  bodies (1996)													[DC96]
 * 3. Weakly compressible SPH for free surface flows (2007)			[BT07]
 * 4. Predictive-Corrective Incompressible SPH (2009)				[SP09]
 * 5. Density contrast SPH interfaces (2008)						[SP08]
 * 6. Versatile Rigid-Fluid Coupling for Incompressible SPH (2012)	[AIS*12]
 * 7. Versatile surface tension and adhesion for SPH fluids (2013)	[AAIT13]
 * 8. SPH Fluids in Computer Graphics (2014)						[IOS*14]
 *
*/
#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1

#include "SphWithBullet.h"

#include "BulletCollision/NarrowPhaseCollision/btMinkowskiPenetrationDepthSolver.h"
#include "BulletCollision/CollisionDispatch/btConvex2dConvex2dAlgorithm.h"


// initialize bullet, should be called in the sph constructor
void SphWithBullet::SphBulletWorld::initialize(const btVector3& gravity)
{
	collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();
	dispatcher = new btCollisionDispatcher(collisionConfiguration);
	dispatcher->setDispatcherFlags(
		dispatcher->getDispatcherFlags()|btCollisionDispatcher::CD_STATIC_STATIC_REPORTED );
	btVector3 worldAabbMin(-1000,-1000,-1000), worldAabbMax(1000,1000,1000);
	int maxProxies = 32766;
	broadphase = new btAxisSweep3(worldAabbMin,worldAabbMax,maxProxies);
	constraintSolver = new btSequentialImpulseConstraintSolver;
	btSoftBodySolver* softBodySolver = 0;
	dynamicsWorld = new btSoftRigidDynamicsWorld(
		dispatcher,broadphase,constraintSolver,collisionConfiguration,softBodySolver);
	dynamicsWorld->setGravity(gravity);

	softBodyWorldInfo = new btSoftBodyWorldInfo();
	softBodyWorldInfo->m_dispatcher = dispatcher;
	softBodyWorldInfo->m_broadphase = broadphase;
	softBodyWorldInfo->m_sparsesdf.Initialize();

	// for 2D
	if( vec_t::dim==2 ){
		btVoronoiSimplexSolver* simplex = new btVoronoiSimplexSolver();
		btMinkowskiPenetrationDepthSolver* pdSolver = new btMinkowskiPenetrationDepthSolver();
		btConvex2dConvex2dAlgorithm::CreateFunc* convexAlgo2d =
			new btConvex2dConvex2dAlgorithm::CreateFunc(simplex,pdSolver);
		dispatcher->registerCollisionCreateFunc(
			CONVEX_2D_SHAPE_PROXYTYPE,CONVEX_2D_SHAPE_PROXYTYPE,convexAlgo2d);
	}

}

// bullet clear, i.e. delete bullet objects, should be called in the sph destructor
void SphWithBullet::SphBulletWorld::destroy()
{
	//remove the btCollisionObject from the dynamics world and delete them
	for(int i=dynamicsWorld->getNumCollisionObjects()-1; i>=0 ; --i) { // must be --i, removeCollisionObject
		btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
			delete body->getMotionState();
		dynamicsWorld->removeCollisionObject( obj );
		delete obj;
	}
	//dynamicsWorld->getCollisionObjectArray().clear(); // removeCollisionObject does the same thing
	//delete collision shapes
	for(size_t i=0; i<collisionShapes.size(); ++i) {
		delete collisionShapes[i]; collisionShapes[i]=0;
	}
	collisionShapes.clear();
	// delete dynamicsworld etc.
	delete dynamicsWorld; dynamicsWorld=0;
	delete constraintSolver; constraintSolver=0;
	delete broadphase; broadphase=0;
	delete dispatcher; dispatcher=0;
	delete collisionConfiguration; collisionConfiguration=0;
	delete softBodyWorldInfo; softBodyWorldInfo=0;

}

// for btCompoundShape, push_back its children(can be compound!) to collisionShapes
void SphWithBullet::SphBulletWorld::pushBackShape(const btCollisionShape* shape)
{
	if( std::find(collisionShapes.begin(),collisionShapes.end(),shape)
	==collisionShapes.end() ){
		collisionShapes.push_back(shape);
	}
	if( shape->isCompound() ){
		btCompoundShape* compound = (btCompoundShape*)shape;
		for(int i=0; i<compound->getNumChildShapes(); ++i){
			pushBackShape( compound->getChildShape(i) );
		}
	}
}

// add solid, set the data of mbt_World, mbt_Solids, m_Solids
void SphWithBullet::addSolid( const Solid& sphsolid, const MBtSolid& mbtsolid )
{
	// btDynamicsWorld
	if(sphsolid.type==Solid::RIGIDBODY){ // rigidbody
		// mbt_World, mbt_Solids
		btRigidBody* prigid = btRigidBody::upcast(mbtsolid.btObject); assert(prigid);
		mbt_World.pushBackShape( prigid->getCollisionShape() );
		mbt_World.dynamicsWorld->addRigidBody( prigid );

	}else if(sphsolid.type==Solid::SOFTBODY){ // softbody
		// mbt_World, mbt_Solids
		btSoftBody* psoft = btSoftBody::upcast(mbtsolid.btObject); assert(psoft);
		psoft->m_cfg.collisions |= btSoftBody::fCollision::VF_SS;
		// need not delete psoft->getCollisionShape(), the deleting causes error
		mbt_World.dynamicsWorld->addSoftBody( psoft );
		
	}else{
		//...
	}

	// m_Solids, mbt_Solids
	m_Solids.push_back( sphsolid ); // not the best efficiency, but ok
	mbt_Solids.push_back( mbtsolid );

	assert( m_Solids.size()==mbt_Solids.size() );
	//m_Solids.back().mbtObject = &(mbt_Solids.back()); // so parameter sphsolid.mbtObject is ignored
	// mbt_Solids may move entirety beacause of push_back !!!!!!!!!!!
	for(size_t i=0; i<m_Solids.size(); ++i){
		m_Solids[i].mbtSolid_ptr = &mbt_Solids[i];
	}
	m_Solids.back().updatePaticles();
	m_Solids.back().adjustSpace(m_TH.spaceMin, m_TH.spaceMax); // modify m_SpaceMin, m_SpaceMax

	prepareSolidPartWeight();

}


// update solids using bullet,
// synchronize boundary particles with rigid's transform or softbody's nodes
void SphWithBullet::updateSolids()
{
	assert(m_Solids.size()==mbt_Solids.size());

	// compute total force and torque from boundaryParticles' forces
	// for rigids, send total force and torque to bullet
	// for softbodys, send boundaryParticles' forces to nodes
	for(int n_r=int(m_Solids.size()),k=0; k<n_r; ++k) 
		if(m_Solids[k].dynamic){
		if(m_Solids[k].type==Solid::RIGIDBODY){ // rigidbody
			btRigidBody* robj = btRigidBody::upcast(mbt_Solids[k].btObject); assert(robj);
			vec_t f; vec3_t t;
			btVector3 o = robj->getCenterOfMassPosition();
			m_Solids[k].rigidForceJoin(btVec3_to_vect(o), f, t);
			if( !(f==vec_t::O && t==vec3_t::O) ){
				robj->activate(true);
				robj->applyCentralForce(vect_to_btVect3(f));
				robj->applyTorque(vec3t_to_btVect3(t));
			}
		}else if(m_Solids[k].type==Solid::SOFTBODY){ // softbody
			btSoftBody* sobj = btSoftBody::upcast(mbt_Solids[k].btObject); assert(sobj);
			std::vector<Solid::Part>& s_ps = m_Solids[k].boundaryParticles;
			int num = int(s_ps.size()); assert(num==sobj->m_nodes.size());
		#pragma omp parallel for
			for(int i=0; i<num; ++i){
				sobj->addForce(vect_to_btVect3(s_ps[i].force), i);
			}
			sobj->applyForces();
		}else{
			// there is no this m_Solids[k].type
		}
	} // each solid

	// bullet simulation
	mbt_World.dynamicsWorld->stepSimulation((btScalar)m_TH.dt, 1, (btScalar)m_TH.dt);

	// synchronize boundary particles with rigid's transform or softbody's nodes
	// for rigids, get transform in bullet, use it to update boundary paritcles
	// for softbody, use nodes to update boundary paritcles
	for(int n_r=int(m_Solids.size()),k=0; k<n_r; ++k) if(m_Solids[k].dynamic){
		m_Solids[k].updatePaticles();
		m_Solids[k].adjustSpace(m_TH.spaceMin, m_TH.spaceMax); // modify m_SpaceMin, m_SpaceMax
	} // each solid

}

// for rigid: set the position and rotation, ignore indices.
// for softbody: set nodes[indices[i]] = trans * innerPositions[indices[i]],
// recall that innerPositions are initial positions of nodes,
// if indices.size()==0, all nodes are processed.
void SphWithBullet::setSolidTransform( int solidIdx,
	const glm::mat4& trans, const std::vector<int>& indices )
{
	assert(0<=solidIdx && solidIdx<int(m_Solids.size()));

	Solid& solid = m_Solids[solidIdx];
	btTransform bt_trans; bt_trans.setFromOpenGLMatrix(&trans[0][0]);
	if(solid.type==Solid::RIGIDBODY){ // rigid body
		btRigidBody* rigid = btRigidBody::upcast(solid.mbtSolid_ptr->btObject); assert(rigid);
		rigid->getMotionState()->setWorldTransform(bt_trans);
	}else if(solid.type==Solid::SOFTBODY){ // soft body
		btSoftBody* soft = btSoftBody::upcast(solid.mbtSolid_ptr->btObject); assert(soft);
		if( indices.size()==0 ){
			for(int i=0; i<soft->m_nodes.size(); ++i)
				soft->m_nodes[i].m_q = soft->m_nodes[i].m_x =
				bt_trans * vect_to_btVect3(solid.innerPositions[i]);
		}else{
			for(size_t i=0; i<indices.size(); ++i)
				soft->m_nodes[indices[i]].m_q = soft->m_nodes[indices[i]].m_x =
				bt_trans * vect_to_btVect3(solid.innerPositions[indices[i]]);
		}
		soft->updateNormals();
	}else{
		//...
	}

	solid.updatePaticles();
	solid.adjustSpace(m_TH.spaceMin, m_TH.spaceMax); // modify m_SpaceMin, m_SpaceMax

}

// for rigid: transform the rigid, WCS-> world coordinate system, ignore indices.
// for softbody: set nodes[indices[i]] = trans * nodes[indices[i]], ignore isWCS,
// if indices.size()==0, all nodes are processed.
void SphWithBullet::addSolidTransform( int solidIdx,
	const glm::mat4& trans, bool isWCS, const std::vector<int>& indices )
{
	assert(0<=solidIdx && solidIdx<int(m_Solids.size()));

	Solid& solid = m_Solids[solidIdx];
	btTransform bt_trans; bt_trans.setFromOpenGLMatrix(&trans[0][0]);
	if(solid.type==Solid::RIGIDBODY){ // rigid body
		btRigidBody* rigid = btRigidBody::upcast(solid.mbtSolid_ptr->btObject); assert(rigid);
		btTransform bt_t; rigid->getMotionState()->getWorldTransform(bt_t);
		if(isWCS) bt_t = bt_trans * bt_t;
		else bt_t = bt_t * bt_trans;
		rigid->getMotionState()->setWorldTransform(bt_t);
	}else if(solid.type==Solid::SOFTBODY){ // soft body
		btSoftBody* soft = btSoftBody::upcast(solid.mbtSolid_ptr->btObject); assert(soft);
		if( indices.size()==0 ){
			/*for(int i=0; i<soft->m_nodes.size(); ++i)
				soft->m_nodes[i].m_x = bt_trans * soft->m_nodes[i].m_x;*/
			soft->transform(bt_trans);
		}else{
			for(size_t i=0; i<indices.size(); ++i){
				soft->m_nodes[indices[i]].m_x = bt_trans * soft->m_nodes[indices[i]].m_x;
				soft->m_nodes[indices[i]].m_q = bt_trans * soft->m_nodes[indices[i]].m_q;
			}
			soft->updateNormals();
		}
	}else{
		//...
	}

	solid.updatePaticles();
	solid.adjustSpace(m_TH.spaceMin, m_TH.spaceMax); // modify m_SpaceMin, m_SpaceMax

}


// for rigidbody, use rigid's transform
// for softbody, use bullet btSoftBody::tNodeArray
void SphSolid::updatePaticles()
{
	boundaryParticles.resize(innerPositions.size());

	if(type==RIGIDBODY){ // rigidbody
		btRigidBody* robj = btRigidBody::upcast(mbtSolid_ptr->btObject); assert(robj);
		btTransform trans; float t[16];
		robj->getMotionState()->getWorldTransform(trans);
		trans.getOpenGLMatrix(t);
		btVector3 v = robj->getLinearVelocity();
		btVector3 w = robj->getAngularVelocity();
		btVector3 o = robj->getCenterOfMassPosition();
		updatePaticles(SphWithBullet::btVec3_to_vect(o), t,
			SphWithBullet::btVec3_to_vect(v), SphWithBullet::btVec3_to_vec3t(w)); // parallel
	}else if(type=SOFTBODY){ // softbody
		btSoftBody* sobj = btSoftBody::upcast(mbtSolid_ptr->btObject); assert(sobj);
		int num = int(boundaryParticles.size()); assert(num==sobj->m_nodes.size());
	#pragma omp parallel for
		for(int i=0; i<num; ++i){
			boundaryParticles[i].position = SphWithBullet::btVec3_to_vect(sobj->m_nodes[i].m_x);
			boundaryParticles[i].velocity = SphWithBullet::btVec3_to_vect(sobj->m_nodes[i].m_v);
		}
	}else{
		// there is no this m_Solids[k].type
	}

}


