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
#include "SphWithIO.h"
#include "BulletCollision/CollisionShapes/btConvex2dShape.h"
#include "vcg_inc.h"


// draw fluid particles or boundary paricles of solids in OpenGL, iteratively
// parameter test use to select particles, i.e., cut down some particles
void SphWithIO::oglDrawFluidParts( void(*draw)(), bool(*test)(const vec_t& p) )  const
{
	assert(draw!=0);

	for(int num=int(m_Fluids.size()),k=0;k<num;++k){
		const std::vector<FluidParticle>& f_pts = m_Fluids[k].fluidParticles;
		for(int n=int(f_pts.size()),i=0; i<n; ++i){
			const vec_t& p = f_pts[i].position;
			if(test!=0 && test(p)==false) continue;
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, f_pts[i].color);
			glPushMatrix();
			glTranslatef( (float)p[0], (float)p[1], vec_t::dim==3?(float)p[2]:0 );
			draw();
			glPopMatrix();
		}
	}

}

void SphWithIO::oglDrawBounParts( void(*draw)(), bool(*test)(const vec_t& p) )  const
{
	assert(draw!=0);

	for(int num=int(m_Solids.size()),k=0; k<num; ++k){
		const std::vector<BoundPart>& r_pts = m_Solids[k].boundaryParticles;
		for(int n=int(r_pts.size()),i=0; i<n; ++i){
			const vec_t& p = r_pts[i].position;
			if(test!=0 && test(p)==false) continue;
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, r_pts[i].color);
			glPushMatrix();
			glTranslatef( (float)p[0], (float)p[1], vec_t::dim==3?(float)p[2]:0 );
			draw();
			glPopMatrix();
		}
	}

}

void SphWithIO::oglDrawSolid( int ridx )  const
{
	assert(m_Solids.size()==mbt_Solids.size());
	assert(ridx>=0 && ridx<int(m_Solids.size()));

	if(m_Solids[ridx].type==Solid::RIGIDBODY){ // rigidbody
		assert(mbt_Solids[ridx].triMesh!=0);
		btRigidBody* prb = btRigidBody::upcast(mbt_Solids[ridx].btObject); assert(prb);
		btTransform trans; float t[16];
		prb->getMotionState()->getWorldTransform(trans);
		trans.getOpenGLMatrix(t);
		glPushMatrix();
		glMultMatrixf(t);
		vcg_drawMesh(*mbt_Solids[ridx].triMesh);
		glPopMatrix();
	}else if(m_Solids[ridx].type==Solid::SOFTBODY){ // softbody
		btSoftBody* psb = btSoftBody::upcast(mbt_Solids[ridx].btObject); assert(psb);
		glBegin(GL_TRIANGLES);
		for(int i=0; i<psb->m_faces.size(); ++i) {
			glNormal3fv(&psb->m_faces[i].m_normal[0]);
			glVertex3fv(&psb->m_faces[i].m_n[0]->m_x[0]);
			glVertex3fv(&psb->m_faces[i].m_n[1]->m_x[0]);
			glVertex3fv(&psb->m_faces[i].m_n[2]->m_x[0]);
		}
		glEnd();
	}else{
		// there is no this m_Solids[k].type
	}

}

void SphWithIO::oglDrawSolid(bool(*test)(int i))  const
{
	for(int i=0;i<getNumSolids();++i){
		if(test!=0 && test(i)==false) continue;
		oglDrawSolid(i);
	}
}


// particles are sampled at regular cubic grid points
// the cube is defined by parameter leftBottom and rightTop(Axis-Aligned box)
// parameter fpart is the initial fluid paricle used to initialize fluid particles
// parameter newFluid = true: add a new fluid, fluidIdx will be ignored
// false: add particles to last fluid and the last three parameters are ignored
void SphWithIO::addFluidCuboid( bool newFluid, int fluidIdx,
	const vec_t& leftBottom, const vec_t& rightTop, const FluidPart& fpart,
	real_t restDensity, real_t viscosity, real_t tension )

{
	if(!newFluid) assert(fluidIdx>=0 && fluidIdx<int(m_Fluids.size()));
// ----------------- need not be parallel ----------------
	m_TH.spaceMin = m_TH.spaceMin.min(leftBottom);
	m_TH.spaceMax = m_TH.spaceMax.max(rightTop);

	veci_t	num = ((rightTop-leftBottom)/m_TH.spacing_r).to<int_t>();
	vec_t posMin = leftBottom + (rightTop-leftBottom-num.to<real_t>()*m_TH.spacing_r)*real_t(0.5);
	
	if(newFluid) m_Fluids.push_back(Fluid(restDensity,viscosity,tension));
	std::vector<FluidPart>& fluidParts =
		newFluid ? m_Fluids.back().fluidParticles
				 : m_Fluids[fluidIdx].fluidParticles;
	//fluidParts.reserve( fluidParts.size() + (num+1).volume()+4 );

	FluidPart part = fpart; veci_t pi(0); //real_t d = m_TH.spacing_r/5;
	while(pi[vec_t::dim-1] <= num[vec_t::dim-1]){
		part.position = posMin + pi.to<real_t>()*m_TH.spacing_r;
		/*int k = 0;
		for(int i=0; i<vec_t::dim; ++i)
			if(pi[i]==0 || pi[i]==num[i]) ++k;
		if(k>=2) for(int i=0; i<vec_t::dim; ++i){
			if(pi[i]==0) part.position[i] += d;
			if(pi[i]==num[i]) part.position[i] -= d;
		}*/
		fluidParts.push_back(part);
		++ pi[0];
		for(int i=0; i<vec_t::dim-1; ++i){
			if(pi[i]>num[i]){
				pi[i] = 0;
				++ pi[i+1];
			}
		}
	}
	
}

// add particles to fluid[fluidIdx]
void SphWithIO::addFluidParts(int fluidIdx, const std::vector<FluidPart>& particles)
{
	assert(fluidIdx>=0 && fluidIdx<int(m_Fluids.size()));

	std::vector<FluidPart>& fluidParts = m_Fluids[fluidIdx].fluidParticles;
	// do not use fluidParts.reserve(), for frequent call, push_back will be fast
	for(size_t i=0; i<particles.size(); ++i) fluidParts.push_back(particles[i]);

}

//see 9
void SphWithIO::removeFluidParts(int fluidIdx, int particle) {
	assert(fluidIdx >= 0 && fluidIdx<int(m_Fluids.size()));
	std::vector<FluidPart>& fluidParts = m_Fluids[fluidIdx].fluidParticles;
	std::vector<FluidPart>::iterator itr = fluidParts.begin();
	fluidParts.erase(itr+particle);
}

//see 9
void SphWithIO::addCandidateFromPLY(const char* meshFileName, const char* sampleFileName, int candidateIdx,
	const CandidatePart& cpart, real_t viscosity, const glm::mat4& transform) {
	return;
}

// the PLY file contain the sampled points
void SphWithIO::addFluidFromPLY( const char* filename,
	bool newFluid, int fluidIdx, const FluidPart& fpart,
	real_t restDensity, real_t viscosity, real_t tension )
{
	if(!newFluid) assert(fluidIdx>=0 && fluidIdx<int(m_Fluids.size()));

	if(newFluid) m_Fluids.push_back(Fluid(restDensity,viscosity,tension));
	std::vector<FluidPart>& fluidParts =
		newFluid ? m_Fluids.back().fluidParticles
		: m_Fluids[fluidIdx].fluidParticles;

	GLTriMesh mesh; vcg_loadFromPLY(mesh, filename, false, false, false);
	std::vector<float> vers; vcg_saveToData(mesh, &vers, 0);
	FluidPart p0 = fpart;
	for(size_t i=0; i<vers.size()/3; ++i){
		p0.position[0] = vers[i*3];
		p0.position[1] = vers[i*3+1];
		if(vec_t::dim==3) p0.position[2] = vers[i*3+2];
		fluidParts.push_back(p0);
	}

}


void SphWithIO::makeCubeSampleParts(
	std::vector<vec_t>& samples, const vec_t& halfWidth, bool haveTop, real_t r) const
{
	vec_t leftBottom=-halfWidth, rightTop=halfWidth;
	vec_t factor = (rightTop-leftBottom)/r;
	for(int i=0; i<vec_t::dim; ++i) factor[i] = std::ceil(factor[i]);
	veci_t num = factor.to<int_t>(); vec_t pro;
	for(int i=0; i<vec_t::dim; ++i){
		if(num[i]>0) pro[i] = (rightTop[i]-leftBottom[i])/num[i];
		else pro[i] = 1;
	}

	//samples.reserve( (num+1).volume()+4 );
	veci_t pi(0);
	while(pi[vec_t::dim-1] <= num[vec_t::dim-1]){
		if( veci_t::O<pi && pi<num )
			goto endAdd;
		if( !haveTop && pi[1]==num[1] ){
			-- pi[1];
			if( veci_t::O<pi && pi<num ){
				++ pi[1];
				goto endAdd;
			}else
				++ pi[1];
		}
		samples.push_back(leftBottom + pi.to<real_t>()*pro);
	endAdd:
		++ pi[0];
		for(int i=0; i<vec_t::dim-1; ++i)
			if( pi[i] > num[i] )
				pi[i] = 0, ++ pi[i+1];
	}

}

void SphWithIO::makeCubeMeshAndBtShape( GLTriMesh& mesh, btCollisionShape** shape,
	const vec_t& halfWidth, bool haveTop, real_t thickness) const
{
#define VERS_CPY(a,b,sMin,sMax,d) \
	for(int si=0; si<sizeof(b)/sizeof(int); ++si) \
	switch(b[si]){ \
	case -2: a.push_back(float(sMin[si%3]-d)); break; \
	case -1: a.push_back(float(sMin[si%3]+d)); break; \
	case  1: a.push_back(float(sMax[si%3]-d)); break; \
	case  2: a.push_back(float(sMax[si%3]+d)); break; \
	case  0: a.push_back(0); \
	}
#define FACS_CPY(a,b) \
	{for(int si=0; si<sizeof(b)/sizeof(int); ++si) \
	a.push_back(b[si]-1);} // reason of -1, indeces of .obj file are from 1, not 0 

	vec_t leftBottom=-halfWidth, rightTop=halfWidth;
	real_t d = thickness/2; float margin = float(thickness/4);
	std::vector<float> vers; std::vector<int> tris;
	if(haveTop){ // have top
		if(vec_t::dim==3){
			int temp_vers[]={-2,-2,2,-2,-2,-2,2,-2,-2,2,-2,2,-2,2,2,-2,2,-2,2,2,-2,2,2,2};
			int temp_facs[]={5,6,2,5,2,1,6,7,3,6,3,2,7,8,4,7,4,3,8,5,1,8,1,4,1,2,3,1,3,4,8,7,6,8,6,5};
			VERS_CPY(vers, temp_vers, leftBottom, rightTop, d);
			FACS_CPY(tris, temp_facs);
			*shape = new btConvexHullShape( &vers[0], int(vers.size())/3, 3*sizeof(float) );
		}else{ // 2D
			int temp_vers[]={-2,2,0,-2,-2,0,2,-2,0,2,2,0}; // without hole
			int temp_facs[]={1,2,3,3,4,1};
			VERS_CPY(vers, temp_vers, leftBottom, rightTop, d);
			FACS_CPY(tris, temp_facs);
			btConvexShape* convexShape =
				new btConvexHullShape(&vers[0], int(vers.size())/3, 3*sizeof(float));
			convexShape->setMargin(margin);
			*shape = new btConvex2dShape(convexShape);
		}
	}else{ // have no top
		std::vector<std::vector<float>> vers_arry;
		std::vector<std::vector<int>> tris_arry;
		if(vec_t::dim==3){
			// from .obj file created by blender
			int temp_vers[]={2,-2,-2,2,-2,2,-2,-2,2,-2,-2,-2,2,2,-2,-2,2,-2,2,2,2,-2,2,2,-1,2,1,-1,2,-1,1,
				2,-1,1,2,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1};
			int temp_facs[]={1,2,3,1,3,4,5,1,4,5,4,6,1,5,7,1,7,2,2,7,8,2,8,3,3,8,6,3,6,4,10,11,5,10,5,6,10,
				6,8,10,8,9,9,8,12,12,8,7,12,7,11,11,7,5,10,9,13,10,13,14,11,10,14,11,14,15,13,16,15,13,15,
				14,9,12,16,9,16,13,12,11,15,12,15,16};
			VERS_CPY(vers, temp_vers, leftBottom, rightTop, d);
			FACS_CPY(tris, temp_facs);
			// convex decomposition, 5 convex elements
			int temp_vers_x0[]={-2,-2,2,-2,-2,-2,-2,2,-2,-2,2,2,-1,2,-1,-1,2,1,-1,-1,1,-1,-1,-1};
			int temp_facs_x0[]={1,4,3,1,3,2,5,3,4,5,4,6,5,6,7,5,7,8,8,2,3,8,3,5,8,7,1,8,1,2,6,4,1,6,1,7};
			int temp_vers_x1[]={2,-2,-2,2,-2,2,2,2,-2,2,2,2,1,2,-1,1,2,1,1,-1,-1,1,-1,1};
			int temp_facs_x1[]={1,3,4,1,4,2,6,4,5,5,4,3,6,8,2,6,2,4,1,7,5,1,5,3,8,7,1,8,1,2,6,5,7,6,7,8};
			int temp_vers_y0[]={2,-2,2,-2,-2,2,2,2,2,-2,2,2,-1,2,1,1,2,1,-1,-1,1,1,-1,1};
			int temp_facs_y0[]={1,3,4,1,4,2,5,4,6,6,4,3,6,3,1,6,1,8,4,5,7,4,7,2,5,6,8,5,8,7,7,8,1,7,1,2};
			int temp_vers_y1[]={2,-2,-2,-2,-2,-2,2,2,-2,-2,2,-2,-1,2,-1,1,2,-1,-1,-1,-1,1,-1,-1};
			int temp_facs_y1[]={3,1,2,3,2,4,5,6,3,5,3,4,6,5,7,6,7,8,2,7,5,2,5,4,1,8,7,1,7,2,1,3,6,1,6,8};
			int temp_vers_z0[]={2,-2,-2,2,-2,2,-2,-2,2,-2,-2,-2,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1};
			int temp_facs_z0[]={1,2,3,1,3,4,8,2,1,8,1,7,5,8,7,5,7,6,6,4,3,6,3,5,5,3,2,5,2,8,6,7,1,6,1,4};
			vers_arry.resize(5); tris_arry.resize(5);
			VERS_CPY(vers_arry[0], temp_vers_x0, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[0], temp_facs_x0);
			VERS_CPY(vers_arry[1], temp_vers_x1, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[1], temp_facs_x1);
			VERS_CPY(vers_arry[2], temp_vers_y0, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[2], temp_facs_y0);
			VERS_CPY(vers_arry[3], temp_vers_y1, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[3], temp_facs_y1);
			VERS_CPY(vers_arry[4], temp_vers_z0, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[4], temp_facs_z0);
			btCompoundShape* compound = new btCompoundShape();
			for(int i=0; i<(int)vers_arry.size(); ++i){
				std::vector<float>& vs = vers_arry[i];
				btConvexHullShape* convexShape = new btConvexHullShape(&vs[0], int(vs.size())/3, 3*sizeof(float));
				convexShape->setMargin(margin);
				btTransform ltrans; ltrans.setIdentity();
				compound->addChildShape(ltrans,convexShape);
			}
			*shape = compound;
		}else{ // 2D
			int temp_vers[]={-2,2,0,-2,-2,0,2,2,0,1,2,0,1,-1,0,-1,-1,0,-1,2,0,2,-2,0};
			int temp_facs[]={5,6,2,5,8,3,2,8,5,6,1,2,6,7,1,4,5,3};
			VERS_CPY(vers, temp_vers, leftBottom, rightTop, d);
			FACS_CPY(tris, temp_facs);
			// convex decomposition, 3 convex elements
			int temp_vers_x0[]={-1,-1,0,-2,-2,0,-2,2,0,-1,2,0};
			int temp_facs_x0[]={1,3,2,1,4,3};
			int temp_vers_x1[]={1,-1,0,2,-2,0,2,2,0,1,2,0};
			int temp_facs_x1[]={1,2,3,4,1,3};
			int temp_vers_y0[]={1,-1,0,-1,-1,0,-2,-2,0,2,-2,0};
			int temp_facs_y0[]={1,2,3,4,1,3};
			vers_arry.resize(3); tris_arry.resize(3);
			VERS_CPY(vers_arry[0], temp_vers_x0, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[0], temp_facs_x0);
			VERS_CPY(vers_arry[1], temp_vers_x1, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[1], temp_facs_x1);
			VERS_CPY(vers_arry[2], temp_vers_y0, leftBottom, rightTop, d);
			FACS_CPY(tris_arry[2], temp_facs_y0);
			// btCompoundShape
			btCompoundShape* compound = new btCompoundShape();
			for(int i=0; i<(int)vers_arry.size(); ++i){
				std::vector<float>& vs = vers_arry[i];
				btConvexHullShape* convexShape = new btConvexHullShape(&vs[0], int(vs.size())/3, 3*sizeof(float));
				convexShape->setMargin(margin);
				btTransform ltrans; ltrans.setIdentity();
				compound->addChildShape(ltrans,new btConvex2dShape(convexShape));
			}
			*shape = compound;
		}
	} // if(haveTop)
	(*shape)->setMargin(margin);

	vcg_loadFromData(mesh, vers, tris);

}


// if mass>0 dynamic=true
void SphWithIO::addRigidCuboid(	const vec_t& halfWidth,
	const BoundPart& bpart, float mass, real_t viscosity,
	bool haveTop, const glm::mat4& transform )
{	
	// btRigidBody
	GLTriMesh* trimesh=new GLTriMesh(); btCollisionShape* cshape;
	makeCubeMeshAndBtShape(*trimesh, &cshape, halfWidth, haveTop, m_TH.spacing_r );
	btRigidBody* prigid; btTransform trans; trans.setFromOpenGLMatrix(&transform[0][0]);
	makeBtRigidBody(&prigid, cshape, mass, trans);

	MBtSolid btsolid(prigid, trimesh);
	Solid sphsolid( mass>0,Solid::RIGIDBODY,viscosity,0 );
	makeCubeSampleParts( sphsolid.innerPositions, halfWidth, haveTop, m_TH.spacing_r );
	sphsolid.boundaryParticles.resize(sphsolid.innerPositions.size(), bpart);

	addSolid( sphsolid, btsolid );

}

// convex hull of the mesh will be the collision shape, i.e., the mesh itself may not
void SphWithIO::addRigidFromPLY( const char* meshFileName, const char* sampleFileName,
	const BoundPart& bpart, float mass, real_t viscosity,
	const glm::mat4& transform )
{
	// btCollisionShape
	GLTriMesh* trimesh = new GLTriMesh(); vcg_loadFromPLY(*trimesh, meshFileName);
	btCollisionShape* cshape;
	makeBtConvexHullShape(&cshape, *trimesh);
	// btRigidBody
	btRigidBody* prigid; btTransform trans; trans.setFromOpenGLMatrix(&transform[0][0]);
	makeBtRigidBody(&prigid, cshape, mass, trans);

	MBtSolid btsolid(prigid, trimesh);
	Solid sphsolid( mass>0,Solid::RIGIDBODY,viscosity,0 );
	// sample points
	GLTriMesh samples; vcg_loadFromPLY(samples, sampleFileName, false, false, false);
	std::vector<float> vers; vcg_saveToData(samples, &vers, 0);
	for(size_t i=0; i<vers.size()/3; ++i){
		sphsolid.innerPositions.push_back( vec_t(vers[i*3],vers[i*3+1],vers[i*3+2]) );
	}
	sphsolid.boundaryParticles.resize(sphsolid.innerPositions.size(), bpart);

	addSolid( sphsolid, btsolid );

}

// "<convexname>0", "<convexname>1", "<convexname>2", ... are the convex parts of the mesh,
// i.e., mainname is used to display in opengl, btCompoundShape of convecnum "<convexname>i"
// is used for collision detection
void SphWithIO::addRigidFromPLY( const char* meshFileName, const char* sampleFileName,
	const char* convexMeshName, int convecnum,
	const BoundPart& bpart, float mass, real_t viscosity,
	const glm::mat4& transform )
{
	// btCollisionShape
	GLTriMesh* trimesh = new GLTriMesh(); vcg_loadFromPLY(*trimesh, meshFileName);
	btCollisionShape* cshape; GLTriMesh m; btTransform btt; btt.setIdentity();
	btCompoundShape* ccompound = new btCompoundShape();
	for(int i=0; i<convecnum; ++i){
		char s[50]; sprintf_s(s,"%d",i);
		vcg_loadFromPLY( m, (convexMeshName+std::string(s)+".ply").c_str(), false, false, false );
		makeBtConvexHullShape(&cshape, m);
		ccompound->addChildShape(btt, cshape);
	}
	ccompound->setMargin(m_TH.spacing_r/4);
	
	// btRigidBody
	btRigidBody* prigid; btTransform trans; trans.setFromOpenGLMatrix(&transform[0][0]);
	makeBtRigidBody(&prigid, ccompound, mass, trans);

	MBtSolid btsolid(prigid, trimesh);
	Solid sphsolid( mass>0,Solid::RIGIDBODY,viscosity,0 );
	// sample points
	GLTriMesh samples; vcg_loadFromPLY(samples, sampleFileName, false, false, false);
	std::vector<float> vers; vcg_saveToData(samples, &vers, 0);
	for(size_t i=0; i<vers.size()/3; ++i){
		sphsolid.innerPositions.push_back( vec_t(vers[i*3],vers[i*3+1],vers[i*3+2]) );
	}
	sphsolid.boundaryParticles.resize(sphsolid.innerPositions.size(), bpart);

	addSolid( sphsolid, btsolid );

}


// mass must >0, dynamic is certainly true
// non-tranformed rope is on x-axis
void SphWithIO::addSoftRopeLine(
	float halfWidth_x, int endFixedBit,
	const BoundPart& bpart, float mass, real_t viscosity,
	const glm::mat4& transform, int twoEndIdx[2] )
{
	;
}

// mass must >0, dynamic is certainly true
// non-tranformed cloth is on XoZ plane
void SphWithIO::addSoftClothPlane(
	float halfWidth_x, float halfWidth_z, float stretchDis_y, int cornerFixedBit,
	const BoundPart& bpart, float mass, real_t viscosity,
	const glm::mat4& transform, int fourCornerIdx[4] )
{
	int resolution_x = int(halfWidth_x*2/m_TH.spacing_r+0.5f)+1;
	int resolution_z = int(halfWidth_z*2/m_TH.spacing_r+0.5f)+1;
	btSoftBody* psb = btSoftBodyHelpers::CreatePatch( *mbt_World.softBodyWorldInfo,
		btVector3(-halfWidth_x,0,-halfWidth_z), btVector3(-halfWidth_x,0,+halfWidth_z), // for corners
		btVector3(+halfWidth_x,0,-halfWidth_z), btVector3(+halfWidth_x,0,+halfWidth_z),
		resolution_x, resolution_z, // x,y resolution
		cornerFixedBit/*1+2+4+8*/, // bit: fix(mass=0) or not for four corners
		true ); // psb->appendLink for diagonal
	psb->getCollisionShape()->setMargin( m_TH.spacing_r/4 );
	btSoftBody::Material* pm=psb->appendMaterial();
	pm->m_kLST = 1.0f; // linear stiffness
	pm->m_kAST = 0.0f; // angular stiffness
	pm->m_kVST = 0.0f; // volume stiffness
	//psb->generateBendingConstraints(2, pm);

	psb->m_nodes[0].m_x[1] =
		psb->m_nodes[resolution_x-1].m_x[1] =
		psb->m_nodes[psb->m_nodes.size()-resolution_x].m_x[1] =
		psb->m_nodes[psb->m_nodes.size()-1].m_x[1] = stretchDis_y;

	if(fourCornerIdx){
		fourCornerIdx[0] = 0;
		fourCornerIdx[1] = resolution_x-1;
		fourCornerIdx[2] = psb->m_nodes.size()-resolution_x;
		fourCornerIdx[3] = psb->m_nodes.size()-1;
	}

	psb->setTotalMass(mass);
	psb->updateNormals();

	psb->m_cfg.piterations = 20;
	psb->m_cfg.citerations = 10;
	psb->m_cfg.diterations = 10;
	//psb->m_cfg.viterations = 10;

	//psb->m_cfg.kSRHR_CL = 1;
	//psb->m_cfg.kSR_SPLT_CL = 0;
	//psb->m_cfg.collisions = btSoftBody::fCollision::CL_SS+btSoftBody::fCollision::CL_RS;
	//psb->m_cfg.collisions |= btSoftBody::fCollision::CL_SELF;
	//psb->generateClusters(0);

	//add the body to the dynamics world
	//m_DynamicsWorld->addSoftBody(psb);

	MBtSolid btsolid(psb, 0);
	Solid sphsolid( true,Solid::SOFTBODY,viscosity,0 );

	btTransform trans; trans.setFromOpenGLMatrix(&transform[0][0]);
	psb->transform(trans);
	// transform nodes, nodes' velocity
	for(int i=0; i<psb->m_nodes.size(); ++i){
		//psb->m_nodes[i].m_x = trans * psb->m_nodes[i].m_x;
		//psb->m_nodes[i].m_v = vect_to_btVect3(bpart.velocity);
		sphsolid.innerPositions.push_back(btVec3_to_vect(psb->m_nodes[i].m_x));
	}
	sphsolid.boundaryParticles.resize(sphsolid.innerPositions.size(), bpart);

	addSolid( sphsolid, btsolid );

}

// tetrahedral mesh
void SphWithIO::addSoftVolumeCuboid(
	const vec_t& halfWidth,
	const BoundPart& bpart, float mass, real_t viscosity,
	const glm::mat4& transform, int eightCornerIdx[8] )
{
	;
}

// each vetex is a node, edge a link
void SphWithIO::addSoftFromPly( const char* meshFileName,
	const BoundPart& bpart, float mass, real_t viscosity,
	const glm::mat4& transform )
{
	GLTriMesh trimesh; vcg_loadFromPLY(trimesh, meshFileName, false, false, false);
	std::vector<float> vers; std::vector<int> tris;
	vcg_saveToData(trimesh, &vers, &tris);
	btSoftBody*	psb=btSoftBodyHelpers::CreateFromTriMesh(
		*mbt_World.softBodyWorldInfo, &vers[0], &tris[0], int(tris.size())/3 );
	psb->generateBendingConstraints(2);
	psb->m_cfg.piterations=2;
	//psb->m_cfg.collisions |= btSoftBody::fCollision::VF_SS;
	psb->randomizeConstraints();
	btTransform trans; trans.setFromOpenGLMatrix(&transform[0][0]);
	psb->transform( trans );
	
	MBtSolid btsolid(psb, 0);
	Solid sphsolid( true,Solid::SOFTBODY,viscosity,0 );

	for(int i=0; i<psb->m_nodes.size(); ++i){
		sphsolid.innerPositions.push_back(btVec3_to_vect(psb->m_nodes[i].m_x));
	}
	sphsolid.boundaryParticles.resize(sphsolid.innerPositions.size(), bpart);

	addSolid( sphsolid, btsolid );

}


void SphWithIO::duplicateLastSolid( const glm::mat4& transform )
{
	if(m_Solids.back().type==Solid::RIGIDBODY){ // rigidbody
		btRigidBody* lastrigid = btRigidBody::upcast(mbt_Solids.back().btObject); assert(lastrigid);
		btTransform trans; trans.setFromOpenGLMatrix(&transform[0][0]);
		btDefaultMotionState* myMotionState = new btDefaultMotionState(trans);
		btVector3 localInertia(0,0,0); btScalar mass = lastrigid->getInvMass()==0 ? 0 : 1/lastrigid->getInvMass();
		if(mass>0) lastrigid->getCollisionShape()->calculateLocalInertia(mass,localInertia);
		btRigidBody* prigid = new btRigidBody( mass,myMotionState,lastrigid->getCollisionShape(),localInertia );

		MBtSolid btsolid = mbt_Solids.back();
		btsolid.btObject = prigid;
		Solid sphsolid = m_Solids.back();

		addSolid( sphsolid, btsolid );

	}else if(m_Solids.back().type==Solid::SOFTBODY){ // softbody
		//...
	}else{
		//...
	}
	
}

void SphWithIO::saveSolidMeshToPLY(
	int solidIdx, const char* filename, bool simulated, bool coordYupToZup) const
{
	assert( 0<=solidIdx && solidIdx<getNumSolids() );

	const Solid& solid = m_Solids[solidIdx];
	if(solid.type==Solid::RIGIDBODY){ // rigid body
		assert(solid.mbtSolid_ptr->triMesh);
		if( simulated ){
			// Opengl(glm) stores matrix in  column-major order. Matrix44 stores matrix in row-major order.
			GLTriMesh mesh; vcg::Matrix44f mat;
			glm::mat4 gm; getRigidBodyTransform(solidIdx, gm);
			for(int i=0; i<16; ++i) mat[i/4][i%4]=gm[i%4][i/4]; // copy gm to mat
			vcg::tri::Append<GLTriMesh,GLTriMesh>::MeshCopy(
				mesh, *solid.mbtSolid_ptr->triMesh, false,false);
			vcg::tri::UpdatePosition<GLTriMesh>::Matrix(mesh, mat, false);
			vcg_saveToPLY(mesh, filename, coordYupToZup, true);
		}else{
			vcg_saveToPLY(*solid.mbtSolid_ptr->triMesh, filename, coordYupToZup, true);
		}
	}else if(solid.type==Solid::SOFTBODY){ // soft body
		if( simulated ){
			std::vector<float> vers; std::vector<int> tris;
			getSoftBodyMeshData(solidIdx, &vers, &tris);
			GLTriMesh mesh; vcg_loadFromData(mesh, vers, tris, false,false,false);
			vcg_saveToPLY(mesh, filename, coordYupToZup, true);
		}else{
			assert( solid.innerPositions.size()==solid.boundaryParticles.size() );
			std::vector<float> vers;
			for(size_t i=0; i<solid.innerPositions.size(); ++i){
				vers.push_back(solid.innerPositions[i][0]);
				vers.push_back(solid.innerPositions[i][1]);
				if(vec_t::dim==3) vers.push_back(solid.innerPositions[i][2]);
			}
			std::vector<int> tris; getSoftBodyMeshData(solidIdx, 0, &tris);
			GLTriMesh mesh; vcg_loadFromData(mesh, vers, tris, false,false,false);
			vcg_saveToPLY(mesh, filename, coordYupToZup, true);
		}
	}else{
		//...
	}

}

// non-transformed rigid mesh, and initial softbody mesh
void SphWithIO::saveAllInitialSolidMeshToPLY(
	const char* namePrefix, bool coordYupToZup, bool staticSolidsTransformed) const
{
	for(int si=0; si<getNumSolids(); ++si){
		string filename(namePrefix); char ss[101];
		sprintf(ss, "_solid_%d_", si); filename += ss;
		const Solid& solid = m_Solids[si];
		if(solid.dynamic) filename += "dynamic_";
		else filename += "static_";
		if(solid.type==Solid::RIGIDBODY)
			if(staticSolidsTransformed && solid.dynamic==false)
				filename += "transformed_rigid.ply";
			else filename += "non-transformed_rigid.ply";
		else if(solid.type==Solid::SOFTBODY) filename += "initial_softbody.ply";
		else {  }
		bool simulated = false;
		if(staticSolidsTransformed && solid.dynamic==false) simulated = true;
		saveSolidMeshToPLY(si, filename.c_str(), simulated, coordYupToZup);
	}

}

// trasnformed rigid mesh, and simulated softbody mesh
void SphWithIO::saveAllSimulatedSolidMeshToPLY(
	const char* namePrefix, bool coordYupToZup, bool staticSolids) const
{
	for(int si=0; si<getNumSolids(); ++si){
		string filename(namePrefix); char ss[101];
		const Solid& solid = m_Solids[si];
		if(staticSolids==true && solid.dynamic==false) continue;
		sprintf(ss, "_solid_%d_", si); filename += ss;
		if(solid.type==Solid::RIGIDBODY) filename += "rigid.ply";
		else if(solid.type==Solid::SOFTBODY) filename += "softbody.ply";
		else { }
		saveSolidMeshToPLY(si, filename.c_str(), true, coordYupToZup);
	}

}


/*
void SphWithIO::setExampleScene(bool rigid)
{
	const real_t vis = 0.05f;
	if(vec_t::dim==3){
		m_TH.spacing_r = 0.1f;
		m_TH.smoothRadius_h = 2*m_TH.spacing_r;
		m_TH.dt = real_t(1.0)*m_TH.spacing_r/m_TH.soundSpeed_cs;

		m_TH.spaceMin.set(-6);
		m_TH.spaceMax.set(6); m_TH.spaceMax[1]=12;

		// container
		BoundPart bp0; bp0.color[0]=bp0.color[1]=bp0.color[2]=0.4f; bp0.color[3]=1;
		bp0.position=bp0.velocity=bp0.force=vec_t::O;
		vec_t rb(4.0f); rb[1]=2.0f;
		addRigidCuboid( rb, 0, vis, bp0, false,
			glm::translate(glm::vec3(0,1.5f+m_TH.spacing_r/2,0)) );

		// cloth
		addSoftCloth( 2, 2, 0.2f, 1+2+4+8, 100, vis, bp0,
			glm::translate(glm::vec3(0,2,0)) );

		// water
		FluidPart fp0; fp0.velocity = vec_t::O; fp0.density=1000;
		fp0.color[1]=fp0.color[2]=0.3f; fp0.color[0]=0.9f; fp0.color[3]=1;
		vec_t lb(-1.5f), rt(1.5f); rt[1]=5; lb[1]+=4.5f; rt[1]+=4.5f;
		addFluidCuboid( true, 0, lb, rt, fp0, 1000, vis, 1 );
		
		if(rigid){
			addRigidCuboid( vec_t(0.5f), 300, vis, bp0, true,
				glm::translate(glm::vec3(-2.1f,5,-2.1f)) );
			duplicateLastSolid(
				glm::translate(glm::vec3(2.1f,5,-2.1f)) );
			duplicateLastSolid(
				glm::translate(glm::vec3(2.1f,5,2.1f)) );
			duplicateLastSolid(
				glm::translate(glm::vec3(-2.1f,5,2.1f)) );

			duplicateLastSolid(
				glm::translate(glm::vec3(-0.8f,10.2f,-0.8f)) );
			duplicateLastSolid(
				glm::translate(glm::vec3(0.8f,10.2f,-0.8f)) );
			duplicateLastSolid(
				glm::translate(glm::vec3(0.8f,10.2f,0.8f)) );
			duplicateLastSolid(
				glm::translate(glm::vec3(-0.8f,10.2f,0.8f)) );
		}

	}else{
		
	}

	logInfoOfScene();

}*/

