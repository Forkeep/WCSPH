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

#include "SphBase.h"


/*
 * ----- fluid and solid --------------------------------------------------------------------------
*/

// Remove particles that are outside of the simulation space defined by spaceMin and spaceMax
// return number of particles removed
int SphFluid::removeOutsideParts(const vec_t& spMin, const vec_t& spMax)
{
    if(!fluidParticles.size()) return 0;

// ----------------- not parallel because of dependency ----------------
    int num = int(fluidParticles.size()),i_front=0, i_back=num-1;
    while(i_front<i_back){
        // find the next particle which is NOT inside space
        while( i_front<i_back && fluidParticles[i_front].position.inside(spMin,spMax) )
            ++ i_front;
        // find the next particle which is inside space
        while( i_back>i_front && !fluidParticles[i_back].position.inside(spMin,spMax) )
            -- i_back;
        if( i_front<i_back ) {
            fluidParticles[i_front] = fluidParticles[i_back];
            ++ i_front; -- i_back;
        }
    }
    if( !fluidParticles[i_front].position.inside(spMin,spMax) )
        -- i_front;
    fluidParticles.resize(i_front+1);
    return num-(i_front+1);
}

// For particles that are outside of the simulation space defined by spaceMin and spaceMax,
// correct their position and velocity, as if spaceMin and spaceMax difine 6 walls
// return number of particles corrected
void SphFluid::correctOutsideParts(const vec_t& spMin, const vec_t& spMax, real_t v0)
{
    int num = int(fluidParticles.size());
#pragma omp parallel for
    for(int i=0; i<num; ++i){
        vec_t& ipos = fluidParticles[i].position;
        vec_t& ivel = fluidParticles[i].velocity;
        for(int k=0; k<vec_t::dim; ++k){ // for each dimension
            if(ipos[k]<=spMin[k]){
                ipos[k] = spMin[k];
                if(ivel[k]<0) ivel[k] = v0;
            }
            if(ipos[k]>=spMax[k]){
                ipos[k] = spMax[k];
                if(ivel[k]>0) ivel[k] = -v0;
            }
        }
    }
}

// For particles that are outside of the simulation space defined by spaceMin and spaceMax,
// change spaceMin and spaceMax to contain them
// return true if spaceMin or spaceMax is changed
bool SphFluid::adjustSpace(vec_t& spMin, vec_t& spMax) const
{
// ----------------- not parallel because of dependency ----------------
    bool re = false;
    vec_t mi(spMin), ma(spMax);
    for(size_t num=fluidParticles.size(),i=0; i<num; ++i){
        if( !fluidParticles[i].position.inside(mi, ma) ){
            re = true;
            mi = mi.min(fluidParticles[i].position);
            ma = ma.max(fluidParticles[i].position);
        }
    }
    spMin = mi, spMax = ma;
    return re;
}


// recompute boundaryParticles from innerPositions, i.e. using rigid's transform
// pragmeter tranform is OpenGL transform matrix (tranform[0-3] is the first column)
void SphSolid::updatePaticles(
    const vec_t& centerOfMass, const float* tranform,
    const vec_t& linearVelocity, const vec3_t& angularVelocity)
{
    boundaryParticles.resize(innerPositions.size());
    int n = int(innerPositions.size());
#pragma omp parallel for
    for(int i=0; i<n; ++i){
        vec_t& p0 = innerPositions[i];
        vec_t& bp = boundaryParticles[i].position;
        if( tranform ){
            if(vec_t::dim==3){
                // note: tranform[0-3] is the first column (not row) of the OpenGL matrix
                bp[0] = tranform[0]*p0[0] + tranform[4]*p0[1] + tranform[8] *p0[2] + tranform[12];
                bp[1] = tranform[1]*p0[0] + tranform[5]*p0[1] + tranform[9] *p0[2] + tranform[13];
                bp[2] = tranform[2]*p0[0] + tranform[6]*p0[1] + tranform[10]*p0[2] + tranform[14];
                real_t alpha =
                    tranform[3]*p0[0] + tranform[7]*p0[1] + tranform[11]*p0[2] + tranform[15];
                bp[0] /= alpha; bp[1] /= alpha; bp[2] /= alpha;
            }else{
                // note: tranform[0-3] is the first column (not row) of the OpenGL matrix
                bp[0] = tranform[0]*p0[0] + tranform[4]*p0[1] + tranform[12];
                bp[1] = tranform[1]*p0[0] + tranform[5]*p0[1] + tranform[13];
                real_t alpha = tranform[3]*p0[0] + tranform[7]*p0[1] + tranform[15];
                bp[0] /= alpha; bp[1] /= alpha;
            }
            
        }else{
            bp = p0;
        }
        boundaryParticles[i].velocity =
            linearVelocity + vec_t(angularVelocity.cross(bp-centerOfMass));
    }
}

// compute total force and torque from boundaryParticles' forces
void SphSolid::rigidForceJoin(
    const vec_t& centerOfMass, vec_t& totalForce, vec3_t& totalTorque) const
{
    vec_t tf(0); vec3_t tt(0);
    size_t n = boundaryParticles.size();
//#pragma omp parallel for reduction(+:tf, tt)
    for(size_t i=0; i<n; ++i){
        const BoundaryParticle& bp = boundaryParticles[i];
        tf += bp.force;
        tt += (bp.position-centerOfMass).cross(bp.force);
    }
    totalForce = tf; totalTorque = tt;
}

// For particles that are outside of the simulation space defined by spaceMin and spaceMax,
// change spaceMin and spaceMax to contain them
// return true if spaceMin or spaceMax is changed
bool SphSolid::adjustSpace(vec_t& spMin, vec_t& spMax) const
{
// ----------------- not parallel because of dependency ----------------
    bool re = false;
    for(size_t num=boundaryParticles.size(),i=0; i<num; ++i){
        if( !boundaryParticles[i].position.inside(spMin, spMax) ){
            re = true;
            spMin = spMin.min(boundaryParticles[i].position);
            spMax = spMax.max(boundaryParticles[i].position);
        }
    }
    return re;
}

// get the index of four corners of softbody
/*void SphSolid::getFourCornerIdxOfCloth(int idx[4]) {
	if(type!=SOFTBODY) return;
	if(innerPositions.size()<9) return; // 3x3

	idx[0] = 0; idx[3] = int(innerPositions.size())-1;
	real_t r = 1.9f*(innerPositions[0]-innerPositions[1]).len_square();
	for(size_t i=0; i<innerPositions.size()-1; ++i){
		if( (innerPositions[i]-innerPositions[i+1]).len_square()>r ){
			idx[1] = int(i); break; }
	}
	for(size_t i=innerPositions.size()-1; i>0; --i){
		if( (innerPositions[i]-innerPositions[i-1]).len_square()>r ){
			idx[2] = int(i); break; }
	}

}*/



/*
 * ----- 2D or 3D SPH simulator base class, implementation ----------------------------------------
*/
#include <ctime>
// the log file used to record the executing of program
// logfile.open is called at SphBase::SphBase()

std::string get_date_time_string(bool year, bool time) {
#ifdef _MSC_VER
    #pragma warning( disable : 4996 )
#endif
    time_t rawtime;
    char buffer[80];
    std::time(&rawtime);
	if(year&&time)
		std::strftime(buffer,80,"%Y-%m-%d %H.%M.%S",std::localtime(&rawtime));
	else if(year)
		std::strftime(buffer,80,"%Y-%m-%d",std::localtime(&rawtime));
	else if(time)
		std::strftime(buffer,80,"%H.%M.%S",std::localtime(&rawtime));
	else
		return std::string();
    return std::string(buffer);
}

// redirect buffer of m_clog to log file "SphBase.log"
// so, one can use m_clog to record something
SphBase::SphBase()
	: m_clog(0), m_ker_h(-1)
{
    // typical value
    m_TH.spacing_r = real_t(0.1);
    m_TH.smoothRadius_h = 2*m_TH.spacing_r;
    m_TH.spaceMin = vec_t(0);
    m_TH.spaceMax = vec_t(1);

    m_TH.dt = real_t(0.001);
    m_TH.gravity_g = vec_t::O; //m_TH.gravity_g[1] = real_t(-9.8);
    m_TH.densityErro_eta = real_t(0.01);
    m_TH.soundSpeed_cs = 100;

    m_TH.frameNumber = 0;
    m_TH.systemTime = 0;

//============================== LOG FILE ====================================
    if(!m_logFile.is_open()){
		m_logFileName += get_date_time_string() + ".log";
        m_logFile.open( m_logFileName );
	}
	if( m_logFile ){
		std::clog.rdbuf( m_logFile.rdbuf() );
		m_clog.rdbuf( m_logFile.rdbuf() );
	}else{ // LOG FILE open failed
		std::cout << "LOG FILE open failed!\n";
	}

    m_clog << "\n****************************** Start at: " << get_date_time_string() << "\n\n";

    // num_threads
    int n = omp_get_num_procs();
    omp_set_num_threads( n );
	#pragma omp parallel
	{
	if( omp_get_thread_num()==0 )
	m_clog << "num_threads: " << omp_get_num_threads() << "\n";
	}
	// dimension
	m_clog << "  dimension: " << vec_t::dim << "\n";
	// 32 or 64 bit, Debug or Release
#ifdef _WIN32
	#ifdef _WIN64
	m_clog << "_WIN64" << "\n";
	#else
	m_clog << "_WIN32" << "\n";
	#endif
#endif
#ifdef _DEBUG
	m_clog << "_DEBUG" << "\n";
#endif
	// macro
#ifdef NEGLECT_NEGATIVE_PRESSURE
	m_clog << "NEGLECT_NEGATIVE_PRESSURE" << "\n";
#endif
#ifdef PREDICTION_PCISPH
	m_clog << "PREDICTION_PCISPH" << "\n";
#endif
#ifdef WC_TIMEADAPTIVE
	m_clog << "WC_TIMEADAPTIVE" << "\n";
#elif WC_ADT
	m_clog << "WC_ADT" << "\n\n";
#endif

}

SphBase::~SphBase()
{
//============================== LOG FILE ====================================
    m_clog << "\n\ntotal sph steps: " << getFrameNumber() << '\n';
    m_clog << " total sph time: " << getSystemTime() << '\n';
    m_clog << "  error count(ZeroDis): " << getErrorCountZeroDis() << '\n';
    m_clog << "error count(ManyNeigb): " << getErrorCountManyNeigb() << '\n';
    m_clog << "      hash code of pos: " << getHashCodeOfPos() << '\n';
    m_clog << "\n****************************** End at: " << get_date_time_string() << "\n\n";

    if(m_logFile.is_open()) m_logFile.close();
	std::clog.rdbuf( std::cout.rdbuf() );
}


// print information of the scene
void SphBase::logInfoOfScene()
{
	int nsoft=0,nstatic=0;
	for(int i=0; i<getNumSolids(); ++i){
		if(m_Solids[i].type==Solid::SOFTBODY) ++nsoft;
		if(m_Solids[i].dynamic==false) ++nstatic;
	} 
	m_clog << "num fluids: " << getNumFluids() << "\n";
	m_clog << "num rigids: " << getNumSolids()
		<< " (rigid:" << getNumSolids()-nsoft << "<" << nstatic << "static>"
		<< " soft:" << nsoft << ")\n";
	m_clog << "   num fluid particles: " << getNumFluidParts() << "\n";
	m_clog << "num boundary particles: " << getNumBounParts() << "\n";
	m_clog << "   num total particles: " << getNumParts() << "\n\n";

	m_clog << "------------------------------ parameters ------------------\n";
	#define PRINT_PARAMETER_LOG(a) m_clog << #a": " << m_TH.a << "\n"
	#define PRINT_PARAMETER_LOGv(a) m_clog << #a": " << m_TH.a[0] << ' ' << m_TH.a[1] \
		<< ' ' << (vec_t::dim==3 ? m_TH.a[2] : 0) << "\n"
	PRINT_PARAMETER_LOG(spacing_r);
	PRINT_PARAMETER_LOG(smoothRadius_h);
	PRINT_PARAMETER_LOGv(spaceMin);
	PRINT_PARAMETER_LOGv(spaceMax);

	PRINT_PARAMETER_LOG(dt);
	PRINT_PARAMETER_LOGv(gravity_g);
	PRINT_PARAMETER_LOG(densityErro_eta);
	PRINT_PARAMETER_LOG(soundSpeed_cs);

	#define PRINT_PARAMETER_LOG_each(a, b) m_clog << #b": " << a[i].b << "\n"
	for(size_t i=0; i<m_Fluids.size(); ++i){
		m_clog << "------------------------------<< fluid " << i << " >>" << "\n";
		PRINT_PARAMETER_LOG_each(m_Fluids, restDensity_rho0);
		PRINT_PARAMETER_LOG_each(m_Fluids, viscosity_alpha);
		PRINT_PARAMETER_LOG_each(m_Fluids, surfaceTension_gamma);
		PRINT_PARAMETER_LOG_each(m_Fluids, fluidParticles.size());
	}
	for(size_t i=0; i<m_Solids.size(); ++i){
		m_clog << "------------------------------<< solid " << i << " >>" << "\n";
		PRINT_PARAMETER_LOG_each(m_Solids, dynamic);
		PRINT_PARAMETER_LOG_each(m_Solids, type);
		PRINT_PARAMETER_LOG_each(m_Solids, viscosity_alpha);
		PRINT_PARAMETER_LOG_each(m_Solids, boundaryParticles.size());
		PRINT_PARAMETER_LOG_each(m_Solids, innerPositions.size());
		PRINT_PARAMETER_LOG_each(m_Solids, mbtSolid_ptr);
	}

	m_clog << "------------------------------ parameters ------------------\n\n";

}


// run one step of sph, i.e. event,advance,++frameNumber,++systemTime
void SphBase::runOneStep()
{
    stepEvent();
    sphStep();
    ++ m_TH.frameNumber;
    m_TH.systemTime += m_TH.dt;
}

// compute hash code of the executing of the program, used for comparison
int SphBase::computeHashCodeOfPos()
{
    // ----------------- need be serial ----------------
    int code = 0;
    for(size_t num=m_Fluids.size(),k=0; k<num; ++k){
        const std::vector<FluidPart>& fps = m_Fluids[k].fluidParticles;
        for(size_t n=fps.size(),i=0; i<n; ++i)
            code += *((int*)&(fps[i].position[1]));
    }
    for(size_t num=m_Solids.size(),k=0; k<num; ++k){
        const std::vector<BoundPart>& bps = m_Solids[k].boundaryParticles;
        for(size_t n=bps.size(),i=0; i<n; ++i)
            code += *((int*)&(bps[i].position[1]));
    }
    m_EC.hashCodeOfPos = code;

    return code;
}

// time integration, Euler-Cromer, [IOS*14]
void SphBase::updateFluids()
{
    // for each fluid
    for(size_t n_f=m_Fluids.size(),k=0; k<n_f; ++k){
        std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
        int num = int(f_parts.size());
        // for each particle
    #pragma omp parallel for
        for(int i=0; i<num; ++i){
            FluidPart& p = f_parts[i];
		#ifdef IISPH
			p.velocity = p.vel_adv + p.acce_presure * m_TH.dt;
		#else
            p.velocity += p.acceleration * m_TH.dt;
		#endif
            p.position += p.velocity * m_TH.dt;
        }
        // as if m_SpaceMin, m_SpaceMax defined 6 walls
        //m_Fluids[k].correctOutsideParts(m_TH.spaceMin, m_TH.spaceMax, 0/*m_TH.spacing_r/30/m_TH.dt*/);
        // Remove Outside Particles
        m_Fluids[k].removeOutsideParts(m_TH.spaceMin, m_TH.spaceMax);
        // modify m_SpaceMin, m_SpaceMax
        //m_Fluids[k].adjustSpace(m_SpaceMin, m_SpaceMax);
    }

}

// compute BoundaryPart::kerSumFromSelf of earch boundary particle
void SphBase::prepareSolidPartWeight()
{
    neighbourSearch(false,true,true);

    // ----------------- need not be parallel ----------------
    real_t w0 = ker_W(0);
    // forearch solid
    for(size_t n_r=m_Solids.size(),k=0; k<n_r; ++k){
        std::vector<BoundPart>& s_parts = m_Solids[k].boundaryParticles;
        const std::vector<NeigbStr>& r_neigbs = mg_NeigbOfSolids[k];
		bool isrigidtype = m_Solids[k].type==Solid::RIGIDBODY;
        // forearch particle
        for(size_t num=s_parts.size(),i=0; i<num; ++i){
            BoundPart& p_a = s_parts[i];
			const Neigb* neigbs = r_neigbs[i].neigs; int n=r_neigbs[i].num;
			if(isrigidtype){
				real_t sum = w0, sum2 = w0;
				for(int j=0; j<n; ++j){
					if(neigbs[j].pidx.toSolidI()==k) sum += ker_W(neigbs[j].dis);
					if(!neigbs[j].pidx.isFluid()) sum2 += ker_W(neigbs[j].dis);
				}
				p_a.kerSumFromSelf = sum;
				p_a.volume = 1/sum2;
			}else{
				real_t sum2 = w0;
				for(int j=0; j<n; ++j){
					if(!neigbs[j].pidx.isFluid()) sum2 += ker_W(neigbs[j].dis);
				}
				p_a.kerSumFromSelf = w0;
				p_a.volume = 1/sum2;
			}
        } // forearch particle
    } // forearch solid

	//updateSolidPartWeight();

}

// compute volume of earch boundary particle
// i.e. using BoundaryPart::kerSumFromSelf to accelerate the computing
void SphBase::updateSolidPartWeight()
{
    // forearch solid
    for(size_t n_r=m_Solids.size(),k=0; k<n_r; ++k){
        std::vector<BoundPart>& s_parts = m_Solids[k].boundaryParticles;
        const std::vector<NeigbStr>& r_neigbs = mg_NeigbOfSolids[k];
        int num = int(s_parts.size());
        // forearch particle
    #pragma omp parallel for
        for(int i=0; i<num; ++i){
            BoundPart& p_a = s_parts[i];
            const Neigb* neigbs = r_neigbs[i].neigs; int n=r_neigbs[i].num;
            real_t volume = p_a.kerSumFromSelf;
            for(int j=0; j<n; ++j){
				if(!neigbs[j].pidx.isFluid()) volume += ker_W(neigbs[j].dis);
            }
            p_a.volume = 1/volume;
        }
    } // forearch solid

}


/*
 * ------------------------------ neighbour search ------------------------------------------------
*/

// find neighbours of each paritcle, parameters are for SphBase::prepareSolidPartWeight
void SphBase::neighbourSearch(
    bool fluidPart, bool solidPart, bool RR_same)
{
	// update internal numbers used to calculate kernel
	ker_update_numbers();

    // insert particles into the grid
    gridSetup();

    int_t n_off = mg_NS.offset.size();
    real_t h2 = m_TH.smoothRadius_h*m_TH.smoothRadius_h;
    // foreach fluid
    if(fluidPart) for(int n_f=int(m_Fluids.size()),k=0; k<n_f; ++k){
        const std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
        std::vector<NeigbStr>& f_neigbs = mg_NeigbOfFluids[k];
        int num = int(f_parts.size());
        PartIdx i, j; i.t=PartIdx::tFromFluidI(k);
        // foreach particle
	#pragma omp parallel for firstprivate(i),private(j),PARALLEL_SCHEDULE
        for(int ii=0; ii<num; ++ii){
		#ifdef WC_TIMEADAPTIVE
			FluidPart& p_a = (FluidPart&)f_parts[ii];
			if( p_a.active ){
				if( p_a.ns_factor<=(char)1 ){
					p_a.ns_factor = (char)4;
		#endif
					NeigbStr& neigbs=f_neigbs[ii]; neigbs.clear();
					const vec_t& pos_a = f_parts[ii].position;
					i.i = ii;
					for(int_t jo=mg_NS.gridOfPos(pos_a),n=0; n<n_off; ++n)
					for(j=mg_NS.gridFirst[jo+mg_NS.offset[n]]; j.valid(); j=getNextPIdx(j)){
						if( i==j ) continue;
						const vec_t& pos_b = getPosOfPIdx(j);
						//real_t len2 = (pos_a-pos_b).len_square();
						real_t len2 = pos_a.dis_square(pos_b);
						if(len2<h2) neigbs.pushBack(j, std::sqrt(len2), m_EC.tooManyNeigb);
					}
		#ifdef WC_TIMEADAPTIVE
				}else{
					--p_a.ns_factor;
					const vec_t& pos_a = f_parts[ii].position;
					Neigb* neigbs=f_neigbs[ii].neigs; int n=int(f_neigbs[ii].num);
					// forearch neighbour
					for(int jj=0; jj<n; ++jj){
						const vec_t& pos_b = getPosOfPIdx(neigbs[jj].pidx);
						real_t len2 = pos_a.dis_square(pos_b);
						neigbs[jj].dis = len2<h2 ? std::sqrt(len2) : m_TH.smoothRadius_h;
					}
				}
			}
		#endif
        }
    } // foreach fluid

    // foreach solid
    if(solidPart) for(int n_r=int(m_Solids.size()),k=0; k<n_r; ++k){
        const std::vector<BoundPart>& s_parts = m_Solids[k].boundaryParticles;
        std::vector<NeigbStr>& r_neigbs = mg_NeigbOfSolids[k];
        int num=int(s_parts.size());
        PartIdx i, j; i.t=PartIdx::tFromSolidI(k);
		bool isrigidtype = m_Solids[k].type==Solid::RIGIDBODY;
        // foreach particle
    #pragma omp parallel for firstprivate(i),private(j),PARALLEL_SCHEDULE
        for(int ii=0; ii<num; ++ii){
            NeigbStr& neigbs=r_neigbs[ii]; neigbs.clear();
            const vec_t& pos_a = s_parts[ii].position;
            i.i = ii;
            for(int_t jo=mg_NS.gridOfPos(pos_a),n=0; n<n_off; ++n)
            for(j=mg_NS.gridFirst[jo+mg_NS.offset[n]]; j.valid(); j=getNextPIdx(j)){
                if( i==j ) continue;
                if( !RR_same && isrigidtype && i.t==j.t ) continue;
                const vec_t& pos_b = getPosOfPIdx(j);
                //real_t len2 = (pos_a-pos_b).len_square();
                real_t len2 = pos_a.dis_square(pos_b);
                if(len2<h2) neigbs.pushBack(j, std::sqrt(len2), m_EC.tooManyNeigb);
            }
        }
    } // foreach solid
    
} // void SphBase::neighbourSearch( bool fluidPart, bool solidPart, bool RR_same )

// initialize the struture of space grid
void SphBase::SpaceGridStruct::initialize(
    const vec_t& spaceMin, const vec_t& spaceMax, real_t radius_h)
{
    const int gridt = 2; // 1: 3x3x3 ,  2: 5x5x5 ,  3: 7x7x7
    // grid information
    min = spaceMin - radius_h*2;
	divisor = gridt/radius_h;
    num = ((spaceMax-spaceMin+radius_h*4)*divisor).to<veci_t::type>() + 2;
    for(int i=0; i<vec_t::dim; ++i){
        factor[i] = 1;
        for(int j=0; j<i; ++j) factor[i] *= num[j];
    }

    // grid data
    gridFirst.resize(num.volume());

    // offset
    veci_t vi(-gridt);
    offset.resize(0);
    while(vi[vec_t::dim-1]<=gridt){
        offset.push_back(vi.dot(factor)); ++ vi[0];
        for(int i=0; i<vec_t::dim-1; ++i)
            if(vi[i]>gridt) vi[i]=-gridt, ++vi[i+1];
    }

} // void SphBase::gridInit(const vec_t& spaceMin, const vec_t& spaceMax, float_t radius_h)

// insert particles into the space grid
void SphBase::gridSetup()
{
    // if the domain defined by m_SpaceMin/m_SpaceMax changed, reinitialize the grid
    //static vec_t smin(1), smax(-1); static real_t r_h(-1);
    if( !(m_TH.spaceMax==mg_NS.lastSpaceMax
		&& m_TH.spaceMin==mg_NS.lastSpaceMin
		&& m_TH.smoothRadius_h==mg_NS.lastSmoothRadius_h) ){
        mg_NS.lastSpaceMax = m_TH.spaceMax;
		mg_NS.lastSpaceMin = m_TH.spaceMin;
		mg_NS.lastSmoothRadius_h = m_TH.smoothRadius_h;
        mg_NS.initialize(m_TH.spaceMin, m_TH.spaceMax, m_TH.smoothRadius_h);
    }

    int n_f = int(m_Fluids.size()), n_r = int(m_Solids.size());
    mg_NS.nextOfFluids.resize(n_f); mg_NeigbOfFluids.resize(n_f);
    for(int k=0; k<n_f; ++k){
        mg_NS.nextOfFluids[k].resize(m_Fluids[k].fluidParticles.size());
        mg_NeigbOfFluids[k].resize(m_Fluids[k].fluidParticles.size());
    }
    mg_NS.nextOfSolids.resize(n_r); mg_NeigbOfSolids.resize(n_r);
    for(int k=0; k<n_r; ++k){
        mg_NS.nextOfSolids[k].resize(m_Solids[k].boundaryParticles.size());
        mg_NeigbOfSolids[k].resize(m_Solids[k].boundaryParticles.size());
    }
    
	// initialize mg_NS.gridFirst, i.e. the tail of the next list
    PartIdx idx0(0,-1); int_t n = mg_NS.gridFirst.size();
#pragma omp parallel for schedule(dynamic,n/160)
	for(int_t i=0; i<n; ++i) mg_NS.gridFirst[i] = idx0;

// ----------------- not parallel because of dependency ----------------
	static std::vector<int_t> grid_tmp;
    // earch solid
    for(int k=0; k<n_r; ++k){
        std::vector<BoundPart>& s_parts = m_Solids[k].boundaryParticles;
        std::vector<PartIdx>& r_nex = mg_NS.nextOfSolids[k];
        PartIdx idx; idx.t=PartIdx::tFromSolidI(k);
		if(grid_tmp.size()<s_parts.size()) grid_tmp.resize(s_parts.size());
		int num = int(s_parts.size());
		// earch particle
	#pragma omp parallel for PARALLEL_SCHEDULE
		for(int i=0; i<num; ++i)
			grid_tmp[i] = mg_NS.gridOfPos(s_parts[i].position);
        // earch particle
        for(int i=0; i<num; ++i){
            int_t grid = grid_tmp[i];
            r_nex[i] = mg_NS.gridFirst[grid];
            idx.i=i; mg_NS.gridFirst[grid]=idx;
        }
    }
    // earch fluid
    for(int k=0; k<n_f; ++k){
        std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
        std::vector<PartIdx>& f_nex = mg_NS.nextOfFluids[k];
        PartIdx idx; idx.t=PartIdx::tFromFluidI(k);
		if(grid_tmp.size()<f_parts.size()) grid_tmp.resize(f_parts.size());
		int num = int(f_parts.size());
		// earch particle
	#pragma omp parallel for PARALLEL_SCHEDULE
		for(int i=0; i<num; ++i)
			grid_tmp[i] = mg_NS.gridOfPos(f_parts[i].position);
        // earch particle
        for(int i=0; i<num; ++i){
            int_t grid = grid_tmp[i];
            f_nex[i] = mg_NS.gridFirst[grid];
            idx.i=i; mg_NS.gridFirst[grid]=idx;
        }
    }

} // void SphBase::gridSetup()

