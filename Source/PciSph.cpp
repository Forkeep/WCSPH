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

#include "PciSph.h"
#include <sstream>
#include <iomanip> 
#include "GL/glew.h"


#define TIME_RECORD
#define ITER_RECORD

#ifdef ITER_RECORD
static std::vector<real_t> errors_max, errors_avr;
#endif

void PciSph::sphStep()
{
#if defined(TIME_RECORD) || defined(ITER_RECORD)
	if( getFrameNumber()==0 ){
		m_clog.setf(std::ios::left);
		m_clog << "Number\tTime\t\t";
	#ifdef TIME_RECORD
		m_clog << "Search\t\tWeight\t\tForce\t\tupFluid\t\tupRigid\t\tnIter\t";
	#endif
	#ifdef ITER_RECORD
		#ifndef TIME_RECORD
			m_clog << "nIter\t";
		#endif
		for(int i=1; i<=10; ++i)
			m_clog << i << "(%):\t\t\t\t\t\t";
		m_clog << "...(%)";
	#endif
		m_clog << "\n";
	}
	m_clog.width(7); m_clog << m_TH.frameNumber << ' ';
	m_clog.width(11); m_clog << m_TH.systemTime << ' ';
#endif

#ifdef TIME_RECORD
	double t1, t2;
	#define CALL_TIME(a) \
		t1=omp_get_wtime(); a(); t2=omp_get_wtime(); \
		m_clog << std::setw(11) << (t2-t1)*1000 << ' '
#else
	#define CALL_TIME(a) a()
#endif

	//neighbourSearch();
	CALL_TIME(neighbourSearch);

	//updateRigidPartWeight();
	CALL_TIME(updateSolidPartWeight);

	//int niter = pci_computeForce();
	CALL_TIME(int niter = pci_computeForce);

	//updateFluids();
	CALL_TIME(updateFluids);

	//updateRigids();
	CALL_TIME(updateSolids);

#ifdef TIME_RECORD
	m_clog.width(7); m_clog << niter << ' ';
#endif
#ifdef ITER_RECORD
	#ifndef TIME_RECORD
		m_clog << std::setw(7) << errors_max.size() << ' ';
	#endif
	for(int n=int(errors_max.size()),i=0; i<n; ++i){
		m_clog << std::setw(11) << errors_max[i]*100 << ' ';
		m_clog << std::setw(15) << errors_avr[i]*100 << ' ';
	}
#endif
#if defined(TIME_RECORD) || defined(ITER_RECORD)
	m_clog << "\n";
#endif

}


// [SP09], Algorithm 2, return iterations
int PciSph::pci_computeForce()
{

	// update the value of m_Pci_delta
	pci_updateDelta();

	// forces except pressure force, Fp=0, p=0
	pci_forceExceptPressure();

#ifdef ITER_RECORD
	errors_max.resize(0);
	errors_avr.resize(0);
#endif

	// prediction-correction iterations
	real_t rho_error, rho_err_avr; int niter=0;
	do{
		pci_predictPosition(); // p.pos_prediction
		pci_predictDensity(rho_error, rho_err_avr); // p.density, p.pressure
		pci_forcePressure(); // p.acce_pressure
		++ niter;
#ifdef ITER_RECORD
		errors_max.push_back(rho_error);
		errors_avr.push_back(rho_err_avr);
#endif
	}while( (rho_error>=m_TH.densityErro_eta || niter<m_minIterations) && niter<m_maxIterations );

	// fp.acceleration += fp.acce_pressure, boundary paritcles' forces
	pci_forceTotal();

	return niter;
}


// update the value of m_Pci_delta
inline void PciSph::pci_updateDelta()
{
	//static real_t pci_dt=0, pci_sp=0, pci_h=0, pci_delta=0;
	if( !(pci_last_dt==m_TH.dt) ){
		pci_last_dt = m_TH.dt;
		if( !(pci_last_sp==m_TH.spacing_r && pci_last_h==m_TH.smoothRadius_h) ){
			pci_last_sp = m_TH.spacing_r; pci_last_h = m_TH.smoothRadius_h;
			int n = (int)std::ceil(m_TH.smoothRadius_h/m_TH.spacing_r);
			// prototype particle from uniform grid
			vec_t grad_sum(0); real_t grad2_sum(0); veci_t vi(-n);
			while(vi[vec_t::dim-1]<=n){
				if(!(vi==veci_t::O)){
					vec_t vf = vi.to<real_t>()*m_TH.spacing_r;
					real_t dis = vf.length();
					if(dis<m_TH.smoothRadius_h){
						real_t grad = -ker_W_grad(dis);
						grad_sum += vf * (grad/dis);
						grad2_sum += grad*grad;
					}
				}
				++vi[0]; for(int i=0; i<vec_t::dim-1; ++i)if(vi[i]>n)vi[i]=-n,++vi[i+1];
			}
			real_t v0 = std::pow(m_TH.spacing_r, vec_t::dim);
			real_t beta = 2 * v0*v0;
			pci_last_delta = -1 / ( beta * ( -grad_sum.dot(grad_sum) - grad2_sum) );
		}
		m_Pci_delta = pci_last_delta / (m_TH.dt*m_TH.dt);
	}

}


// paritlce-particle interaction, [Mon92], [BT07]
inline void PciSph::pci_fluidPartForceExceptPressure_fsame(
	FluidPart& fa, const FluidPart& fb, const real_t& dis,
	const real_t& fm0, const real_t& alpha, const real_t& gamma)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position-fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	real_t acce = 0;
	// viscosity
	real_t pro = (fa.velocity-fb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fa.density + fb.density );
		real_t pi = -nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		acce += grad * ( - fm0 * pi );
	}
	if(acce){ xab *= acce; fa.acceleration += xab; }
	// surface tension
	// ...

}

// multiphase fluid, [SP08]
inline void PciSph::pci_fluidPartForceExceptPressure_fdiff(
	FluidPart& fa, const FluidPart& fb, const real_t& dis,
	const real_t& fma, const real_t& fmb, const real_t& alpha)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position-fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	real_t acce = 0;
	// viscosity
	real_t pro = (fa.velocity-fb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fa.density + fb.density );
		real_t pi = -nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		acce += grad * ( - (fma+fmb)/2 * pi );
	}
	if(acce){ xab *= acce; fa.acceleration += xab; }

}

// fluid-rigid coupling, [AIS*12]
inline void PciSph::pci_fluidPartForceExceptPressure_bound(
	FluidPart& fa, const BoundPart& rb, const real_t& dis,
	const real_t& frho0, const real_t& r_alpha)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position - rb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	real_t acce = 0;
	// viscosity
	real_t pro = (fa.velocity-rb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*r_alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fa.density*2 );
		real_t pi = -nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		acce += grad * (- frho0*rb.volume * pi );
	}
	if(acce){ xab *= acce; fa.acceleration += xab; }

}

// forces except pressure force, Fp=0, p=0
void PciSph::pci_forceExceptPressure()
{
	real_t v0 = std::pow(m_TH.spacing_r, vec_t::dim);
	// foreach fluid
	for(int n_f=int(m_Fluids.size()),k=0; k<n_f; ++k){
		std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
		const std::vector<NeigbStr>& f_neigbs = mg_NeigbOfFluids[k];
		real_t rho0 = m_Fluids[k].restDensity_rho0;
		real_t fm0 = rho0 * v0;
		real_t alpha = m_Fluids[k].viscosity_alpha;
		real_t gamma = m_Fluids[k].surfaceTension_gamma;
		int num = int(f_parts.size());
		// foreach particle
	#pragma omp parallel for
		for(int i=0; i<num; ++i){
			FluidPart& p_a = f_parts[i];
			p_a.presure = 0; p_a.acce_pressure = vec_t::O;
			p_a.acceleration = m_TH.gravity_g;
			const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
			// forearch neighbour
			for(int j=0; j<n; ++j){
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
					const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
					int idx_b = neigbs[j].pidx.toFluidI();
					if(idx_b==k){
						// the same fluid
						pci_fluidPartForceExceptPressure_fsame(
							p_a, p_b, neigbs[j].dis, fm0, alpha, gamma);
					}else{
						// different fluid
						real_t fmb = v0 * m_Fluids[idx_b].restDensity_rho0;
						real_t b_alpha = m_Fluids[idx_b].viscosity_alpha;
						pci_fluidPartForceExceptPressure_fdiff(
							p_a, p_b, neigbs[j].dis, fm0, fmb, (alpha+b_alpha)/2);
					}
				}else{ // boundary neighbour
					const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
					int idx_b = neigbs[j].pidx.toSolidI();
					real_t r_alpha = m_Solids[idx_b].viscosity_alpha;
					pci_fluidPartForceExceptPressure_bound(
						p_a, p_b, neigbs[j].dis, rho0, r_alpha);
				}
			}//neighbour
		}//particle
	}//fluid

}


// predict velocity and position from v(t) and X(t) using a_vg and a_p
void PciSph::pci_predictPosition()
{
	// for each fluid
	for(int n_f=int(m_Fluids.size()),k=0; k<n_f; ++k){
		std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
		int num = int(f_parts.size()); vec_t pos_pre;
		// for each particle
	#pragma omp parallel for private(pos_pre)
		for(int i=0; i<num; ++i){
			FluidPart& p = f_parts[i];
			//p.pos_prediction = p.position + (p.velocity + (p.acceleration+p.acce_pressure)*m_TH.dt) * m_TH.dt;
			pos_pre = p.acceleration; pos_pre += p.acce_pressure;
			pos_pre *= m_TH.dt; pos_pre += p.velocity; pos_pre *= m_TH.dt;
			pos_pre += p.position;
			p.pos_prediction = pos_pre;
		}
	}

}


// predict density using x*(t+1)
void PciSph::pci_predictDensity(real_t& rho_erro_max, real_t& rho_erro_avr)
{
	real_t v0 = std::pow( m_TH.spacing_r, vec_t::dim );
	real_t error_max=0; double error_avr=0; int error_avr_n=0;
	// foreach fluid
	for(int n_f=int(m_Fluids.size()),k=0; k<n_f; ++k){
		std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
		std::vector<NeigbStr>& f_neigbs = mg_NeigbOfFluids[k];
		real_t rho0 = m_Fluids[k].restDensity_rho0;
		real_t fm0 = rho0 * v0;
		real_t wm0 = ker_W(0) * fm0;
		real_t den_error = 0;
		int num = int(f_parts.size());
		// foreach particle
	#pragma omp parallel for reduction(+:error_avr,error_avr_n)
		for(int i=0; i<num; ++i){
			FluidPart& p_a = f_parts[i];
			real_t density = wm0;
			Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
			// forearch neighbour
			for(int j=0; j<n; ++j){
				// distance need be recomputed
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
					const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
					//int idx_b = neigbs[j].pidx.toFluidI();
					real_t dis = p_a.pos_prediction.distance(p_b.pos_prediction);
					if(dis<m_TH.smoothRadius_h) density += fm0 * ker_W(dis);
				}else{ // boundary neighbour
					const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
					//int idx_b = neigbs[j].pidx.toRigidI();
					real_t dis = p_a.pos_prediction.distance(p_b.position);
					if(dis<m_TH.smoothRadius_h) density += rho0*p_b.volume * ker_W(dis);
				}
			}//neighbour
			p_a.density = density;
			real_t thoerr = density-rho0;
			// density error, clamp negative pressure to 0
		#ifdef NEGLECT_NEGATIVE_PRESSURE
			if(thoerr<0) thoerr=0;
		#endif
			p_a.presure += pci_errorFunc(thoerr);
			thoerr = std::abs(thoerr);
		#pragma omp critical
			den_error = std::max(den_error, thoerr);
			error_avr += thoerr/rho0; ++error_avr_n;
		}//particle
		error_max = std::max(error_max, den_error/rho0);
	}//fluid

	rho_erro_max = error_max;
	rho_erro_avr = real_t(error_avr/error_avr_n);

}


// paritlce-particle interaction, [Mon92], [BT07]
inline void PciSph::pci_fluidPartForcePressure_fsame(
	FluidPart& fa, const FluidPart& fb, const real_t& dis,
	const real_t& fm0)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position-fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	// momentum
	real_t acce = grad * ( -fm0 * (fa.presure/(fa.density*fa.density) + fb.presure/(fb.density*fb.density)) );
	xab *= acce; fa.acce_pressure += xab;

}

// multiphase fluid, [SP08]
inline void PciSph::pci_fluidPartForcePressure_fdiff(
	FluidPart& fa, const FluidPart& fb, const real_t& dis,
	const real_t& fma, const real_t& fmb)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position-fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	real_t dalta_a=fa.density/fma, dalta_b=fb.density/fmb;
	// momentum
	real_t acce = grad * ( - 1/fma * (fa.presure/(dalta_a*dalta_a) + fb.presure/(dalta_b*dalta_b)) );
	xab *= acce; fa.acce_pressure += xab;

}

// fluid-rigid coupling, [AIS*12]
inline void PciSph::pci_fluidPartForcePressure_bound(
	FluidPart& fa, const BoundPart& rb, const real_t& dis,
	const real_t& frho0)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position - rb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	real_t acce = 0;
	// momentum
	if(fa.presure>0)
		acce += grad * ( - frho0*rb.volume * (2*fa.presure/(fa.density*fa.density)) );
	if(acce){ xab *= acce; fa.acce_pressure += xab; }

}

// compute pressure force using cumulative pressure
void PciSph::pci_forcePressure()
{
	real_t v0 = std::pow( m_TH.spacing_r, vec_t::dim );
	// foreach fluid
	for(int n_f=int(m_Fluids.size()),k=0; k<n_f; ++k){
		std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
		const std::vector<NeigbStr>& f_neigbs = mg_NeigbOfFluids[k];
		real_t rho0 = m_Fluids[k].restDensity_rho0;
		real_t fm0 = rho0 * v0;
		int num = int(f_parts.size());
		// foreach particle
	#pragma omp parallel for
		for(int i=0; i<num; ++i){
			FluidPart& p_a = f_parts[i];
			p_a.acce_pressure = vec_t::O;
			const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
			// forearch neighbour
			for(int j=0; j<n; ++j){
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
					const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
					int idx_b = neigbs[j].pidx.toFluidI();
					if(idx_b==k){
						// the same fluid
						pci_fluidPartForcePressure_fsame(p_a, p_b, neigbs[j].dis, fm0);
					}else{
						// different fluid
						real_t fmb = v0 * m_Fluids[idx_b].restDensity_rho0;
						pci_fluidPartForcePressure_fdiff(p_a, p_b, neigbs[j].dis, fm0, fmb);
					}
				}else{ // boundary neighbour
					const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
					//int idx_b = neigbs[j].pidx.toRigidI;
					pci_fluidPartForcePressure_bound(p_a, p_b, neigbs[j].dis, rho0);
				}
			}//neighbour
		}//particle
	}//fluid

}


// fluid-rigid coupling, [AIS*12]
inline void PciSph::pci_boundPartForce_f(
	BoundPart& ra, const FluidPart& fb, const real_t& dis,
	const real_t& frho0, const real_t& fm0, const real_t& r_alpha )
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = ra.position - fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	real_t force=0;
	// momentum
	if(fb.presure>0)
		force += grad * ( - fm0*frho0*ra.volume * (2*fb.presure/(fb.density*fb.density)) );
	// viscosity
	real_t pro = (ra.velocity-fb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*r_alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fb.density*2 );
		real_t pi = -nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		force += grad * (- fm0*frho0*ra.volume * pi );
	}
	if(force){ xab *= force; ra.force += xab; }

}

// fp.acceleration += fp.acce_pressure, boundary paritcles' forces
void PciSph::pci_forceTotal()
{
	// foreach fluid
	for(int n_f=int(m_Fluids.size()),k=0; k<n_f; ++k){
		std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
		int num = int(f_parts.size());
		// foreach particle
	#pragma omp parallel for
		for(int i=0; i<num; ++i){
			FluidPart& p = f_parts[i];
			p.acceleration += p.acce_pressure;
		}//particle
	}//fluid

	real_t v0 = std::pow( m_TH.spacing_r, vec_t::dim );
	// forearch rigid
	for(int n_r=int(m_Solids.size()),k=0; k<n_r; ++k) if(m_Solids[k].dynamic){
		std::vector<BoundPart>& r_parts = m_Solids[k].boundaryParticles;
		const std::vector<NeigbStr>& r_neigbs = mg_NeigbOfSolids[k];
		real_t alpha = m_Solids[k].viscosity_alpha;
		int num = int(r_parts.size());
		// foreach particle
	#pragma omp parallel for
		for(int i=0; i<num; ++i){
			BoundPart& p_a = r_parts[i];
			p_a.force = vec_t::O;
			const Neigb* neigbs=r_neigbs[i].neigs; int n=r_neigbs[i].num;
			// forearch neighbour
			for(int j=0; j<n; ++j){
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
					const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
					int idx_b = neigbs[j].pidx.toFluidI();
					real_t frho0 = m_Fluids[idx_b].restDensity_rho0;
					pci_boundPartForce_f(p_a, p_b, neigbs[j].dis, frho0, frho0*v0, alpha);
				}
			}//neighbour
		}//particle
	}//rigid

}


// draw fluid particles or boundary paricles of rigids in OpenGL, iteratively
void PciSph::oglDrawFluidParts_pre( void(*draw)() )
{
	for(int num=int(m_Fluids.size()),k=0;k<num;++k){
		const std::vector<FluidParticle>& f_pts = m_Fluids[k].fluidParticles;
		for(int n=int(f_pts.size()),i=0; i<n; ++i){
			const vec_t& p = f_pts[i].pos_prediction;
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, f_pts[i].color);
			glPushMatrix();
			if(vec_t::dim==3)
				glTranslatef( (float)p[0], (float)p[1], (float)p[2] );
			else glTranslatef( (float)p[0], (float)p[1], 0 );
			draw();
			glPopMatrix();
		}
	}
}



