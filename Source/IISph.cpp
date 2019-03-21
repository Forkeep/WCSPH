
#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1

#include "IISph.h"
#include <iomanip>
#include "iostream"
using namespace std;

// override SphBase::runOneStep()
void IISph::sphStep()
{
#define TIME_RECORD

#ifdef TIME_RECORD
	if( getFrameNumber()==0 ){
		m_clog.setf(std::ios::left);
		//m_clog << "Number\tSearch\t\tWeight\t\tPressure\tForce\t\tconstrainDt\tadaptiveDt\tupFluid\t\tupRigid\t\tDT\t\t\tpercentOfActive\n";
	#ifdef II_TIMEADAPTIVE
		m_clog << "Number\tTime\t\tselectActiv\tSearch\t\tWeight\t\tPressure\tForce\t\tconstrainDt\tadaptiveDt\tupFluid\t\tupRigid\t\tDT\t\t\tpercentOfActive\n";
	#else
		#ifdef II_ADT
		m_clog << "Number\tTime\t\tSearch\t\tWeight\t\tPressure\tForce\t\tadaptiveDt\tupFluid\t\tupRigid\t\tDT\n";
		#else
		m_clog << "Number\tTime\t\tSearch\t\tWeight\t\tPressure\tForce\t\tupFluid\t\tupRigid\t\tDT\n";
		#endif
	#endif
	}
	m_clog.width(7); m_clog << m_TH.frameNumber << ' ';
	m_clog.width(11); m_clog << m_TH.systemTime << ' ';

	//#define CALL_TIME(a) this->a()
	double t1, t2;
	#define CALL_TIME(a) t1=omp_get_wtime(); a; t2=omp_get_wtime(); \
						 m_clog.width(11); m_clog <<  (t2-t1)*1000 << ' '
#else
	#define CALL_TIME(a) a
#endif

#ifdef II_TIMEADAPTIVE
	CALL_TIME(selectActiveii());
#endif

	//neighbourSearch();
	CALL_TIME(neighbourSearch());

	//updateRigidPartWeight()
	CALL_TIME(updateSolidPartWeight());

	//computePressure();
	CALL_TIME(ii_computePressure());
		//computeForce();
	CALL_TIME(ii_computeForcePressure());

	//CALL_TIME(ii_adaptiveDt());

#ifdef II_TIMEADAPTIVE
	CALL_TIME(wc_updateFluids());
#else
	//updateFluids();
	CALL_TIME(updateFluids());
#endif

	//updateRigids();
	CALL_TIME(updateSolids());

	m_clog.width(11); m_clog << m_TH.dt << ' ';
#ifdef II_TIMEADAPTIVE
	m_clog.width(11); m_clog << wc_percentOfActive << ' ';
#endif

#ifdef TIME_RECORD
	m_clog << '\n';
#endif

}

// compute fluid particles' pressure, WCSPH or PCISPH
void IISph::ii_computePressure()
{
	//compute acce_adv;
	ii_forceExceptPressure();

	real_t v0 = std::pow( m_TH.spacing_r, vec_t::dim );
	// foreach fluid
	for(int n_f=int(m_Fluids.size()),k=0; k<n_f; ++k){
		std::vector<FluidPart>& f_parts = m_Fluids[k].fluidParticles;
		const std::vector<NeigbStr>& f_neigbs = mg_NeigbOfFluids[k];
		real_t rho0 = m_Fluids[k].restDensity_rho0;
		real_t fm0 = rho0 * v0;
		real_t wm0 = ker_W(0) * fm0;
		real_t h = 0.2;
		int num = int(f_parts.size());
		// first loop for foreach particle 
	#pragma omp parallel for
		for(int i=0; i<num; ++i){
			vec_t tempValue1 = vec_t::O;
			vec_t grad = vec_t::O;
			vec_t ni = vec_t::O;

			FluidPart& p_a = f_parts[i];
			real_t density = wm0;
			const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
			// forearch neighbour for computing rho
			for(int j=0; j<n; ++j){
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
						density += fm0 * ker_W(neigbs[j].dis);
					}else{ // boundary neighbour
						const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
						density += rho0*p_b.volume * ker_W(neigbs[j].dis);
					}
			}
			p_a.density = density;

			// forearch neighbour for computing dii
			for(int j=0; j<n; ++j){
				grad = vec_t::O;
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
						const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
						grad = (p_a.position-p_b.position)*(-ker_W_grad(neigbs[j].dis)/neigbs[j].dis);
						tempValue1 += grad * (-fm0/(p_a.density*p_a.density));

						//compute ni(��������)
						ni += grad*(fm0/p_b.density);
					}else{ // boundary neighbour
						const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
						grad = (p_a.position-p_b.position)*(-ker_W_grad(neigbs[j].dis)/neigbs[j].dis);
						tempValue1 += grad * (-(rho0*p_b.volume)/(p_a.density*p_a.density));
					}
			}
			p_a.dii = tempValue1 * m_TH.dt * m_TH.dt;
			ni = ni * m_TH.h;
			p_a.n = ni;

			//compute vel_adv & position_adv
			p_a.vel_adv = p_a.velocity + p_a.acce_adv * m_TH.dt;
			//p_a.linear_advec = p_a.acce_adv * m_TH.dt;
			p_a.position_adv = p_a.position + p_a.velocity * m_TH.dt;
			p_a.dis_adv = (p_a.position - p_a.position_adv);
		}

		if (m_TH.enable_vortex) {
			//1.5th loop forearch neighbour for computing vortex1 & 2
#pragma omp parallel for
			for (int i = 0; i < num; ++i) {
				vec_t grad = vec_t::O;
				FluidPart& p_a = f_parts[i];
				p_a.vortex = vec_t(0, 0, 0);
				p_a.vortex2 = vec_t(0, 0, 0);
				p_a.vortex2_adv = vec_t(0, 0, 0);
				p_a.vortex_velocity_refinement = vec_t(0, 0, 0);
				p_a.stream = vec_t(0, 0, 0);
				p_a.vortex_refinement_parameter = vec_t(0, 0, 0);

				const Neigb* neigbs = f_neigbs[i].neigs; int n = f_neigbs[i].num;
				// forearch neighbour for computing vortex 1 & 2
				for (int j = 0; j < n; ++j) {
					if (neigbs[j].pidx.isFluid()) { // fluid neighbour
						const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
						grad = (p_a.position - p_b.position)*(-ker_W_grad(neigbs[j].dis) / neigbs[j].dis);
						p_a.vortex += ((p_a.vel_adv - p_b.vel_adv)*(1 / p_a.density * fm0)).cross(grad);
						p_a.vortex2 += ((p_a.velocity - p_b.velocity)*(1 / p_a.density * fm0)).cross(grad);
					}
					else { // boundary neighbour
						const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
						grad = (p_a.position - p_b.position)*(-ker_W_grad(neigbs[j].dis) / neigbs[j].dis);
						p_a.vortex += ((p_a.vel_adv - p_a.vel_adv)*(1 / p_a.density * fm0)).cross(grad);
						p_a.vortex2 += ((p_a.velocity - p_a.velocity)*(1 / p_a.density * fm0)).cross(grad);
					}
				}
				p_a.vortex2_adv = p_a.vortex2;
			}

			//1.55th loop for refinement of vortex2
#pragma omp parallel for
			for (int i = 0; i < num; ++i) {
				vec_t grad = vec_t::O;
				FluidPart& p_a = f_parts[i];
				double radius2xy = 0;
				double radius2xz = 0;
				double radius2yz = 0;
				double threshold = 0.0000001;
				if(abs(p_a.vortex2[2])>threshold)
					radius2xy = (pow(p_a.velocity[0], 2) + pow(p_a.velocity[1], 2)) / pow(p_a.vortex2[2], 2) / 4;
				if (abs(p_a.vortex2[1])>threshold)
					radius2xz = (pow(p_a.velocity[0], 2) + pow(p_a.velocity[2], 2)) / pow(p_a.vortex2[1], 2) / 4;
				if (abs(p_a.vortex2[0])>threshold)
					radius2yz = (pow(p_a.velocity[2], 2) + pow(p_a.velocity[1], 2)) / pow(p_a.vortex2[0], 2) / 4;
				

				double xy2 = pow(p_a.dis_adv[0], 2) + pow(p_a.dis_adv[1], 2);
				double yz2 = pow(p_a.dis_adv[1], 2) + pow(p_a.dis_adv[2], 2);
				double xz2 = pow(p_a.dis_adv[0], 2) + pow(p_a.dis_adv[2], 2);

				//printf("%f %f %f %f %f %f\n", radius2xy, radius2xz, radius2yz, xy2, yz2, xz2);

				if (yz2 > 0 && abs(p_a.vortex2[0])>threshold)
					p_a.vortex_refinement_parameter[0] = sqrt(yz2 / (yz2 + radius2yz));
				if (xz2 > 0 && abs(p_a.vortex2[1])>threshold)
					p_a.vortex_refinement_parameter[1] = sqrt(xz2 / (xz2 + radius2xz));
				if (xy2 > 0 && abs(p_a.vortex2[2])>threshold)
					p_a.vortex_refinement_parameter[2] = sqrt(xy2 / (xy2 + radius2xy));

				vec_t vis_vel = p_a.vis_acce*m_TH.dt;

				/*real_t tmp2 = p_a.velocity.dot(p_a.velocity);
				real_t tmp3 = (p_a.velocity - vis_vel).dot(p_a.velocity - vis_vel);
				real_t energy_loss;
				if (tmp2 > 0.0001 && tmp3<tmp2) {
					energy_loss = tmp3 / tmp2;
					energy_loss = 1 - energy_loss;
				}
				else {
					energy_loss = 0;
				}
				p_a.energy_loss = sqrt(energy_loss);
				p_a.vis_acce = vec_t(0, 0, 0);*/

				real_t tmp2 = p_a.velocity.dot(p_a.velocity);
				real_t tmp3 = (p_a.velocity + vis_vel).dot(p_a.velocity + vis_vel);
				real_t energy_add;
				if (tmp2 > 0.0001 && tmp3 > tmp2) {
					energy_add = tmp3 / tmp2;
					energy_add = energy_add - 1;
				}
				else {
					energy_add = 0;
				}
				p_a.energy_loss = sqrt(energy_add);
				p_a.vis_acce = vec_t(0, 0, 0);

			}

			//1.6th loop forearch particle compute stream function
#pragma omp parallel for
			for (int i = 0; i < num; ++i) {
				//vec_t grad = vec_t::O;
				FluidPart& p_a = f_parts[i];
				

				const Neigb* neigbs = f_neigbs[i].neigs; int n = f_neigbs[i].num;
				for (int j = 0; j < n; ++j) {
					if (neigbs[j].pidx.isFluid()) { // fluid neighbour
						const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
						/*p_a.stream += ((p_b.vortex2_adv - p_b.vortex)*m_TH.particle_volume) / (4 * 3.1415926*neigbs[j].dis);*/
						//p_a.stream += ((p_b.vortex2*0.016)*m_TH.particle_volume) / (4 * 3.1415926*neigbs[j].dis);
						real_t alpha = 0.6;
						p_a.stream += (
							(p_b.vortex2*p_b.energy_loss*alpha)
							*m_TH.particle_volume) / (4 * 3.1415926*neigbs[j].dis);
						


						//refinement = 0;

						/*real_t refinement = p_a.velocity.length() / (m_TH.spacing_r / m_TH.dt / 100);
						p_a.stream += (
							(p_b.vortex2*p_b.vortex_refinement_parameter / (1 + refinement) *1)
							*m_TH.particle_volume) / (4 * 3.1415926*neigbs[j].dis);*/


						/*if (refinement < 1 ) {
							p_a.stream += (
								(p_b.vortex2*p_b.vortex_refinement_parameter /(1+ refinement) * 0.1)
								*m_TH.particle_volume) / (4 * 3.1415926*neigbs[j].dis);
						}*/
					}
				}
			}

			//1.7th loop forearch particle compute vortex_velocity_refinement
#pragma omp parallel for
			for (int i = 0; i < num; ++i) {
				vec_t grad = vec_t::O;
				FluidPart& p_a = f_parts[i];
				const Neigb* neigbs = f_neigbs[i].neigs; int n = f_neigbs[i].num;
				for (int j = 0; j < n; ++j) {
					if (neigbs[j].pidx.isFluid()) { // fluid neighbour
						const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
						grad = (p_a.position - p_b.position)*(-ker_W_grad(neigbs[j].dis) / neigbs[j].dis);
						p_a.vortex_velocity_refinement += ((p_a.stream - p_b.stream)*(1 / p_a.density * fm0)).cross(grad);
					}
					else {
						/*p_a.vortex_velocity_refinement = vec_t(0, 0, 0);
						break;*/
					}
				}
				if (p_a.density > 200)
					p_a.vel_adv += p_a.vortex_velocity_refinement;
					;
			}


		}

		// second loop for foreach particle 
	#pragma omp parallel for
		for(int i=0; i<num; ++i){
			real_t tempValue1 = 0;
			real_t tempValue2 = 0;
			vec_t grad = vec_t::O;
			vec_t dji = vec_t::O;

			FluidPart& p_a = f_parts[i];
			const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
			// forearch neighbour for computing rho_adv,aii
			for(int j=0; j<n; ++j){
				grad = vec_t::O;
				dji = vec_t::O;
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
						const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
						grad = (p_a.position-p_b.position)*(-ker_W_grad(neigbs[j].dis)/neigbs[j].dis);
						//compute rho_adv
						tempValue1 += grad.dot((p_a.vel_adv -p_b.vel_adv))*fm0*m_TH.dt;
						//compute aii
						dji = grad * (fm0*m_TH.dt*m_TH.dt/(p_a.density*p_a.density));
						tempValue2 += (p_a.dii - dji).dot(grad) * fm0;
						
					}else{ // boundary neighbour
						const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
						grad = (p_a.position-p_b.position)*(-ker_W_grad(neigbs[j].dis)/neigbs[j].dis);
						//compute rho_adv
						tempValue1 += grad.dot((p_a.velocity-p_b.velocity))*(rho0*p_b.volume)*m_TH.dt;
						//compute aii
						tempValue2 += (p_a.dii).dot(grad) * (rho0*p_b.volume);
					}
			}
			if(tempValue2 == 0){tempValue2 = 1;}
			p_a.aii = tempValue2;
			p_a.rho_adv = p_a.density + tempValue1;
			p_a.p_l = (0.5*p_a.presure);
		}

		//pressure slover
		int l = 0;
		real_t averageDensityError = 1e10;
		real_t sumDensities = 0;
		long nbParticlesInSummation = 0;
		double eta = 0.001*rho0;
		while((averageDensityError > eta) || l<2){
			
			// third loop for foreach particle
		#pragma omp parallel for
			for(int i=0; i<num; ++i){
				vec_t tempValue1 = vec_t::O;
				vec_t grad = vec_t::O;

				FluidPart& p_a = f_parts[i];
				const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
				// forearch neighbour 
				for(int j=0; j<n; ++j){
					grad = vec_t::O;
					if(neigbs[j].pidx.isFluid()){ // fluid neighbour
							const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
							grad = (p_a.position-p_b.position)*(-ker_W_grad(neigbs[j].dis)/neigbs[j].dis);
							tempValue1 += grad * (-(fm0*p_b.p_l)/(p_b.density*p_b.density));
						}else{ // boundary neighbour
							const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
							//??????????
						}
				}
				p_a.sum_dijpj = tempValue1 * (m_TH.dt*m_TH.dt);
			}

			// fourth loop for foreach particle
		#pragma omp parallel for
			for(int i=0; i<num; ++i){
				vec_t grad = vec_t::O;
				vec_t tempValue1 = vec_t::O;
				real_t finalterm = 0;
				real_t t = 0;

				FluidPart& p_a = f_parts[i];
				const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
				// forearch neighbour 
				for(int j=0; j<n; ++j){
					vec_t dji = vec_t::O;
					grad = vec_t::O;
					if(neigbs[j].pidx.isFluid()){ // fluid neighbour
							const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
							//compute dji
							grad = (p_a.position-p_b.position)*(-ker_W_grad(neigbs[j].dis)/neigbs[j].dis);
							dji = grad*m_TH.dt*m_TH.dt*fm0/(p_a.density*p_a.density);

							 // compute sum(djk*pk)-dji*pi
							tempValue1 = (p_a.sum_dijpj - p_b.dii*p_b.p_l)-(p_b.sum_dijpj-dji*p_a.p_l);
							finalterm =finalterm+(fm0 * tempValue1.dot(grad));

						}else{ // boundary neighbour
							const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
							grad = (p_a.position-p_b.position)*(-ker_W_grad(neigbs[j].dis)/neigbs[j].dis);
							finalterm=finalterm + (rho0*p_b.volume * (p_a.sum_dijpj).dot(grad));
						}
				}
				t = 0.5*p_a.p_l + (0.5/p_a.aii) * (rho0 - p_a.rho_adv - finalterm);
				if(t < 0){
					t = 0;
				}else{
					sumDensities += p_a.rho_adv + p_a.aii*p_a.p_l + finalterm;
					++nbParticlesInSummation;
				}
				p_a.p_l = t;
				p_a.presure = t;
				l = l+1;
			}
			if(nbParticlesInSummation > 0){
				averageDensityError = (sumDensities / num) -rho0;
			}else{
				 averageDensityError = 0.0;
			}
			
		}
	}
}


// paritlce-particle interaction, [Mon92], [BT07]
inline void IISph::ii_fluidPartForceExceptPressure_fsame(
	FluidPart& fa, const FluidPart& fb, const real_t& dis,
	const real_t& fm0, const real_t& alpha, const real_t& gamma)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position-fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	vec_t grad2 = (fa.position - fb.position)*(-ker_W_grad(dis) / dis);
	real_t acce = 0;
	// viscosity
	real_t pro = (fa.velocity-fb.velocity).dot( xab );
	vec_t tmp = vec_t(0, 0, 0);

	tmp = grad2 * pro / (dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h) * fm0 / fb.density * 10 * alpha;

	fa.vis_acce += tmp;
	fa.acce_adv += tmp;

	// surface tension
	// ...

}

// multiphase fluid, [SP08]
inline void IISph::ii_fluidPartForceExceptPressure_fdiff(
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
	if(acce){ xab *= acce; fa.acce_adv += xab; }

}

// fluid-rigid coupling, [AIS*12]
inline void IISph::ii_fluidPartForceExceptPressure_bound(
	FluidPart& fa, const BoundPart& rb, const real_t& dis,
	const real_t& frho0, const real_t& r_alpha)
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position - rb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	vec_t grad2 = (fa.position - rb.position)*(-ker_W_grad(dis) / dis);
	real_t acce = 0;
	// viscosity
	real_t pro = (fa.velocity-rb.velocity).dot( xab );
	/*if( pro<0 ){
		real_t nu = 2*r_alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fa.density*2 );
		real_t pi = -nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		acce += grad * (- frho0*rb.volume * pi );
	}
	if(acce){ xab *= acce; fa.acce_adv += xab; }*/
	vec_t tmp = grad2 * pro / (dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h) * (frho0*rb.volume) / fa.density * 10 * r_alpha;
	fa.vis_acce += tmp;
	fa.acce_adv += tmp;
}

void IISph::ii_forceExceptPressure(){
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
			p_a.presure = 0; 
			p_a.acce_adv = m_TH.gravity_g;
			//p_a.acce_adv = vec_t::O;
			const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
			// forearch neighbour
			for(int j=0; j<n; ++j){
				if(neigbs[j].pidx.isFluid()){ // fluid neighbour
					const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
					int idx_b = neigbs[j].pidx.toFluidI();
					if(idx_b==k){
						// the same fluid
						ii_fluidPartForceExceptPressure_fsame(
							p_a, p_b, neigbs[j].dis, fm0, alpha, gamma);
					}else{
						// different fluid
						real_t fmb = v0 * m_Fluids[idx_b].restDensity_rho0;
						real_t b_alpha = m_Fluids[idx_b].viscosity_alpha;
						ii_fluidPartForceExceptPressure_fdiff(
							p_a, p_b, neigbs[j].dis, fm0, fmb, (alpha+b_alpha)/2);

					}
				}else{ // boundary neighbour
					const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
					int idx_b = neigbs[j].pidx.toSolidI();
					real_t r_alpha = m_Solids[idx_b].viscosity_alpha;
					ii_fluidPartForceExceptPressure_bound(
						p_a, p_b, neigbs[j].dis, rho0, r_alpha);
				}
			}
		}
	}
}

// paritlce-particle interaction, [Mon92], [BT07]
inline void IISph::ii_fluidPartForcePressure_fsame(
	FluidPart& fa, const FluidPart& fb, const real_t& dis,
	const real_t& fm0, const real_t& alpha, const real_t& gamma )
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position-fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	// momentum
	real_t acce = grad
		* ( -fm0 * (fa.presure/(fa.density*fa.density) + fb.presure/(fb.density*fb.density)) );
	// viscosity
	/*real_t pro = (fa.velocity-fb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fa.density + fb.density );
		real_t pi =
			-nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		acce += grad * ( - fm0 * pi );
	}*/
	xab *= acce; fa.acce_presure += xab;
	
	// surface tension
	//real_t h = 0.2;
	//real_t r = 1;
	////compute cohesion
	//vec_t xab_coh = fa.position-fb.position;
	//double pai = 3.14;
	//real_t Kij = 2*1000/(fa.density+fb.density);
	//real_t C = 0;
	//if(2*dis>m_TH.h && dis<=m_TH.h){
	//	C = (32/(pai*pow(m_TH.h,9)))*pow(m_TH.h-dis,3)*pow(dis,3);
	//}else if(dis>0 && 2*dis<=m_TH.h){
	//	C = (32/(pai*pow(m_TH.h,9)))*(2*pow(m_TH.h-dis,3)*pow(dis,3)-pow(m_TH.h,6)/64);
	//}else{
	//	C = 0;
	//}
	//real_t acce_coh = -r*fm0*C/dis;
	//xab_coh *= acce_coh;

	////compute curvature
	//real_t acce_cur = -m_TH.r;
	//vec_t xab_cur = fa.n-fb.n;
	//xab_cur *= acce_cur;

	//vec_t acce_st = (xab_coh+xab_cur)*Kij;
	//fa.acce_presure += acce_st;
}
// multiphase fluid, [SP08]
inline void IISph::ii_fluidPartForcePressure_fdiff(
	FluidPart& fa, const FluidPart& fb, const real_t& dis,
	const real_t& fma, const real_t& fmb, const real_t& alpha )
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position-fb.position;
	real_t grad = -ker_W_grad(dis)/dis;
	real_t dalta_a=fa.density/fma, dalta_b=fb.density/fmb;
	// momentum
	real_t acce = grad
		* ( - 1/fma * (fa.presure/(dalta_a*dalta_a) + fb.presure/(dalta_b*dalta_b)) );
	// viscosity
	real_t pro = (fa.velocity-fb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fa.density + fb.density );
		real_t pi =
			-nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		acce += grad * ( - (fma+fmb)/2 * pi );
	}
	xab *= acce; fa.acce_presure += xab;


}
// fluid-rigid coupling, [AIS*12]
inline void IISph::ii_fluidPartForcePressure_bound(
	FluidPart& fa, const BoundPart& rb, const real_t& dis,
	const real_t& frho0, const real_t& r_alpha )
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }
	vec_t xab = fa.position - rb.position;
	real_t grad = -ker_W_grad(dis)/dis, acce=0;
	// momentum
	if(fa.presure>0)
		acce += grad * ( - frho0*rb.volume * (fa.presure/(fa.density*fa.density)*2) );
	// viscosity
	/*real_t pro = (fa.velocity-rb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*r_alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fa.density*2 );
		real_t pi =
			-nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		acce += grad * (- frho0*rb.volume * pi );
	}*/
	xab *= acce; fa.acce_presure += xab;
	
	//// surface tension & adhesion
	//real_t bt = 1;
	//real_t h = 0.2;
	//real_t A = 0;
	//vec_t xab_adh = fa.position - rb.position;

	//if(2*dis>m_TH.h && dis<=m_TH.h){
	//	A = 0.007/(pow(m_TH.h,13/4))*pow((-4*dis*dis/m_TH.h+6*dis-2*m_TH.h),1/4);
	//}else{
	//	A = 0;
	//}
	//real_t acce_adh = -m_TH.bt*frho0*rb.volume*A/dis;
	//xab_adh *= acce_adh;
	//fa.acce_presure += xab_adh;
}
// fluid-rigid coupling, [AIS*12]
inline void IISph::ii_boundPartForcePressure_f(
	BoundPart& ra, const FluidPart& fb, const real_t& dis,
	const real_t& frho0, const real_t& fm0, const real_t& r_alpha )
{
	if( dis==0 ){ ++m_EC.zeroDis; return; }

	vec_t xab = ra.position - fb.position;
	real_t grad = -ker_W_grad(dis)/dis, force=0;
	// momentum
	if(fb.presure>0)
		force += grad * ( - fm0*frho0*ra.volume * (fb.presure/(fb.density*fb.density)*2) );
	// viscosity
	real_t pro = (ra.velocity-fb.velocity).dot( xab );
	if( pro<0 ){
		real_t nu = 2*r_alpha*m_TH.smoothRadius_h*m_TH.soundSpeed_cs / ( fb.density*2 );
		real_t pi =
			-nu * pro / ( dis*dis + real_t(0.01)*m_TH.smoothRadius_h*m_TH.smoothRadius_h );
		force += grad * (- fm0*frho0*ra.volume * pi );
	}
	xab *= force; ra.force += xab;

}
// compute gravity, pressure and friction force
void IISph::ii_computeForcePressure()
{
	real_t v0 = std::pow( m_TH.spacing_r, vec_t::dim );
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
				p_a.acce_presure = vec_t::O;
				const Neigb* neigbs=f_neigbs[i].neigs; int n=f_neigbs[i].num;
				// forearch neighbour
				for(int j=0; j<n; ++j){
					if(neigbs[j].pidx.isFluid()){ // fluid neighbour
						const FluidPart& p_b = getFluidPartOfIdx(neigbs[j].pidx);
						int idx_b = neigbs[j].pidx.toFluidI();
						// the same fluid
						if(idx_b==k){
							ii_fluidPartForcePressure_fsame(p_a, p_b, neigbs[j].dis, fm0, alpha, gamma);
						// different fluid
						}else{
							real_t fmb = v0 * m_Fluids[idx_b].restDensity_rho0;
							real_t b_alpha = m_Fluids[idx_b].viscosity_alpha;
							ii_fluidPartForcePressure_fdiff(
								p_a, p_b, neigbs[j].dis, fm0, fmb, (alpha+b_alpha)/2);
						}
					}else{ // boundary neighbour
						const BoundPart& p_b = getBoundPartOfIdx(neigbs[j].pidx);
						int idx_b = neigbs[j].pidx.toSolidI();
						real_t r_alpha = m_Solids[idx_b].viscosity_alpha;
						ii_fluidPartForcePressure_bound(p_a, p_b, neigbs[j].dis, rho0, r_alpha);
					}
				}//nneighbour
			#ifdef II_ADT
				// the h in the paper is r, and the smoothing_h is 2h
				real_t v = !(p_a.velocity==vec_t::O)
					? m_Lambda_v*(m_TH.spacing_r/p_a.velocity.length()) : 1;
				real_t f = !(p_a.acceleration==vec_t::O)
					? m_Lambda_f*std::sqrt(m_TH.spacing_r/p_a.acceleration.length()) : 1;
				p_a.dt = std::min(v, f);
			#endif
		}//particle
	}//fluid

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
					ii_boundPartForcePressure_f(p_a, p_b, neigbs[j].dis, frho0, frho0*v0, alpha);
				}
			}//neighbour
		}//particle
	}//rigid

}

inline void IISph::fluidToSolid(int fluidIdx, const char* meshFileName,
	const char* sampleFileName, real_t viscosity, const glm::mat4& transform) {

}