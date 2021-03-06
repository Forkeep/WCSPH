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

#ifndef PCI_SPH_H_
#define PCI_SPH_H_


#include "SphWithIO.h"


// for PCISPH, need to define macro PREDICTION_PCISPH,
// i.e. project's properties > configuration properties > C/C++ > preprocessor
#ifndef PREDICTION_PCISPH
#error Do not define macro: "PREDICTION_PCISPH"
#endif


class PciSph : public SphWithIO{

public:

	PciSph() { m_minIterations=3; m_maxIterations=10; }
	virtual ~PciSph() {}

	// override SphBase::runOneStep()
	virtual void sphStep();

	void oglDrawFluidParts_pre( void(*draw)() );

	// overwrite
	virtual void logInfoOfScene()
	{
		SphBase::logInfoOfScene();

		m_clog << "--------------------- PCISPH ---------\n";
		m_clog << "m_Pci_delta: " << m_Pci_delta << "\n";
		m_clog << "m_minIterations: " << m_minIterations << "\n";
		m_clog << "m_maxIterations: " << m_maxIterations << "\n";
		m_clog << "--------------------- PCISPH ---------\n\n";
	}

protected:
	
	real_t m_Pci_delta; // [SP09, eq.(8)]
	int m_minIterations, m_maxIterations;

	// [SP09], Algorithm 2, return iterations
	int pci_computeForce();

	real_t pci_last_dt, pci_last_sp, pci_last_h, pci_last_delta;
	// update the value of m_Pci_delta
	inline void pci_updateDelta();

	// correction value of pressure from tho_error*(t+1)
	inline real_t pci_errorFunc(const real_t& rho_error)
	{ return m_Pci_delta*rho_error; }


	// paritlce-particle interaction, [Mon92], [BT07]
	inline void pci_fluidPartForceExceptPressure_fsame(
		FluidPart& fa, const FluidPart& fb, const real_t& dis,
		const real_t& fm0, const real_t& alpha, const real_t& gamma);
	// multiphase fluid, [SP08]
	inline void pci_fluidPartForceExceptPressure_fdiff(
		FluidPart& fa, const FluidPart& fb, const real_t& dis,
		const real_t& fma, const real_t& fmb, const real_t& alpha);
	// fluid-rigid coupling, [AIS*12]
	inline void pci_fluidPartForceExceptPressure_bound(
		FluidPart& fa, const BoundPart& rb, const real_t& dis,
		const real_t& frho0, const real_t& r_alpha);
	// forces except pressure force, Fp=0, p=0
	void pci_forceExceptPressure();


	// predict velocity and position from v(t) and X(t) using a_vg and a_p
	void pci_predictPosition();


	// predict density using x*(t+1)
	void pci_predictDensity(real_t& rho_erro_max, real_t& rho_erro_avr);


	// paritlce-particle interaction, [Mon92], [BT07]
	inline void pci_fluidPartForcePressure_fsame(
		FluidPart& fa, const FluidPart& fb, const real_t& dis,
		const real_t& fm0);
	// multiphase fluid, [SP08]
	inline void pci_fluidPartForcePressure_fdiff(
		FluidPart& fa, const FluidPart& fb, const real_t& dis,
		const real_t& fma, const real_t& fmb);
	// fluid-rigid coupling, [AIS*12]
	inline void pci_fluidPartForcePressure_bound(
		FluidPart& fa, const BoundPart& rb, const real_t& dis,
		const real_t& frho0);
	// compute pressure force using cumulative pressure
	void pci_forcePressure();


	// fluid-rigid coupling, [AIS*12]
	inline void pci_boundPartForce_f(
		BoundPart& ra, const FluidPart& fb, const real_t& dis,
		const real_t& frho0, const real_t& fm0, const real_t& r_alpha);
	// fp.acceleration += fp.acce_pressure, boundary paritcles' forces
	void pci_forceTotal();

	
};


#endif // #ifndef PCI_SPH_H_


