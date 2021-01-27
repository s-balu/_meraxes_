#include <assert.h>
#include <math.h>

/* Local includes */
#include "blackhole_feedback.h"
#include "cooling.h"
#include "core/cooling.h"
#include "core/misc_tools.h"
#include "meraxes.h"

/**
 * @brief  The mass of the gas from the HotGas that cools down to the ColdGas of the central galaxy. Section 2.3 of DRAGONS3.
 *
 * The cooling happens only if the T_vir of the FoF group is above 10^4 K. All infalling mass is assumed to have been shocked to
 * the virial temperature of the FoF group, hence T_vir = 35.9 * V_vir(in kms-1)^2. A cooling radius r_cool is defined from the
 * T_vir value. Depending on whether r_cool is less or greater than the R_vir, different cooling rate prescriptions are employed.
 *
 * @param gal The galaxy onto which the gas is cooling.
 */
double gas_cooling(galaxy_t* gal)
{
  double cooling_mass = 0.0;

  /*! We only need to do cooling if there is anything to cool! */
  if (gal->HotGas > 1e-10) {
    fof_group_t* fof_group = gal->Halo->FOFGroup;

    /*! Calculate the halo virial temperature T_vir. N.B. This assumes ionised gas with mu=0.59... */
    double Tvir = 35.9 * fof_group->Vvir * fof_group->Vvir; // internal units (Kelvin)

    /*! There is no cooling if we are below 10^4 K */
    if (Tvir >= 1e4) {
      double t_cool, max_cooling_mass;
      double logZ, lambda, x, rho_r_cool, r_cool, isothermal_norm;
      run_units_t* units = &(run_globals.units);
      double max_cooling_mass_factor = run_globals.params.physics.MaxCoolingMassFactor;

      /*! Following Croton+ 2006, the maximum cooling time t_cool is set to be the dynamical time of the host dark matter halo */
      t_cool = fof_group->Rvir / fof_group->Vvir; // internal units

      /*! Compute the log10(metallicity) value */
      if (gal->MetalsHotGas > 0)
        logZ = log10(calc_metallicity(gal->HotGas, gal->MetalsHotGas));
      else
        logZ = -10.0;

      /*! Interpolate the temperature and metallicity dependant cooling rate funtion lambda(T,Z).
      Based on Sutherland & Dopita 1993 */
      lambda = interpolate_cooling_rate(log10(Tvir), logZ);

      /*! Following equation (3) of Croton+ 2006, calculate the hot gas density at
      the radius r_cool (i.e. where the cooling time is equal to `t_cool` above).
      (From the eq. 2 of DRAGONS3. Instead of solving for t_cool, we have found the rho_hot value) */ 
      x = PROTONMASS * BOLTZMANN * Tvir / lambda;              // now this has units sec g/cm^3
      x /= (units->UnitDensity_in_cgs * units->UnitTime_in_s); // now in internal units
      rho_r_cool = x / t_cool * 0.885;                         // 0.885 = 3/2 * mu, mu=0.59 for a fully ionized gas

      // TODO: We can actually get mu from the cooling tables of Sutherland &
      // Dopita for T>=1e4K.  We should do this rather than assuming the value.
      // THIS WILL BE ESPECIALLY IMPORTANT FOR TINY TIAMAT.

      /*! Assuming an isothermal density profile extending to Rvir, the cooling radius r_cool is calculated. */
      assert(rho_r_cool > 0);
      isothermal_norm = gal->HotGas / (4. * M_PI * fof_group->Rvir);
      r_cool = sqrt(isothermal_norm / rho_r_cool);
      gal->Rcool = r_cool;

      /*! The maximum amount of gas we can possibly cool is limited by the amount of mass within the free fall radius */
      
      /*! dt gives the time difference between when the halo (and hence the galaxy) was identified. Multiply by this factor
      as a kind of linear approximation to find the cooling mass.*/
      max_cooling_mass = max_cooling_mass_factor * gal->HotGas / t_cool * gal->dt;
      
      /*! Cases i & ii in Section 2.3 of DRAGONS3. */
      if (r_cool > fof_group->Rvir)
        /*! Here we are in the rapid cooling regime and we accrete all gas within the free-fall radius */
        cooling_mass = max_cooling_mass;
      // cooling_mass = gal->HotGas;
      else {
        /*! Here we are in the hot halo regime (but still limited by what's inside the free-fall radius). Eq 4 of DRAGONS3 */
        cooling_mass = max_cooling_mass / fof_group->Rvir * r_cool;
        if (cooling_mass > max_cooling_mass)
          cooling_mass = max_cooling_mass;
      }

      /*! Do one last sanity check to ensure we aren't cooling more gas than is available etc. */
      if (cooling_mass > gal->HotGas)
        cooling_mass = gal->HotGas;

      if (run_globals.params.physics.Flag_BHFeedback)
        cooling_mass -= radio_mode_BH_heating(gal, cooling_mass, x);

      if (cooling_mass < 0)
        cooling_mass = 0.0;
    }
  }
  return cooling_mass;
}

/**
 * @brief  Adds the cooling_mass calculated by gas_cooling() to the galaxy.
 * 
 * The cooling_mass is stripped from the HotGas and is added into the ColdGas. At the same time, the metallicity of the
 * HotGas and ColdGas are also updated accordingly.
 * 
 * @param gal The galaxy onto which the gas is cooling.
 * @param cooling_mass The mass of the cooled gas that gets added.
 */
void cool_gas_onto_galaxy(galaxy_t* gal, double cooling_mass)
{
  double cooling_metals;

  if (cooling_mass > gal->HotGas)
    cooling_mass = gal->HotGas;

  // what mass of metals is coming along with this cooling gas?
  cooling_metals = cooling_mass * calc_metallicity(gal->HotGas, gal->MetalsHotGas);

  // save the cooling mass
  gal->Mcool = cooling_mass;

  // update the galaxy reservoirs
  gal->HotGas -= cooling_mass;
  gal->MetalsHotGas -= cooling_metals;
  gal->ColdGas += cooling_mass;
  gal->MetalsColdGas += cooling_metals;
}
