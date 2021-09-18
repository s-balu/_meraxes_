#include <assert.h>
#include <math.h>

#include "blackhole_feedback.h"
#include "cooling.h"
#include "core/cooling.h"
#include "core/misc_tools.h"
#include "meraxes.h"

double gas_cooling(galaxy_t* gal)
{
  double cooling_mass = 0.0;

  // we only need to do cooling if there is anything to cool!
  if (gal->HotGas > 1e-10) {
    fof_group_t* fof_group = gal->Halo->FOFGroup;

    // calculate the halo virial temperature
    // N.B. This assumes ionised gas with mu=0.59...
    double Tvir = 35.9 * fof_group->Vvir * fof_group->Vvir; // internal units (Kelvin)

    // If we are below 10^4 K then no cooling either
    if (Tvir >= 1e4) {
      double t_cool, max_cooling_mass;
      double logZ, lambda, x, rho_r_cool, r_cool, isothermal_norm;
      run_units_t* units = &(run_globals.units);
      double max_cooling_mass_factor = run_globals.params.physics.MaxCoolingMassFactor;

      // following Croton+ 2006, we set the maximum cooling time to be the
      // dynamical time of the host dark matter halo
      t_cool = fof_group->Rvir / fof_group->Vvir; // internal units

      // get the log10(metallicity) value
      if (gal->MetalsHotGas > 0)
        logZ = log10(calc_metallicity(gal->HotGas, gal->MetalsHotGas));
      else
        logZ = -10.0;

      // interpolate the temperature and metallicity dependant cooling rate (lambda)
      lambda = interpolate_cooling_rate(log10(Tvir), logZ);

      // following equation (3) of Croton+ 2006, calculate the hot gas density at
      // the radius r_cool (i.e. where the cooling time is equal to `t_cool`
      // above)
      x = PROTONMASS * BOLTZMANN * Tvir / lambda;              // now this has units sec g/cm^3
      x /= (units->UnitDensity_in_cgs * units->UnitTime_in_s); // now in internal units
      rho_r_cool = x / t_cool * 0.885;                         // 0.885 = 3/2 * mu, mu=0.59 for a fully ionized gas

      // TODO: We can actually get mu from the cooling tables of Sutherland &
      // Dopita for T>=1e4K.  We should do this rather than assuming the value.
      // THIS WILL BE ESPECIALLY IMPORTANT FOR TINY TIAMAT.

      // under the assumption of an isothermal density profile extending to Rvir,
      // now calculate the cooling radius
      assert(rho_r_cool > 0);
      isothermal_norm = gal->HotGas / (4. * M_PI * fof_group->Rvir);
      r_cool = sqrt(isothermal_norm / rho_r_cool);
      //      gal->Rcool = r_cool;

      // the maximum amount of gas we can possibly cool is limited by the amount
      // of mass within the free fall radius
      max_cooling_mass = max_cooling_mass_factor * gal->HotGas / t_cool * gal->dt;

      if (r_cool > fof_group->Rvir)
        // here we are in the rapid cooling regime and we accrete all gas within
        // the free-fall radius
        cooling_mass = max_cooling_mass;
      // cooling_mass = gal->HotGas;
      else {
        // here we are in the hot halo regime (but still limited by what's inside the free-fall radius)
        cooling_mass = max_cooling_mass / fof_group->Rvir * r_cool;
        if (cooling_mass > max_cooling_mass)
          cooling_mass = max_cooling_mass;
      }

      // do one last sanity check to ensure we aren't cooling more gas than is available etc.
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

void cool_gas_onto_galaxy(galaxy_t* gal, double cooling_mass)
{
  double cooling_metals;

  if (cooling_mass > gal->HotGas)
    cooling_mass = gal->HotGas;

  // what mass of metals is coming along with this cooling gas?
  cooling_metals = cooling_mass * calc_metallicity(gal->HotGas, gal->MetalsHotGas);

  // save the cooling mass
  // gal->Mcool = cooling_mass;

  // update the galaxy reservoirs
  gal->HotGas -= cooling_mass;
  gal->MetalsHotGas -= cooling_metals;
  gal->ColdGas += cooling_mass;
  gal->MetalsColdGas += cooling_metals;
}
