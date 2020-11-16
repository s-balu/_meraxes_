/* Local includes */
#include "reincorporation.h"
#include "core/misc_tools.h"
#include "meraxes.h"

/**
 * @brief  Add the reincorporated mass computed by the reincorporate_ejected_gas() to the HotGas component of the central galaxy.
 * 
 * The reincorporated gas is stripped from the EjectedGas and added to the HotGas component. The metallicities of the HotGas
 * and the EjectedGas are modified accordingly.
 * 
 * @param gal The central galaxy of the FoF group onto which the reincorporated mass is added.
 * @param reincorporated Mass of the reincorporated gas.
 */
static void update_reservoirs_from_reincorporation(galaxy_t* gal, double reincorporated)
{
  double metals = reincorporated * calc_metallicity(gal->EjectedGas, gal->MetalsEjectedGas);

  gal->EjectedGas -= reincorporated;
  gal->MetalsEjectedGas -= metals;
  gal->HotGas += reincorporated;
  gal->MetalsHotGas += metals;
}

/**
 * @brief  Compute the reincorporated mass of the central galaxy.
 *
 *A portion of the EjectedGas of the central gas can cool down and fall back onto the HotGas halo. Two prescriptions
 *are employed to find this reincorporated mass. This process is assumed to happen over a dynamical time scale of the halo.
 * 
 * @param gal The central galaxy of the FoF group.
 */
void reincorporate_ejected_gas(galaxy_t* gal)
{ 
  double ReincorporationEff = run_globals.params.physics.ReincorporationEff;

  if (gal->EjectedGas > 0 && ReincorporationEff > 0.) {
    int ReincorporationModel = run_globals.params.physics.ReincorporationModel;
    fof_group_t* fof_group = gal->Halo->FOFGroup;
    double reincorporated = 0.;
    double t_dyn = fof_group->Rvir / fof_group->Vvir;
    double t_rein;

    switch (ReincorporationModel) {
      case 1:
        /*! Allow some of the ejected gas associated with the central to be reincorporated
        following the prescription of Guo 2010 (which is actually almost identical to SAGE). */
        reincorporated = ReincorporationEff * gal->EjectedGas * (gal->dt / t_dyn);
        break;
      case 2:
        /*! Following the prescription of Henriques et al. 2013 */
        t_rein = ReincorporationEff / (fof_group->Mvir / run_globals.params.Hubble_h); // Unit: Myr
        // Convert to Meraxes intrinsic unit
        t_rein /= run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h;
        if (t_rein < t_dyn)
          t_rein = t_dyn;
        reincorporated = gal->EjectedGas * (gal->dt / t_rein);
        break;
      default:
        mlog_error("Unknown ReincorporationModel!");
        ABORT(EXIT_FAILURE);
        break;
    }

    // ensure consistency
    if (reincorporated > gal->EjectedGas)
      reincorporated = gal->EjectedGas;

    /*! Update the baryonic reservoirs */
    update_reservoirs_from_reincorporation(gal, reincorporated);
  }
}
