/* Local includes */
#include "evolve.h"
#include "blackhole_feedback.h"
#include "cooling.h"
#include "core/stellar_feedback.h"
#include "infall.h"
#include "meraxes.h"
#include "mergers.h"
#include "reincorporation.h"
#include "star_formation.h"
#include "supernova_feedback.h"

/**
 * @brief  Evolves the galaxies forward in time
 *
 * @param fof_group The Friends-of-Friends struct that contains all the FoF groups in the simulation
 * @param snapshot The snapshot value at which the galaxies's evolution are to be computed
 * @param NGal Total number of galaxies in the simulation
 * @param Nfof Number of FoF groups
 */
//! Evolve existing galaxies forward in time
int evolve_galaxies(fof_group_t* fof_group, int snapshot, int NGal, int NFof)
{
  galaxy_t* gal = NULL;
  halo_t* halo = NULL;
  
  /*! Number of galaxies in the simulation */
  int gal_counter = 0;

  /*! Number of dead galaxies */
  int dead_gals = 0;

  /*! Mass of the infalling gas */
  double infalling_gas = 0;
  
  /*! The mass that falls from the FoF group down to the central halo/galaxy ?? */ 
  double cooling_mass = 0;
  
  /*! The number of steps that you need to hop from one snapshot to the next. Currently  NSteps = 1 ALWAYS */
  int NSteps = run_globals.params.NSteps;
  
  /*! Whether Instantaneous Recycling Approximation is enabled or not. 
  If so, the supernova feedback is instantaneous instead of spread over a number of snapshots */
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);

  mlog("Doing physics...", MLOG_OPEN | MLOG_TIMERSTART);
  // pre-calculate feedback tables for each lookback snapshot
  compute_stellar_feedback_tables(snapshot);

  /*! Loop over all the FoF groups in the present snapshot */
  for (int i_fof = 0; i_fof < NFof; i_fof++) {
    /*! Skip to the next FoF group if this one is empty i.e. if the halo is not occupied*/
    // First check to see if this FOF group is empty.  If it is then skip it.
    if (fof_group[i_fof].FirstOccupiedHalo == NULL)
      continue;

    /*! Compute the amount of infalling gas that is to be added to the FoF group */
    infalling_gas = gas_infall(&(fof_group[i_fof]), snapshot);
    
    for (int i_step = 0; i_step < NSteps; i_step++) {
      halo = fof_group[i_fof].FirstHalo;
      /* Start with the central halo and go "down" in the halo structure of an FoF group */
      while (halo != NULL) {
        gal = halo->Galaxy;

        while (gal != NULL) {
          if (gal->Type == 0) {
            cooling_mass = gas_cooling(gal);

            add_infall_to_hot(gal, infalling_gas / ((double)NSteps));

            reincorporate_ejected_gas(gal);

            cool_gas_onto_galaxy(gal, cooling_mass);
          }

          if (gal->Type < 3) {
            if (!Flag_IRA)
              delayed_supernova_feedback(gal, snapshot);

            if (gal->BlackHoleAccretingColdMass > 0)
              previous_merger_driven_BH_growth(gal);

            insitu_star_formation(gal, snapshot);

            // If this is a type 2 then decrement the merger clock
            if (gal->Type == 2)
              gal->MergTime -= gal->dt;
          }

          if (i_step == NSteps - 1)
            gal_counter++;

          gal = gal->NextGalInHalo;
        }

        halo = halo->NextHaloInFOFGroup;
      }

      // Check for mergers
      halo = fof_group[i_fof].FirstHalo;
      while (halo != NULL) {
        gal = halo->Galaxy;
        while (gal != NULL) {
          if (gal->Type == 2)
            // If the merger clock has run out or our target halo has already
            // merged then process a merger event.
            if ((gal->MergTime < 0) || (gal->MergerTarget->Type == 3))
              merge_with_target(gal, &dead_gals, snapshot);

          gal = gal->NextGalInHalo;
        }
        halo = halo->NextHaloInFOFGroup;
      }
    }
  }

  if (gal_counter + (run_globals.NGhosts) != NGal) {
    mlog_error("We have not processed the expected number of galaxies...");
    mlog("gal_counter = %d but NGal = %d", MLOG_MESG, gal_counter, NGal);
    ABORT(EXIT_FAILURE);
  }

  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

  return gal_counter - dead_gals;
}

void passively_evolve_ghost(galaxy_t* gal, int snapshot)
{
  // Passively evolve ghosts.
  // Currently, this just means evolving their stellar pops...

  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);

  if (!Flag_IRA)
    delayed_supernova_feedback(gal, snapshot);
}
