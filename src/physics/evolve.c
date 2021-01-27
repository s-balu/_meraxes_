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
 * @brief  Evolves the existing galaxies forward in time.
 *
 * @param fof_group The fof struct that contains all the FoF groups in the simulation.
 * @param snapshot The snapshot value at which the galaxies's evolution is to be computed.
 * @param NGal Total number of galaxies currently exist in the simulation.
 * @param Nfof Number of FoF groups.
 */
int evolve_galaxies(fof_group_t* fof_group, int snapshot, int NGal, int NFof)
{
  /*! galaxy struct to store the galaxy properties. */
  galaxy_t* gal = NULL;

  /*! halo_t struct to store the halo properties. */
  halo_t* halo = NULL;
  
  /*! Number of galaxies in the simulation. */
  int gal_counter = 0;

  /*! Number of dead galaxies. */
  int dead_gals = 0;

  /*! Mass of the infalling gas that is to be added to the FoF group. This mass is assumed to be shocked to the virial
  temperature of the host FoF group. */
  double infalling_gas = 0;
  
  /*! The mass that cools down from the HotGas to the ColdGas. */ 
  double cooling_mass = 0;
  
  /*! The number of steps that you need to hop from one snapshot to the next. Currently NSteps = 1 ALWAYS. */
  int NSteps = run_globals.params.NSteps;
  
  /*! Whether Instantaneous Recycling Approximation is enabled or not. 
  If so, the supernova feedback is instantaneous instead of spread over a number of snapshots */
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);

  mlog("Doing physics...", MLOG_OPEN | MLOG_TIMERSTART);
  // pre-calculate feedback tables for each lookback snapshot
  compute_stellar_feedback_tables(snapshot);

  /*! Loop over all the FoF groups in the snapshot. Loops over all the structures in the present snapshot.
  The hierarchy is FoF Group -> Haloes -> Galaxies */
  for (int i_fof = 0; i_fof < NFof; i_fof++) {

    /*! Skip to the next FoF group if there is no halo present */
    // First check to see if this FOF group is empty.  If it is then skip it.
    if (fof_group[i_fof].FirstOccupiedHalo == NULL)
      continue;

    infalling_gas = gas_infall(&(fof_group[i_fof]), snapshot);

    /*! This loop essentially runs only once since NSteps = 1. !!! */
    for (int i_step = 0; i_step < NSteps; i_step++) {
      halo = fof_group[i_fof].FirstHalo;
      /* Start with the central halo and go "down" in the halo structure of an FoF group */
      while (halo != NULL) {
        gal = halo->Galaxy;
        
        /* Cycle through all the galaxies in the present halo. Depending on the Type of the galaxy, different
        physics are implemented. There are currently four types of galaxies that are tracked:
        Type 0 : Central galaxy
        Type 1 : Satellite galaxy
        Type 2 : Halo-less galaxy: The halo of this galaxy has merged with a bigger halo. Also known
                            as Ghost galaxies. Becomes a Type3 galaxy once the mereger clock runs out.
        Type 3 : Dead galaxy: The galaxy has merged to the central galaxy. Type2 evolves into these.
        */
        while (gal != NULL) {
          
          /*! Only do this if the galaxy is the central galaxy (Type 0). */
          if (gal->Type == 0) {

            /*! The cooling mass that cools down from the HotGas to the ColdGas. */
            cooling_mass = gas_cooling(gal);

            /* The infalling gas is added to the HotGas. */
            add_infall_to_hot(gal, infalling_gas / ((double)NSteps));

            /*! A portion of the EjectedGas is added to the HotGas. */
            reincorporate_ejected_gas(gal);

            /*! The cooling_mass is added to the ColdGas from the HotGas. */
            cool_gas_onto_galaxy(gal, cooling_mass);
          }

          /*! Do this for all galaxies except the Type 3 ones which are dead galaxies. */
          if (gal->Type < 3) {

            /*! Apply the SN feedback as spread across the snapshot if the Flag_IRA is not set ON. */
            if (!Flag_IRA)
              delayed_supernova_feedback(gal, snapshot);

            if (gal->BlackHoleAccretingColdMass > 0)
              previous_merger_driven_BH_growth(gal);

            insitu_star_formation(gal, snapshot);

            /*! If this is a Type 2 galaxy then decrement the merger clock */
            if (gal->Type == 2)
              gal->MergTime -= gal->dt;
          }

          if (i_step == NSteps - 1)
            gal_counter++;

          gal = gal->NextGalInHalo;
        }

        halo = halo->NextHaloInFOFGroup;
      }

      
      /*! Check for mergers */
      halo = fof_group[i_fof].FirstHalo;
      
      while (halo != NULL) {
        gal = halo->Galaxy;
        while (gal != NULL) {
          /*! Type 2 galaxies are the ones that are about to merge */
          if (gal->Type == 2)
            
            /*! If the merger clock has run out OR our target halo has already merged then process a merger event. */
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

/**
 * @brief Passively evolve the ghost galaxies.
 * 
 * @param gal The Type 3 (dead) galaxy.
 * @param snapshot The snapshot value at which the galaxies's evolution are to be computed.
 */
void passively_evolve_ghost(galaxy_t* gal, int snapshot)
{
  // Passively evolve ghosts.
  // Currently, this just means evolving their stellar pops...

  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);

  if (!Flag_IRA)
    delayed_supernova_feedback(gal, snapshot);
}