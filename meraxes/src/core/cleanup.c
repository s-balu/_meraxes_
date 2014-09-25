#include "meraxes.h"
#include <hdf5.h>

void cleanup(run_globals_t *run_globals)
{
  SID_log("Running cleanup...", SID_LOG_OPEN);

  free_halo_storage(run_globals);

  cleanup_mags(run_globals);

  if (run_globals->RequestedForestId)
    SID_free(SID_FARG run_globals->RequestedForestId);

#ifdef USE_TOCF
  if (run_globals->params.TOCF_Flag)
    free_reionization_grids(run_globals);
#endif

  SID_log("Freeing hdf5 related stuff...", SID_LOG_OPEN);
  SID_free(SID_FARG run_globals->hdf5props.field_types);
  SID_free(SID_FARG run_globals->hdf5props.field_names);
  SID_free(SID_FARG run_globals->hdf5props.dst_field_sizes);
  SID_free(SID_FARG run_globals->hdf5props.dst_offsets);
  H5Tclose(run_globals->hdf5props.array3f_tid);
  SID_log(" ...done", SID_LOG_CLOSE);

  gsl_rng_free(run_globals->random_generator);

  SID_free(SID_FARG run_globals->LTTime);
  SID_free(SID_FARG run_globals->ZZ);
  SID_free(SID_FARG run_globals->AA);

  SID_log(" ...done", SID_LOG_CLOSE);

#ifdef DEBUG
  // close the debug file
  fclose(meraxes_debug_file);
#endif

  // close the log file
  // if(SID.n_proc > 1)
  // {
  //   fflush(SID.fp_log);
  //   fclose(SID.fp_log);
  // }
}
