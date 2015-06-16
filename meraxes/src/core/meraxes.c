#define _MAIN
#include "meraxes.h"
#include <sys/stat.h>

#ifdef USE_TOCF
#include <21cmfast.h>
#endif


int main(int argc, char **argv)
{
  // init SID
  SID_init(&argc, &argv, NULL, NULL);

  struct stat filestatus;
  run_globals_t run_globals;

  // char log_fname[50];
  // FILE *log_file = NULL;
  // if(SID.n_proc > 1)
  // {
  //   SID.flag_log_allranks = 1;
  //   sprintf(log_fname, "rank_%d.log", SID.My_rank);
  //   log_file = fopen(log_fname, "w");
  //   SID.fp_log = log_file;
  // }

  // deal with any input arguments
  if (argc != 2)
  {
    SID_log("\n  usage: %s <parameterfile>\n\n", SID_LOG_COMMENT, argv[0]);
    ABORT(1);
  }

#ifdef DEBUG
  // open the debug file for this core
  char debug_fname[50];
  sprintf(debug_fname, "debug_%d.txt", SID.My_rank);
  meraxes_debug_file = fopen(debug_fname, "w");
  // if(SID.My_rank==0)
  //   mpi_debug_here();
#endif

#ifdef USE_TOCF
  // Note that this must be done *before* we read the parameter file as we may
  // want to overwrite some of the set defaults.
  init_default_tocf_params();
#endif

  // read the input parameter file
  read_parameter_file(&run_globals, argv[1], 0);

  // Check to see if the output directory exists and if not, create it
  if (stat(run_globals.params.OutputDir, &filestatus) != 0)
    mkdir(run_globals.params.OutputDir, 02755);

#ifdef USE_TOCF
    // Check to see if the tocf_logs directory exists and if not, create it
    if (stat(tocf_params.logfile_dir, &filestatus) != 0)
        mkdir(tocf_params.logfile_dir, 02755);
#endif

  // initiate meraxes
  init_meraxes(&run_globals);

  // calculate the output hdf5 file properties for later use
  calc_hdf5_props(&run_globals);

  // Run the model!
  if (!run_globals.params.FlagInteractive)
    dracarys(&run_globals);
  else
  {
    while (run_globals.params.FlagInteractive)
    {
      dracarys(&run_globals);
      continue_prompt(&run_globals, argv[1]);
    }
  }

  // cleanup
  cleanup(&run_globals);

  SID_exit(EXIT_SUCCESS);
}
