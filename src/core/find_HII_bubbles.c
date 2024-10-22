#include "meraxes.h"
#include <fftw3-mpi.h>
#include <math.h>

#include "recombinations.c"

/*
 * This code is a re-write of the modified version of 21cmFAST used in Mutch et
 * al. (2016; Meraxes paper).  The original code was written by Andrei Mesinger
 * with additions as detailed in Sobacchi & Mesinger (2013abc).  Updates were
 * subsequently made by Simon Mutch & Paul Geil.
 *
 * Inclusion of electron fraction (X-ray heating) and inhomogeneous recombinations
 * by Bradley Greig. Relevant functions taken from public version of 21cmFAST.
 */

double RtoM(double R)
{
    // All in internal units
    int filter = run_globals.params.ReionRtoMFilterType;
    double OmegaM = run_globals.params.OmegaM;
    double RhoCrit = run_globals.RhoCrit;

    switch (filter) {
        case 0: //top hat M = (4/3) PI <rho> R^3
            return (4.0 / 3.0) * M_PI * pow(R, 3) * (OmegaM * RhoCrit);
        case 1: //gaussian: M = (2PI)^1.5 <rho> R^3
            return pow(2 * M_PI, 1.5) * OmegaM * RhoCrit * pow(R, 3);
        default: // filter not defined
            mlog_error("Unrecognised filter (%d). Aborting...", filter);
            ABORT(EXIT_FAILURE);
            break;
    }

    return -1;
}

void _find_HII_bubbles(int snapshot)
{
    // TODO: TAKE A VERY VERY CLOSE LOOK AT UNITS!!!!

    double box_size = run_globals.params.BoxSize; // Mpc/h
    int ReionGridDim = run_globals.params.ReionGridDim;
    double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
    double cell_length_factor = L_FACTOR;
    double total_n_cells = pow((double)ReionGridDim, 3);
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
    int slab_n_real = local_nix * ReionGridDim * ReionGridDim;
    int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
    int flag_ReionUVBFlag = run_globals.params.ReionUVBFlag;
    double ReionEfficiency = run_globals.params.physics.ReionEfficiency;
    double ReionNionPhotPerBary = run_globals.params.physics.ReionNionPhotPerBary;
    run_units_t* units = &(run_globals.units);
    float J_21_aux = 0;
    double J_21_aux_constant;
    double density_over_mean;
    double sfr_density;
    double f_coll_stars;
    double electron_fraction;
    double Gamma_R_prefactor;

    double dNrec, rec, Gamma_R;
    float fabs_dtdz, ZSTEP, z_eff;

    double redshift = run_globals.ZZ[snapshot];
    double prev_redshift;
    if(snapshot==0) {
        prev_redshift = run_globals.ZZ[snapshot];
    }
    else {
        prev_redshift = run_globals.ZZ[snapshot-1];
    }

    ZSTEP = (float)(prev_redshift - redshift);
    fabs_dtdz = (float)fabs(dtdz((float)redshift) / run_globals.params.Hubble_h);

    // Initialise interpolation tables for inhomogeneous recombinations
    // TODO(merge): Does this need to be done every call?
    if(run_globals.params.Flag_IncludeRecombinations) {
        init_MHR();
    }

    int i_real;
    int i_padded;

    // This parameter choice is sensitive to noise on the cell size, at least for the typical
    // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
    if ((box_size / (double)ReionGridDim) < 1.0) // Fairly arbitrary length based on 2 runs Sobacchi did
        cell_length_factor = 1.0;

    // Init J_21
    float* J_21 = run_globals.reion_grids.J_21;
    if (flag_ReionUVBFlag)
        for (int ii = 0; ii < slab_n_real; ii++)
            J_21[ii] = 0.0;

    // Init xH
    float* xH = run_globals.reion_grids.xH;
    for (int ii = 0; ii < slab_n_real; ii++)
        xH[ii] = 1.0;

    // Init r_bubble
    float* r_bubble = run_globals.reion_grids.r_bubble;
    for (int ii = 0; ii < slab_n_real; ii++)
        r_bubble[ii] = 0.0;

    // #ifdef DEBUG
    //   {
    //     char fname_debug[STRLEN];
    //     sprintf(fname_debug, "pre_stars-%d.dat", run_globals.mpi_rank);
    //     FILE *fd_debug = fopen(fname_debug, "wb");
    //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
    //       for(int iy = 0; iy < ReionGridDim; iy++)
    //         for(int iz = 0; iz < ReionGridDim; iz++)
    //           fwrite(&(run_globals.reion_grids.stars[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
    //     fclose(fd_debug);
    //   }
    // #endif

    // #ifdef DEBUG
    //   {
    //     char fname_debug[STRLEN];
    //     sprintf(fname_debug, "pre_sfr-%d.dat", run_globals.mpi_rank);
    //     FILE *fd_debug = fopen(fname_debug, "wb");
    //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
    //       for(int iy = 0; iy < ReionGridDim; iy++)
    //         for(int iz = 0; iz < ReionGridDim; iz++)
    //           fwrite(&(run_globals.reion_grids.sfr[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
    //     fclose(fd_debug);
    //   }
    // #endif

    // Forward fourier transform to obtain k-space fields
    // TODO: Ensure that fftwf_mpi_init has been called and fftwf_mpi_cleanup will be called
    // TODO: Don't use estimate and calculate plan in code init
    float* deltax = run_globals.reion_grids.deltax;
    float* deltax_temp = run_globals.reion_grids.deltax_temp;

    // Make a copy of the box for FFT'ing
    memcpy(deltax_temp, deltax, sizeof(fftwf_complex) * slab_n_complex);

    fftwf_complex* deltax_unfiltered = (fftwf_complex*)deltax_temp; // WATCH OUT!
    fftwf_complex* deltax_filtered = run_globals.reion_grids.deltax_filtered;
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, deltax_temp, deltax_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    float* stars = run_globals.reion_grids.stars;
    float* stars_temp = run_globals.reion_grids.stars_temp;

    // Make a copy of the box for FFT'ing
    memcpy(stars_temp, stars, sizeof(fftwf_complex) * slab_n_complex);

    fftwf_complex* stars_unfiltered = (fftwf_complex*)stars_temp; // WATCH OUT!
    fftwf_complex* stars_filtered = run_globals.reion_grids.stars_filtered;
    plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, stars_temp, stars_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    float* sfr = run_globals.reion_grids.sfr;
    float* sfr_temp = run_globals.reion_grids.sfr_temp;

    // Make a copy of the box for FFT'ing
    memcpy(sfr_temp, sfr, sizeof(fftwf_complex) * slab_n_complex);

    fftwf_complex* sfr_unfiltered = (fftwf_complex*)sfr_temp; // WATCH OUT!
    fftwf_complex* sfr_filtered = run_globals.reion_grids.sfr_filtered;
    plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr_temp, sfr_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // The free electron fraction from X-rays
    float* x_e_box;
    fftwf_complex* x_e_unfiltered;
    fftwf_complex* x_e_filtered;
    if(run_globals.params.Flag_IncludeSpinTemp) {
        x_e_box = run_globals.reion_grids.x_e_box;
        x_e_unfiltered = (fftwf_complex*)x_e_box; // WATCH OUT!
        x_e_filtered = run_globals.reion_grids.x_e_filtered;
        plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, x_e_box, x_e_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }

    // Fields relevant for computing the inhomogeneous recombinations
    float *z_re, *Gamma12, *N_rec_prev, *N_rec;
    fftwf_complex* N_rec_unfiltered;
    fftwf_complex* N_rec_filtered;
    if(run_globals.params.Flag_IncludeRecombinations) {
        z_re = run_globals.reion_grids.z_re;
        Gamma12 = run_globals.reion_grids.Gamma12;

        N_rec = run_globals.reion_grids.N_rec;
        N_rec_prev = run_globals.reion_grids.N_rec_prev;

        // Make a copy of the box for FFT'ing
        memcpy(N_rec_prev, N_rec, sizeof(fftwf_complex) * slab_n_complex);

        N_rec_unfiltered = (fftwf_complex*)N_rec_prev; // WATCH OUT!
        N_rec_filtered = run_globals.reion_grids.N_rec_filtered;
        plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, N_rec_prev, N_rec_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    }

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    // TODO: Double check that looping over correct number of elements here
    for (int ii = 0; ii < slab_n_complex; ii++) {
        deltax_unfiltered[ii] /= total_n_cells;
        stars_unfiltered[ii] /= total_n_cells;
        sfr_unfiltered[ii] /= total_n_cells;
        if(run_globals.params.Flag_IncludeRecombinations) {
            N_rec_unfiltered[ii] /= total_n_cells;
        }
        if(run_globals.params.Flag_IncludeSpinTemp) {
            x_e_unfiltered[ii] /= total_n_cells;
        }
    }

    // Loop through filter radii
    double ReionRBubbleMax;
    if(run_globals.params.Flag_IncludeRecombinations) {
        ReionRBubbleMax = run_globals.params.physics.ReionRBubbleMaxRecomb; // Mpc/h
    }
    else {
        ReionRBubbleMax = run_globals.params.physics.ReionRBubbleMax; // Mpc/h
    }
    double ReionRBubbleMin = run_globals.params.physics.ReionRBubbleMin; // Mpc/h
    double R = fmin(ReionRBubbleMax, L_FACTOR * box_size); // Mpc/h
    double ReionDeltaRFactor = run_globals.params.ReionDeltaRFactor;
    double ReionGammaHaloBias = run_globals.params.physics.ReionGammaHaloBias;

    bool flag_last_filter_step = false;

    // set recombinations to zero (for case when recombinations are not used)
    rec = 0.0;

    while (!flag_last_filter_step) {
        // check to see if this is our last filtering step
        if (((R / ReionDeltaRFactor) <= (cell_length_factor * box_size / (double)ReionGridDim))
                || ((R / ReionDeltaRFactor) <= ReionRBubbleMin)) {
            flag_last_filter_step = true;
            R = cell_length_factor * box_size / (double)ReionGridDim;
        }

        // DEBUG
        // mlog("R = %.2e (h=0.678 -> %.2e)", MLOG_MESG, R, R/0.678);
        mlog(".", MLOG_CONT);

        // copy the k-space grids
        memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
        memcpy(stars_filtered, stars_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
        memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex);

        if(run_globals.params.Flag_IncludeRecombinations) {
            memcpy(N_rec_filtered, N_rec_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
        }
        if(run_globals.params.Flag_IncludeSpinTemp) {
            memcpy(x_e_filtered, x_e_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
        }

        // do the filtering unless this is the last filter step
        int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);
        if (!flag_last_filter_step) {
            filter(deltax_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
            filter(stars_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
            filter(sfr_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);

            if(run_globals.params.Flag_IncludeRecombinations) {
                filter(N_rec_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
            }
            if(run_globals.params.Flag_IncludeSpinTemp) {
                filter(x_e_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
            }
        }

        // inverse fourier transform back to real space
        plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, deltax_filtered, (float*)deltax_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, stars_filtered, (float*)stars_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr_filtered, (float*)sfr_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        if(run_globals.params.Flag_IncludeRecombinations) {
            plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, N_rec_filtered, (float*)N_rec_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        }

        if(run_globals.params.Flag_IncludeSpinTemp) {
            plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, x_e_filtered, (float*)x_e_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        }

        // Perform sanity checks to account for aliasing effects
        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                    ((float*)deltax_filtered)[i_padded] = fmaxf(((float*)deltax_filtered)[i_padded], -1 + REL_TOL);
                    ((float*)stars_filtered)[i_padded] = fmaxf(((float*)stars_filtered)[i_padded], 0.0);
                    ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

                    if(run_globals.params.Flag_IncludeRecombinations) {
                        ((float*)N_rec_filtered)[i_padded] = fmaxf(((float*)N_rec_filtered)[i_padded], 0.0);
                    }
                    if(run_globals.params.Flag_IncludeSpinTemp) {
                        ((float*)x_e_filtered)[i_padded] = fmaxf(((float*)x_e_filtered)[i_padded], 0.0);
                        ((float*)x_e_filtered)[i_padded] = fminf(((float*)x_e_filtered)[i_padded], 0.999);
                    }
                }

        // #ifdef DEBUG
        //   {
        //     char fname_debug[STRLEN];
        //     sprintf(fname_debug, "post_stars-%d_R%.3f.dat", run_globals.mpi_rank, R);
        //     FILE *fd_debug = fopen(fname_debug, "wb");
        //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
        //       for(int iy = 0; iy < ReionGridDim; iy++)
        //         for(int iz = 0; iz < ReionGridDim; iz++)
        //           fwrite(&(run_globals.reion_grids.stars_filtered[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
        //     fclose(fd_debug);
        //   }
        //   {
        //     char fname_debug[STRLEN];
        //     sprintf(fname_debug, "post_sfr-%d_R%.3f.dat", run_globals.mpi_rank, R);
        //     FILE *fd_debug = fopen(fname_debug, "wb");
        //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
        //       for(int iy = 0; iy < ReionGridDim; iy++)
        //         for(int iz = 0; iz < ReionGridDim; iz++)
        //           fwrite(&(run_globals.reion_grids.sfr_filtered[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
        //     fclose(fd_debug);
        //   }
        // #endif

        /*
         * Main loop through the box...
         */

        J_21_aux_constant = (1.0 + redshift) * (1.0 + redshift) / (4.0 * M_PI)
            * run_globals.params.physics.ReionAlphaUV * PLANCK
            * 1e21 // * run_globals.params.physics.ReionEscapeFrac
            * R * units->UnitLength_in_cm * ReionNionPhotPerBary / PROTONMASS
            * units->UnitMass_in_g / pow(units->UnitLength_in_cm, 3) / units->UnitTime_in_s;

        if(run_globals.params.Flag_IncludeRecombinations) {
            Gamma_R_prefactor = (1.0 + redshift) * (1.0 + redshift) * R * (units->UnitLength_in_cm / run_globals.params.Hubble_h ) * SIGMA_HI *
                run_globals.params.physics.ReionAlphaUV / (run_globals.params.physics.ReionAlphaUV+2.75) / 1.0e-12;  // Converting R h^-1 to R.
        }

        // DEBUG
        // for (int ix=0; ix<local_nix; ix++)
        //   for (int iy=0; iy<ReionGridDim; iy++)
        //     for (int iz=0; iz<ReionGridDim; iz++)
        //       if (((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] > 0)
        //       {
        //         mlog("J_21_aux ==========", MLOG_OPEN);
        //         mlog("redshift = %.2f", MLOG_MESG, redshift);
        //         mlog("alpha = %.2e", MLOG_MESG, run_globals.params.physics.ReionAlphaUV);
        //         mlog("PLANCK = %.2e", MLOG_MESG, PLANCK);
        //         mlog("ReionEscapeFrac = %.2f", MLOG_MESG, ReionEscapeFrac);
        //         mlog("ReionNionPhotPerBary = %.1f", MLOG_MESG, ReionNionPhotPerBary);
        //         mlog("UnitMass_in_g = %.2e", MLOG_MESG, units->UnitMass_in_g);
        //         mlog("UnitLength_in_cm = %.2e", MLOG_MESG, units->UnitLength_in_cm);
        //         mlog("PROTONMASS = %.2e", MLOG_MESG, PROTONMASS);
        //         mlog("-> J_21_aux_constant = %.2e", MLOG_MESG, J_21_aux_constant);
        //         mlog("pixel_volume = %.2e", MLOG_MESG, pixel_volume);
        //         mlog("sfr_density[%d, %d, %d] = %.2e", MLOG_MESG, ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] / pixel_volume);
        //         mlog("-> J_21_aux = %.2e", MLOG_MESG, ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] / pixel_volume * J_21_aux_constant);
        //         mlog("==========", MLOG_CLOSE);
        //         if (run_globals.mpi_rank == 0)
        //           ABORT(EXIT_SUCCESS);
        //       }

        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

                    density_over_mean = 1.0 + (double)((float*)deltax_filtered)[i_padded];

                    f_coll_stars = (double)((float*)stars_filtered)[i_padded] / (RtoM(R) * density_over_mean)
                        * (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;

                    sfr_density = (double)((float*)sfr_filtered)[i_padded] / pixel_volume; // In internal units

                    // #ifdef DEBUG
                    //           if(abs(redshift - 23.074) < 0.01)
                    //             debug("%d, %d, %d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
                    //                 run_globals.mpi_rank, ix, iy, iz,
                    //                 R, density_over_mean, f_coll_stars,
                    //                 ((float *)stars_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)],
                    //                 RtoM(R), (4.0/3.0)*M_PI*pow(R,3.0), pixel_volume, 1.0/ReionEfficiency, sfr_density, ((float*)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)], (float)(sfr_density * J_21_aux_constant));
                    // #endif

                    // Calculate the recombinations within the cell
                    if(run_globals.params.Flag_IncludeRecombinations) {
                        Gamma_R = Gamma_R_prefactor * sfr_density * (units->UnitMass_in_g / units->UnitTime_in_s) * pow( units->UnitLength_in_cm / run_globals.params.Hubble_h, -3. )
                            *  ReionNionPhotPerBary / PROTONMASS; // Convert pixel volume (Mpc/h)^3 -> (cm)^3
                        rec = (double)((float*)N_rec_filtered)[i_padded] / density_over_mean;
                    }

                    // Account for the partial ionisation of the cell from X-rays
                    if(run_globals.params.Flag_IncludeSpinTemp) {
                        electron_fraction = 1.0 - ((float*)x_e_filtered)[i_padded];
                    }
                    else {
                        electron_fraction = 1.0;
                    }

                    if (flag_ReionUVBFlag)
                        J_21_aux = (float)(sfr_density * J_21_aux_constant);

                    // Modified reionisation condition, including recombinations and partial ionisations from X-rays
                    // Check if ionised!
                    if (f_coll_stars > ( electron_fraction / ReionEfficiency ) * (1. + rec) ) // IONISED!!!!
                    {
                        // If it is the first crossing of the ionisation barrier for this cell (largest R), let's record J_21
                        if (xH[i_real] > REL_TOL) {
                            if (flag_ReionUVBFlag)
                                J_21[i_real] = J_21_aux;

                            // Store the ionisation background and the reionisation redshift for each cell
                            if(run_globals.params.Flag_IncludeRecombinations) {
                                Gamma12[i_real] = (float)Gamma_R;
                                if(z_re[i_real] < 0) {
                                    z_re[i_real] = (float)redshift;
                                }
                            }
                        }

                        // Mark as ionised
                        xH[i_real] = 0;

                        // Record radius
                        r_bubble[i_real] = (float)R;
                    }
                    // Check if this is the last filtering step.
                    // If so, assign partial ionisations to those cells which aren't fully ionised
                    else if (flag_last_filter_step && (xH[i_real] > REL_TOL))
                    {
                        xH[i_real] = (float)(electron_fraction - f_coll_stars * ReionEfficiency);
                        if(xH[i_real] < 0.) {
                            xH[i_real] = (float)0.;
                        }
                        else if (xH[i_real] > 1.0) {
                            xH[i_real] = (float)1.;
                        }
                    }

                    // Check if new ionisation
                    float* z_in = run_globals.reion_grids.z_at_ionization;
                    if ((xH[i_real] < REL_TOL) && (z_in[i_real] < 0)) // New ionisation!
                    {
                        z_in[i_real] = (float)redshift;
                        if (flag_ReionUVBFlag)
                            run_globals.reion_grids.J_21_at_ionization[i_real] = J_21_aux * (float)ReionGammaHaloBias;
                    }
                }
        // iz
        R /= ReionDeltaRFactor;
    }

    // Find the volume and mass weighted neutral fractions
    // TODO: The deltax grid will have rounding errors from forward and reverse
    //       FFT. Should cache deltax slabs prior to ffts and reuse here.
    double volume_weighted_global_xH = 0.0;
    double mass_weighted_global_xH = 0.0;
    double mass_weight = 0.0;

    for (int ix = 0; ix < local_nix; ix++)
        for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
                i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                volume_weighted_global_xH += (double)xH[i_real];
                density_over_mean = 1.0 + (double)((float*)deltax_filtered)[i_padded];
                mass_weighted_global_xH += (double)(xH[i_real]) * density_over_mean;
                mass_weight += density_over_mean;
            }

    MPI_Allreduce(MPI_IN_PLACE, &volume_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &mass_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &mass_weight, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

    volume_weighted_global_xH /= total_n_cells;
    mass_weighted_global_xH /= mass_weight;
    run_globals.reion_grids.volume_weighted_global_xH = volume_weighted_global_xH;
    run_globals.reion_grids.mass_weighted_global_xH = mass_weighted_global_xH;

    if(run_globals.params.Flag_IncludeRecombinations) {
        // Store the resultant recombination grid
        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                    i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

                    density_over_mean = 1.0 + deltax[i_padded];
                    z_eff = (float)((1. + redshift) * pow(density_over_mean, 1.0/3.0) - 1);
                    dNrec = splined_recombination_rate(z_eff, Gamma12[i_real]) * fabs_dtdz * ZSTEP * (1. - xH[i_real]);
                    N_rec[i_padded] += dNrec;

                }

        // TODO(merge): Does this need to be done every call?
        free_MHR();
    }

}

// This function makes sure that the right version of find_HII_bubbles() gets called.
void find_HII_bubbles(int snapshot, timer_info* timer_total)
{
    // Call the version of find_HII_bubbles we've been passed (and time it)
    double redshift = run_globals.ZZ[snapshot];
    timer_info timer;
#ifdef USE_CUDA
#ifdef USE_CUFFT
    mlog("Calling pure-GPU version of find_HII_bubbles() for snap=%d/z=%.2lf...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, redshift);
#else
    mlog("Calling hybrid-GPU/FFTW version of find_HII_bubbles() for snap=%d/z=%.2lf...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, redshift);
#endif
    // Run the GPU version of _find_HII_bubbles()
    timer_start(&timer);

    int flag_write_validation_data = false;
    _find_HII_bubbles_gpu(redshift, flag_write_validation_data);
#else
    // Run the Meraxes version of _find_HII_bubbles()
    mlog("Calling pure-CPU version of find_HII_bubbles() for snap=%d/z=%.2lf...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, redshift);
    timer_start(&timer);
    _find_HII_bubbles(snapshot);
#endif
    timer_stop(&timer);
    timer_stop(timer_total);
    timer_gpu += timer_delta(timer);
    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
    mlog("Total time spent in find_HII_bubbles vs. total run time (snapshot %d ): %.2f of %.2f s", MLOG_MESG, snapshot, timer_gpu, timer_delta(*timer_total));
}
