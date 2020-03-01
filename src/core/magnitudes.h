#ifndef MAGNITUDES_H
#define MAGNITUDES_H

#include "meraxes.h"

#ifdef CALC_MAGS
#include <sector.h>
#endif

#define TOL 1e-30 // Minimum Flux

typedef struct mag_params_t {
    int targetSnap[MAGS_N_SNAPS];
    int nBeta;
    int nRest;
    int minZ;
    int maxZ;
    int nMaxZ;
    double tBC;
    int iAgeBC[MAGS_N_SNAPS];
    size_t totalSize;
    double *working;
    double *inBC;
    double *outBC;
    double *centreWaves;
    double *logWaves;
} mag_params_t;

enum core {MASTER};

#ifdef __cplusplus
extern "C" {
#endif

void init_luminosities(struct galaxy_t *gal);
void add_luminosities(mag_params_t *miniSpectra, struct galaxy_t *gal, int snapshot, double metals, double sfr);
void merge_luminosities(struct galaxy_t *target, struct galaxy_t *gal);
void init_templates_mini(mag_params_t *miniSpectra, char *fName, double *LTTime, int *targetSnaps, double *redshifts, double *betaBands, int nBeta, double *restBands, int nRest, double tBC);
void init_magnitudes(void);
void cleanup_mags(void);
void get_output_magnitudes(float *target, struct galaxy_t *gal, int snapshot);

#ifdef __cplusplus
}
#endif

#endif