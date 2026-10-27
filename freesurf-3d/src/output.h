#ifndef FS3_OUTPUT_H
#define FS3_OUTPUT_H

#include <stddef.h>

#include "diagnostics.h"
#include "sim.h"

/* IEEE binary16 storage type (a GNU C extension before C23) */
__extension__ typedef _Float16 fs_half;

#define OUT_NSTORED 5  /* stored per frame: C, p_dyn, u, v, w */
#define OUT_NDISPLAY 5 /* shown by the viewer: C, p, p_dyn, speed, vorticity */

/*
 * Output for the web viewer, in <output_dir>/:
 *   frames/frame_NNNNNN.bin      OUT_NSTORED float16 fields (C, p_dyn, u, v, w), cell-centred,
 *                                field-major, x fastest then y then z
 *   surface/surface_NNNNNN.bin   nx*ny float32 surface height eta(x, y), x fastest
 *   manifest.json                run metadata, display ranges, frame times
 * plus the ParaView files of vtk_io.h when vtk_output is set.
 */
typedef struct {
    size_t n, cap;
    double *times;
    double fmin[OUT_NDISPLAY], fmax[OUT_NDISPLAY];
    double eta_min, eta_max; /* surface elevation eta - depth */
    double *p_hydro;         /* flat-surface hydrostatic pressure per layer */
    float *cu, *cv, *cw;     /* cell-centred velocity scratch */
    fs_half *half;
    float *surface;
    double wall_start;
} Output;

int output_init(Output *o, const Sim *s);
int output_frame(Output *o, Sim *s);
int output_manifest(const Output *o, const Sim *s, const char *status, const char *message, const Summary *sum);
void output_free(Output *o);

#endif
