#ifndef FS_OUTPUT_H
#define FS_OUTPUT_H

#include <stddef.h>

#include "diagnostics.h"
#include "sim.h"

#define OUT_NFIELDS 5

/*
 * Output for the web viewer, in <output_dir>/:
 *   frames/frame_NNNNNN.bin  OUT_NFIELDS float32 fields (C, p, p_dyn, speed, vorticity),
 *                            field-major, each row-major with x fastest and the bottom row first
 *   manifest.json            run metadata, global field ranges and frame times
 * plus the ParaView files of vtk_io.h when vtk_output is set.
 */
typedef struct {
    size_t n, cap;
    double *times;
    double fmin[OUT_NFIELDS], fmax[OUT_NFIELDS];
    double *p_hydro; /* flat-surface hydrostatic pressure per row */
    float *buf;
    double wall_start;
} Output;

int output_init(Output *o, const Sim *s);
int output_frame(Output *o, Sim *s);
int output_manifest(const Output *o, const Sim *s, const char *status, const char *message, const Summary *sum);
void output_free(Output *o);

#endif
