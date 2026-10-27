#ifndef FS3_VTK_IO_H
#define FS3_VTK_IO_H

#include <stddef.h>

#include "sim.h"

/*
 * ParaView output in <output_dir>/vtk/:
 *   fields_NNNNNN.vti  ImageData, cell data C, p, velocity (appended raw Float32)
 *   fields.pvd         time series (rewritten after every frame)
 * Use a Contour filter at C = 0.5 to see the free surface.
 */
int vtk_write_frame(Sim *s, int index, const double *times, size_t nframes);

#endif
