#ifndef FS_VTK_IO_H
#define FS_VTK_IO_H

#include <stddef.h>

#include "sim.h"

/*
 * ParaView output in <output_dir>/vtk/:
 *   fields_NNNNNN.vti     ImageData, cell data C, p, rho, vorticity, velocity (appended raw Float32)
 *   interface_NNNNNN.vtp  PLIC interface segments
 *   fields.pvd, interface.pvd  time series (rewritten after every frame)
 * The tank's vertical axis z is written as VTK's y so ParaView's 2D view shows it upright.
 */
int vtk_write_frame(Sim *s, int index, const double *times, size_t nframes);

#endif
