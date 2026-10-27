#include "vtk_io.h"
#include "util.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int write_block(FILE *f, const float *data, size_t count)
{
    uint64_t nbytes = (uint64_t)count * sizeof(float);
    if (fwrite(&nbytes, sizeof nbytes, 1, f) != 1)
        return -1;
    return fwrite(data, sizeof(float), count, f) == count ? 0 : -1;
}

static int write_vti(const Sim *s, const char *path, double t)
{
    const int nx = s->nx, ny = s->ny, nz = s->nz;
    const size_t n = (size_t)nx * ny * nz;
    float *buf = malloc(5 * n * sizeof(float));
    if (!buf)
        return -1;
    float *C = buf, *p = buf + n, *vel = buf + 2 * n;
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t m = ((size_t)k * ny + j) * nx + i;
                C[m] = (float)s->C[IX(s, i, j, k)];
                p[m] = (float)sim_cell_pressure(s, i, j, k);
                vel[3 * m] = (float)sim_cell_u(s, i, j, k);
                vel[3 * m + 1] = (float)sim_cell_v(s, i, j, k);
                vel[3 * m + 2] = (float)sim_cell_w(s, i, j, k);
            }

    char tmp[FS_PATH_MAX + 80];
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f) {
        free(buf);
        return -1;
    }
    const uint64_t block = 8 + n * sizeof(float);
    fprintf(f, "<?xml version=\"1.0\"?>\n");
    fprintf(f, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"%s\" header_type=\"UInt64\">\n",
            host_is_little_endian() ? "LittleEndian" : "BigEndian");
    fprintf(f, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"%.10g %.10g %.10g\">\n", nx, ny,
            nz, s->dx, s->dy, s->dz);
    fprintf(f, "    <FieldData>\n      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">%.10g</DataArray>\n    </FieldData>\n", t);
    fprintf(f, "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n", nx, ny, nz);
    fprintf(f, "      <CellData Scalars=\"C\" Vectors=\"velocity\">\n");
    fprintf(f, "        <DataArray type=\"Float32\" Name=\"C\" format=\"appended\" offset=\"0\"/>\n");
    fprintf(f, "        <DataArray type=\"Float32\" Name=\"p\" format=\"appended\" offset=\"%llu\"/>\n", (unsigned long long)block);
    fprintf(f, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%llu\"/>\n",
            (unsigned long long)(2 * block));
    fprintf(f, "      </CellData>\n    </Piece>\n  </ImageData>\n  <AppendedData encoding=\"raw\">\n_");
    int rc = write_block(f, C, n);
    if (rc == 0)
        rc = write_block(f, p, n);
    if (rc == 0)
        rc = write_block(f, vel, 3 * n);
    fprintf(f, "\n  </AppendedData>\n</VTKFile>\n");
    free(buf);
    if (rc != 0) {
        fclose(f);
        remove(tmp);
        return -1;
    }
    return atomic_commit(f, tmp, path);
}

static int write_pvd(const char *path, const double *times, size_t n)
{
    char tmp[FS_PATH_MAX + 80];
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f)
        return -1;
    fprintf(f, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n  <Collection>\n");
    for (size_t m = 0; m < n; m++)
        fprintf(f, "    <DataSet timestep=\"%.10g\" file=\"fields_%06zu.vti\"/>\n", times[m], m);
    fprintf(f, "  </Collection>\n</VTKFile>\n");
    return atomic_commit(f, tmp, path);
}

int vtk_write_frame(Sim *s, int index, const double *times, size_t nframes)
{
    char dir[FS_PATH_MAX + 8], path[FS_PATH_MAX + 72];
    snprintf(dir, sizeof dir, "%s/vtk", s->prm.output_dir);
    if (mkdir_p(dir) != 0)
        return -1;
    snprintf(path, sizeof path, "%s/fields_%06d.vti", dir, index);
    if (write_vti(s, path, s->t) != 0)
        return -1;
    snprintf(path, sizeof path, "%s/fields.pvd", dir);
    return write_pvd(path, times, nframes);
}
