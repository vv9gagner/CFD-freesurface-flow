#include "vtk_io.h"
#include "util.h"
#include "vof.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define VOF_EPS 1e-6

static int write_block(FILE *f, const float *data, size_t count)
{
    uint64_t nbytes = (uint64_t)count * sizeof(float);
    if (fwrite(&nbytes, sizeof nbytes, 1, f) != 1)
        return -1;
    return fwrite(data, sizeof(float), count, f) == count ? 0 : -1;
}

static int write_vti(const Sim *s, const char *path, double t)
{
    const int nx = s->nx, nz = s->nz;
    const size_t n = (size_t)nx * (size_t)nz;
    float *buf = malloc(7 * n * sizeof(float));
    if (!buf)
        return -1;
    float *C = buf, *p = buf + n, *rho = buf + 2 * n, *vort = buf + 3 * n, *vel = buf + 4 * n;
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t m = (size_t)j * nx + i, k = IX(s, i, j);
            C[m] = (float)s->C[k];
            p[m] = (float)sim_cell_pressure(s, i, j);
            rho[m] = (float)s->rho[k];
            vort[m] = (float)sim_cell_vorticity(s, i, j);
            vel[3 * m] = (float)sim_cell_u(s, i, j);
            vel[3 * m + 1] = (float)sim_cell_w(s, i, j);
            vel[3 * m + 2] = 0.0f;
        }

    char tmp[FS_PATH_MAX + 8];
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f) {
        free(buf);
        return -1;
    }
    const char *names[] = { "C", "p", "rho", "vorticity" };
    const uint64_t scalar_block = 8 + n * sizeof(float);
    fprintf(f, "<?xml version=\"1.0\"?>\n");
    fprintf(f, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"%s\" header_type=\"UInt64\">\n",
            host_is_little_endian() ? "LittleEndian" : "BigEndian");
    fprintf(f, "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\"%.10g %.10g 1\">\n",
            nx, nz, s->dx, s->dz);
    fprintf(f, "    <FieldData>\n      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">%.10g</DataArray>\n    </FieldData>\n", t);
    fprintf(f, "    <Piece Extent=\"0 %d 0 %d 0 0\">\n", nx, nz);
    fprintf(f, "      <CellData Scalars=\"C\" Vectors=\"velocity\">\n");
    for (int a = 0; a < 4; a++)
        fprintf(f, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"%llu\"/>\n",
                names[a], (unsigned long long)(a * scalar_block));
    fprintf(f, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%llu\"/>\n",
            (unsigned long long)(4 * scalar_block));
    fprintf(f, "      </CellData>\n    </Piece>\n  </ImageData>\n  <AppendedData encoding=\"raw\">\n_");
    int rc = 0;
    for (int a = 0; a < 4 && rc == 0; a++)
        rc = write_block(f, buf + (size_t)a * n, n);
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

static int write_vtp(Sim *s, const char *path)
{
    vof_reconstruct(s, s->C);
    size_t cap = 256, nseg = 0;
    double *seg = malloc(cap * 4 * sizeof(double));
    if (!seg)
        return -1;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++) {
            size_t k = IX(s, i, j);
            double c = s->C[k], pts[4];
            if (c <= VOF_EPS || c >= 1.0 - VOF_EPS)
                continue;
            if (!plic_segment(s->pn_x[k], s->pn_z[k], s->palpha[k], pts))
                continue;
            if (nseg == cap) {
                double *grown = realloc(seg, 2 * cap * 4 * sizeof(double));
                if (!grown) {
                    free(seg);
                    return -1;
                }
                seg = grown;
                cap *= 2;
            }
            double *q = seg + 4 * nseg++;
            q[0] = (i + pts[0]) * s->dx;
            q[1] = (j + pts[1]) * s->dz;
            q[2] = (i + pts[2]) * s->dx;
            q[3] = (j + pts[3]) * s->dz;
        }

    char tmp[FS_PATH_MAX + 8];
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f) {
        free(seg);
        return -1;
    }
    fprintf(f, "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
    fprintf(f, "  <PolyData>\n    <Piece NumberOfPoints=\"%zu\" NumberOfLines=\"%zu\">\n", 2 * nseg, nseg);
    fprintf(f, "      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (size_t m = 0; m < nseg; m++)
        fprintf(f, "          %.9g %.9g 0 %.9g %.9g 0\n", seg[4 * m], seg[4 * m + 1], seg[4 * m + 2], seg[4 * m + 3]);
    fprintf(f, "        </DataArray>\n      </Points>\n      <Lines>\n");
    fprintf(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for (size_t m = 0; m < nseg; m++)
        fprintf(f, "          %zu %zu\n", 2 * m, 2 * m + 1);
    fprintf(f, "        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for (size_t m = 0; m < nseg; m++)
        fprintf(f, "          %zu\n", 2 * m + 2);
    fprintf(f, "        </DataArray>\n      </Lines>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
    free(seg);
    return atomic_commit(f, tmp, path);
}

static int write_pvd(const char *path, const char *stem, const char *ext, const double *times, size_t n)
{
    char tmp[FS_PATH_MAX + 8];
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f)
        return -1;
    fprintf(f, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n  <Collection>\n");
    for (size_t m = 0; m < n; m++)
        fprintf(f, "    <DataSet timestep=\"%.10g\" file=\"%s_%06zu.%s\"/>\n", times[m], stem, m, ext);
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
    snprintf(path, sizeof path, "%s/interface_%06d.vtp", dir, index);
    if (write_vtp(s, path) != 0)
        return -1;
    snprintf(path, sizeof path, "%s/fields.pvd", dir);
    if (write_pvd(path, "fields", "vti", times, nframes) != 0)
        return -1;
    snprintf(path, sizeof path, "%s/interface.pvd", dir);
    return write_pvd(path, "interface", "vtp", times, nframes);
}
