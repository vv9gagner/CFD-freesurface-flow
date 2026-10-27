#ifndef FS_UTIL_H
#define FS_UTIL_H

#include <stddef.h>
#include <stdio.h>

#define FS_PI 3.14159265358979323846

static inline double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double sq(double x) { return x * x; }

/* Create a directory and all missing parents (like `mkdir -p`). Returns 0 on success. */
int mkdir_p(const char *path);

/* Monotonic wall-clock time in seconds. */
double wall_seconds(void);

/* 1 if the host stores multi-byte numbers little-endian. */
int host_is_little_endian(void);

/* JSON helpers: strings are escaped, non-finite numbers become null. */
void json_write_string(FILE *f, const char *s);
void json_write_number(FILE *f, double v);

/* Atomic file replacement: write to "<path>.tmp", then rename over <path>.
 * Readers (the web server) never see a half-written file. */
FILE *atomic_open(const char *path, char *tmp, size_t tmplen);
int atomic_commit(FILE *f, const char *tmp, const char *path);

#endif
