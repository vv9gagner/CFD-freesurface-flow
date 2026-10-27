#define _POSIX_C_SOURCE 200809L
#include "util.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

int mkdir_p(const char *path)
{
    char buf[4096];
    size_t len = strlen(path);
    if (len == 0 || len >= sizeof buf)
        return -1;
    memcpy(buf, path, len + 1);
    for (char *q = buf + 1; *q; q++) {
        if (*q == '/') {
            *q = '\0';
            if (mkdir(buf, 0755) != 0 && errno != EEXIST)
                return -1;
            *q = '/';
        }
    }
    if (mkdir(buf, 0755) != 0 && errno != EEXIST)
        return -1;
    return 0;
}

double wall_seconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int host_is_little_endian(void)
{
    const uint16_t one = 1;
    unsigned char b;
    memcpy(&b, &one, 1);
    return b == 1;
}

void json_write_string(FILE *f, const char *s)
{
    fputc('"', f);
    for (const unsigned char *c = (const unsigned char *)s; *c; c++) {
        switch (*c) {
        case '"': fputs("\\\"", f); break;
        case '\\': fputs("\\\\", f); break;
        case '\n': fputs("\\n", f); break;
        case '\r': fputs("\\r", f); break;
        case '\t': fputs("\\t", f); break;
        default:
            if (*c < 0x20)
                fprintf(f, "\\u%04x", *c);
            else
                fputc(*c, f);
        }
    }
    fputc('"', f);
}

void json_write_number(FILE *f, double v)
{
    if (isfinite(v))
        fprintf(f, "%.10g", v);
    else
        fputs("null", f);
}

FILE *atomic_open(const char *path, char *tmp, size_t tmplen)
{
    int n = snprintf(tmp, tmplen, "%s.tmp", path);
    if (n < 0 || (size_t)n >= tmplen)
        return NULL;
    return fopen(tmp, "wb");
}

int atomic_commit(FILE *f, const char *tmp, const char *path)
{
    int bad = ferror(f);
    if (fclose(f) != 0 || bad) {
        remove(tmp);
        return -1;
    }
    if (rename(tmp, path) != 0) {
        remove(tmp);
        return -1;
    }
    return 0;
}
