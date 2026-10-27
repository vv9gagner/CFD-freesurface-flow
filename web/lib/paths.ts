import path from "node:path";

import type { Dim } from "./params";

// Server-only: where the solvers, their case files and run outputs live.
// `next dev` runs from web/, so the project root defaults to its parent.
export const ROOT = path.resolve(process.env.FREESURF_ROOT ?? path.join(process.cwd(), ".."));
export const RUNS_DIR = path.join(ROOT, "runs");

export const SOLVERS: Record<Dim, { app: string; bin: string; caseFile: string }> = {
  2: {
    app: path.join(ROOT, "freesurf-2d"),
    bin: path.join(ROOT, "freesurf-2d", "build", "freesurf"),
    caseFile: path.join(ROOT, "freesurf-2d", "cases", "standing_wave.cfg"),
  },
  3: {
    app: path.join(ROOT, "freesurf-3d"),
    bin: path.join(ROOT, "freesurf-3d", "build", "freesurf3d"),
    caseFile: path.join(ROOT, "freesurf-3d", "cases", "standing_wave_3d.cfg"),
  },
};

const RUN_ID_RE = /^[a-z0-9-]{1,64}$/;

/** Validates a run id before it touches the filesystem (no path traversal). */
export function runDir(id: string): string | null {
  return RUN_ID_RE.test(id) ? path.join(RUNS_DIR, id) : null;
}
