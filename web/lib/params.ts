// Run parameters editable from the web form. Shared by the client form and the API
// route, so validation is identical on both sides. Fluids are water and air.

export type Dim = 2 | 3;

export type RunParams = {
  dim: Dim;
  nx: number;
  ny: number; // 3D only
  nz: number;
  Lx: number;
  Ly: number; // 3D only
  Lz: number;
  depth: number;
  amplitude: number;
  mode: number; // mode along x (m)
  mode_y: number; // 3D only (n)
  t_end: number;
  output_dt: number;
  wall_bc: "slip" | "noslip";
  vtk_output: boolean;
};

export type NumericKey = Exclude<keyof RunParams, "dim" | "wall_bc" | "vtk_output">;

export type NumericField = {
  key: NumericKey;
  label: string;
  unit: string;
  step: number;
  integer?: boolean;
  /** [min, max] per dimension; a field without limits for a dimension is not used there */
  limits: Partial<Record<Dim, [number, number]>>;
};

export const DEFAULT_PARAMS: Record<Dim, RunParams> = {
  2: {
    dim: 2,
    nx: 128,
    ny: 1,
    nz: 64,
    Lx: 1.0,
    Ly: 0.5,
    Lz: 0.5,
    depth: 0.25,
    amplitude: 0.01,
    mode: 1,
    mode_y: 0,
    t_end: 10,
    output_dt: 0.025,
    wall_bc: "slip",
    vtk_output: false,
  },
  3: {
    dim: 3,
    nx: 64,
    ny: 32,
    nz: 32,
    Lx: 1.0,
    Ly: 0.5,
    Lz: 0.5,
    depth: 0.25,
    amplitude: 0.01,
    mode: 1,
    mode_y: 1,
    t_end: 10,
    output_dt: 0.05,
    wall_bc: "slip",
    vtk_output: false,
  },
};

export const NUMERIC_FIELDS: NumericField[] = [
  { key: "Lx", label: "Tank length", unit: "m", step: 0.1, limits: { 2: [0.1, 10], 3: [0.1, 10] } },
  { key: "Ly", label: "Tank width", unit: "m", step: 0.1, limits: { 3: [0.05, 10] } },
  { key: "Lz", label: "Tank height", unit: "m", step: 0.05, limits: { 2: [0.05, 5], 3: [0.05, 5] } },
  { key: "nx", label: "Cells in x", unit: "", step: 16, integer: true, limits: { 2: [16, 512], 3: [8, 128] } },
  { key: "ny", label: "Cells in y", unit: "", step: 8, integer: true, limits: { 3: [4, 128] } },
  { key: "nz", label: "Cells in z", unit: "", step: 8, integer: true, limits: { 2: [8, 256], 3: [8, 64] } },
  { key: "depth", label: "Water depth", unit: "m", step: 0.01, limits: { 2: [0.01, 5], 3: [0.01, 5] } },
  { key: "amplitude", label: "Amplitude", unit: "m", step: 0.005, limits: { 2: [0, 1], 3: [0, 1] } },
  { key: "mode", label: "Mode m (x)", unit: "", step: 1, integer: true, limits: { 2: [1, 8], 3: [0, 8] } },
  { key: "mode_y", label: "Mode n (y)", unit: "", step: 1, integer: true, limits: { 3: [0, 8] } },
  { key: "t_end", label: "End time", unit: "s", step: 0.5, limits: { 2: [0.1, 60], 3: [0.1, 30] } },
  { key: "output_dt", label: "Frame interval", unit: "s", step: 0.005, limits: { 2: [0.005, 1], 3: [0.01, 1] } },
];

export function fieldsFor(dim: Dim): NumericField[] {
  return NUMERIC_FIELDS.filter((f) => f.limits[dim]);
}

export const MAX_FRAMES: Record<Dim, number> = { 2: 2000, 3: 600 };
const MAX_CELLS: Record<Dim, number> = { 2: 512 * 256, 3: 128 * 64 * 64 };

const RHO_L = 1000;
const RHO_G = 1.2;
const MU_L = 1e-3;
const G = 9.81;

/** Linear standing-wave period for water under an air layer in a closed tank (same formula as the solvers). */
export function analyticPeriod(p: Pick<RunParams, "dim" | "Lx" | "Ly" | "Lz" | "depth" | "mode" | "mode_y">): number {
  const kx = (p.mode * Math.PI) / p.Lx;
  const ky = p.dim === 3 ? (p.mode_y * Math.PI) / p.Ly : 0;
  const k = Math.hypot(kx, ky);
  const den = RHO_L / Math.tanh(k * p.depth) + RHO_G / Math.tanh(k * (p.Lz - p.depth));
  return (2 * Math.PI) / Math.sqrt((G * k * (RHO_L - RHO_G)) / den);
}

export function cellCount(p: Pick<RunParams, "dim" | "nx" | "ny" | "nz">): number {
  return p.dim === 3 ? p.nx * p.ny * p.nz : p.nx * p.nz;
}

export function frameCount(p: Pick<RunParams, "t_end" | "output_dt">): number {
  return Math.ceil(p.t_end / p.output_dt - 1e-9) + 1;
}

/** Disk bytes per frame: 5 float32 fields in 2D; 5 float16 fields plus a float32 surface map in 3D. */
export function frameBytes(p: Pick<RunParams, "dim" | "nx" | "ny" | "nz">): number {
  return p.dim === 3 ? 5 * cellCount(p) * 2 + p.nx * p.ny * 4 : 5 * cellCount(p) * 4;
}

/** Rough wall time from measured throughput (Apple M-series; 3D on 4 threads). */
export function estimateSeconds(p: RunParams): number {
  const h = Math.min(p.Lx / p.nx, p.Lz / p.nz, p.dim === 3 ? p.Ly / p.ny : Infinity);
  const nuMax = MU_L / RHO_G;
  const dt = Math.min(0.005, (p.dim === 3 ? 0.1 : 0.125) * (h * h) / nuMax, 0.5 * Math.sqrt(h / G));
  const perCellStep = p.dim === 3 ? 6.6e-7 : 1.2e-6;
  return cellCount(p) * (p.t_end / dt) * perCellStep;
}

export type ParamErrors = Partial<Record<keyof RunParams | "form", string>>;

export function validateParams(input: unknown): { params: RunParams | null; errors: ParamErrors } {
  if (typeof input !== "object" || input === null) {
    return { params: null, errors: { form: "expected a JSON object" } };
  }
  const raw = input as Record<string, unknown>;
  const dimRaw = raw.dim === undefined ? 2 : Number(raw.dim);
  if (dimRaw !== 2 && dimRaw !== 3) return { params: null, errors: { dim: "must be 2 or 3" } };
  const dim = dimRaw as Dim;
  const out: RunParams = { ...DEFAULT_PARAMS[dim] };
  const errors: ParamErrors = {};

  for (const f of fieldsFor(dim)) {
    if (raw[f.key] === undefined) continue;
    const v = typeof raw[f.key] === "string" ? Number(raw[f.key]) : raw[f.key];
    const [min, max] = f.limits[dim]!;
    if (typeof v !== "number" || !Number.isFinite(v)) errors[f.key] = "must be a number";
    else if (f.integer && !Number.isInteger(v)) errors[f.key] = "must be an integer";
    else if (v < min || v > max) errors[f.key] = `must be between ${min} and ${max}`;
    else out[f.key] = v;
  }
  if (raw.wall_bc !== undefined) {
    if (raw.wall_bc === "slip" || raw.wall_bc === "noslip") out.wall_bc = raw.wall_bc;
    else errors.wall_bc = "must be slip or noslip";
  }
  if (raw.vtk_output !== undefined) {
    if (typeof raw.vtk_output === "boolean") out.vtk_output = raw.vtk_output;
    else errors.vtk_output = "must be true or false";
  }
  if (Object.keys(errors).length) return { params: null, errors };

  if (out.depth - out.amplitude <= 0) errors.amplitude = "must be smaller than the depth";
  else if (out.depth + out.amplitude >= out.Lz) errors.depth = "depth + amplitude must stay below the tank height";
  const dz = out.Lz / out.nz;
  const aspectX = out.Lx / out.nx / dz;
  if (aspectX < 0.5 || aspectX > 2) errors.nz = "cells must be roughly cubic (dx/dz between 0.5 and 2)";
  if (dim === 3) {
    const aspectY = out.Ly / out.ny / dz;
    if (aspectY < 0.5 || aspectY > 2) errors.ny = "cells must be roughly cubic (dy/dz between 0.5 and 2)";
    if (out.mode === 0 && out.mode_y === 0) errors.mode_y = "modes m and n cannot both be 0";
  }
  if (cellCount(out) > MAX_CELLS[dim]) errors.nx = `grid too large (max ${MAX_CELLS[dim].toLocaleString("en")} cells)`;
  if (frameCount(out) > MAX_FRAMES[dim]) errors.output_dt = `too many frames (max ${MAX_FRAMES[dim]})`;

  return Object.keys(errors).length ? { params: null, errors } : { params: out, errors };
}

/** Command-line overrides for the solver of the run's dimension. */
export function solverArgs(p: RunParams): string[] {
  const common = [
    `nx=${p.nx}`,
    `nz=${p.nz}`,
    `Lx=${p.Lx}`,
    `Lz=${p.Lz}`,
    `depth=${p.depth}`,
    `amplitude=${p.amplitude}`,
    `t_end=${p.t_end}`,
    `output_dt=${p.output_dt}`,
    `wall_bc=${p.wall_bc}`,
    `vtk_output=${p.vtk_output ? 1 : 0}`,
  ];
  if (p.dim === 2) return [...common, `mode=${p.mode}`];
  return [...common, `ny=${p.ny}`, `Ly=${p.Ly}`, `mode_x=${p.mode}`, `mode_y=${p.mode_y}`];
}

/** Older runs were saved before 3D existed: fill in what they lack. */
export function normalizeParams(p: Partial<RunParams> | undefined): RunParams {
  const dim: Dim = p?.dim === 3 ? 3 : 2;
  return { ...DEFAULT_PARAMS[dim], ...p, dim };
}
