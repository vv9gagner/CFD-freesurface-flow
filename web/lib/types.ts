import type { Dim, RunParams } from "./params";

export type RunStatus = "starting" | "running" | "done" | "failed" | "cancelled" | "interrupted";

export const ACTIVE_STATUSES: RunStatus[] = ["starting", "running"];

export type FieldInfo = {
  name: string;
  label: string;
  unit: string;
  min: number | null;
  max: number | null;
};

export type Summary = {
  period_measured: number | null;
  period_analytic: number | null;
  period_error: number | null;
  crossings: number;
  volume_drift: number | null;
  energy_ratio: number | null;
  max_div: number | null;
};

/** manifest.json written by the solvers ("freesurf-web-1" in 2D, "freesurf3d-web-1" in 3D). */
export type Manifest = {
  format: string;
  dim?: Dim; // missing in manifests written before 3D existed (= 2)
  status: RunStatus;
  message: string | null;
  case: string;
  nx: number;
  nz: number;
  Lx: number;
  Lz: number;
  t: number;
  t_end: number;
  step: number;
  progress: number;
  wall_time: number;
  analytic_period: number | null;
  params: Record<string, number | string | boolean>;
  fields: FieldInfo[];
  frame_bytes: number;
  frame_times: number[];
  summary: Summary | null;
  // 3D only
  ny?: number;
  Ly?: number;
  depth?: number;
  p_hydro?: number[];
  surface_range?: { min: number | null; max: number | null };
  surface_bytes?: number;
};

export type RunInfo = {
  id: string;
  status: RunStatus;
  createdAt: string;
  params: RunParams;
  progress: number;
  frames: number;
  message: string | null;
};

export type RunDetail = RunInfo & { manifest: Manifest | null };
