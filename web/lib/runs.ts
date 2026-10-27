import fs from "node:fs/promises";
import path from "node:path";

import { normalizeParams, type RunParams } from "./params";
import { RUNS_DIR, runDir } from "./paths";
import { isAlive, stopRun } from "./runner";
import type { Manifest, RunDetail, RunInfo, RunStatus } from "./types";

// Server-only: reads run directories written by the solver (manifest.json) and
// by the runner (run.json, status.json).

/** A run that has shown no sign of life for this long is reported as interrupted. */
const STALE_MS = 60_000;

async function readJson<T>(file: string): Promise<T | null> {
  try {
    return JSON.parse(await fs.readFile(file, "utf8")) as T;
  } catch {
    return null;
  }
}

async function mtimeMs(file: string): Promise<number | null> {
  try {
    return (await fs.stat(file)).mtimeMs;
  } catch {
    return null;
  }
}

type RunMeta = { id: string; createdAt: string; params: RunParams };
type StatusFile = { status: RunStatus; message: string };

async function loadRun(id: string): Promise<RunDetail | null> {
  const dir = runDir(id);
  if (!dir) return null;
  const meta = await readJson<RunMeta>(path.join(dir, "run.json"));
  const manifest = await readJson<Manifest>(path.join(dir, "manifest.json"));
  if (!meta && !manifest) return null;
  const statusFile = await readJson<StatusFile>(path.join(dir, "status.json"));
  const alive = isAlive(id);

  let status: RunStatus;
  let message = manifest?.message ?? null;
  if (statusFile && (!manifest || manifest.status === "running")) {
    status = statusFile.status;
    message = statusFile.message;
  } else if (manifest) {
    status = manifest.status;
    if (status === "running" && !alive) {
      const last = await mtimeMs(path.join(dir, "manifest.json"));
      if (last !== null && Date.now() - last > STALE_MS) {
        status = "interrupted";
        message = "the solver stopped without finishing (server restarted?)";
      }
    }
  } else {
    const created = meta ? Date.parse(meta.createdAt) : 0;
    status = alive || Date.now() - created < STALE_MS ? "starting" : "interrupted";
  }

  return {
    id,
    status,
    message,
    createdAt: meta?.createdAt ?? new Date(0).toISOString(),
    params: normalizeParams(meta?.params ?? { dim: manifest?.dim === 3 ? 3 : 2 }),
    progress: status === "done" ? 1 : manifest?.progress ?? 0,
    frames: manifest?.frame_times.length ?? 0,
    manifest,
  };
}

export async function listRuns(): Promise<RunInfo[]> {
  let names: string[];
  try {
    names = await fs.readdir(RUNS_DIR);
  } catch {
    return [];
  }
  const runs = await Promise.all(names.map((n) => loadRun(n)));
  return runs
    .filter((r): r is RunDetail => r !== null)
    .map(
      (r): RunInfo => ({
        id: r.id,
        status: r.status,
        message: r.message,
        createdAt: r.createdAt,
        params: r.params,
        progress: r.progress,
        frames: r.frames,
      }),
    )
    .sort((a, b) => b.createdAt.localeCompare(a.createdAt));
}

export async function getRun(id: string): Promise<RunDetail | null> {
  return loadRun(id);
}

export async function deleteRun(id: string): Promise<boolean> {
  const dir = runDir(id);
  if (!dir) return false;
  await stopRun(id);
  try {
    await fs.access(dir);
  } catch {
    return false;
  }
  await fs.rm(dir, { recursive: true, force: true });
  return true;
}

/** Path of a frame or surface file, or null for an invalid id / index. */
export function framePath(id: string, index: string, kind: "frame" | "surface" = "frame"): string | null {
  const dir = runDir(id);
  if (!dir || !/^\d{1,6}$/.test(index)) return null;
  const folder = kind === "frame" ? "frames" : "surface";
  return path.join(dir, folder, `${kind}_${index.padStart(6, "0")}.bin`);
}

/** The run's manifest alone (cheaper than getRun when only the grid is needed). */
export async function readManifest(id: string): Promise<Manifest | null> {
  const dir = runDir(id);
  return dir ? readJson<Manifest>(path.join(dir, "manifest.json")) : null;
}
