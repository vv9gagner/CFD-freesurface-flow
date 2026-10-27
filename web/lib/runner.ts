import { spawn, type ChildProcess } from "node:child_process";
import { randomBytes } from "node:crypto";
import fs from "node:fs";
import path from "node:path";

import { solverArgs, type Dim, type RunParams } from "./params";
import { ROOT, RUNS_DIR, SOLVERS, runDir } from "./paths";

// Server-only: launches the C solvers as child processes and tracks the live ones.
// Route handlers can be bundled separately, so the registry lives on globalThis.

export const MAX_CONCURRENT = 2;
/** A 3D run already uses several cores; running two at once only slows both down. */
export const MAX_CONCURRENT_3D = 1;

type Registry = Map<string, { child: ChildProcess; dim: Dim }>;
const g = globalThis as typeof globalThis & { __freesurfRuns?: Registry };
const registry: Registry = (g.__freesurfRuns ??= new Map());

export function isAlive(id: string): boolean {
  return registry.has(id);
}

function makeId(): string {
  const d = new Date();
  const p = (n: number) => String(n).padStart(2, "0");
  const stamp = `${d.getFullYear()}${p(d.getMonth() + 1)}${p(d.getDate())}-${p(d.getHours())}${p(d.getMinutes())}${p(d.getSeconds())}`;
  return `${stamp}-${randomBytes(2).toString("hex")}`;
}

function logTail(dir: string, lines = 12): string {
  try {
    return fs.readFileSync(path.join(dir, "log.txt"), "utf8").trimEnd().split("\n").slice(-lines).join("\n");
  } catch {
    return "";
  }
}

function writeStatus(dir: string, status: string, message: string) {
  try {
    fs.writeFileSync(path.join(dir, "status.json"), JSON.stringify({ status, message }, null, 2));
  } catch {
    // the run directory may already have been deleted
  }
}

export type StartResult = { id: string } | { error: string; httpStatus: number };

export function startRun(params: RunParams): StartResult {
  const solver = SOLVERS[params.dim];
  try {
    fs.accessSync(solver.bin, fs.constants.X_OK);
  } catch {
    const app = path.relative(ROOT, solver.app);
    return { error: `Solver not found at ${solver.bin}. Run \`make -C ${app}\` in ${ROOT} first.`, httpStatus: 500 };
  }
  if (registry.size >= MAX_CONCURRENT) {
    return { error: `At most ${MAX_CONCURRENT} simulations can run at once.`, httpStatus: 409 };
  }
  const running3d = [...registry.values()].filter((r) => r.dim === 3).length;
  if (params.dim === 3 && running3d >= MAX_CONCURRENT_3D) {
    return { error: "A 3D simulation is already running; wait for it to finish.", httpStatus: 409 };
  }

  const id = makeId();
  const dir = runDir(id)!;
  fs.mkdirSync(RUNS_DIR, { recursive: true });
  fs.mkdirSync(dir);
  fs.writeFileSync(path.join(dir, "run.json"), JSON.stringify({ id, createdAt: new Date().toISOString(), params }, null, 2));

  const log = fs.openSync(path.join(dir, "log.txt"), "a");
  let child: ChildProcess;
  try {
    child = spawn(solver.bin, [solver.caseFile, ...solverArgs(params), `output_dir=${dir}`, "print_every=200"], {
      cwd: ROOT,
      stdio: ["ignore", log, log],
    });
  } finally {
    fs.closeSync(log);
  }

  registry.set(id, { child, dim: params.dim });
  child.on("error", (err) => {
    registry.delete(id);
    writeStatus(dir, "failed", `could not start solver: ${err.message}`);
  });
  child.on("exit", (code, signal) => {
    registry.delete(id);
    // exit 0 = done and 3 = cancelled both leave a final manifest; anything else gets a status file
    if (signal) writeStatus(dir, "cancelled", `terminated by ${signal}`);
    else if (code !== 0 && code !== 3) writeStatus(dir, "failed", logTail(dir) || `solver exited with code ${code}`);
  });
  return { id };
}

/** Ask a running solver to stop (it writes a final "cancelled" manifest); resolves once it has exited. */
export function stopRun(id: string, timeoutMs = 5000): Promise<void> {
  const entry = registry.get(id);
  if (!entry) return Promise.resolve();
  const { child } = entry;
  return new Promise((resolve) => {
    const timer = setTimeout(() => {
      child.kill("SIGKILL");
      resolve();
    }, timeoutMs);
    child.once("exit", () => {
      clearTimeout(timer);
      resolve();
    });
    child.kill("SIGTERM");
  });
}
