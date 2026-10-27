import assert from "node:assert/strict";
import { test } from "node:test";

import { extractSlice, sliceBytes, sliceShape, type Grid3 } from "./slices.ts";

const g: Grid3 = { nx: 6, ny: 5, nz: 4, Lx: 1.2, Ly: 1.0, Lz: 0.8, pHydro: [30, 20, 10, 0] };
const N = g.nx * g.ny * g.nz;
const dx = g.Lx / g.nx, dy = g.Ly / g.ny, dz = g.Lz / g.nz;

/** Stored frame (C, p_dyn, u, v, w) built from functions of the cell-centre coordinates. */
function frameOf(f: (x: number, y: number, z: number) => [number, number, number, number, number]) {
  const frame = new Float32Array(5 * N);
  for (let k = 0; k < g.nz; k++)
    for (let j = 0; j < g.ny; j++)
      for (let i = 0; i < g.nx; i++) {
        const c = (k * g.ny + j) * g.nx + i;
        const vals = f((i + 0.5) * dx, (j + 0.5) * dy, (k + 0.5) * dz);
        vals.forEach((v, comp) => (frame[comp * N + c] = v));
      }
  return frame;
}

test("slice shapes and byte sizes", () => {
  assert.deepEqual(sliceShape(g, "xz").axes, ["x", "z"]);
  assert.equal(sliceShape(g, "yz").positions, g.nx);
  assert.equal(sliceBytes(g, "xy"), 5 * g.nx * g.ny * 4);
});

test("scalar fields land at the right cells, p adds the hydrostatic profile", () => {
  const frame = frameOf((x, y, z) => [x + 10 * y + 100 * z, z, 0, 0, 0]);
  const at = 2;
  const s = extractSlice(g, frame, "xz", at);
  const n1 = g.nx, M = g.nx * g.nz;
  for (let k = 0; k < g.nz; k++)
    for (let i = 0; i < g.nx; i++) {
      const x = (i + 0.5) * dx, y = (at + 0.5) * dy, z = (k + 0.5) * dz;
      assert.ok(Math.abs(s[k * n1 + i] - (x + 10 * y + 100 * z)) < 1e-4);
      assert.ok(Math.abs(s[M + k * n1 + i] - (z + g.pHydro[k])) < 1e-4);
    }
  const xy = extractSlice(g, frame, "xy", 1);
  assert.ok(Math.abs(xy[g.nx * g.ny + 0] - (dz * 1.5 + g.pHydro[1])) < 1e-4);
});

test("speed and the normal vorticity of solid-body rotations (vorticity = 2 omega)", () => {
  const om = 3;
  // rotation about y: u = -om (z - zc), w = om (x - xc); as drawn in xz (x right, z up) it is counter-clockwise
  const rotY = frameOf((x, _y, z) => [1, 0, -om * (z - 0.4), 0, om * (x - 0.6)]);
  const xz = extractSlice(g, rotY, "xz", 2);
  const M = g.nx * g.nz;
  for (let m = 0; m < M; m++) assert.ok(Math.abs(xz[4 * M + m] - 2 * om) < 1e-4);
  const x0 = 0.5 * dx, z0 = 0.5 * dz;
  assert.ok(Math.abs(xz[3 * M] - Math.hypot(om * (z0 - 0.4), om * (x0 - 0.6))) < 1e-4);

  // rotation about z: u = -om (y - yc), v = om (x - xc), seen in xy
  const rotZ = frameOf((x, y) => [1, 0, -om * (y - 0.5), om * (x - 0.6), 0]);
  const xy = extractSlice(g, rotZ, "xy", 1);
  const Mxy = g.nx * g.ny;
  for (let m = 0; m < Mxy; m++) assert.ok(Math.abs(xy[4 * Mxy + m] - 2 * om) < 1e-4);

  // rotation about x: v = -om (z - zc), w = om (y - yc), seen in yz (y right, z up)
  const rotX = frameOf((_x, y, z) => [1, 0, 0, -om * (z - 0.4), om * (y - 0.5)]);
  const yz = extractSlice(g, rotX, "yz", 3);
  const Myz = g.ny * g.nz;
  for (let m = 0; m < Myz; m++) assert.ok(Math.abs(yz[4 * Myz + m] - 2 * om) < 1e-4);
});
