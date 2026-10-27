// Slices of a 3D solver frame, returned in the same layout as a 2D frame so the viewer can
// draw them with the 2D heatmap: 5 float32 fields (C, p, p_dyn, speed, vorticity), each
// n1*n2 values with the slice's horizontal axis fastest and its bottom row first.
// Pure functions (no Node or Next imports) so they can be unit-tested with `node --test`.

export type SliceView = "xz" | "yz" | "xy";

export const SLICE_VIEWS: SliceView[] = ["xz", "yz", "xy"];

export type Grid3 = {
  nx: number;
  ny: number;
  nz: number;
  Lx: number;
  Ly: number;
  Lz: number;
  pHydro: ArrayLike<number>; // flat-surface hydrostatic pressure per z layer
};

export type SliceShape = {
  n1: number; // cells along the slice's horizontal axis
  n2: number; // cells along its vertical axis
  L1: number;
  L2: number;
  axes: [string, string]; // axis names, e.g. ["x", "z"]
  normal: "x" | "y" | "z"; // axis the slice position runs along
  positions: number; // number of slice positions along the normal
};

export function sliceShape(g: Grid3, view: SliceView): SliceShape {
  switch (view) {
    case "xz":
      return { n1: g.nx, n2: g.nz, L1: g.Lx, L2: g.Lz, axes: ["x", "z"], normal: "y", positions: g.ny };
    case "yz":
      return { n1: g.ny, n2: g.nz, L1: g.Ly, L2: g.Lz, axes: ["y", "z"], normal: "x", positions: g.nx };
    case "xy":
      return { n1: g.nx, n2: g.ny, L1: g.Lx, L2: g.Ly, axes: ["x", "y"], normal: "z", positions: g.nz };
  }
}

/** Cell-centre coordinate of slice position `at` along the view's normal axis. */
export function slicePosition(g: Grid3, view: SliceView, at: number): number {
  const s = sliceShape(g, view);
  const L = s.normal === "x" ? g.Lx : s.normal === "y" ? g.Ly : g.Lz;
  return ((at + 0.5) * L) / s.positions;
}

/** Bytes of one slice frame as returned to the viewer. */
export function sliceBytes(g: Grid3, view: SliceView): number {
  const s = sliceShape(g, view);
  return 5 * s.n1 * s.n2 * 4;
}

/**
 * Extract slice `at` from a stored frame: 5 fields (C, p_dyn, u, v, w) of nx*ny*nz values,
 * field-major, x fastest then y then z. Vorticity is the component normal to the slice,
 * signed so that positive means counter-clockwise as drawn (horizontal axis right, vertical up):
 * xz shows dw/dx - du/dz, yz shows dw/dy - dv/dz, xy shows dv/dx - du/dy.
 */
export function extractSlice(g: Grid3, frame: ArrayLike<number>, view: SliceView, at: number): Float32Array {
  const { nx, ny, nz } = g;
  const N = nx * ny * nz;
  const s = sliceShape(g, view);
  const { n1, n2 } = s;
  const h1 = s.L1 / n1;
  const h2 = s.L2 / n2;
  const out = new Float32Array(5 * n1 * n2);
  const M = n1 * n2;

  // 3D index of slice cell (a, b)
  const cell = (a: number, b: number): number => {
    if (view === "xz") return (b * ny + at) * nx + a;
    if (view === "yz") return (b * ny + a) * nx + at;
    return (at * ny + b) * nx + a;
  };
  const layer = (b: number): number => (view === "xy" ? at : b);
  // velocity components along the slice's horizontal (a) and vertical (b) axes
  const compA = view === "yz" ? 3 : 2; // u, except v for yz
  const compB = view === "xy" ? 3 : 4; // w, except v for xy
  const vel = (comp: number, a: number, b: number) => frame[comp * N + cell(a, b)];
  // central differences, one-sided at the walls
  const dA = (comp: number, a: number, b: number) =>
    a === 0 ? (vel(comp, 1, b) - vel(comp, 0, b)) / h1
    : a === n1 - 1 ? (vel(comp, a, b) - vel(comp, a - 1, b)) / h1
    : (vel(comp, a + 1, b) - vel(comp, a - 1, b)) / (2 * h1);
  const dB = (comp: number, a: number, b: number) =>
    b === 0 ? (vel(comp, a, 1) - vel(comp, a, 0)) / h2
    : b === n2 - 1 ? (vel(comp, a, b) - vel(comp, a, b - 1)) / h2
    : (vel(comp, a, b + 1) - vel(comp, a, b - 1)) / (2 * h2);

  for (let b = 0; b < n2; b++)
    for (let a = 0; a < n1; a++) {
      const c = cell(a, b);
      const m = b * n1 + a;
      const pDyn = frame[N + c];
      const u = frame[2 * N + c], v = frame[3 * N + c], w = frame[4 * N + c];
      out[m] = frame[c];
      out[M + m] = pDyn + g.pHydro[layer(b)];
      out[2 * M + m] = pDyn;
      out[3 * M + m] = Math.sqrt(u * u + v * v + w * w);
      out[4 * M + m] = dA(compB, a, b) - dB(compA, a, b);
    }
  return out;
}
