// Colour scales for the field heatmap.
//   sequential (magnitude: C, p, speed): one blue hue, evenly stepped in OKLab lightness
//   diverging  (sign: p_dyn, vorticity): blue <-> red with a neutral gray midpoint
// Dark mode uses its own steps: near-zero recedes toward the dark surface.

export type ScaleKind = "sequential" | "diverging";

const SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"];
const DIVERGING_LIGHT = ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#f0efec", "#f8d2ce", "#ec8d88", "#d03f3e", "#7d1a19"];
const DIVERGING_DARK = ["#cde2fb", "#6da7ec", "#2a78d6", "#1c4f8f", "#383835", "#8a2a28", "#d03f3e", "#ec8d88", "#f8d2ce"];

export function scaleKind(field: string): ScaleKind {
  return field === "p_dyn" || field === "vorticity" ? "diverging" : "sequential";
}

export function rampStops(kind: ScaleKind, dark: boolean): string[] {
  if (kind === "diverging") return dark ? DIVERGING_DARK : DIVERGING_LIGHT;
  return dark ? [...SEQUENTIAL].reverse() : SEQUENTIAL;
}

function hexToRgb(hex: string): [number, number, number] {
  return [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16)) as [number, number, number];
}

const lutCache = new Map<string, Uint8ClampedArray>();

/** 256-entry RGBA lookup table, piecewise linear between evenly spaced stops. */
export function colorLut(kind: ScaleKind, dark: boolean): Uint8ClampedArray {
  const key = `${kind}-${dark}`;
  const cached = lutCache.get(key);
  if (cached) return cached;
  const stops = rampStops(kind, dark).map(hexToRgb);
  const lut = new Uint8ClampedArray(256 * 4);
  for (let n = 0; n < 256; n++) {
    const x = (n / 255) * (stops.length - 1);
    const k = Math.min(Math.floor(x), stops.length - 2);
    const f = x - k;
    for (let c = 0; c < 3; c++) lut[4 * n + c] = Math.round(stops[k][c] * (1 - f) + stops[k + 1][c] * f);
    lut[4 * n + 3] = 255;
  }
  lutCache.set(key, lut);
  return lut;
}

export function cssGradient(kind: ScaleKind, dark: boolean): string {
  return `linear-gradient(to right, ${rampStops(kind, dark).join(", ")})`;
}
