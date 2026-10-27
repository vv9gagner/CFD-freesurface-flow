/** Compact number formatting for field values, axis ticks and readouts. */
export function formatValue(v: number): string {
  const a = Math.abs(v);
  if (!Number.isFinite(v)) return "–";
  if (a === 0) return "0";
  if (a >= 1e4 || a < 1e-3) return v.toExponential(2);
  if (a >= 100) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  if (a >= 1) return v.toFixed(2);
  return v.toFixed(3);
}

/** "Nice" tick positions from 0 to max (steps of 1, 2 or 5 x 10^k). */
export function niceTicks(max: number, target = 5): number[] {
  const raw = max / target;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const norm = raw / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const ticks: number[] = [];
  for (let t = 0; t <= max + step * 1e-6; t += step) ticks.push(Number(t.toPrecision(10)));
  return ticks;
}
