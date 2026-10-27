import { cssGradient, type ScaleKind } from "@/lib/colormaps";
import { formatValue } from "@/lib/format";

type Props = { kind: ScaleKind; dark: boolean; vmin: number; vmax: number; label: string; unit: string };

export default function ColorLegend({ kind, dark, vmin, vmax, label, unit }: Props) {
  return (
    <div className="w-full max-w-xs">
      <div className="mb-1 text-xs text-ink-2">
        {label}
        {unit && <span className="text-muted"> [{unit}]</span>}
      </div>
      <div className="h-2.5 rounded-sm ring-1 ring-line" style={{ background: cssGradient(kind, dark) }} />
      <div className="mt-0.5 flex justify-between text-[11px] tabular-nums text-muted">
        <span>{formatValue(vmin)}</span>
        {kind === "diverging" && <span>0</span>}
        <span>{formatValue(vmax)}</span>
      </div>
    </div>
  );
}
