"use client";

import { useEffect, useRef, useState } from "react";

import { formatValue, niceTicks } from "@/lib/format";

// A plane of cells: nx along the horizontal axis, nz along the vertical axis (for 3D slices
// these are the slice's own axes, named by `axes`, e.g. ["y", "z"] or ["x", "y"]).
type Props = {
  values: Float32Array; // one field, nx*nz, horizontal axis fastest, bottom row first
  liquid: Float32Array; // liquid fraction C, for the surface line
  nx: number;
  nz: number;
  Lx: number;
  Lz: number;
  axes?: [string, string];
  lut: Uint8ClampedArray;
  vmin: number;
  vmax: number;
  smooth: boolean;
  showSurface: boolean; // only meaningful when the vertical axis is z
  dark: boolean;
  label: string;
  unit: string;
};

type Hover = { i: number; j: number; px: number; py: number };

export default function FieldCanvas(props: Props) {
  const { values, liquid, nx, nz, Lx, Lz, lut, vmin, vmax, smooth, showSurface, dark, label, unit } = props;
  const [axH, axV] = props.axes ?? ["x", "z"];
  const wrapRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const offRef = useRef<HTMLCanvasElement | null>(null);
  const [size, setSize] = useState({ w: 0, h: 0 });
  const [hover, setHover] = useState<Hover | null>(null);

  // fit the tank's aspect ratio into the available width (and at most 70% of the viewport height)
  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const ro = new ResizeObserver(([entry]) => {
      const maxH = window.innerHeight * 0.7;
      let w = entry.contentRect.width;
      let h = (w * Lz) / Lx;
      if (h > maxH) {
        h = maxH;
        w = (h * Lx) / Lz;
      }
      setSize({ w: Math.floor(w), h: Math.floor(h) });
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, [Lx, Lz]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || size.w === 0) return;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(size.w * dpr);
    canvas.height = Math.round(size.h * dpr);

    const off = (offRef.current ??= document.createElement("canvas"));
    off.width = nx;
    off.height = nz;
    const octx = off.getContext("2d");
    const ctx = canvas.getContext("2d");
    if (!octx || !ctx) return;

    const img = octx.createImageData(nx, nz);
    const span = vmax - vmin || 1;
    for (let j = 0; j < nz; j++) {
      const row = (nz - 1 - j) * nx; // image rows run top-down, the tank bottom-up
      for (let i = 0; i < nx; i++) {
        let t = (values[j * nx + i] - vmin) / span;
        t = t < 0 ? 0 : t > 1 ? 1 : t;
        const src = Math.round(t * 255) * 4;
        const dst = (row + i) * 4;
        img.data[dst] = lut[src];
        img.data[dst + 1] = lut[src + 1];
        img.data[dst + 2] = lut[src + 2];
        img.data[dst + 3] = 255;
      }
    }
    octx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = smooth;
    ctx.drawImage(off, 0, 0, canvas.width, canvas.height);

    if (showSurface) {
      // surface height per column, eta_i = dz * sum_j C(i, j)
      const W = canvas.width;
      const H = canvas.height;
      ctx.beginPath();
      for (let i = 0; i < nx; i++) {
        let col = 0;
        for (let j = 0; j < nz; j++) col += liquid[j * nx + i];
        const x = ((i + 0.5) / nx) * W;
        const y = H - (col / nz) * H;
        if (i === 0) ctx.moveTo(0, y);
        ctx.lineTo(x, y);
        if (i === nx - 1) ctx.lineTo(W, y);
      }
      ctx.lineJoin = "round";
      ctx.strokeStyle = dark ? "#1a1a19" : "#fcfcfb";
      ctx.lineWidth = 4 * dpr;
      ctx.stroke();
      ctx.strokeStyle = dark ? "#f0efec" : "#0b0b0b";
      ctx.lineWidth = 1.5 * dpr;
      ctx.stroke();
    }
  }, [values, liquid, nx, nz, lut, vmin, vmax, smooth, showSurface, dark, size]);

  function onPointerMove(e: React.PointerEvent<HTMLCanvasElement>) {
    const r = e.currentTarget.getBoundingClientRect();
    const px = e.clientX - r.left;
    const py = e.clientY - r.top;
    const i = Math.min(nx - 1, Math.max(0, Math.floor((px / r.width) * nx)));
    const j = Math.min(nz - 1, Math.max(0, nz - 1 - Math.floor((py / r.height) * nz)));
    setHover({ i, j, px, py });
  }

  const xTicks = niceTicks(Lx);
  const zTicks = niceTicks(Lz, 4);
  const cellW = size.w / nx;
  const cellH = size.h / nz;
  const hoverValue = hover ? values[hover.j * nx + hover.i] : 0;

  return (
    <div className="grid grid-cols-[2.75rem_1fr] grid-rows-[auto_auto]">
      <div className="relative" style={{ height: size.h }}>
        {zTicks.map((t) => (
          <span
            key={t}
            className="absolute right-1.5 -translate-y-1/2 text-[11px] tabular-nums text-muted"
            style={{ top: size.h - (t / Lz) * size.h }}
          >
            {t}
          </span>
        ))}
        <span className="absolute top-1/2 -left-1 -translate-y-1/2 -rotate-90 text-[11px] text-muted">{axV} [m]</span>
      </div>
      <div ref={wrapRef} className="relative min-w-0">
        <canvas
          ref={canvasRef}
          style={{ width: size.w, height: size.h }}
          className="block cursor-crosshair rounded-sm ring-1 ring-line"
          onPointerMove={onPointerMove}
          onPointerLeave={() => setHover(null)}
          role="img"
          aria-label={`${label} field on a ${nx} by ${nz} grid`}
        />
        {hover && (
          <>
            <div
              className="pointer-events-none absolute ring-1 ring-ink"
              style={{
                left: hover.i * cellW,
                top: (nz - 1 - hover.j) * cellH,
                width: Math.max(cellW, 2),
                height: Math.max(cellH, 2),
              }}
            />
            <div
              className="pointer-events-none absolute z-10 rounded-md border border-line bg-surface px-2 py-1 text-xs whitespace-nowrap shadow-sm"
              style={{
                left: hover.px + (hover.px > size.w - 180 ? -12 : 12),
                top: hover.py + (hover.py > size.h - 60 ? -48 : 12),
                transform: hover.px > size.w - 180 ? "translateX(-100%)" : undefined,
              }}
            >
              <div className="font-medium tabular-nums text-ink">
                {label}: {formatValue(hoverValue)} {unit}
              </div>
              <div className="tabular-nums text-ink-2">
                {axH} {(((hover.i + 0.5) * Lx) / nx).toFixed(3)} m · {axV} {(((hover.j + 0.5) * Lz) / nz).toFixed(3)} m
              </div>
            </div>
          </>
        )}
      </div>
      <div />
      <div className="relative h-8" style={{ width: size.w }}>
        {xTicks.map((t) => (
          <span
            key={t}
            className="absolute top-1 -translate-x-1/2 text-[11px] tabular-nums text-muted"
            style={{ left: (t / Lx) * size.w }}
          >
            {t}
          </span>
        ))}
        <span className="absolute right-0 bottom-0 text-[11px] text-muted">{axH} [m]</span>
      </div>
    </div>
  );
}
