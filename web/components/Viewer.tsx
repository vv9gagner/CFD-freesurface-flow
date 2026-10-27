"use client";

import dynamic from "next/dynamic";
import { useEffect, useMemo, useRef, useState, useSyncExternalStore } from "react";

import { colorLut, scaleKind } from "@/lib/colormaps";
import { formatValue } from "@/lib/format";
import { FrameStore, fieldSlice, sliceRange } from "@/lib/frames";
import { SLICE_VIEWS, sliceShape, slicePosition, type Grid3, type SliceView } from "@/lib/slices";
import { ACTIVE_STATUSES, type FieldInfo, type Manifest, type RunDetail } from "@/lib/types";
import ColorLegend from "./ColorLegend";
import FieldCanvas from "./FieldCanvas";
import StatusBadge from "./StatusBadge";
import TimeBar from "./TimeBar";

// WebGL only exists in the browser
const Surface3D = dynamic(() => import("./Surface3D"), {
  ssr: false,
  loading: () => <p className="py-16 text-center text-sm text-muted">Loading 3D view…</p>,
});

const FIELD_HELP: Record<string, string> = {
  C: "Liquid volume fraction: 1 in water, 0 in air",
  p: "Pressure relative to a top corner of the tank",
  p_dyn: "Pressure minus the hydrostatic pressure of the flat surface",
  speed: "Magnitude of the velocity",
  vorticity: "Vorticity normal to the plane shown (positive = counter-clockwise as drawn)",
};

type View = "surface" | SliceView;
const VIEW_LABEL: Record<View, string> = { surface: "3D surface", xz: "Slice x–z", yz: "Slice y–z", xy: "Slice x–y" };

function useDarkMode(): boolean {
  return useSyncExternalStore(
    (cb) => {
      const mq = window.matchMedia("(prefers-color-scheme: dark)");
      mq.addEventListener("change", cb);
      return () => mq.removeEventListener("change", cb);
    },
    () => window.matchMedia("(prefers-color-scheme: dark)").matches,
    () => false,
  );
}

const noSubscribe = () => () => {};
const zero = () => 0;

/** One frame cache per run and per view (each slice position has its own frames). */
function useFrameStore(runId: string, frameBytes: number, query: string): FrameStore | null {
  return useMemo(() => (frameBytes ? new FrameStore(runId, frameBytes, query) : null), [runId, frameBytes, query]);
}

type Plane = { n1: number; n2: number; L1: number; L2: number; axes: [string, string] };

/** The plane drawn on screen (none for the 3D surface) and how its frames are fetched. */
function displayOf(
  manifest: Manifest | null,
  grid: Grid3 | null,
  sliceView: SliceView | null,
  at: number,
): { plane: Plane | null; query: string; frameBytes: number } {
  if (manifest && !grid)
    return {
      plane: { n1: manifest.nx, n2: manifest.nz, L1: manifest.Lx, L2: manifest.Lz, axes: ["x", "z"] },
      query: "",
      frameBytes: manifest.frame_bytes,
    };
  if (grid && sliceView) {
    const s = sliceShape(grid, sliceView);
    return {
      plane: { n1: s.n1, n2: s.n2, L1: s.L1, L2: s.L2, axes: s.axes },
      query: `?view=${sliceView}&at=${at}`,
      frameBytes: 5 * s.n1 * s.n2 * 4,
    };
  }
  if (grid) return { plane: null, query: "?view=surface", frameBytes: manifest?.surface_bytes ?? grid.nx * grid.ny * 4 };
  return { plane: null, query: "", frameBytes: 0 };
}

function downloadCsv(name: string, values: Float32Array, pl: Plane) {
  const lines = [`${pl.axes[0]},${pl.axes[1]},${name}`];
  for (let j = 0; j < pl.n2; j++)
    for (let i = 0; i < pl.n1; i++)
      lines.push(
        `${(((i + 0.5) * pl.L1) / pl.n1).toFixed(6)},${(((j + 0.5) * pl.L2) / pl.n2).toFixed(6)},${values[j * pl.n1 + i]}`,
      );
  const url = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/csv" }));
  const a = document.createElement("a");
  a.href = url;
  a.download = `${name}.csv`;
  a.click();
  URL.revokeObjectURL(url);
}

export default function Viewer({ runId }: { runId: string }) {
  const [detail, setDetail] = useState<RunDetail | null>(null);
  const [loadError, setLoadError] = useState<string | null>(null);
  const [frameError, setFrameError] = useState<string | null>(null);
  const [index, setIndex] = useState(0);
  const [follow, setFollow] = useState(true);
  const [playing, setPlaying] = useState(false);
  const [speed, setSpeed] = useState(1);
  const [fieldName, setFieldName] = useState("C");
  const [rangeMode, setRangeMode] = useState<"run" | "frame">("run");
  const [smooth, setSmooth] = useState(true);
  const [showSurface, setShowSurface] = useState(true);
  const [view, setView] = useState<View>("surface");
  const [sliceAt, setSliceAt] = useState<Partial<Record<SliceView, number>>>({});
  const [exaggeration, setExaggeration] = useState<number | null>(null); // null = automatic
  const dark = useDarkMode();

  const followRef = useRef(follow);
  useEffect(() => {
    followRef.current = follow;
  }, [follow]);

  // run metadata; polled while the solver is still working
  useEffect(() => {
    let stopped = false;
    let timer: ReturnType<typeof setTimeout>;
    const tick = async () => {
      try {
        const res = await fetch(`/api/runs/${runId}`, { cache: "no-store" });
        if (res.status === 404) {
          if (!stopped) setLoadError("This run does not exist (it may have been deleted).");
          return;
        }
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const d = (await res.json()) as RunDetail;
        if (stopped) return;
        setDetail(d);
        setLoadError(null);
        const active = ACTIVE_STATUSES.includes(d.status);
        if (active && followRef.current && d.frames > 0) setIndex(d.frames - 1);
        if (active) timer = setTimeout(tick, 1000);
      } catch (e) {
        if (stopped) return;
        setLoadError(`Could not load run: ${(e as Error).message}`);
        timer = setTimeout(tick, 3000);
      }
    };
    tick();
    return () => {
      stopped = true;
      clearTimeout(timer);
    };
  }, [runId]);

  const manifest = detail?.manifest ?? null;
  const active = detail ? ACTIVE_STATUSES.includes(detail.status) : false;
  const nFrames = manifest?.frame_times.length ?? 0;
  const last = Math.max(0, nFrames - 1);
  const shown = Math.min(Math.max(index, 0), last);

  // what is on screen: the 2D field, a 3D slice, or the 3D surface
  const is3d = manifest?.dim === 3 && manifest.ny !== undefined && manifest.Ly !== undefined;
  const grid: Grid3 | null =
    is3d && manifest
      ? { nx: manifest.nx, ny: manifest.ny!, nz: manifest.nz, Lx: manifest.Lx, Ly: manifest.Ly!, Lz: manifest.Lz, pHydro: [] }
      : null;
  const depth = manifest?.depth ?? detail?.params.depth ?? 0;
  const defaultAt = (v: SliceView): number => {
    if (!grid) return 0;
    if (v === "xz") return Math.floor(grid.ny / 2);
    if (v === "yz") return Math.floor(grid.nx / 2);
    return Math.min(grid.nz - 1, Math.max(0, Math.floor(depth / (grid.Lz / grid.nz)) - 1)); // just below the surface
  };
  const sliceView: SliceView | null = grid && view !== "surface" ? view : null;
  const at = sliceView ? (sliceAt[sliceView] ?? defaultAt(sliceView)) : 0;

  const { plane, query, frameBytes } = displayOf(manifest, grid, sliceView, at);
  const store = useFrameStore(runId, frameBytes, query);
  // re-render whenever a frame lands in the cache
  useSyncExternalStore(store?.subscribe ?? noSubscribe, store?.getVersion ?? zero, zero);
  const current = store ? store.peekNearest(shown) : null;

  // the frame on screen
  useEffect(() => {
    if (!store || nFrames === 0) return;
    store.load(shown).catch((e: Error) => setFrameError(e.message));
  }, [store, shown, nFrames]);

  // background prefetch of a window of frames around the current one (all of them if they fit)
  const anchor = Math.floor(shown / 32) * 32;
  useEffect(() => {
    if (!store || nFrames === 0) return;
    let stopped = false;
    const cap = store.capacity;
    const lo = nFrames <= cap ? 0 : Math.max(0, Math.min(anchor - Math.floor(cap / 4), nFrames - cap));
    const hi = Math.min(nFrames, lo + cap);
    const queue: number[] = [];
    for (let i = anchor; i < hi; i++) if (!store.has(i)) queue.push(i);
    for (let i = Math.min(anchor, hi) - 1; i >= lo; i--) if (!store.has(i)) queue.push(i);
    const worker = async () => {
      while (!stopped && queue.length) {
        const i = queue.shift()!;
        try {
          await store.load(i);
        } catch {
          // a missing frame shows up when it is displayed
        }
      }
    };
    for (let w = 0; w < 4; w++) worker();
    return () => {
      stopped = true;
    };
  }, [store, nFrames, anchor]);

  // playback in simulation time: `speed` simulated seconds per wall-clock second
  const indexRef = useRef(shown);
  useEffect(() => {
    indexRef.current = shown;
  }, [shown]);
  const times = manifest?.frame_times;
  useEffect(() => {
    if (!playing || !store || !times) return;
    let raf = 0;
    let prev = performance.now();
    let acc = 0;
    const step = (now: number) => {
      acc += ((now - prev) / 1000) * speed;
      prev = now;
      const i = indexRef.current;
      if (i >= times.length - 1) {
        if (!active) {
          setPlaying(false);
          return;
        }
        acc = 0; // live run: wait for the next frame
      } else {
        const dt = times[i + 1] - times[i];
        if (acc >= dt) {
          if (store.has(i + 1)) {
            acc = Math.min(acc - dt, dt);
            indexRef.current = i + 1;
            setIndex(i + 1);
          } else {
            store.load(i + 1).catch(() => {});
          }
        }
      }
      raf = requestAnimationFrame(step);
    };
    raf = requestAnimationFrame(step);
    return () => cancelAnimationFrame(raf);
  }, [playing, speed, store, times, active]);

  function seek(i: number) {
    setFollow(false);
    setIndex(Math.min(Math.max(i, 0), last));
  }

  function togglePlay() {
    if (playing) {
      setPlaying(false);
      return;
    }
    setFollow(false);
    if (shown >= last && !active) setIndex(0);
    setPlaying(true);
  }

  // keyboard: ← → step, space play/pause (ignored while typing in form fields)
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const el = e.target as HTMLElement | null;
      const tag = el?.tagName;
      const isRange = tag === "INPUT" && (el as HTMLInputElement).type === "range";
      if ((tag === "INPUT" && !isRange) || tag === "SELECT" || tag === "TEXTAREA" || tag === "BUTTON") return;
      if (e.key === " ") {
        e.preventDefault();
        togglePlay();
      } else if (!isRange && (e.key === "ArrowLeft" || e.key === "ArrowRight")) {
        e.preventDefault();
        seek(shown + (e.key === "ArrowRight" ? 1 : -1));
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  });

  const fieldIdx = manifest ? Math.max(0, manifest.fields.findIndex((f) => f.name === fieldName)) : 0;
  const field: FieldInfo | undefined = manifest?.fields[fieldIdx];
  const values = current && plane ? fieldSlice(current.data, fieldIdx, plane.n1, plane.n2) : null;
  const liquid = current && plane ? fieldSlice(current.data, 0, plane.n1, plane.n2) : null;

  const kind = scaleKind(field?.name ?? "C");
  let vmin = 0;
  let vmax = 1;
  if (field && field.name !== "C") {
    [vmin, vmax] = rangeMode === "frame" && values ? sliceRange(values) : [field.min ?? 0, field.max ?? 1];
    if (kind === "diverging") {
      const m = Math.max(Math.abs(vmin), Math.abs(vmax));
      [vmin, vmax] = [-m, m];
    }
    if (!(vmax > vmin)) vmax = vmin + 1e-9;
  }
  const lut = colorLut(kind, dark);

  // 3D surface: symmetric elevation range for the colours, and an automatic vertical exaggeration
  const surfaceRange =
    Math.max(Math.abs(manifest?.surface_range?.min ?? 0), Math.abs(manifest?.surface_range?.max ?? 0)) ||
    detail?.params.amplitude ||
    0.01;
  const autoExaggeration = grid ? Math.min(100, Math.max(1, Math.round((0.15 * grid.Lz) / surfaceRange))) : 1;
  const exag = exaggeration ?? autoExaggeration;

  if (loadError && !detail) return <p className="text-critical">{loadError}</p>;
  if (!detail) return <p className="text-muted">Loading run…</p>;

  const p = detail.params;
  const sum = manifest?.summary;
  const pct = (v: number | null | undefined, digits = 2) => (v == null ? "–" : `${(v * 100).toFixed(digits)} %`);
  const showingSurface3d = !!grid && view === "surface";
  const verticalPlane = !sliceView || sliceView !== "xy";

  return (
    <div className="mx-auto flex max-w-6xl flex-col gap-4">
      <header className="flex flex-wrap items-baseline gap-x-4 gap-y-1">
        <h2 className="text-lg font-semibold tracking-tight">
          {p.dim === 3 ? `3D standing wave, mode (${p.mode}, ${p.mode_y})` : `2D standing wave, mode ${p.mode}`}
        </h2>
        <StatusBadge status={detail.status} progress={detail.progress} />
        <span className="text-xs tabular-nums text-muted">{runId}</span>
      </header>
      {active && (
        <div className="h-1.5 overflow-hidden rounded-full bg-gridline">
          <div className="h-full rounded-full bg-accent transition-[width]" style={{ width: `${Math.max(1, detail.progress * 100)}%` }} />
        </div>
      )}
      {detail.message && detail.status !== "done" && (
        <pre className="overflow-x-auto rounded-md border border-line bg-surface p-3 text-xs whitespace-pre-wrap text-critical">
          {detail.message}
        </pre>
      )}

      {!manifest ? (
        <p className="text-muted">
          {active ? "Waiting for the solver to write its first frame…" : "This run produced no frames."}
        </p>
      ) : (
        <>
          {grid && (
            <div className="flex flex-wrap items-end gap-x-4 gap-y-2">
              <div role="tablist" aria-label="View" className="flex overflow-hidden rounded-md border border-line text-sm">
                {(["surface", ...SLICE_VIEWS] as View[]).map((v) => (
                  <button
                    key={v}
                    type="button"
                    role="tab"
                    aria-selected={view === v}
                    onClick={() => setView(v)}
                    className={`px-3 py-1 ${view === v ? "bg-accent text-accent-ink" : "bg-surface text-ink hover:bg-surface-2"}`}
                  >
                    {VIEW_LABEL[v]}
                  </button>
                ))}
              </div>
              {sliceView && (
                <label className="flex min-w-56 flex-1 items-center gap-2 text-xs text-ink-2">
                  <span className="whitespace-nowrap">Slice at</span>
                  <input
                    type="range"
                    min={0}
                    max={sliceShape(grid, sliceView).positions - 1}
                    value={at}
                    onChange={(e) => setSliceAt({ ...sliceAt, [sliceView]: Number(e.target.value) })}
                    aria-label="Slice position"
                    className="min-w-0 flex-1"
                  />
                  <span className="w-24 tabular-nums text-ink">
                    {sliceShape(grid, sliceView).normal} = {slicePosition(grid, sliceView, at).toFixed(3)} m
                  </span>
                </label>
              )}
              {showingSurface3d && (
                <label className="flex min-w-56 flex-1 items-center gap-2 text-xs text-ink-2">
                  <span className="whitespace-nowrap">Vertical exaggeration</span>
                  <input
                    type="range"
                    min={1}
                    max={100}
                    value={exag}
                    onChange={(e) => setExaggeration(Number(e.target.value))}
                    aria-label="Vertical exaggeration"
                    className="min-w-0 flex-1"
                  />
                  <span className="w-10 tabular-nums text-ink">{exag}×</span>
                </label>
              )}
            </div>
          )}

          {showingSurface3d ? (
            <div className="flex flex-wrap items-end justify-between gap-2">
              <p className="text-xs text-muted">
                Drag to rotate, scroll to zoom. Heights are exaggerated {exag}× around the mean depth h = {depth} m.
              </p>
              <ColorLegend kind="diverging" dark={dark} vmin={-surfaceRange} vmax={surfaceRange} label="Surface elevation η − h" unit="m" />
            </div>
          ) : (
            <>
              <div className="flex flex-wrap items-end gap-x-4 gap-y-2">
                <label className="flex flex-col gap-0.5 text-xs text-ink-2">
                  Field
                  <select
                    value={field?.name}
                    onChange={(e) => setFieldName(e.target.value)}
                    className="rounded-md border border-line bg-surface px-2 py-1 text-sm text-ink"
                  >
                    {manifest.fields.map((f) => (
                      <option key={f.name} value={f.name}>
                        {f.label}
                        {f.unit ? ` [${f.unit}]` : ""}
                      </option>
                    ))}
                  </select>
                </label>
                <fieldset className="flex flex-col gap-0.5 text-xs text-ink-2" disabled={field?.name === "C"}>
                  <legend className="mb-0.5">Colour range</legend>
                  <div className="flex overflow-hidden rounded-md border border-line text-sm">
                    {(["run", "frame"] as const).map((m) => (
                      <button
                        key={m}
                        type="button"
                        onClick={() => setRangeMode(m)}
                        aria-pressed={rangeMode === m}
                        className={`px-2 py-1 disabled:opacity-50 ${rangeMode === m ? "bg-accent text-accent-ink" : "bg-surface text-ink hover:bg-surface-2"}`}
                      >
                        {m === "run" ? "Whole run" : "This frame"}
                      </button>
                    ))}
                  </div>
                </fieldset>
                {verticalPlane && (
                  <label className="flex items-center gap-1.5 pb-1 text-xs text-ink-2">
                    <input type="checkbox" checked={showSurface} onChange={(e) => setShowSurface(e.target.checked)} />
                    Surface line
                  </label>
                )}
                <label className="flex items-center gap-1.5 pb-1 text-xs text-ink-2">
                  <input type="checkbox" checked={smooth} onChange={(e) => setSmooth(e.target.checked)} />
                  Smooth
                </label>
                <button
                  type="button"
                  disabled={!values}
                  onClick={() =>
                    values &&
                    field &&
                    plane &&
                    downloadCsv(
                      `${field.name}${sliceView ? `_${sliceView}${at}` : ""}_t${(manifest.frame_times[current!.index] ?? 0).toFixed(3)}`,
                      values,
                      plane,
                    )
                  }
                  className="mb-0.5 rounded-md border border-line px-2 py-1 text-xs text-ink-2 hover:bg-surface-2 disabled:opacity-40"
                >
                  Download CSV
                </button>
                <div className="ml-auto">
                  {field && <ColorLegend kind={kind} dark={dark} vmin={vmin} vmax={vmax} label={field.label} unit={field.unit} />}
                </div>
              </div>
              <p className="-mt-2 text-xs text-muted">{field && FIELD_HELP[field.name]}</p>
            </>
          )}

          <div className="rounded-lg border border-line bg-surface p-3">
            {showingSurface3d && current && grid ? (
              <Surface3D
                eta={current.data}
                nx={grid.nx}
                ny={grid.ny}
                Lx={grid.Lx}
                Ly={grid.Ly}
                Lz={grid.Lz}
                depth={depth}
                range={surfaceRange}
                exaggeration={exag}
                dark={dark}
              />
            ) : values && liquid && field && plane ? (
              <FieldCanvas
                values={values}
                liquid={liquid}
                nx={plane.n1}
                nz={plane.n2}
                Lx={plane.L1}
                Lz={plane.L2}
                axes={plane.axes}
                lut={lut}
                vmin={vmin}
                vmax={vmax}
                smooth={smooth}
                showSurface={showSurface && verticalPlane}
                dark={dark}
                label={field.label}
                unit={field.unit}
              />
            ) : (
              <p className="py-16 text-center text-sm text-muted">{frameError ?? "Loading frame…"}</p>
            )}
          </div>

          <TimeBar
            index={shown}
            times={manifest.frame_times}
            tEnd={manifest.t_end}
            playing={playing}
            speed={speed}
            live={active}
            follow={follow}
            loaded={store?.size ?? 0}
            onIndex={seek}
            onTogglePlay={togglePlay}
            onSpeed={setSpeed}
            onFollow={(f) => {
              setFollow(f);
              if (f) setIndex(last);
            }}
          />

          <div className="grid gap-4 md:grid-cols-2">
            <section className="rounded-lg border border-line bg-surface p-3">
              <h3 className="mb-2 text-sm font-semibold">Setup</h3>
              <dl className="grid grid-cols-[auto_1fr] gap-x-4 gap-y-1 text-sm tabular-nums">
                <dt className="text-muted">Tank</dt>
                <dd>
                  {grid ? `${grid.Lx} × ${grid.Ly} × ${grid.Lz}` : `${p.Lx} × ${p.Lz}`} m,{" "}
                  {p.wall_bc === "slip" ? "free-slip" : "no-slip"} walls
                </dd>
                <dt className="text-muted">Grid</dt>
                <dd>
                  {grid ? `${grid.nx} × ${grid.ny} × ${grid.nz}` : `${manifest.nx} × ${manifest.nz}`} cells (Δx ={" "}
                  {formatValue(manifest.Lx / manifest.nx)} m)
                </dd>
                <dt className="text-muted">Surface</dt>
                <dd>
                  h = {p.depth} m, A = {p.amplitude} m, mode {p.dim === 3 ? `(${p.mode}, ${p.mode_y})` : p.mode}
                </dd>
                <dt className="text-muted">Fluids</dt>
                <dd>water / air (ρ 1000 / 1.2 kg/m³)</dd>
                <dt className="text-muted">Progress</dt>
                <dd>
                  t = {manifest.t.toFixed(3)} / {manifest.t_end} s, {manifest.step} steps, {manifest.wall_time.toFixed(1)} s wall time
                </dd>
              </dl>
            </section>
            <section className="rounded-lg border border-line bg-surface p-3">
              <h3 className="mb-2 text-sm font-semibold">Check against linear theory</h3>
              <dl className="grid grid-cols-[auto_1fr] gap-x-4 gap-y-1 text-sm tabular-nums">
                <dt className="text-muted">Analytic period</dt>
                <dd>{manifest.analytic_period ? `${manifest.analytic_period.toFixed(4)} s` : "–"}</dd>
                <dt className="text-muted">Measured period</dt>
                <dd>
                  {sum?.period_measured != null
                    ? `${sum.period_measured.toFixed(4)} s (${sum.period_error! >= 0 ? "+" : ""}${pct(sum.period_error, 2)}, ${sum.crossings} zero crossings)`
                    : sum
                      ? "– (run too short)"
                      : "available when the run finishes"}
                </dd>
                <dt className="text-muted">Wave energy kept</dt>
                <dd>{sum?.energy_ratio != null ? pct(sum.energy_ratio, 1) + " (last vs first period)" : "–"}</dd>
                <dt className="text-muted">Volume drift</dt>
                <dd>{sum?.volume_drift != null ? sum.volume_drift.toExponential(1) : "–"}</dd>
                <dt className="text-muted">max |∇·u|</dt>
                <dd>{sum?.max_div != null ? `${sum.max_div.toExponential(1)} 1/s` : "–"}</dd>
              </dl>
            </section>
          </div>
        </>
      )}
    </div>
  );
}
