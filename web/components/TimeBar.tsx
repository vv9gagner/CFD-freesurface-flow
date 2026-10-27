"use client";

type Props = {
  index: number;
  times: number[];
  tEnd: number;
  playing: boolean;
  speed: number;
  live: boolean;
  follow: boolean;
  loaded: number;
  onIndex: (i: number) => void;
  onTogglePlay: () => void;
  onSpeed: (s: number) => void;
  onFollow: (f: boolean) => void;
};

const SPEEDS = [0.25, 0.5, 1, 2, 4];

export default function TimeBar(p: Props) {
  const last = Math.max(0, p.times.length - 1);
  const t = p.times[p.index] ?? 0;
  const btn = "rounded-md border border-line px-2 py-1 text-sm text-ink hover:bg-surface-2 disabled:opacity-40";

  return (
    <div className="flex flex-col gap-2 rounded-lg border border-line bg-surface p-3">
      <div className="flex items-center gap-3">
        <span className="w-28 shrink-0 text-sm font-medium tabular-nums text-ink">t = {t.toFixed(3)} s</span>
        <input
          type="range"
          min={0}
          max={last}
          step={1}
          value={p.index}
          onChange={(e) => p.onIndex(Number(e.target.value))}
          aria-label="Simulation time"
          aria-valuetext={`t = ${t.toFixed(3)} seconds`}
          className="h-2 min-w-0 flex-1 cursor-pointer"
          disabled={p.times.length < 2}
        />
        <span className="w-20 shrink-0 text-right text-xs tabular-nums text-muted">/ {p.tEnd} s</span>
      </div>
      <div className="flex flex-wrap items-center gap-2">
        <button type="button" className={btn} onClick={() => p.onIndex(Math.max(0, p.index - 1))} aria-label="Previous frame" disabled={p.index <= 0}>
          ◀︎
        </button>
        <button type="button" className={`${btn} w-20`} onClick={p.onTogglePlay} disabled={p.times.length < 2}>
          {p.playing ? "❚❚ Pause" : "▶ Play"}
        </button>
        <button type="button" className={btn} onClick={() => p.onIndex(Math.min(last, p.index + 1))} aria-label="Next frame" disabled={p.index >= last}>
          ▶︎
        </button>
        <label className="ml-1 flex items-center gap-1 text-xs text-ink-2">
          Speed
          <select
            value={p.speed}
            onChange={(e) => p.onSpeed(Number(e.target.value))}
            className="rounded-md border border-line bg-page px-1 py-0.5 text-xs text-ink"
          >
            {SPEEDS.map((s) => (
              <option key={s} value={s}>
                {s}×
              </option>
            ))}
          </select>
        </label>
        {p.live && (
          <label className="flex items-center gap-1 text-xs text-ink-2">
            <input type="checkbox" checked={p.follow} onChange={(e) => p.onFollow(e.target.checked)} />
            Follow live
          </label>
        )}
        <span className="ml-auto text-xs tabular-nums text-muted">
          frame {p.index + 1}/{p.times.length} · cached {p.loaded}
        </span>
      </div>
    </div>
  );
}
