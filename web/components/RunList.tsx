"use client";

import Link from "next/link";
import { usePathname, useRouter } from "next/navigation";
import { useState } from "react";

import { ACTIVE_STATUSES, type RunInfo } from "@/lib/types";
import { useRuns } from "./RunsProvider";
import StatusBadge from "./StatusBadge";

function describe(r: RunInfo): string {
  const p = r.params;
  if (p.dim === 3)
    return `${p.nx}×${p.ny}×${p.nz} · h=${p.depth} · A=${p.amplitude} · (${p.mode},${p.mode_y}) · ${p.t_end}s`;
  return `${p.nx}×${p.nz} · h=${p.depth} · A=${p.amplitude} · n=${p.mode} · ${p.t_end}s`;
}

function formatCreated(iso: string): string {
  const d = new Date(iso);
  return d.toLocaleString(undefined, { month: "short", day: "numeric", hour: "2-digit", minute: "2-digit", second: "2-digit" });
}

export default function RunList() {
  const { runs, error, refresh } = useRuns();
  const pathname = usePathname();
  const router = useRouter();
  const [deleting, setDeleting] = useState<string | null>(null);

  async function remove(r: RunInfo) {
    const active = ACTIVE_STATUSES.includes(r.status);
    const what = active ? "Stop and delete this running simulation?" : "Delete this run and all its frames?";
    if (!window.confirm(`${what}\n\n${describe(r)}`)) return;
    setDeleting(r.id);
    try {
      await fetch(`/api/runs/${r.id}`, { method: "DELETE" });
      await refresh();
      if (pathname === `/runs/${r.id}`) router.push("/");
    } finally {
      setDeleting(null);
    }
  }

  return (
    <section className="flex flex-col gap-2">
      <h2 className="text-sm font-semibold text-ink">Runs</h2>
      {error && <p className="text-xs text-critical">{error}</p>}
      {runs === null && !error && <p className="text-xs text-muted">Loading…</p>}
      {runs?.length === 0 && <p className="text-xs text-muted">No runs yet. Start one above.</p>}
      <ul className="flex flex-col gap-1.5">
        {runs?.map((r) => {
          const selected = pathname === `/runs/${r.id}`;
          const active = ACTIVE_STATUSES.includes(r.status);
          return (
            <li key={r.id} className="group relative">
              <Link
                href={`/runs/${r.id}`}
                aria-current={selected ? "page" : undefined}
                className={`block rounded-md border px-3 py-2 pr-9 transition-colors ${
                  selected ? "border-accent bg-surface-2" : "border-line hover:bg-surface-2"
                }`}
              >
                <div className="flex items-center justify-between gap-2">
                  <span className="flex items-center gap-1.5 text-xs tabular-nums text-ink-2">
                    <span className="rounded border border-line px-1 text-[10px] font-semibold text-ink">{r.params.dim}D</span>
                    {formatCreated(r.createdAt)}
                  </span>
                  <StatusBadge status={r.status} progress={r.progress} />
                </div>
                <div className="mt-0.5 truncate text-xs tabular-nums text-ink">{describe(r)}</div>
                {active && (
                  <div className="mt-1.5 h-1 overflow-hidden rounded-full bg-gridline">
                    <div className="h-full rounded-full bg-accent" style={{ width: `${Math.max(2, r.progress * 100)}%` }} />
                  </div>
                )}
              </Link>
              <button
                type="button"
                onClick={() => remove(r)}
                disabled={deleting === r.id}
                aria-label={active ? "Stop and delete run" : "Delete run"}
                title={active ? "Stop and delete" : "Delete"}
                className="absolute top-1.5 right-1.5 rounded p-1 text-xs text-muted opacity-60 hover:bg-gridline hover:text-critical hover:opacity-100 focus:opacity-100 disabled:opacity-30"
              >
                ✕
              </button>
            </li>
          );
        })}
      </ul>
    </section>
  );
}
