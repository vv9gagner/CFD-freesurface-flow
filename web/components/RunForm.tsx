"use client";

import { useRouter } from "next/navigation";
import { useState } from "react";

import {
  DEFAULT_PARAMS,
  analyticPeriod,
  estimateSeconds,
  fieldsFor,
  frameBytes,
  frameCount,
  validateParams,
  type Dim,
  type NumericKey,
  type ParamErrors,
  type RunParams,
} from "@/lib/params";
import { useRuns } from "./RunsProvider";

type FormState = Record<NumericKey, string> & Pick<RunParams, "wall_bc" | "vtk_output">;

function initialState(dim: Dim): FormState {
  const d = DEFAULT_PARAMS[dim];
  const s = { wall_bc: d.wall_bc, vtk_output: d.vtk_output } as FormState;
  for (const f of fieldsFor(dim)) s[f.key] = String(d[f.key]);
  return s;
}

function formatBytes(n: number): string {
  if (n >= 1e9) return `${(n / 1e9).toFixed(1)} GB`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(0)} MB`;
  return `${(n / 1e3).toFixed(0)} kB`;
}

function formatDuration(sec: number): string {
  if (sec < 90) return `${Math.max(1, Math.round(sec))} s`;
  if (sec < 5400) return `${Math.round(sec / 60)} min`;
  return `${(sec / 3600).toFixed(1)} h`;
}

export default function RunForm() {
  const router = useRouter();
  const { refresh } = useRuns();
  const [dim, setDim] = useState<Dim>(2);
  const [forms, setForms] = useState<Record<Dim, FormState>>(() => ({ 2: initialState(2), 3: initialState(3) }));
  const [errors, setErrors] = useState<ParamErrors>({});
  const [submitting, setSubmitting] = useState(false);
  const [serverError, setServerError] = useState<string | null>(null);

  const form = forms[dim];
  const setForm = (next: FormState) => setForms({ ...forms, [dim]: next });
  const { params } = validateParams({ ...form, dim });
  const period = params ? analyticPeriod(params) : null;

  async function submit(e: React.FormEvent) {
    e.preventDefault();
    const check = validateParams({ ...form, dim });
    setErrors(check.errors);
    setServerError(null);
    if (!check.params) return;
    setSubmitting(true);
    try {
      const res = await fetch("/api/runs", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(check.params),
      });
      const body = await res.json();
      if (!res.ok) {
        setErrors(body.errors ?? {});
        setServerError(body.error ?? `HTTP ${res.status}`);
        return;
      }
      await refresh();
      router.push(`/runs/${body.id}`);
    } catch (err) {
      setServerError((err as Error).message);
    } finally {
      setSubmitting(false);
    }
  }

  const inputClass =
    "w-full rounded-md border border-line bg-page px-2 py-1 text-sm tabular-nums text-ink focus:outline-2 focus:outline-accent";

  return (
    <form onSubmit={submit} className="flex flex-col gap-3" noValidate>
      <div className="flex items-center justify-between gap-2">
        <h2 className="text-sm font-semibold text-ink">New simulation</h2>
        <div role="radiogroup" aria-label="Dimension" className="flex overflow-hidden rounded-md border border-line text-xs">
          {([2, 3] as Dim[]).map((d) => (
            <button
              key={d}
              type="button"
              role="radio"
              aria-checked={dim === d}
              onClick={() => {
                setDim(d);
                setErrors({});
                setServerError(null);
              }}
              className={`px-3 py-1 font-medium ${dim === d ? "bg-accent text-accent-ink" : "bg-page text-ink-2 hover:bg-surface-2"}`}
            >
              {d}D
            </button>
          ))}
        </div>
      </div>
      <div className="grid grid-cols-2 gap-x-3 gap-y-2">
        {fieldsFor(dim).map((f) => {
          const [min, max] = f.limits[dim]!;
          return (
            <label key={f.key} className="flex flex-col gap-0.5 text-xs text-ink-2">
              <span>
                {f.label}
                {f.unit && <span className="text-muted"> [{f.unit}]</span>}
              </span>
              <input
                type="number"
                inputMode="decimal"
                min={min}
                max={max}
                step={f.step}
                value={form[f.key]}
                onChange={(e) => setForm({ ...form, [f.key]: e.target.value })}
                aria-invalid={!!errors[f.key]}
                className={inputClass}
              />
              {errors[f.key] && <span className="text-critical">{errors[f.key]}</span>}
            </label>
          );
        })}
        <label className="flex flex-col gap-0.5 text-xs text-ink-2">
          <span>Walls</span>
          <select
            value={form.wall_bc}
            onChange={(e) => setForm({ ...form, wall_bc: e.target.value as RunParams["wall_bc"] })}
            className={inputClass}
          >
            <option value="slip">free-slip</option>
            <option value="noslip">no-slip</option>
          </select>
        </label>
      </div>
      <label className="flex items-center gap-2 text-xs text-ink-2">
        <input
          type="checkbox"
          checked={form.vtk_output}
          onChange={(e) => setForm({ ...form, vtk_output: e.target.checked })}
        />
        Also write VTK files for ParaView
      </label>

      <dl className="grid grid-cols-[auto_1fr] gap-x-3 gap-y-0.5 rounded-md bg-surface-2 px-3 py-2 text-xs">
        <dt className="text-muted">Linear period</dt>
        <dd className="tabular-nums text-ink">{period ? `${period.toFixed(4)} s` : "–"}</dd>
        <dt className="text-muted">Frames</dt>
        <dd className="tabular-nums text-ink">
          {params ? `${frameCount(params)} (${formatBytes(frameCount(params) * frameBytes(params))})` : "–"}
        </dd>
        <dt className="text-muted">Run time</dt>
        <dd className="tabular-nums text-ink">{params ? `≈ ${formatDuration(estimateSeconds(params))}` : "–"}</dd>
      </dl>

      {(serverError || errors.form) && (
        <p role="alert" className="text-xs text-critical">
          {serverError ?? errors.form}
        </p>
      )}
      <div className="flex gap-2">
        <button
          type="submit"
          disabled={submitting}
          className="flex-1 rounded-md bg-accent px-3 py-1.5 text-sm font-medium text-accent-ink disabled:opacity-60"
        >
          {submitting ? "Starting…" : `Run ${dim}D simulation`}
        </button>
        <button
          type="button"
          onClick={() => {
            setForm(initialState(dim));
            setErrors({});
            setServerError(null);
          }}
          className="rounded-md border border-line px-3 py-1.5 text-sm text-ink-2 hover:bg-surface-2"
        >
          Reset
        </button>
      </div>
    </form>
  );
}
