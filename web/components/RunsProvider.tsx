"use client";

import { createContext, useCallback, useContext, useEffect, useState } from "react";

import { ACTIVE_STATUSES, type RunInfo } from "@/lib/types";

type RunsContextValue = {
  runs: RunInfo[] | null;
  error: string | null;
  refresh: () => Promise<void>;
};

const RunsContext = createContext<RunsContextValue | null>(null);

export function useRuns(): RunsContextValue {
  const ctx = useContext(RunsContext);
  if (!ctx) throw new Error("useRuns must be used inside <RunsProvider>");
  return ctx;
}

/** Keeps the run list fresh: polls quickly while any run is active, slowly otherwise. */
export function RunsProvider({ children }: { children: React.ReactNode }) {
  const [runs, setRuns] = useState<RunInfo[] | null>(null);
  const [error, setError] = useState<string | null>(null);

  const refresh = useCallback(async () => {
    try {
      const res = await fetch("/api/runs", { cache: "no-store" });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const body = (await res.json()) as { runs: RunInfo[] };
      setRuns(body.runs);
      setError(null);
    } catch (e) {
      setError(`Could not load runs (${(e as Error).message})`);
    }
  }, []);

  const anyActive = runs?.some((r) => ACTIVE_STATUSES.includes(r.status)) ?? false;

  useEffect(() => {
    let timer: ReturnType<typeof setTimeout>;
    let stopped = false;
    const tick = async () => {
      await refresh();
      if (!stopped) timer = setTimeout(tick, anyActive ? 1500 : 10000);
    };
    timer = setTimeout(tick, 0);
    return () => {
      stopped = true;
      clearTimeout(timer);
    };
  }, [refresh, anyActive]);

  return <RunsContext.Provider value={{ runs, error, refresh }}>{children}</RunsContext.Provider>;
}
