import type { RunStatus } from "@/lib/types";

// Status is always shown as icon + label, never colour alone.
const STYLES: Record<RunStatus, { icon: string; label: string; className: string }> = {
  starting: { icon: "◌", label: "Starting", className: "text-accent" },
  running: { icon: "●", label: "Running", className: "text-accent" },
  done: { icon: "✓", label: "Done", className: "text-good" },
  failed: { icon: "✕", label: "Failed", className: "text-critical" },
  cancelled: { icon: "■", label: "Cancelled", className: "text-muted" },
  interrupted: { icon: "!", label: "Interrupted", className: "text-warning" },
};

export default function StatusBadge({ status, progress }: { status: RunStatus; progress?: number }) {
  const s = STYLES[status];
  return (
    <span className="inline-flex items-center gap-1 text-xs font-medium whitespace-nowrap">
      <span aria-hidden className={`${s.className} ${status === "running" ? "animate-pulse" : ""}`}>
        {s.icon}
      </span>
      <span className="text-ink-2">
        {s.label}
        {status === "running" && progress !== undefined ? ` ${Math.floor(progress * 100)}%` : ""}
      </span>
    </span>
  );
}
