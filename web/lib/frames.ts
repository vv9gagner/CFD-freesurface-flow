// Client-side cache of frames. A frame is float32 data: a 2D frame or 3D slice holds
// 5 fields of n1*n2 values (field-major, horizontal axis fastest, bottom row first);
// a 3D surface frame holds the height map eta(x, y).
// Subscribers are notified whenever a frame arrives, so React can read the
// cache through useSyncExternalStore.

export class FrameStore {
  private cache = new Map<number, Float32Array>(); // Map order doubles as LRU order
  private inflight = new Map<number, Promise<Float32Array>>();
  private listeners = new Set<() => void>();
  private version = 0;

  /** `query` is appended to the frame URL, e.g. "?view=xz&at=16" for a 3D slice. */
  constructor(
    readonly runId: string,
    readonly frameBytes: number,
    readonly query = "",
    readonly budgetBytes = 300 * 1024 * 1024,
  ) {}

  /** How many frames fit in the memory budget. */
  get capacity(): number {
    return Math.max(8, Math.floor(this.budgetBytes / this.frameBytes));
  }

  get size(): number {
    return this.cache.size;
  }

  subscribe = (listener: () => void) => {
    this.listeners.add(listener);
    return () => this.listeners.delete(listener);
  };

  getVersion = () => this.version;

  has(index: number): boolean {
    return this.cache.has(index);
  }

  /** Cached frame `index`, or the closest cached frame before it (to avoid flicker while loading). */
  peekNearest(index: number): { index: number; data: Float32Array } | null {
    for (let i = index; i >= 0 && i > index - 64; i--) {
      const data = this.cache.get(i);
      if (data) return { index: i, data };
    }
    return null;
  }

  /** Fetch a frame (deduplicated; requests are shared, so they are never aborted). */
  load(index: number): Promise<Float32Array> {
    const hit = this.cache.get(index);
    if (hit) {
      this.cache.delete(index);
      this.cache.set(index, hit);
      return Promise.resolve(hit);
    }
    const pending = this.inflight.get(index);
    if (pending) return pending;

    const p = fetch(`/api/runs/${this.runId}/frames/${index}${this.query}`)
      .then(async (res) => {
        if (!res.ok) throw new Error(`frame ${index}: HTTP ${res.status}`);
        const buf = await res.arrayBuffer();
        if (buf.byteLength !== this.frameBytes) {
          throw new Error(`frame ${index}: ${buf.byteLength} bytes, expected ${this.frameBytes}`);
        }
        const data = new Float32Array(buf);
        this.cache.set(index, data);
        while (this.cache.size > this.capacity) {
          const oldest = this.cache.keys().next().value as number;
          this.cache.delete(oldest);
        }
        this.version++;
        this.listeners.forEach((l) => l());
        return data;
      })
      .finally(() => this.inflight.delete(index));
    this.inflight.set(index, p);
    return p;
  }
}

/** One field of a frame as a view (no copy). */
export function fieldSlice(frame: Float32Array, field: number, nx: number, nz: number): Float32Array {
  const n = nx * nz;
  return frame.subarray(field * n, (field + 1) * n);
}

export function sliceRange(data: Float32Array): [number, number] {
  let lo = Infinity;
  let hi = -Infinity;
  for (let k = 0; k < data.length; k++) {
    const v = data[k];
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return [lo, hi];
}
