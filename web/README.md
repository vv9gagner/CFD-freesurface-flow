# Freesurf web viewer

This is a Next.js 16 front end for both solver apps (`../freesurf-2d`, `../freesurf-3d`). With it you can:
- start 2D or 3D simulations,
- watch them fill in while they run,
- scrub through any run with a time bar.

2D runs and 3D slices are drawn as heatmaps with a choice of field. 3D runs also have an interactive three.js view of the free surface.

```bash
make -C ..          # build both solvers first (GNU GCC)
npm install
npm run dev         # http://localhost:3000
npm test            # unit tests for the 3D slice extraction (node --test)
```

The project root defaults to `..`; set `FREESURF_ROOT` to use another one. Runs are written to `<root>/runs/<id>/`.

| Path | Purpose |
|---|---|
| `app/api/runs/` | `GET` lists runs; `POST` starts one (body: `dim` plus the run parameters as JSON) |
| `app/api/runs/[id]/` | `GET` returns the run with its manifest; `DELETE` stops it and removes it |
| `app/api/runs/[id]/frames/[n]/` | `GET` returns a 2D frame, a 3D surface (`?view=surface`) or a 3D slice (`?view=xz\|yz\|xy&at=k`) |
| `lib/params.ts` | Parameter schema per dimension, validation, time and disk estimates (shared by the form and the API) |
| `lib/slices.ts` | Turns stored float16 3D frames into slices in the 2D frame layout (C, p, p_dyn, speed, normal vorticity) |
| `lib/runner.ts` | Spawns the solvers (at most 2 at once, at most 1 of them 3D) and tracks live processes |
| `lib/frames.ts` | Client frame cache (LRU, about 300 MB) and prefetch, one per view |
| `components/Viewer.tsx` | View tabs, field selector, heatmap, legend, time bar, results |
| `components/Surface3D.tsx` | three.js free-surface view (client only) |
