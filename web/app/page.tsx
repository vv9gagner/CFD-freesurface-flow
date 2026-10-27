export default function Home() {
  return (
    <div className="mx-auto max-w-2xl py-10">
      <h2 className="text-2xl font-semibold tracking-tight">Sloshing waves in a closed tank</h2>
      <p className="mt-3 text-ink-2">
        Each run solves the incompressible Navier–Stokes equations for water and air together. The free surface is
        tracked with a geometric volume-of-fluid method. It starts as a cosine and oscillates as a standing wave.
      </p>
      <ul className="mt-4 space-y-1 text-ink-2">
        <li>
          <b className="text-ink">2D</b>: a vertical slice of the tank,{" "}
          <span className="font-mono text-sm">η(x) = h + A·cos(mπx/L)</span>.
        </li>
        <li>
          <b className="text-ink">3D</b>: the whole rectangular tank,{" "}
          <span className="font-mono text-sm">η(x,y) = h + A·cos(mπx/Lx)·cos(nπy/Ly)</span>.
        </li>
      </ul>
      <ol className="mt-6 list-decimal space-y-2 pl-5 text-ink-2">
        <li>
          Choose 2D or 3D in the panel, set the tank, grid and wave parameters, then press <b className="text-ink">Run</b>.
        </li>
        <li>The run opens immediately and new frames appear while the solver is working.</li>
        <li>
          Drag the time bar to scrub through the run. You can also press play, or step with the ← → keys. Switch the
          field between liquid fraction, pressure, dynamic pressure, speed and vorticity.
        </li>
        <li>
          3D runs open on an interactive surface view; drag to rotate it. Slices through the tank (x–z, y–z, x–y) show
          the fields inside, and a slider moves the slice.
        </li>
      </ol>
      <p className="mt-6 text-sm text-muted">Pick a previous run from the list to explore it again.</p>
    </div>
  );
}
