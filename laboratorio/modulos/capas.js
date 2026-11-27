// Solver compartido por los módulos de dos capas (olas.js y mezcla.js)
// ===================== Two-dimensional Boussinesq flow over a sill =====================
// vorticity–streamfunction, x periodic with sponges, rigid lid, barotropic tide through psi_top
class InternalWaves {
  constructor(nx = 300, nz = 75, Lx = 4, H = 1) {
    this.nx = nx; this.nz = nz; this.Lx = Lx; this.H = H;
    this.dx = Lx / nx; this.dz = H / nz;
    const n = nx * (nz + 1); this.n = n;
    this.solid = new Uint8Array(n);
    this.hx = new Float64Array(nx);
    this.kb = new Int32Array(nx);
    this.psibt = new Float64Array(n);
    this.psip = new Float64Array(n);
    this.psi = new Float64Array(n);
    this.om = new Float64Array(n);
    this.b = new Float64Array(n);
    this.b0 = new Float64Array(n);
    this.u = new Float64Array(n);
    this.w = new Float64Array(n);
    this.t1 = new Float64Array(n);
    this.t2 = new Float64Array(n);
    this.db = 1.0; this.U0 = 0.2; this.hs = 0.5; this.T = 24;
    this.zi = 0.66; this.delta = 0.03; this.ws = 0.28;
    this.nu = 3e-4; this.kappa = 3e-4;
    this.nonhyd = true; this.sor = 1.9; this.iters = 18;
    this.t = 0; this.G = 0; this.Q = 0; this.Ls = 0.45; // sponge width at the x ends (0 = fully periodic)
    this.build(true);
  }
  idx(i, k) { return k * this.nx + i; }
  setB0() {
    const { nx, nz, dz } = this;
    for (let k = 0; k <= nz; k++) {
      const v = 0.5 * this.db * Math.tanh((k * dz - this.zi) / this.delta);
      for (let i = 0; i < nx; i++) this.b0[k * nx + i] = v;
    }
  }
  build(resetState) {
    const { nx, nz, dx, dz } = this;
    const x0 = this.Lx / 2;
    for (let i = 0; i < nx; i++) {
      const s = (i * dx - x0) / this.ws;
      this.hx[i] = this.hs * Math.exp(-s * s);
      let kb = 1;
      while (kb * dz <= this.hx[i]) kb++;
      this.kb[i] = kb;
      for (let k = 0; k <= nz; k++) this.solid[k * nx + i] = k < kb ? 1 : 0;
    }
    this.setB0();
    if (resetState) {
      this.b.set(this.b0); this.om.fill(0); this.psip.fill(0); this.t = 0;
    } else {
      for (let j = 0; j < this.n; j++) if (this.solid[j]) { this.b[j] = this.b0[j]; this.om[j] = 0; this.psip[j] = 0; }
    }
    this.buildBarotropic();
    if (!resetState) this.solvePsi(this.nonhyd ? 400 : 0);
  }
  buildBarotropic() {
    const { nx, nz, kb, psibt } = this;
    // hydrostatic (column-linear) profile, also the first guess for the potential flow
    for (let i = 0; i < nx; i++) {
      const k0 = kb[i] - 1;
      for (let k = 0; k <= nz; k++) psibt[k * nx + i] = k <= k0 ? 0 : (k - k0) / (nz - k0);
    }
    if (!this.nonhyd) return;
    const solid = this.solid, w = 1.92;
    for (let it = 0; it < 2500; it++) {
      for (let k = 1; k < nz; k++) {
        const r = k * nx;
        for (let i = 0; i < nx; i++) {
          const j = r + i; if (solid[j]) continue;
          const e = i + 1 < nx ? j + 1 : r, wv = i > 0 ? j - 1 : r + nx - 1;
          const nv = (psibt[e] + psibt[wv] + psibt[j + nx] + psibt[j - nx]) * 0.25;
          psibt[j] += w * (nv - psibt[j]);
        }
      }
    }
  }
  solvePsi(iters) {
    const { nx, nz, dz, psip, om, solid, kb } = this;
    const h2 = dz * dz;
    if (this.nonhyd) {
      const w = this.sor;
      for (let it = 0; it < iters; it++) {
        for (let k = 1; k < nz; k++) {
          const r = k * nx;
          for (let i = 0; i < nx; i++) {
            const j = r + i; if (solid[j]) continue;
            const e = i + 1 < nx ? j + 1 : r, wv = i > 0 ? j - 1 : r + nx - 1;
            const nv = (psip[e] + psip[wv] + psip[j + nx] + psip[j - nx] - h2 * om[j]) * 0.25;
            psip[j] += w * (nv - psip[j]);
          }
        }
      }
    } else {
      // psi_zz = omega, column by column (Thomas algorithm)
      const cp = this.t1, dp = this.t2;
      for (let i = 0; i < nx; i++) {
        const k0 = kb[i];
        // unknowns k0..nz-1 ; a=1,b=-2,c=1
        let prevc = 0, prevd = 0;
        for (let k = k0; k < nz; k++) {
          const j = k * nx + i;
          const den = -2 - (k > k0 ? prevc : 0);
          const c = 1 / den;
          const d = (h2 * om[j] - (k > k0 ? prevd : 0)) / den;
          cp[j] = c; dp[j] = d; prevc = c; prevd = d;
        }
        let next = 0;
        for (let k = nz - 1; k >= k0; k--) {
          const j = k * nx + i;
          const v = dp[j] - cp[j] * next;
          psip[j] = v; next = v;
        }
        for (let k = 0; k < k0; k++) psip[k * nx + i] = 0;
        psip[nz * nx + i] = 0;
      }
    }
  }
  velocities() {
    const { nx, nz, dx, dz, psi, psip, psibt, u, w, solid } = this;
    const Q = this.Q;
    for (let j = 0; j < this.n; j++) psi[j] = Q * psibt[j] + psip[j];
    for (let k = 0; k <= nz; k++) {
      const r = k * nx;
      for (let i = 0; i < nx; i++) {
        const j = r + i;
        if (solid[j]) { u[j] = 0; w[j] = 0; continue; }
        if (k === nz) { u[j] = (psi[j] - psi[j - nx]) / dz; w[j] = 0; continue; }
        u[j] = (psi[j + nx] - psi[j - nx]) / (2 * dz);
        const e = i + 1 < nx ? j + 1 : r, wv = i > 0 ? j - 1 : r + nx - 1;
        w[j] = -(psi[e] - psi[wv]) / (2 * dx);
      }
    }
  }
  interp(a, x, z) {
    const { nx, nz, dx, dz } = this;
    let fx = x / dx; let i0 = Math.floor(fx); const ax = fx - i0;
    i0 = ((i0 % nx) + nx) % nx; const i1 = i0 + 1 < nx ? i0 + 1 : 0;
    let fz = z / dz; if (fz < 0) fz = 0; if (fz > nz - 1e-9) fz = nz - 1e-9;
    const k0 = Math.floor(fz), az = fz - k0;
    const r0 = k0 * nx, r1 = r0 + nx;
    return (1 - az) * ((1 - ax) * a[r0 + i0] + ax * a[r0 + i1]) + az * ((1 - ax) * a[r1 + i0] + ax * a[r1 + i1]);
  }
  // Catmull-Rom interpolation, clipped to the local min/max (keeps the interface sharp without overshoots)
  interpC(a, x, z) {
    const { nx, nz, dx, dz } = this;
    let fx = x / dx; let i1 = Math.floor(fx); const tx = fx - i1;
    let fz = z / dz; if (fz < 0) fz = 0; if (fz > nz - 1e-9) fz = nz - 1e-9;
    const k1 = Math.floor(fz), tz = fz - k1;
    i1 = ((i1 % nx) + nx) % nx;
    const i0 = i1 > 0 ? i1 - 1 : nx - 1, i2 = i1 + 1 < nx ? i1 + 1 : 0, i3 = i2 + 1 < nx ? i2 + 1 : 0;
    const ks = [k1 > 0 ? k1 - 1 : 0, k1, k1 + 1, k1 + 2 <= nz ? k1 + 2 : nz];
    const cr = (p0, p1, p2, p3, t) => p1 + 0.5 * t * (p2 - p0 + t * (2 * p0 - 5 * p1 + 4 * p2 - p3 + t * (3 * (p1 - p2) + p3 - p0)));
    const row = new Array(4);
    for (let m = 0; m < 4; m++) {
      const r = ks[m] * nx;
      row[m] = cr(a[r + i0], a[r + i1], a[r + i2], a[r + i3], tx);
    }
    let v = cr(row[0], row[1], row[2], row[3], tz);
    const r1 = k1 * nx, r2 = (k1 + 1) * nx;
    const p = a[r1 + i1], q = a[r1 + i2], s2 = a[r2 + i1], t2 = a[r2 + i2];
    const lo = Math.min(p, q, s2, t2), hi = Math.max(p, q, s2, t2);
    return v < lo ? lo : v > hi ? hi : v;
  }
  // convective adjustment: what hydrostatic ocean models do when heavy water ends up above light water
  convect() {
    const { nx, nz, b, kb } = this;
    for (let pass = 0; pass < 3; pass++) {
      for (let i = 0; i < nx; i++) {
        for (let k = kb[i]; k < nz; k++) {
          const j = k * nx + i, ju = j + nx;
          if (b[ju] < b[j]) { const m = 0.5 * (b[ju] + b[j]); b[ju] = m; b[j] = m; }
        }
      }
    }
  }
  step(dt) {
    const { nx, nz, dx, dz, u, w, om, b, b0, solid } = this;
    this.Q = this.U0 * this.H * Math.sin(2 * Math.PI * this.t / this.T);
    this.velocities();
    const omN = this.t1, bN = this.t2;
    // semi-Lagrangian advection (RK2 backtrace)
    for (let k = 0; k <= nz; k++) {
      const z = k * dz, r = k * nx;
      for (let i = 0; i < nx; i++) {
        const j = r + i;
        if (solid[j]) { omN[j] = 0; bN[j] = b0[j]; continue; }
        const x = i * dx;
        const xm = x - 0.5 * dt * u[j], zm = z - 0.5 * dt * w[j];
        const um = this.interp(u, xm, zm), wm = this.interp(w, xm, zm);
        const xd = x - dt * um, zd = z - dt * wm;
        omN[j] = k === nz ? 0 : this.interp(om, xd, zd);
        bN[j] = this.interpC(b, xd, zd);
      }
    }
    // baroclinic torque + diffusion + sponge
    const cx = dt / (2 * dx), nuc = this.nu * dt / (dx * dx), kac = this.kappa * dt / (dx * dx);
    const Ls = this.Ls, rmax = 1.2, Lx = this.Lx;
    for (let k = 1; k <= nz; k++) {
      const r = k * nx;
      for (let i = 0; i < nx; i++) {
        const j = r + i;
        if (solid[j]) { om[j] = 0; b[j] = b0[j]; continue; }
        const e = i + 1 < nx ? j + 1 : r, wv = i > 0 ? j - 1 : r + nx - 1;
        if (k === nz) {
          om[j] = 0;
          b[j] = bN[j] + kac * (bN[e] + bN[wv] + 2 * bN[j - nx] - 4 * bN[j]);
        } else {
          let o = omN[j] - cx * (bN[e] - bN[wv]);
          o += nuc * (omN[e] + omN[wv] + omN[j + nx] + omN[j - nx] - 4 * omN[j]);
          om[j] = o;
          b[j] = bN[j] + kac * (bN[e] + bN[wv] + bN[j + nx] + bN[j - nx] - 4 * bN[j]);
        }
        const x = i * dx; const d = Math.min(x, Lx - x);
        if (d < Ls) {
          const s = 1 - d / Ls; const rr = rmax * s * s * dt;
          om[j] /= 1 + rr;
          b[j] = (b[j] + rr * b0[j]) / (1 + rr);
        }
      }
    }
    if (!this.nonhyd) this.convect();
    this.t += dt;
    this.Q = this.U0 * this.H * Math.sin(2 * Math.PI * this.t / this.T);
    this.solvePsi(this.iters);
  }
  // composite internal Froude number over the sill crest
  froude() {
    const { nx, nz, dz, b, u } = this;
    this.velocities();
    const ic = Math.round(nx / 2), k0 = this.kb[ic];
    let zint = null;
    for (let k = k0; k < nz; k++) {
      const a = b[k * nx + ic], c = b[(k + 1) * nx + ic];
      if (a <= 0 && c > 0) { zint = (k + a / (a - c)) * dz; break; }
    }
    const hb = this.hx[ic];
    if (zint === null) { const bb = b[k0 * nx + ic]; zint = bb > 0 ? hb : this.H; }
    const h2 = Math.max(zint - hb, 0.01), h1 = Math.max(this.H - zint, 0.01);
    let s1 = 0, n1 = 0, s2 = 0, n2 = 0;
    for (let k = k0; k <= nz; k++) {
      const z = k * dz, v = u[k * nx + ic];
      if (z > zint) { s1 += v; n1++; } else { s2 += v; n2++; }
    }
    const u1 = n1 ? s1 / n1 : 0, u2 = n2 ? s2 / n2 : 0;
    const G2 = (u1 * u1) / (this.db * h1) + (u2 * u2) / (this.db * h2);
    return { G: Math.sqrt(G2), zint, h1, h2 };
  }
  healthy() {
    for (let j = 0; j < this.n; j += 53) { const v = this.om[j]; if (!(Math.abs(v) < 1e4)) return false; }
    return true;
  }
}

// colour table for density, shared by both modules
const LUT_RHO = new Uint32Array(256);
(() => {
  const dense = hexToRgb('#182c66'), midc = hexToRgb('#2a7d9b'), fresh = hexToRgb('#9be8cf');
  for (let i = 0; i < 256; i++) { const t = i / 255; LUT_RHO[i] = pack(t < 0.5 ? mix(dense, midc, t / 0.5) : mix(midc, fresh, (t - 0.5) / 0.5)); }
})();
