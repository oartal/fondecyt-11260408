// ===================== Lattice Boltzmann D2Q9 (BGK + Smagorinsky) =====================
const LB_CX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const LB_CY = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const LB_W = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
const LB_OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];

class LBM {
  constructor(nx, ny) {
    this.nx = nx; this.ny = ny; this.N = nx * ny;
    this.f = new Float32Array(9 * this.N);
    this.g = new Float32Array(9 * this.N);
    this.solid = new Uint8Array(this.N);
    this.ux = new Float32Array(this.N);
    this.uy = new Float32Array(this.N);
    this.Uset = 0.1; this.U = 0.02; this.nu = 0.02; this.Cs = 0.12;
    this.sponge = 36; this.kick = 0.15;
    this.reset();
  }
  reset() {
    const N = this.N, f = this.f, U = this.U;
    for (let k = 0; k < 9; k++) {
      const cu = LB_CX[k] * U;
      const v = LB_W[k] * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * U * U);
      f.fill(v, k * N, (k + 1) * N);
    }
    // tiny random seed so symmetric shapes can still shed vortices
    // plus a gentle sideways sway, so the wake does not take ages to start wobbling
    const ph = Math.random() * 6.283, nx = this.nx;
    for (let i = 0; i < N; i++) {
      const e = (Math.random() - 0.5) * 2e-4 + this.kick * U * Math.sin(ph + (i % nx) * 0.07);
      f[2 * N + i] += e; f[4 * N + i] -= e;
    }
    this.ux.fill(U); this.uy.fill(0);
    for (let i = 0; i < N; i++) if (this.solid[i]) { this.ux[i] = 0; }
    this.steps = 0;
  }
  step() {
    const nx = this.nx, ny = this.ny, N = this.N, f = this.f, g = this.g, s = this.solid;
    const uxA = this.ux, uyA = this.uy;
    this.U += (this.Uset - this.U) * 0.004;
    const U = this.U;
    const tau0 = 3 * this.nu + 0.5;
    const smag = 18 * Math.SQRT2 * this.Cs * this.Cs;
    const spw = this.sponge, sp0 = nx - spw;
    const N2 = 2 * N, N3 = 3 * N, N4 = 4 * N, N5 = 5 * N, N6 = 6 * N, N7 = 7 * N, N8 = 8 * N;
    const usqIn = 1.5 * U * U;
    const q0 = (4 / 9) * (1 - usqIn), q1 = (1 / 9) * (1 + 3 * U + 4.5 * U * U - usqIn), q2 = (1 / 9) * (1 - usqIn),
      q3 = (1 / 9) * (1 - 3 * U + 4.5 * U * U - usqIn), q5 = (1 / 36) * (1 + 3 * U + 4.5 * U * U - usqIn),
      q6 = (1 / 36) * (1 - 3 * U + 4.5 * U * U - usqIn);
    for (let y = 0; y < ny; y++) {
      const ru = (y + 1 < ny ? y + 1 : 0) * nx;
      const rd = (y > 0 ? y - 1 : ny - 1) * nx;
      const r = y * nx;
      // inlet column: equilibrium at U
      const i0 = r;
      if (!s[i0]) {
        g[i0] = q0; g[N + i0] = q1; g[N2 + i0] = q2; g[N3 + i0] = q3; g[N4 + i0] = q2;
        g[N5 + i0] = q5; g[N6 + i0] = q6; g[N7 + i0] = q6; g[N8 + i0] = q5;
        uxA[i0] = U; uyA[i0] = 0;
      }
      for (let x = 1; x < nx - 1; x++) {
        const idx = r + x;
        if (s[idx]) { uxA[idx] = 0; uyA[idx] = 0; continue; }
        let j;
        const f0 = f[idx];
        j = idx - 1;       const f1 = s[j] ? f[N3 + idx] : f[N + j];
        j = rd + x;        const f2 = s[j] ? f[N4 + idx] : f[N2 + j];
        j = idx + 1;       const f3 = s[j] ? f[N + idx] : f[N3 + j];
        j = ru + x;        const f4 = s[j] ? f[N2 + idx] : f[N4 + j];
        j = rd + x - 1;    const f5 = s[j] ? f[N7 + idx] : f[N5 + j];
        j = rd + x + 1;    const f6 = s[j] ? f[N8 + idx] : f[N6 + j];
        j = ru + x + 1;    const f7 = s[j] ? f[N5 + idx] : f[N7 + j];
        j = ru + x - 1;    const f8 = s[j] ? f[N6 + idx] : f[N8 + j];
        const rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
        const ir = 1 / rho;
        const ux = (f1 - f3 + f5 - f6 - f7 + f8) * ir;
        const uy = (f2 - f4 + f5 + f6 - f7 - f8) * ir;
        const usq = 1.5 * (ux * ux + uy * uy);
        const r9 = rho / 9, r36 = rho / 36;
        const fe0 = 4 * r9 * (1 - usq);
        const fe1 = r9 * (1 + 3 * ux + 4.5 * ux * ux - usq);
        const fe2 = r9 * (1 + 3 * uy + 4.5 * uy * uy - usq);
        const fe3 = r9 * (1 - 3 * ux + 4.5 * ux * ux - usq);
        const fe4 = r9 * (1 - 3 * uy + 4.5 * uy * uy - usq);
        const a5 = ux + uy, a6 = uy - ux;
        const fe5 = r36 * (1 + 3 * a5 + 4.5 * a5 * a5 - usq);
        const fe6 = r36 * (1 + 3 * a6 + 4.5 * a6 * a6 - usq);
        const fe7 = r36 * (1 - 3 * a5 + 4.5 * a5 * a5 - usq);
        const fe8 = r36 * (1 - 3 * a6 + 4.5 * a6 * a6 - usq);
        const n1 = f1 - fe1, n2 = f2 - fe2, n3 = f3 - fe3, n4 = f4 - fe4;
        const n5 = f5 - fe5, n6 = f6 - fe6, n7 = f7 - fe7, n8 = f8 - fe8;
        const pxx = n1 + n3 + n5 + n6 + n7 + n8;
        const pyy = n2 + n4 + n5 + n6 + n7 + n8;
        const pxy = n5 - n6 + n7 - n8;
        const Q = Math.sqrt(pxx * pxx + pyy * pyy + 2 * pxy * pxy);
        let tau = 0.5 * (tau0 + Math.sqrt(tau0 * tau0 + smag * Q * ir));
        if (x > sp0) { const q = (x - sp0) / spw; tau += 1.5 * q * q; }
        const om = 1 / tau;
        g[idx] = f0 + om * (fe0 - f0);
        g[N + idx] = f1 - om * n1;
        g[N2 + idx] = f2 - om * n2;
        g[N3 + idx] = f3 - om * n3;
        g[N4 + idx] = f4 - om * n4;
        g[N5 + idx] = f5 - om * n5;
        g[N6 + idx] = f6 - om * n6;
        g[N7 + idx] = f7 - om * n7;
        g[N8 + idx] = f8 - om * n8;
        uxA[idx] = ux; uyA[idx] = uy;
      }
      // outlet: zero gradient
      const io = r + nx - 1, ip = io - 1;
      for (let k = 0; k < 9; k++) g[k * N + io] = g[k * N + ip];
      uxA[io] = uxA[ip]; uyA[io] = uyA[ip];
    }
    this.f = g; this.g = f;
    this.steps++;
  }
  healthy() {
    const ux = this.ux, n = this.N;
    for (let i = 0; i < n; i += 97) { const v = ux[i]; if (!(v > -1 && v < 1)) return false; }
    return true;
  }
  // characteristic length: vertical extent (frontal height) of solid
  frontalHeight() {
    const nx = this.nx, ny = this.ny, s = this.solid;
    let cnt = 0;
    for (let y = 0; y < ny; y++) { for (let x = 0; x < nx; x++) if (s[y * nx + x]) { cnt++; break; } }
    return cnt;
  }
}

Lab.register('re', {
  html: `
    <section class="module" id="mod-re" role="tabpanel" aria-labelledby="tab-re">
      <div class="stage" id="stage-re">
        <canvas id="cv-re" aria-label="Simulación del flujo alrededor de un obstáculo"></canvas>
        <div class="tag tl"><span class="arrow">→ → →</span> el fluido entra por la izquierda</div>
        <div class="tag tr">Simulación en 2D</div>
        <div class="tag bl">Clic izquierdo: dibuja una roca · Clic derecho: bórrala</div>
        <div class="demo-chip" id="demo-re" hidden>Modo demostración · mueve el mouse para jugar</div>
        <div class="toast" id="toast-re" hidden></div>
      </div>
      <div class="panel">
        <div class="card controls">
          <div class="ctl wide">
            <div class="label">Fluido <span class="fine">↑ ↓</span></div>
            <div class="seg" id="re-fluids" style="--n:3"></div>
          </div>
          <div class="ctl wide">
            <label class="label" for="re-speed">Velocidad <span class="fine">← →</span></label>
            <input type="range" id="re-speed" min="0" max="100" step="1" value="60">
            <div class="ends"><span>lento</span><span>rápido</span></div>
          </div>
          <div class="ctl">
            <label class="label" for="re-size">Tamaño</label>
            <input type="range" id="re-size" min="5" max="14" step="1" value="9">
            <div class="ends"><span>chico</span><span>grande</span></div>
          </div>
          <div class="ctl">
            <div class="label">Vista <span class="fine">V</span></div>
            <div class="seg" id="re-views" style="--n:2"></div>
          </div>
          <div class="ctl wide">
            <div class="label">Objeto</div>
            <div class="seg" id="re-shapes" style="--n:4"></div>
          </div>
        </div>
        <div class="card gauge">
          <div class="label">¿Cómo pasa el fluido?</div>
          <div class="lights" id="re-lights" role="status" aria-live="polite">
            <div class="light" style="--c: var(--calm)"><i></i>Ordenado</div>
            <div class="light" style="--c: var(--warn)"><i></i>Remolinos</div>
            <div class="light" style="--c: var(--hot)"><i></i>Revuelto</div>
          </div>
          <p class="sub" id="re-regime-sub">—</p>
        </div>
        <div class="card explain">
          <h2>¿Qué pasa detrás de la roca?</h2>
          <p>Si el fluido va lento, o es espeso como la miel, pasa ordenado.</p>
          <p>Si va rápido, deja remolinos que se sueltan a un lado y al otro. En el mar pasa lo mismo detrás de las islas y las rocas.</p>
          <p class="fine">Los científicos miden esto con el número de Reynolds.</p>
        </div>
      </div>
    </section>`,
  init() {
    'use strict';
    // ==================================================================
    //                          MÓDULO 1: REYNOLDS
    // ==================================================================
    const RE = {
      nx: 300, ny: 100, cm: 0.5, // each cell = 0.5 cm
      fluids: [
        // K = Reynolds number on screen at full speed for the default rock (the slider spans K/12 … K)
        { id: 'miel', name: 'Miel', hint: 'muy espesa', K: 14 },
        { id: 'aceite', name: 'Aceite', hint: 'espeso', K: 110 },
        { id: 'agua', name: 'Agua', hint: 'nada espesa', K: 800 },
      ],
      shapes: [
        { id: 'circulo', name: 'Roca' },
        { id: 'placa', name: 'Placa' },
        { id: 'ala', name: 'Ala' },
        { id: 'islas', name: 'Islas' },
      ],
      views: [{ id: 'remolinos', name: 'Remolinos' }, { id: 'tinta', name: 'Tinta' }],
      fluid: 2, shape: 0, view: 0, size: 9, speedT: 60,
      ReSim: 150, ReMax: 900, D: 18, msStep: 1.5, state: -2,
    };
    const lbm = new LBM(RE.nx, RE.ny);
    const cvRe = $('#cv-re'), ctxRe = cvRe.getContext('2d');
    const offRe = document.createElement('canvas'); offRe.width = RE.nx; offRe.height = RE.ny;
    const offReCtx = offRe.getContext('2d');
    const imgRe = offReCtx.createImageData(RE.nx, RE.ny);
    const pixRe = new Uint32Array(imgRe.data.buffer);

    // colour tables
    const LUT_VORT = new Uint32Array(512), LUT_SPEED = new Uint32Array(256);
    let ROCK_RE = 0;
    function buildReLuts() {
      const bg = hexToRgb('#0f2331'), neg = hexToRgb('#46c2ff'), pos = hexToRgb('#ff6b5a');
      const negHi = hexToRgb('#d8f3ff'), posHi = hexToRgb('#ffe1c2');
      for (let i = 0; i < 512; i++) {
        const v = (i - 255.5) / 255.5; const a = Math.min(1, Math.abs(v));
        const t = Math.pow(a, 0.7);
        let c = mix(bg, v < 0 ? neg : pos, Math.min(1, t * 1.4));
        if (t > 0.72) c = mix(c, v < 0 ? negHi : posHi, (t - 0.72) / 0.28 * 0.6);
        LUT_VORT[i] = pack(c);
      }
      const s0 = hexToRgb('#0b1c29'), s1 = hexToRgb('#1d4a66'), s2 = hexToRgb('#3d86a6');
      for (let i = 0; i < 256; i++) { const t = i / 255; LUT_SPEED[i] = pack(t < 0.6 ? mix(s0, s1, t / 0.6) : mix(s1, s2, (t - 0.6) / 0.4)); }
      ROCK_RE = pack(hexToRgb('#c8b28a'));
    }
    buildReLuts();

    function buildShape() {
      const { nx, ny } = RE; const s = lbm.solid; s.fill(0);
      const r = RE.size, cx = Math.round(nx * 0.22), cy = ny / 2 + 0.3;
      const disk = (x0, y0, rr) => {
        for (let y = Math.floor(y0 - rr - 1); y <= y0 + rr + 1; y++) for (let x = Math.floor(x0 - rr - 1); x <= x0 + rr + 1; x++) {
          if (x < 3 || x > nx - 4 || y < 0 || y >= ny) continue;
          const dx = x - x0, dy = y - y0; if (dx * dx + dy * dy <= rr * rr) s[y * nx + x] = 1;
        }
      };
      const id = RE.shapes[RE.shape].id;
      if (id === 'circulo') disk(cx, cy, r);
      else if (id === 'placa') {
        for (let y = Math.round(cy - r); y <= Math.round(cy + r); y++) for (let x = cx - 1; x <= cx + 1; x++) s[y * nx + x] = 1;
      } else if (id === 'ala') {
        const c = 4.2 * r, ang = -16 * Math.PI / 180, ca = Math.cos(ang), sa = Math.sin(ang);
        for (let y = 0; y < ny; y++) for (let x = 3; x < nx - 3; x++) {
          const X = (x - cx) * ca + (y - cy) * sa, Y = -(x - cx) * sa + (y - cy) * ca;
          const xi = (X + c / 2) / c; if (xi < 0 || xi > 1) continue;
          const t = 5 * 0.16 * c * (0.2969 * Math.sqrt(xi) - 0.126 * xi - 0.3516 * xi * xi + 0.2843 * xi ** 3 - 0.1015 * xi ** 4);
          if (Math.abs(Y) <= t) s[y * nx + x] = 1;
        }
      } else if (id === 'islas') {
        disk(cx - 0.6 * r, cy - 1.5 * r, 0.6 * r);
        disk(cx + 1.4 * r, cy + 0.3 * r, 0.75 * r);
        disk(cx - 0.4 * r, cy + 1.8 * r, 0.5 * r);
      }
      lbm.reset();
      particlesRe.n = 0; resetProbeRe();
      updateReParams();
    }

    // Every position of the controls maps to a Reynolds number the screen really computes (1 … 900).
    function updateReParams() {
      const f = RE.fluids[RE.fluid];
      const Dc = lbm.frontalHeight(); RE.D = Dc;
      if (Dc < 1) { RE.ReSim = NaN; lbm.Uset = 0.08; lbm.nu = 0.02; return; }
      const t = RE.speedT / 100;
      const Re = clamp(f.K * Math.pow(12, t - 1) * Dc / 18, 1, RE.ReMax);
      RE.ReSim = Re;
      let U = 0.03 + 0.09 * t, nu = U * Dc / Re; // faster slider = visibly faster flow
      if (nu > 0.5) { nu = 0.5; U = Re * nu / Dc; }
      lbm.Uset = U; lbm.nu = nu;
    }

    // how much the wake wobbles in time, measured on the simulation itself
    const probeRe = (() => {
      const idx = [];
      for (let x = Math.round(RE.nx * 0.36); x < RE.nx * 0.88; x += 12) for (let y = 3; y < RE.ny; y += 6) idx.push(y * RE.nx + x);
      return { idx, m: new Float32Array(idx.length), v: new Float32Array(idx.length), level: 0 };
    })();
    function resetProbeRe() { probeRe.m.fill(0); probeRe.v.fill(0); probeRe.level = 0; }
    function sampleProbeRe(nsteps) {
      const P = probeRe, a = Math.min(0.2, nsteps / 1000), uy = lbm.uy; let s = 0;
      for (let p = 0; p < P.idx.length; p++) {
        const d = uy[P.idx[p]] - P.m[p];
        P.m[p] += a * d; P.v[p] += a * (d * d - P.v[p]); s += P.v[p];
      }
      P.level = Math.sqrt(s / P.idx.length) / Math.max(lbm.U, 1e-4);
    }
    const RE_SUBS = {
      '-1': 'Dibuja una roca con el mouse sobre el agua.',
      0: 'El fluido rodea el objeto suave y parejo, sin soltar remolinos.',
      1: 'Detrás del objeto se sueltan remolinos, uno a cada lado.',
      2: 'Los remolinos salen grandes y fuertes, y revuelven toda el agua detrás del objeto.',
    };
    function updateRegimeRe() {
      let st;
      if (!(RE.D >= 1)) st = -1;
      else if (probeRe.level < (RE.state >= 1 ? 0.06 : 0.08)) st = 0;
      else st = RE.ReSim < 350 ? 1 : 2;
      if (st === RE.state) return;
      RE.state = st;
      [...$('#re-lights').children].forEach((el, i) => el.classList.toggle('on', i === st));
      $('#re-regime-sub').textContent = RE_SUBS[st];
    }

    // particles (dye)
    const particlesRe = { x: new Float32Array(7000), y: new Float32Array(7000), c: new Uint8Array(7000), n: 0 };
    function sampleRe(a, x, y) {
      const { nx, ny } = RE;
      const x0 = clamp(Math.floor(x), 0, nx - 2), ax = clamp(x - x0, 0, 1);
      let y0 = Math.floor(y); const ay = y - y0; y0 = ((y0 % ny) + ny) % ny; const y1 = y0 + 1 < ny ? y0 + 1 : 0;
      return (1 - ay) * ((1 - ax) * a[y0 * nx + x0] + ax * a[y0 * nx + x0 + 1]) + ay * ((1 - ax) * a[y1 * nx + x0] + ax * a[y1 * nx + x0 + 1]);
    }
    function advectParticlesRe(nsteps) {
      const P = particlesRe, { nx, ny } = RE;
      // emit
      for (let row = 0; row < 17 && P.n < P.x.length; row++) {
        P.x[P.n] = 1.2; P.y[P.n] = 3 + row * 5.75 + (Math.random() - 0.5) * 0.4; P.c[P.n] = row & 1; P.n++;
      }
      let w = 0;
      for (let i = 0; i < P.n; i++) {
        let x = P.x[i], y = P.y[i];
        x += sampleRe(lbm.ux, x, y) * nsteps; y += sampleRe(lbm.uy, x, y) * nsteps;
        if (y < 0) y += ny; if (y >= ny) y -= ny;
        if (x >= nx - 2) continue;
        const xi = Math.round(x), yi = Math.round(y) % ny;
        if (lbm.solid[yi * nx + clamp(xi, 0, nx - 1)]) continue;
        P.x[w] = x; P.y[w] = y; P.c[w] = P.c[i]; w++;
      }
      P.n = w;
    }

    function drawRe() {
      const dpr = fitCanvas(cvRe);
      const W = cvRe.width, H = cvRe.height, { nx, ny } = RE;
      const R = letterbox(W, H, nx / ny); RE.rect = R;
      const ux = lbm.ux, uy = lbm.uy, s = lbm.solid;
      const tinta = RE.views[RE.view].id === 'tinta';
      const U = Math.max(lbm.U, 1e-4);
      const vs = RE.D > 0 ? 255.5 / (3.2 * U / Math.max(RE.D, 4)) : 2000;
      const ss = 255 / (1.8 * U);
      for (let y = 0; y < ny; y++) {
        const yu = (y + 1 < ny ? y + 1 : 0) * nx, yd = (y > 0 ? y - 1 : ny - 1) * nx, r = y * nx;
        for (let x = 0; x < nx; x++) {
          const i = r + x;
          if (s[i]) { pixRe[i] = ROCK_RE; continue; }
          if (tinta) {
            const sp = Math.sqrt(ux[i] * ux[i] + uy[i] * uy[i]);
            pixRe[i] = LUT_SPEED[Math.min(255, (sp * ss) | 0)];
          } else {
            const xe = x + 1 < nx ? i + 1 : i, xw = x > 0 ? i - 1 : i;
            const om = 0.5 * ((uy[xe] - uy[xw]) - (ux[yu + x] - ux[yd + x]));
            let k = (om * vs + 255.5) | 0; k = k < 0 ? 0 : k > 511 ? 511 : k;
            pixRe[i] = LUT_VORT[k];
          }
        }
      }
      offReCtx.putImageData(imgRe, 0, 0);
      ctxRe.fillStyle = css('--deep'); ctxRe.fillRect(0, 0, W, H);
      ctxRe.imageSmoothingEnabled = true; ctxRe.imageSmoothingQuality = 'high';
      ctxRe.drawImage(offRe, R.x, R.y, R.w, R.h);
      if (tinta) {
        const P = particlesRe, sx = R.w / nx, sy = R.h / ny, sz = Math.max(1.6, 1.3 * dpr);
        const cols = [css('--buoy'), css('--fresh')];
        for (let c = 0; c < 2; c++) {
          ctxRe.fillStyle = cols[c];
          for (let i = 0; i < P.n; i++) if (P.c[i] === c) ctxRe.fillRect(R.x + P.x[i] * sx - sz / 2, R.y + P.y[i] * sy - sz / 2, sz, sz);
        }
      }
    }

    function stepRe(budget) {
      const n = clamp(Math.round(budget / RE.msStep), 1, 18);
      const t0 = performance.now();
      for (let i = 0; i < n; i++) lbm.step();
      const ms = (performance.now() - t0) / n;
      RE.msStep += (ms - RE.msStep) * 0.1;
      if (RE.views[RE.view].id === 'tinta') advectParticlesRe(n);
      sampleProbeRe(n);
      if (!lbm.healthy()) { buildShape(); toast('#toast-re', 'La simulación se reinició'); }
    }

    // painting with the mouse
    let painting = null, lastCell = null;
    function cellFromEvent(e) {
      const R = RE.rect; if (!R) return null;
      const b = cvRe.getBoundingClientRect(), dpr = cvRe.width / b.width;
      const px = (e.clientX - b.left) * dpr, py = (e.clientY - b.top) * dpr;
      return { x: (px - R.x) / R.w * RE.nx, y: (py - R.y) / R.h * RE.ny };
    }
    function paintAt(c, solid) {
      const { nx, ny } = RE, N = lbm.N, f = lbm.f, br = 2.6;
      for (let y = Math.floor(c.y - br); y <= c.y + br; y++) for (let x = Math.floor(c.x - br); x <= c.x + br; x++) {
        if (x < 3 || x > nx - 4 || y < 0 || y >= ny) continue;
        const dx = x - c.x, dy = y - c.y; if (dx * dx + dy * dy > br * br) continue;
        const i = y * nx + x;
        if (solid) lbm.solid[i] = 1;
        else if (lbm.solid[i]) { lbm.solid[i] = 0; for (let k = 0; k < 9; k++) f[k * N + i] = LB_W[k]; }
      }
    }
    cvRe.addEventListener('contextmenu', (e) => e.preventDefault());
    cvRe.addEventListener('pointerdown', (e) => {
      painting = !(e.button === 2 || e.shiftKey); lastCell = cellFromEvent(e);
      cvRe.setPointerCapture(e.pointerId); if (lastCell) paintAt(lastCell, painting);
    });
    cvRe.addEventListener('pointermove', (e) => {
      if (painting === null) return;
      const c = cellFromEvent(e); if (!c || !lastCell) return;
      const d = Math.hypot(c.x - lastCell.x, c.y - lastCell.y), n = Math.max(1, Math.ceil(d / 1.2));
      for (let i = 1; i <= n; i++) paintAt({ x: lastCell.x + (c.x - lastCell.x) * i / n, y: lastCell.y + (c.y - lastCell.y) * i / n }, painting);
      lastCell = c;
    });
    const endPaint = () => { if (painting !== null) { painting = null; updateReParams(); } };
    cvRe.addEventListener('pointerup', endPaint); cvRe.addEventListener('pointercancel', endPaint);

    // controls
    const setFluidUI = segButtons($('#re-fluids'), RE.fluids, (i) => { RE.fluid = i; setFluidUI(i); updateReParams(); });
    const setShapeUI = segButtons($('#re-shapes'), RE.shapes, (i) => { RE.shape = i; setShapeUI(i); buildShape(); });
    const setViewUI = segButtons($('#re-views'), RE.views, (i) => { RE.view = i; setViewUI(i); particlesRe.n = 0; });
    setFluidUI(RE.fluid); setShapeUI(RE.shape); setViewUI(RE.view);
    RE.speedT = +$('#re-speed').value;
    $('#re-speed').addEventListener('input', (e) => { RE.speedT = +e.target.value; updateReParams(); });
    $('#re-size').addEventListener('input', (e) => { RE.size = +e.target.value; buildShape(); });
    buildShape();
    return {
      step: stepRe, draw: drawRe, readout: updateRegimeRe, every: 6, reset: buildShape, lr: '#re-speed',
      ud(up) { const i = clamp(RE.fluid + (up ? -1 : 1), 0, RE.fluids.length - 1); $('#re-fluids').children[i].click(); },
      key(k) { if (k === 'v' || k === 'V') { $('#re-views').children[(RE.view + 1) % 2].click(); return true; } return false; },
      demo: [80, 100, 45, 100].map((speed, p) => () => {
        if (p === 0) { RE.shape = 0; setShapeUI(0); RE.view = 0; setViewUI(0); buildShape(); }
        RE.fluid = [0, 1, 2, 2][p]; setFluidUI(RE.fluid); setRange('#re-speed', speed);
      }),
    };
  },
});
