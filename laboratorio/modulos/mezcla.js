// ===================== Stratified shear layer (Kelvin–Helmholtz) =====================
// same solver, flat bottom, no tide, fully periodic in x. The horizontal-mean shear and
// density profiles are nudged back to the chosen targets, so billows keep forming.
class ShearLayer extends InternalWaves {
  constructor(nx = 240, nz = 60, Lx = 4, H = 1) {
    super(nx, nz, Lx, H);
    this.hs = 0; this.U0 = 0; this.Ls = 0; this.nonhyd = true;
    this.zi = 0.5 * (H + this.dz); this.delta = 0.06;
    this.nu = 1.2e-4; this.kappa = 1.2e-4; this.iters = 20;
    this.dU = 0.6; this.db = 0.5; this.tau = 7; this.noise = 0.06;
    this.omT = new Float64Array(nz + 1);
    this.restart();
  }
  buildBarotropic() { /* no tide in this module */ }
  // gradient Richardson number at the interface for the chosen profiles: N² / S²
  richardson() { return 2 * this.db * this.delta / (this.dU * this.dU); }
  setTargets() {
    this.setB0();
    for (let k = 0; k <= this.nz; k++) {
      const s = 1 / Math.cosh((k * this.dz - this.zi) / this.delta);
      this.omT[k] = 0.5 * this.dU / this.delta * s * s;
    }
  }
  restart() {
    this.build(true); this.setTargets();
    const { nx, nz } = this;
    for (let k = 1; k < nz; k++) for (let i = 0; i < nx; i++) this.om[k * nx + i] = this.omT[k];
    this.solvePsi(600);
  }
  step(dt) {
    const { nx, nz, om, b, b0, omT } = this;
    const r = dt / this.tau, amp = this.noise * this.dU / this.delta;
    for (let k = 1; k < nz; k++) {
      const row = k * nx; let so = 0, sb = 0;
      for (let i = 0; i < nx; i++) { so += om[row + i]; sb += b[row + i]; }
      const co = r * (omT[k] - so / nx), cb = r * (b0[row] - sb / nx);
      const seed = Math.abs(k * this.dz - this.zi) < 2 * this.delta;
      for (let i = 0; i < nx; i++) {
        om[row + i] += co + (seed ? amp * (Math.random() - 0.5) * dt : 0);
        b[row + i] += cb;
      }
    }
    super.step(dt);
  }
  // thickness of the blended zone, from the horizontal-mean density profile (fraction of the depth)
  mixing() {
    const { nx, nz, b } = this; const q = 2 / this.db; let s = 0;
    for (let k = 1; k <= nz; k++) {
      const row = k * nx; let m = 0;
      for (let i = 0; i < nx; i++) m += b[row + i];
      m = m / nx * q; s += 1 - m * m;
    }
    return s / nz;
  }
  // rms vertical displacement of the interface (in units of H): how rolled-up it is
  wiggle() {
    const { nx, nz, dz, w } = this; let s = 0;
    const k = Math.round(this.zi / dz);
    for (let i = 0; i < nx; i++) { const v = w[k * nx + i]; s += v * v; }
    return Math.sqrt(s / nx) / Math.max(this.dU, 1e-6);
  }
}

Lab.register('ri', {
  html: `
    <section class="module" id="mod-ri" role="tabpanel" aria-labelledby="tab-ri">
      <div class="stage" id="stage-ri">
        <canvas id="cv-ri" aria-label="Simulación de dos capas de agua que se deslizan una sobre otra"></canvas>
        <div class="tag tl">Dos capas que se mueven en sentido contrario</div>
        <div class="tag tr">Simulación en 2D · prototipo</div>
        <div class="demo-chip" id="demo-ri" hidden>Modo demostración · mueve el mouse para jugar</div>
        <div class="toast" id="toast-ri" hidden></div>
      </div>
      <div class="panel">
        <div class="card controls">
          <div class="ctl wide">
            <label class="label" for="ri-shear">Diferencia de velocidad <span class="fine">← →</span></label>
            <input type="range" id="ri-shear" min="0" max="100" step="1" value="50">
            <div class="ends"><span>poca</span><span>mucha</span></div>
          </div>
          <div class="ctl wide">
            <label class="label" for="ri-dens">Diferencia de peso entre las aguas <span class="fine">↑ ↓</span></label>
            <input type="range" id="ri-dens" min="0" max="100" step="1" value="12">
            <div class="ends"><span>poca · verano seco</span><span>mucha · invierno lluvioso</span></div>
          </div>
          <div class="ctl wide">
            <button class="pill" id="ri-reset" type="button">Empezar de nuevo <span class="fine">R</span></button>
          </div>
        </div>
        <div class="card gauge">
          <div class="label">Número de Richardson · capas contra corriente</div>
          <div class="bignum"><span class="sym">Ri =</span><span class="val" id="ri-val">—</span></div>
          <div class="regime" id="ri-regime"><i></i><span>—</span></div>
          <p class="sub" id="ri-regime-sub">—</p>
          <div class="gaugebox" id="ri-gauge"></div>
          <div class="mixrow"><span class="label">Mezcla</span><div class="mixbar"><i id="ri-mix"></i></div><output id="ri-mix-out">0 %</output></div>
        </div>
        <div class="card explain">
          <h2>¿Se mezcla o no se mezcla?</h2>
          <p>Arriba va el agua dulce y abajo el agua salada, que es más pesada. Cada una se mueve para un lado distinto.</p>
          <p>Si la diferencia de velocidad es grande, la frontera se enrolla en remolinos y las aguas se mezclan. Si la diferencia de peso es grande, las capas aguantan separadas.</p>
          <p class="fine">Los científicos miden esta pelea con el número de Richardson. Cuando es chico, gana la corriente.</p>
        </div>
      </div>
    </section>`,
  init() {
    'use strict';
    // ==================================================================
    //                    MÓDULO 3: RICHARDSON
    // ==================================================================
    const RI = { nx: 240, nz: 60, dt: 0.03, msStep: 5, mix: 0, regime: '' };
    const sl = new ShearLayer(RI.nx, RI.nz, 4, 1);
    const cvRi = $('#cv-ri'), ctxRi = cvRi.getContext('2d');
    const offRi = document.createElement('canvas'); offRi.width = RI.nx; offRi.height = RI.nz + 1;
    const offRiCtx = offRi.getContext('2d');
    const imgRi = offRiCtx.createImageData(RI.nx, RI.nz + 1);
    const pixRi = new Uint32Array(imgRi.data.buffer);
    const riMap = {
      dU: (t) => 0.4 + 0.6 * t / 100,
      db: (t) => 0.05 + 1.45 * t / 100,
    };
    const PRI = 420;
    const particlesRi = { x: new Float32Array(PRI), z: new Float32Array(PRI), life: new Float32Array(PRI) };
    function spawnRi(i) {
      particlesRi.x[i] = Math.random() * sl.Lx; particlesRi.z[i] = 0.04 + Math.random() * (sl.H - 0.08);
      particlesRi.life[i] = 4 + Math.random() * 10;
    }
    function resetParticlesRi() { for (let i = 0; i < PRI; i++) spawnRi(i); }
    function restartRi() { sl.restart(); sl.velocities(); resetParticlesRi(); RI.mix = 0; }
    function setShear(t) {
      const nd = riMap.dU(t), r = nd / sl.dU; sl.dU = nd;
      for (let j = 0; j < sl.n; j++) sl.om[j] *= r;
      sl.setTargets(); sl.solvePsi(60);
    }
    function setDens(t) {
      const nd = riMap.db(t), r = nd / sl.db; sl.db = nd;
      for (let j = 0; j < sl.n; j++) sl.b[j] *= r;
      sl.setTargets();
    }
    sl.dU = riMap.dU(+$('#ri-shear').value); sl.db = riMap.db(+$('#ri-dens').value);
    restartRi();
    $('#ri-shear').addEventListener('input', (e) => setShear(+e.target.value));
    $('#ri-dens').addEventListener('input', (e) => setDens(+e.target.value));
    $('#ri-reset').addEventListener('click', (e) => { restartRi(); if (e.detail > 0) e.currentTarget.blur(); });

    // In this small simulation the layers start to roll up near Ri ≈ 0.15; the theoretical limit is 0.25.
    const needleRi = buildGauge($('#ri-gauge'), {
      aria: 'Escala del número de Richardson', min: 0.005, max: 2, map: Math.log10,
      zones: [
        { a: 0.005, b: 0.15, color: css('--hot') },
        { a: 0.15, b: 0.25, color: css('--warn') },
        { a: 0.25, b: 2, color: css('--calm') },
      ],
      ticks: [{ v: 0.01, label: '0.01' }, { v: 0.1, label: '0.1' }, { v: 0.25, label: '0.25' }, { v: 1, label: '1' }],
      marks: [{ v: 0.027, label: 'gana la corriente' }, { v: 0.19, label: 'límite' }, { v: 0.7, label: 'ganan las capas' }],
    });
    function setRegimeRi(Ri) {
      let c, t, s;
      if (Ri < 0.15) { c = '--hot'; t = 'Gana la corriente'; s = 'La frontera se enrolla en remolinos y el agua dulce se mezcla con la salada.'; }
      else if (Ri < 0.25) { c = '--warn'; t = 'Al límite'; s = 'Las capas aguantan apenas. Con un poco más de velocidad se enrollan.'; }
      else { c = '--calm'; t = 'Ganan las capas'; s = 'La diferencia de peso mantiene las dos aguas separadas, cada una en su capa.'; }
      if (RI.regime === t) return;
      RI.regime = t;
      const box = $('#ri-regime');
      box.querySelector('i').style.background = `var(${c})`; box.querySelector('span').textContent = t; $('#ri-regime-sub').textContent = s;
    }
    function updateRiReadout() {
      const Ri = sl.richardson();
      $('#ri-val').textContent = Ri.toFixed(2); needleRi(Ri); setRegimeRi(Ri);
      // 0 % = two clean layers, 100 % = blended zone as thick as 60 % of the depth
      const m = clamp((sl.mixing() - 0.15) / 0.45, 0, 1);
      RI.mix += (m - RI.mix) * 0.12;
      const pct = Math.round(RI.mix * 100);
      $('#ri-mix').style.width = pct + '%'; $('#ri-mix-out').textContent = pct + ' %';
    }
    function stepRi(budget) {
      const n = clamp(Math.round(budget / RI.msStep), 1, 6);
      const t0 = performance.now();
      for (let i = 0; i < n; i++) sl.step(RI.dt);
      const ms = (performance.now() - t0) / n;
      RI.msStep += (ms - RI.msStep) * 0.1;
      sl.velocities();
      const P = particlesRi, tsum = n * RI.dt;
      for (let i = 0; i < PRI; i++) {
        let x = P.x[i], z = P.z[i];
        x += sl.interp(sl.u, x, z) * tsum; z += sl.interp(sl.w, x, z) * tsum;
        if (x < 0) x += sl.Lx; if (x >= sl.Lx) x -= sl.Lx;
        P.life[i] -= tsum;
        if (P.life[i] <= 0 || z <= 0.02 || z >= sl.H - 0.003) { spawnRi(i); continue; }
        P.x[i] = x; P.z[i] = z;
      }
      if (!sl.healthy()) { restartRi(); toast('#toast-ri', 'La simulación se reinició'); }
    }
    function drawRi() {
      const dpr = fitCanvas(cvRi);
      const W = cvRi.width, H = cvRi.height, { nx, nz } = RI;
      const R = letterbox(W, H, 3);
      const b = sl.b, inv = 1 / sl.db;
      for (let k = 0; k <= nz; k++) {
        const row = (nz - k) * nx, r = Math.max(k, 1) * nx;
        for (let i = 0; i < nx; i++) {
          let t = b[r + i] * inv + 0.5;
          t = t < 0 ? 0 : t > 1 ? 1 : t;
          pixRi[row + i] = LUT_RHO[(t * 255) | 0];
        }
      }
      offRiCtx.putImageData(imgRi, 0, 0);
      const c = ctxRi;
      c.fillStyle = css('--deep'); c.fillRect(0, 0, W, H);
      c.imageSmoothingEnabled = true; c.imageSmoothingQuality = 'high';
      c.drawImage(offRi, R.x, R.y, R.w, R.h);
      const X = (x) => R.x + (x / sl.dx + 0.5) / nx * R.w;
      const Y = (z) => R.y + ((nz - z / sl.dz) + 0.5) / (nz + 1) * R.h;
      c.fillStyle = 'rgba(233, 242, 244, 0.55)';
      const sz = Math.max(1.6, 1.4 * dpr), P = particlesRi;
      for (let i = 0; i < PRI; i++) c.fillRect(X(P.x[i]) - sz / 2, Y(P.z[i]) - sz / 2, sz, sz);
      // labels and layer arrows (arrow length follows the chosen velocity difference)
      const fs = Math.round(clamp(R.h * 0.05, 12 * dpr, 22 * dpr));
      c.font = `700 ${fs}px ${css('--f-body')}`; c.textBaseline = 'middle'; c.textAlign = 'left';
      const len = R.w * (0.05 + 0.11 * (sl.dU - 0.4) / 0.6), x0 = R.x + 16 * dpr, hd = fs * 0.5;
      const arrow = (y, dir, col) => {
        const a = dir > 0 ? x0 : x0 + len, bb = dir > 0 ? x0 + len : x0;
        c.strokeStyle = col; c.fillStyle = col; c.lineWidth = 3 * dpr; c.lineCap = 'round';
        c.beginPath(); c.moveTo(a, y); c.lineTo(bb - dir * hd, y); c.stroke();
        c.beginPath(); c.moveTo(bb, y); c.lineTo(bb - dir * hd * 1.4, y - hd * 0.8); c.lineTo(bb - dir * hd * 1.4, y + hd * 0.8); c.closePath(); c.fill();
      };
      const yTop = R.y + R.h * 0.14, yBot = R.y + R.h * 0.86;
      arrow(yTop, 1, '#0a1822'); arrow(yBot, -1, 'rgba(233, 242, 244, 0.9)');
      c.fillStyle = '#0a1822'; c.fillText('agua dulce', x0 + len + 12 * dpr, yTop);
      c.fillStyle = 'rgba(233, 242, 244, 0.9)'; c.fillText('agua salada', x0 + len + 12 * dpr, yBot);
    }
    return {
      step: stepRi, draw: drawRi, readout: updateRiReadout, every: 3, reset: restartRi, lr: '#ri-shear', ud: '#ri-dens',
      demo: [
        () => { setRange('#ri-shear', 80); setRange('#ri-dens', 5); restartRi(); },
        null,
        () => setRange('#ri-dens', 75),
      ],
    };
  },
});
