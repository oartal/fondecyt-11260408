Lab.register('iw', {
  html: `
    <section class="module" id="mod-iw" role="tabpanel" aria-labelledby="tab-iw">
      <div class="stage" id="stage-iw">
        <canvas id="cv-iw" aria-label="Simulación de ondas internas sobre un umbral submarino"></canvas>
        <div class="tag tl" id="iw-tide">Marea</div>
        <div class="tag tr" id="iw-modetag">Cálculo completo</div>
        <div class="demo-chip" id="demo-iw" hidden>Modo demostración · mueve el mouse para jugar</div>
        <div class="toast" id="toast-iw" hidden></div>
      </div>
      <div class="panel">
        <div class="card controls">
          <div class="ctl wide">
            <label class="label" for="iw-tidein">Fuerza de la marea <span class="fine">← →</span></label>
            <input type="range" id="iw-tidein" min="0" max="100" step="1" value="45">
            <div class="ends"><span>marea muerta</span><span>marea viva</span></div>
          </div>
          <div class="ctl wide">
            <label class="label" for="iw-fresh">Agua dulce de ríos y lluvia <span class="fine">↑ ↓</span></label>
            <input type="range" id="iw-fresh" min="0" max="100" step="1" value="50">
            <div class="ends"><span>verano seco</span><span>invierno lluvioso</span></div>
          </div>
          <div class="ctl wide">
            <label class="label" for="iw-sill">Altura del umbral <output id="iw-sill-out">—</output></label>
            <input type="range" id="iw-sill" min="0" max="100" step="1" value="55">
          </div>
          <div class="ctl wide">
            <div class="label">Forma de calcular <span class="fine">H</span></div>
            <div class="seg" id="iw-modes" style="--n:2"></div>
          </div>
        </div>
        <div class="card gauge">
          <div class="label">Número de Froude · corriente contra olas</div>
          <div class="bignum"><span class="sym">Fr =</span><span class="val" id="fr-val">—</span></div>
          <div class="regime" id="fr-regime"><i></i><span>—</span></div>
          <p class="sub" id="fr-regime-sub">—</p>
          <div class="gaugebox" id="fr-gauge"></div>
        </div>
        <div class="card explain">
          <h2>¿Puede haber olas debajo del mar?</h2>
          <p>El agua dulce de ríos y lluvia es más liviana y flota sobre el agua salada. La frontera entre las dos puede ondularse como una ola, aunque arriba el mar se vea tranquilo.</p>
          <p>Esas olas escondidas nacen cuando la marea pasa sobre un cerro submarino, que se llama <b>umbral</b>.</p>
          <p class="fine" id="iw-modeexplain"></p>
        </div>
      </div>
    </section>`,
  init() {
    'use strict';
    // ==================================================================
    //                    MÓDULO 2: ONDAS INTERNAS
    // ==================================================================
    const IW = { nx: 240, nz: 60, dt: 0.02, msStep: 4, Hm: 100, G: 0, frame: 0 };
    const iw = new InternalWaves(IW.nx, IW.nz, 4, 1);
    iw.T = 16;
    const cvIw = $('#cv-iw'), ctxIw = cvIw.getContext('2d');
    const offIw = document.createElement('canvas'); offIw.width = IW.nx; offIw.height = IW.nz + 1;
    const offIwCtx = offIw.getContext('2d');
    const imgIw = offIwCtx.createImageData(IW.nx, IW.nz + 1);
    const pixIw = new Uint32Array(imgIw.data.buffer);

    const iwMap = {
      U0: (t) => 0.04 + 0.46 * t / 100,
      db: (t) => 0.35 + 1.65 * t / 100,
      hs: (t) => 0.25 + 0.37 * t / 100,
    };
    function iwLabels() {
      $('#iw-sill-out').textContent = Math.round(iw.hs * IW.Hm) + ' m de alto';
    }
    let sillTimer = null;
    $('#iw-tidein').addEventListener('input', (e) => { iw.U0 = iwMap.U0(+e.target.value); iwLabels(); });
    $('#iw-fresh').addEventListener('input', (e) => {
      const nd = iwMap.db(+e.target.value), r = nd / iw.db; iw.db = nd;
      for (let j = 0; j < iw.n; j++) iw.b[j] *= r;
      iw.setB0(); iwLabels();
    });
    $('#iw-sill').addEventListener('input', (e) => {
      iw.hs = iwMap.hs(+e.target.value); iwLabels();
      clearTimeout(sillTimer); sillTimer = setTimeout(() => { iw.build(false); resetParticlesIw(); }, 180);
    });
    const modes = [{ id: 'hid', name: 'Simplificado', hint: 'hidrostático' }, { id: 'nohid', name: 'Completo', hint: 'no hidrostático' }];
    const setModeUI = segButtons($('#iw-modes'), modes, (i) => setMode(i === 1));
    function setMode(nonhyd) {
      iw.nonhyd = nonhyd;
      // hydrostatic models need more mixing/diffusion to stay smooth
      iw.nu = iw.kappa = nonhyd ? 3e-4 : 2e-3;
      iw.buildBarotropic(); iw.solvePsi(nonhyd ? 400 : 0);
      setModeUI(nonhyd ? 1 : 0);
      $('#iw-modetag').textContent = nonhyd ? 'Cálculo completo' : 'Cálculo simplificado';
      $('#iw-modeexplain').textContent = nonhyd
        ? 'Cálculo completo. Incluye el sube y baja del agua, y por eso aparecen olas cortas. Es el que estudiamos en este proyecto.'
        : 'Cálculo simplificado. Es más rápido y lo usa la mayoría de los modelos del mar, pero pierde las olas cortas. ¿Notas la diferencia?';
    }
    iw.U0 = iwMap.U0(+$('#iw-tidein').value);
    iw.db = iwMap.db(+$('#iw-fresh').value);
    iw.hs = iwMap.hs(+$('#iw-sill').value);
    iw.build(true);
    setMode(true); iwLabels();

    const needleFr = buildGauge($('#fr-gauge'), {
      aria: 'Escala del número de Froude interno', min: 0, max: 2.5, map: (v) => v,
      zones: [
        { a: 0, b: 0.8, color: css('--calm') },
        { a: 0.8, b: 1.2, color: css('--warn') },
        { a: 1.2, b: 2.5, color: css('--hot') },
      ],
      ticks: [{ v: 0, label: '0', anchor: 'start' }, { v: 0.5, label: '0.5' }, { v: 1, label: '1' }, { v: 1.5, label: '1.5' }, { v: 2, label: '2' }, { v: 2.5, label: '2.5', anchor: 'end' }],
      marks: [{ v: 0.4, label: 'ganan las olas' }, { v: 1, label: 'empate' }, { v: 1.85, label: 'gana la corriente' }],
    });
    function setRegimeFr(G) {
      const box = $('#fr-regime'), dot = box.querySelector('i'), txt = box.querySelector('span'), sub = $('#fr-regime-sub');
      let c, t, s;
      if (G < 0.8) { c = '--calm'; t = 'Ganan las olas'; s = 'Las olas escondidas viajan más rápido que la corriente y se alejan del umbral.'; }
      else if (G < 1.2) { c = '--warn'; t = 'Empate'; s = 'La corriente va tan rápido como las olas, que quedan atrapadas sobre el umbral.'; }
      else { c = '--hot'; t = 'Gana la corriente'; s = 'La corriente arrastra las olas. Detrás del umbral el agua cae como una cascada bajo el mar.'; }
      if (txt.textContent !== t) { dot.style.background = `var(${c})`; txt.textContent = t; sub.textContent = s; }
    }

    // tracer particles
    const PIW = 520;
    const particlesIw = { x: new Float32Array(PIW), z: new Float32Array(PIW), life: new Float32Array(PIW) };
    function hAt(x) { const s = (x - iw.Lx / 2) / iw.ws; return iw.hs * Math.exp(-s * s); }
    function spawnIw(i) {
      const x = Math.random() * iw.Lx, h = hAt(x);
      particlesIw.x[i] = x; particlesIw.z[i] = h + 0.02 + Math.random() * (iw.H - h - 0.04);
      particlesIw.life[i] = 4 + Math.random() * 10;
    }
    function resetParticlesIw() { for (let i = 0; i < PIW; i++) spawnIw(i); }
    resetParticlesIw();
    function advectParticlesIw(tsum) {
      const P = particlesIw;
      for (let i = 0; i < PIW; i++) {
        let x = P.x[i], z = P.z[i];
        x += iw.interp(iw.u, x, z) * tsum; z += iw.interp(iw.w, x, z) * tsum;
        if (x < 0) x += iw.Lx; if (x >= iw.Lx) x -= iw.Lx;
        P.life[i] -= tsum;
        if (P.life[i] <= 0 || z <= hAt(x) + 0.005 || z >= iw.H - 0.003) { spawnIw(i); continue; }
        P.x[i] = x; P.z[i] = z;
      }
    }

    function stepIw(budget) {
      const n = clamp(Math.round(budget / IW.msStep), 1, 8);
      const t0 = performance.now();
      for (let i = 0; i < n; i++) iw.step(IW.dt);
      const ms = (performance.now() - t0) / n;
      IW.msStep += (ms - IW.msStep) * 0.1;
      iw.velocities();
      advectParticlesIw(n * IW.dt);
      if (!iw.healthy()) { iw.build(true); resetParticlesIw(); toast('#toast-iw', 'La simulación se reinició'); }
    }

    function drawIw() {
      const dpr = fitCanvas(cvIw);
      const W = cvIw.width, H = cvIw.height, { nx, nz } = IW;
      const R = letterbox(W, H, 3);
      const b = iw.b, solid = iw.solid, inv = 1 / iw.db;
      for (let k = 0; k <= nz; k++) {
        const row = (nz - k) * nx, r = k * nx;
        for (let i = 0; i < nx; i++) {
          const j = r + i;
          let t = solid[j] ? 0 : (b[j] * inv + 0.5);
          t = t < 0 ? 0 : t > 1 ? 1 : t;
          pixIw[row + i] = LUT_RHO[(t * 255) | 0];
        }
      }
      offIwCtx.putImageData(imgIw, 0, 0);
      const c = ctxIw;
      c.fillStyle = css('--deep'); c.fillRect(0, 0, W, H);
      c.imageSmoothingEnabled = true; c.imageSmoothingQuality = 'high';
      c.drawImage(offIw, R.x, R.y, R.w, R.h);
      const X = (x) => R.x + (x / iw.dx + 0.5) / nx * R.w;
      const Y = (z) => R.y + ((nz - z / iw.dz) + 0.5) / (nz + 1) * R.h;
      // sill (crisp)
      c.beginPath(); c.moveTo(R.x, R.y + R.h);
      for (let i = 0; i <= 120; i++) { const x = i / 120 * iw.Lx; c.lineTo(R.x + i / 120 * R.w, Y(hAt(x))); }
      c.lineTo(R.x + R.w, R.y + R.h); c.closePath();
      const grd = c.createLinearGradient(0, Y(iw.hs), 0, R.y + R.h);
      grd.addColorStop(0, css('--rock')); grd.addColorStop(1, '#6f5f45');
      c.fillStyle = grd; c.fill();
      // particles
      c.fillStyle = 'rgba(233, 242, 244, 0.55)';
      const sz = Math.max(1.6, 1.4 * dpr), P = particlesIw;
      for (let i = 0; i < PIW; i++) c.fillRect(X(P.x[i]) - sz / 2, Y(P.z[i]) - sz / 2, sz, sz);
      // surface
      c.strokeStyle = 'rgba(233, 242, 244, 0.8)'; c.lineWidth = 2 * dpr;
      c.beginPath(); c.moveTo(R.x, R.y + dpr); c.lineTo(R.x + R.w, R.y + dpr); c.stroke();
      // labels
      const fs = Math.round(clamp(R.h * 0.045, 12 * dpr, 20 * dpr));
      c.font = `700 ${fs}px ${css('--f-body')}`; c.textBaseline = 'middle';
      c.fillStyle = '#0a1822'; c.textAlign = 'left';
      c.fillText('agua dulce', R.x + 14 * dpr, R.y + R.h * 0.2);
      c.fillStyle = 'rgba(233, 242, 244, 0.9)';
      c.fillText('agua salada', R.x + 14 * dpr, R.y + R.h * 0.82);
      c.fillStyle = '#2b2416'; c.textAlign = 'center';
      c.fillText('umbral', X(iw.Lx / 2), Y(iw.hs) + 0.55 * (R.y + R.h - Y(iw.hs)));
      c.fillStyle = 'rgba(10, 24, 34, 0.8)'; c.textAlign = 'center';
      c.font = `700 ${Math.round(fs * 0.8)}px ${css('--f-body')}`;
      c.fillText('arriba, la superficie se ve tranquila', R.x + R.w / 2, R.y + fs * 0.9);
      // depth axis
      c.textAlign = 'right'; c.fillStyle = 'rgba(233, 242, 244, 0.75)'; c.strokeStyle = 'rgba(233, 242, 244, 0.5)'; c.lineWidth = dpr;
      for (const d of [50, 100]) {
        const yy = Y(1 - d / IW.Hm) + (d === 100 ? -fs * 0.7 : 0);
        c.beginPath(); c.moveTo(R.x + R.w - 6 * dpr, Y(1 - d / IW.Hm)); c.lineTo(R.x + R.w, Y(1 - d / IW.Hm)); c.stroke();
        c.lineWidth = 3 * dpr; c.strokeStyle = 'rgba(10, 24, 34, 0.75)'; c.strokeText(d + ' m', R.x + R.w - 10 * dpr, yy);
        c.fillText(d + ' m', R.x + R.w - 10 * dpr, yy); c.lineWidth = dpr; c.strokeStyle = 'rgba(233, 242, 244, 0.5)';
      }
    }
    function updateTideTag() {
      const q = Math.sin(2 * Math.PI * iw.t / iw.T);
      const n = Math.round(Math.abs(q) * 4);
      const bars = '▮'.repeat(n) + '▯'.repeat(4 - n);
      const hr = ((iw.t / iw.T) % 1) * 12.4;
      const dir = Math.abs(q) < 0.12 ? 'Estoa (cambio de marea)' : q > 0 ? 'Corriente de marea →' : '← Corriente de marea';
      $('#iw-tide').innerHTML = `<span class="arrow">${bars}</span> ${dir} · hora ${hr.toFixed(1)} de 12.4`;
    }
    return {
      step: stepIw, draw: drawIw, every: 3, lr: '#iw-tidein', ud: '#iw-fresh',
      readout() {
        const F = iw.froude(); IW.G += (F.G - IW.G) * 0.25;
        $('#fr-val').textContent = IW.G.toFixed(2); needleFr(IW.G); setRegimeFr(IW.G);
        updateTideTag();
      },
      reset() { iw.build(true); resetParticlesIw(); },
      key(k) { if (k === 'h' || k === 'H') { setMode(!iw.nonhyd); return true; } return false; },
      demo: [
        () => { setRange('#iw-tidein', 70); iw.build(true); resetParticlesIw(); setMode(true); },
        null,
        () => setMode(false),
      ],
    };
  },
});
