// ==================================================================
//                    MÓDULO 4: ROSSBY
// ==================================================================
// Ajuste geostrófico de una mancha de agua liviana, vista desde arriba.
// Aguas someras lineales de una capa con gravedad reducida, en plano f:
//   eta_t = -H (u_x + v_y) ;  u_t = +f v - g' eta_x ;  v_t = -f u - g' eta_y
// Malla C, avance forward–backward, paredes con esponja en los bordes.
// Número en pantalla: Ro = c / (f L) = radio de Rossby / radio de la mancha, con c = sqrt(g' H).
class RotatingSW {
  constructor(nx = 240, ny = 80, dx = 1250) {
    this.nx = nx; this.ny = ny; this.dx = dx; this.n = nx * ny;
    this.eta = new Float64Array(this.n); this.u = new Float64Array(this.n); this.v = new Float64Array(this.n);
    this.H = 20; this.gp = 0.02; this.f = -1e-4; this.dt = 300; this.t = 0; this.sp = 14;
    this.damp = new Float64Array(this.n);
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const d = Math.min(i, nx - 1 - i, j, ny - 1 - j);
      const s = d < this.sp ? 1 - d / this.sp : 0;
      this.damp[j * nx + i] = 1 / (1 + 0.12 * s * s);
    }
  }
  get c() { return Math.sqrt(this.gp * this.H); }
  clear() { this.eta.fill(0); this.u.fill(0); this.v.fill(0); this.t = 0; }
  blob(x0, y0, L, A) { // x0, y0, L in metres; A in metres of interface displacement
    const { nx, ny, dx, eta } = this;
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const rx = (i + 0.5) * dx - x0, ry = (j + 0.5) * dx - y0;
      eta[j * nx + i] += A * Math.exp(-(rx * rx + ry * ry) / (L * L));
    }
  }
  step() {
    const { nx, ny, dx, eta, u, v, damp, dt, f } = this;
    const ch = dt * this.H / dx, cg = dt * this.gp / dx, fd = f * dt;
    // u[j,i] lives on the east face of cell i, v[j,i] on the north face of cell j
    for (let j = 0; j < ny; j++) {
      const r = j * nx, rs = r - nx;
      for (let i = 0; i < nx; i++) {
        const k = r + i;
        const uw = i > 0 ? u[k - 1] : 0, vs = j > 0 ? v[rs + i] : 0;
        eta[k] -= ch * (u[k] - uw + v[k] - vs);
      }
    }
    const flip = (Math.round(this.t / dt) & 1) === 1;
    const doU = () => {
      for (let j = 0; j < ny; j++) {
        const r = j * nx, rs = r - nx;
        for (let i = 0; i < nx - 1; i++) {
          const k = r + i;
          const vb = 0.25 * (v[k] + v[k + 1] + (j > 0 ? v[rs + i] + v[rs + i + 1] : 0));
          u[k] += fd * vb - cg * (eta[k + 1] - eta[k]);
        }
        u[r + nx - 1] = 0;
      }
    };
    const doV = () => {
      for (let j = 0; j < ny - 1; j++) {
        const r = j * nx, rn = r + nx;
        for (let i = 0; i < nx; i++) {
          const k = r + i;
          const ub = 0.25 * (u[k] + u[rn + i] + (i > 0 ? u[k - 1] + u[rn + i - 1] : 0));
          v[k] += -fd * ub - cg * (eta[rn + i] - eta[k]);
        }
      }
      v.fill(0, (ny - 1) * nx, ny * nx);
    };
    if (flip) { doV(); doU(); } else { doU(); doV(); }
    for (let k = 0; k < this.n; k++) { const d = damp[k]; if (d < 1) { eta[k] *= d; u[k] *= d; v[k] *= d; } }
    this.t += dt;
  }
  // velocity at cell centre (i, j may be fractional cell coordinates)
  vel(x, y) {
    const { nx, ny, u, v } = this;
    const i = Math.min(nx - 1, Math.max(0, Math.floor(x))), j = Math.min(ny - 1, Math.max(0, Math.floor(y)));
    const k = j * nx + i;
    return [0.5 * (u[k] + (i > 0 ? u[k - 1] : 0)), 0.5 * (v[k] + (j > 0 ? v[k - nx] : 0))];
  }
  volume(x0, y0, R) { // volume of light water (positive eta) within R of (x0, y0)
    const { nx, ny, dx, eta } = this; let s = 0;
    const i0 = Math.max(0, Math.floor((x0 - R) / dx)), i1 = Math.min(nx - 1, Math.ceil((x0 + R) / dx));
    const j0 = Math.max(0, Math.floor((y0 - R) / dx)), j1 = Math.min(ny - 1, Math.ceil((y0 + R) / dx));
    for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
      const rx = (i + 0.5) * dx - x0, ry = (j + 0.5) * dx - y0;
      if (rx * rx + ry * ry <= R * R) s += eta[j * nx + i];
    }
    return s * dx * dx;
  }
  healthy() { for (let k = 0; k < this.n; k += 37) if (!(Math.abs(this.eta[k]) < 1e3)) return false; return true; }
}

Lab.register('ro', {
  html: `
    <section class="module" id="mod-ro" role="tabpanel" aria-labelledby="tab-ro">
      <div class="stage" id="stage-ro">
        <canvas id="cv-ro" aria-label="Simulación, vista desde arriba, de una mancha de agua dulce sobre un mar que gira con la Tierra"></canvas>
        <div class="tag tl" id="ro-time">Tiempo</div>
        <div class="tag tr">Vista desde arriba · prototipo</div>
        <div class="tag bl">Clic: suelta otra mancha de agua dulce</div>
        <div class="demo-chip" id="demo-ro" hidden>Modo demostración · mueve el mouse para jugar</div>
        <div class="toast" id="toast-ro" hidden></div>
      </div>
      <div class="panel">
        <div class="card controls">
          <div class="ctl wide">
            <label class="label" for="ro-lat">Latitud <output id="ro-lat-out">—</output></label>
            <input type="range" id="ro-lat" min="0" max="70" step="1" value="42">
            <div class="ends"><span>ecuador <span class="fine">←</span></span><span><span class="fine">→</span> cerca del polo</span></div>
          </div>
          <div class="ctl wide">
            <label class="label" for="ro-size">Tamaño de la mancha <output id="ro-size-out">—</output></label>
            <input type="range" id="ro-size" min="4" max="40" step="1" value="20">
            <div class="ends"><span>chica <span class="fine">↓</span></span><span><span class="fine">↑</span> grande</span></div>
          </div>
          <div class="ctl">
            <div class="label">Hemisferio <span class="fine">N</span></div>
            <div class="seg" id="ro-hemis" style="--n:2"></div>
          </div>
          <div class="ctl">
            <div class="label">&nbsp;</div>
            <button class="pill" id="ro-reset" type="button">Soltar de nuevo <span class="fine">R</span></button>
          </div>
        </div>
        <div class="card gauge">
          <div class="label">Número de Rossby · olas contra el giro de la Tierra</div>
          <div class="bignum"><span class="sym">Ro =</span><span class="val" id="ro-val">—</span></div>
          <div class="regime" id="ro-regime"><i></i><span>—</span></div>
          <p class="sub" id="ro-regime-sub">—</p>
          <div class="gaugebox" id="ro-gauge"></div>
          <div class="mixrow"><span class="label">Atrapada</span><div class="mixbar"><i id="ro-trap"></i></div><output id="ro-trap-out">— %</output></div>
        </div>
        <div class="card explain">
          <h2>¿Por qué los remolinos del mar duran tanto?</h2>
          <p>La Tierra gira, y todo lo que se mueve sobre ella se desvía: hacia la izquierda en el hemisferio sur y hacia la derecha en el norte.</p>
          <p>Una mancha chica de agua dulce se desparrama en olas antes de notar el giro. Una mancha grande no alcanza a escapar: el giro la atrapa y queda dando vueltas como un remolino.</p>
          <p class="fine" id="ro-fine">Los científicos comparan esas dos fuerzas con el número de Rossby.</p>
        </div>
      </div>
    </section>`,
  init() {
    'use strict';
    const RO = { nx: 240, ny: 80, dxm: 1250, A: 10, stepsPerFrame: 2, lat: 42, Lkm: 20, south: true, x0: 0, y0: 0, V0: 1, trap: 1, regime: '' };
    const sw = new RotatingSW(RO.nx, RO.ny, RO.dxm);
    const OMEGA = 7.2921e-5;
    const cv = $('#cv-ro'), ctx = cv.getContext('2d');
    const off = document.createElement('canvas'); off.width = RO.nx; off.height = RO.ny;
    const offCtx = off.getContext('2d');
    const img = offCtx.createImageData(RO.nx, RO.ny);
    const pix = new Uint32Array(img.data.buffer);
    const LUT = new Uint32Array(256);
    (() => {
      const low = hexToRgb('#0b1838'), sea = hexToRgb('#1d5878'), fresh = hexToRgb('#9be8cf');
      for (let i = 0; i < 256; i++) { const t = i / 255; LUT[i] = pack(t < 0.5 ? mix(low, sea, t / 0.5) : mix(sea, fresh, (t - 0.5) / 0.5)); }
    })();

    const rossbyRadius = () => { const f = Math.abs(sw.f); return f > 1e-9 ? sw.c / f : Infinity; }; // metres
    const rossby = () => rossbyRadius() / (RO.Lkm * 1000);
    function setF() { sw.f = 2 * OMEGA * Math.sin(RO.lat * Math.PI / 180) * (RO.south ? -1 : 1); }
    function labels() {
      $('#ro-lat-out').textContent = RO.lat === 0 ? '0° · ecuador' : RO.lat + '° ' + (RO.south ? 'S' : 'N');
      $('#ro-size-out').textContent = RO.Lkm + ' km de radio';
      const Rd = rossbyRadius();
      $('#ro-fine').textContent = 'Los científicos comparan esas dos fuerzas con el número de Rossby. El tamaño donde empatan se llama radio de Rossby' +
        (isFinite(Rd) ? ': aquí mide ' + fmtNum(Rd / 1000) + ' km.' : '. En el ecuador el giro no se siente.');
    }

    const NP = 700;
    const P = { x: new Float32Array(NP), y: new Float32Array(NP), life: new Float32Array(NP) };
    function spawn(i) { P.x[i] = 2 + Math.random() * (RO.nx - 4); P.y[i] = 2 + Math.random() * (RO.ny - 4); P.life[i] = 20 + Math.random() * 60; }
    function drop(xm, ym, keep) {
      if (!keep) sw.clear();
      RO.x0 = xm; RO.y0 = ym;
      sw.blob(xm, ym, RO.Lkm * 1000, RO.A);
      RO.V0 = Math.max(1e-9, sw.volume(xm, ym, 2 * RO.Lkm * 1000)); RO.trap = 1; RO.tDrop = sw.t;
    }
    function restart() { setF(); drop(RO.nx * RO.dxm / 2, RO.ny * RO.dxm / 2, false); for (let i = 0; i < NP; i++) spawn(i); labels(); }

    let timer = null;
    const later = () => { clearTimeout(timer); timer = setTimeout(restart, 200); };
    $('#ro-lat').addEventListener('input', (e) => { RO.lat = +e.target.value; setF(); labels(); later(); });
    $('#ro-size').addEventListener('input', (e) => { RO.Lkm = +e.target.value; labels(); later(); });
    const setHemiUI = segButtons($('#ro-hemis'), [{ id: 'sur', name: 'Sur' }, { id: 'norte', name: 'Norte' }], (i) => { RO.south = i === 0; setHemiUI(i); restart(); });
    setHemiUI(0);
    $('#ro-reset').addEventListener('click', (e) => { restart(); if (e.detail > 0) e.currentTarget.blur(); });
    cv.addEventListener('pointerdown', (e) => {
      const R = RO.rect; if (!R) return;
      const b = cv.getBoundingClientRect(), dpr = cv.width / b.width;
      const x = ((e.clientX - b.left) * dpr - R.x) / R.w, y = 1 - ((e.clientY - b.top) * dpr - R.y) / R.h;
      if (x < 0.06 || x > 0.94 || y < 0.1 || y > 0.9) return;
      drop(x * RO.nx * RO.dxm, y * RO.ny * RO.dxm, true);
    });

    const needle = buildGauge($('#ro-gauge'), {
      aria: 'Escala del número de Rossby', min: 0.05, max: 20, map: Math.log10,
      zones: [
        { a: 0.05, b: 0.5, color: css('--calm') },
        { a: 0.5, b: 2, color: css('--warn') },
        { a: 2, b: 20, color: css('--hot') },
      ],
      ticks: [{ v: 0.1, label: '0.1' }, { v: 0.5, label: '0.5' }, { v: 1, label: '1' }, { v: 2, label: '2' }, { v: 10, label: '10' }],
      marks: [{ v: 0.16, label: 'gana el giro' }, { v: 1, label: 'empate' }, { v: 6.3, label: 'ganan las olas' }],
    });
    function setRegime(Ro) {
      let c, t, s;
      if (Ro < 0.5) { c = '--calm'; t = 'Gana el giro'; s = 'La mancha casi no se desparrama: queda atrapada y gira como un remolino.'; }
      else if (Ro < 2) { c = '--warn'; t = 'Empate'; s = 'Una parte de la mancha escapa en olas y la otra queda girando.'; }
      else { c = '--hot'; t = 'Ganan las olas'; s = 'La mancha se deshace en olas que se alejan, antes de notar que la Tierra gira.'; }
      if (RO.regime === t) return;
      RO.regime = t;
      const box = $('#ro-regime');
      box.querySelector('i').style.background = `var(${c})`; box.querySelector('span').textContent = t; $('#ro-regime-sub').textContent = s;
    }

    function step() {
      for (let s = 0; s < RO.stepsPerFrame; s++) sw.step();
      const k = RO.stepsPerFrame * sw.dt / RO.dxm;
      for (let i = 0; i < NP; i++) {
        const w = sw.vel(P.x[i], P.y[i]);
        P.x[i] += w[0] * k; P.y[i] += w[1] * k; P.life[i] -= 1 / 60;
        if (P.life[i] <= 0 || P.x[i] < 1 || P.x[i] > RO.nx - 1 || P.y[i] < 1 || P.y[i] > RO.ny - 1) spawn(i);
      }
      if (!sw.healthy()) { restart(); toast('#toast-ro', 'La simulación se reinició'); }
    }
    function draw() {
      const dpr = fitCanvas(cv);
      const W = cv.width, H = cv.height, { nx, ny } = RO;
      const R = letterbox(W, H, nx / ny); RO.rect = R;
      const eta = sw.eta, g = 0.5 / RO.A;
      for (let j = 0; j < ny; j++) {
        const row = (ny - 1 - j) * nx, r = j * nx;
        for (let i = 0; i < nx; i++) {
          const a = eta[r + i] * g * 2;                       // -1 … 1
          let t = 0.5 + 0.5 * Math.sign(a) * Math.pow(Math.min(1, Math.abs(a)), 0.55);
          pix[row + i] = LUT[(t * 255) | 0];
        }
      }
      offCtx.putImageData(img, 0, 0);
      const c = ctx;
      c.fillStyle = css('--deep'); c.fillRect(0, 0, W, H);
      c.imageSmoothingEnabled = true; c.imageSmoothingQuality = 'high';
      c.drawImage(off, R.x, R.y, R.w, R.h);
      const X = (x) => R.x + x / nx * R.w, Y = (y) => R.y + (1 - y / ny) * R.h;
      // particles with a tail that shows where they come from
      const tail = 12 * 3600 / RO.dxm; // cells travelled in 12 h at 1 m/s
      c.lineCap = 'round'; c.lineWidth = Math.max(1.4, 1.3 * dpr);
      c.strokeStyle = 'rgba(233, 242, 244, 0.5)'; c.fillStyle = 'rgba(233, 242, 244, 0.9)';
      c.beginPath();
      for (let i = 0; i < NP; i++) {
        const w = sw.vel(P.x[i], P.y[i]);
        c.moveTo(X(P.x[i]), Y(P.y[i])); c.lineTo(X(P.x[i] - w[0] * tail), Y(P.y[i] - w[1] * tail));
      }
      c.stroke();
      const sz = Math.max(1.8, 1.6 * dpr);
      for (let i = 0; i < NP; i++) c.fillRect(X(P.x[i]) - sz / 2, Y(P.y[i]) - sz / 2, sz, sz);
      // Rossby radius around the release point
      const fs = Math.round(clamp(R.h * 0.045, 12 * dpr, 20 * dpr));
      c.font = `700 ${fs}px ${css('--f-body')}`; c.textBaseline = 'middle';
      const Rd = rossbyRadius();
      if (isFinite(Rd) && Rd < 0.6 * ny * RO.dxm) {
        const cx = X(RO.x0 / RO.dxm), cy = Y(RO.y0 / RO.dxm), rr = Rd / RO.dxm / nx * R.w;
        c.setLineDash([6 * dpr, 6 * dpr]); c.strokeStyle = css('--buoy'); c.lineWidth = 2 * dpr;
        c.beginPath(); c.arc(cx, cy, rr, 0, 2 * Math.PI); c.stroke(); c.setLineDash([]);
        c.fillStyle = css('--buoy'); c.textAlign = 'center';
        c.fillText('radio de Rossby', cx, cy + rr + fs * 0.9 < R.y + R.h - fs ? cy + rr + fs * 0.9 : cy - rr - fs * 0.9);
      }
      // scale bar
      const barKm = 50, bw = barKm * 1000 / RO.dxm / nx * R.w, bx = R.x + R.w - bw - 16 * dpr, by = R.y + R.h - 18 * dpr;
      c.strokeStyle = 'rgba(233, 242, 244, 0.85)'; c.lineWidth = 2 * dpr;
      c.beginPath(); c.moveTo(bx, by - 5 * dpr); c.lineTo(bx, by); c.lineTo(bx + bw, by); c.lineTo(bx + bw, by - 5 * dpr); c.stroke();
      c.fillStyle = 'rgba(233, 242, 244, 0.9)'; c.textAlign = 'center';
      c.font = `700 ${Math.round(fs * 0.8)}px ${css('--f-body')}`;
      c.fillText(barKm + ' km', bx + bw / 2, by - fs * 0.75);
    }
    function readout() {
      const Ro = rossby();
      $('#ro-val').textContent = isFinite(Ro) ? (Ro >= 10 ? Ro.toFixed(0) : Ro.toFixed(2)) : '∞';
      needle(isFinite(Ro) ? Ro : 20); setRegime(isFinite(Ro) ? Ro : 99);
      const m = clamp(sw.volume(RO.x0, RO.y0, 2 * RO.Lkm * 1000) / RO.V0, 0, 1);
      RO.trap += (m - RO.trap) * 0.15;
      const pct = Math.round(RO.trap * 100);
      $('#ro-trap').style.width = pct + '%'; $('#ro-trap-out').textContent = pct + ' %';
      const days = (sw.t - (RO.tDrop || 0)) / 86400;
      const turn = RO.lat === 0 ? 'sin desvío' : RO.south ? 'desvío a la izquierda ↺' : 'desvío a la derecha ↻';
      $('#ro-time').innerHTML = `<span class="arrow">día ${days.toFixed(1)}</span> · ${turn}`;
    }

    restart();
    return {
      step, draw, readout, every: 3, reset: restart, lr: '#ro-lat', ud: '#ro-size', udStep: 3,
      key(k) { if (k === 'n' || k === 'N') { $('#ro-hemis').children[RO.south ? 1 : 0].click(); return true; } return false; },
      demo: [
        () => { setRange('#ro-lat', 42); setRange('#ro-size', 30); },
        () => { setRange('#ro-size', 6); },
        () => { setRange('#ro-lat', 4); setRange('#ro-size', 25); },
      ],
      _dbg: { RO, sw, rossby },
    };
  },
});
