'use strict';
// ---------- utilidades comunes ----------
const $ = (s) => document.querySelector(s);
const cssCache = {};
const css = (name) => cssCache[name] || (cssCache[name] = getComputedStyle(document.documentElement).getPropertyValue(name).trim());
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const THIN = ' ';

// ---------- number formatting (SI: thin space for 5+ digits, decimal point) ----------
function groupInt(n) {
  const s = String(Math.round(n));
  if (s.length <= 4) return s;
  return s.replace(/\B(?=(\d{3})+(?!\d))/g, THIN);
}
function fmtNum(v) {
  if (!isFinite(v)) return '—';
  if (v >= 100) { const p = Math.pow(10, Math.floor(Math.log10(v)) - 2); return groupInt(Math.round(v / p) * p); }
  if (v >= 1) return String(Number(v.toPrecision(3)));
  if (v >= 0.001) return String(Number(v.toPrecision(2)));
  const e = Math.floor(Math.log10(v)); const m = v / Math.pow(10, e);
  const sup = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  return m.toFixed(1) + '×10' + String(e).split('').map((c) => sup[c]).join('');
}
function hexToRgb(h) { const n = parseInt(h.replace('#', ''), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
function mix(a, b, t) { return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]; }
const LITTLE = new Uint8Array(new Uint32Array([1]).buffer)[0] === 1;
function pack(c, a = 255) {
  const r = c[0] | 0, g = c[1] | 0, b = c[2] | 0;
  return LITTLE ? ((a << 24) | (b << 16) | (g << 8) | r) >>> 0 : ((r << 24) | (g << 16) | (b << 8) | a) >>> 0;
}

// ---------- canvas helpers ----------
function fitCanvas(cv) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(1, Math.round(cv.clientWidth * dpr)), h = Math.max(1, Math.round(cv.clientHeight * dpr));
  if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
  return dpr;
}
function letterbox(W, H, A) {
  if (W / H > A) { const w = H * A; return { x: (W - w) / 2, y: 0, w, h: H }; }
  const h = W / A; return { x: 0, y: (H - h) / 2, w: W, h };
}
function toast(id, msg) {
  const el = $(id); el.textContent = msg; el.hidden = false;
  clearTimeout(el._t); el._t = setTimeout(() => { el.hidden = true; }, 2600);
}
function segButtons(container, items, onPick) {
  container.innerHTML = '';
  items.forEach((it, i) => {
    const b = document.createElement('button');
    b.type = 'button'; b.id = container.id + '-' + it.id;
    b.innerHTML = it.name + (it.hint ? `<small>${it.hint}</small>` : '');
    b.addEventListener('click', (e) => { onPick(i); if (e.detail > 0) b.blur(); });
    container.appendChild(b);
  });
  return (sel) => [...container.children].forEach((b, i) => b.setAttribute('aria-pressed', String(i === sel)));
}

// ---------- gauge (SVG) ----------
const SVGNS = 'http://www.w3.org/2000/svg';
function el(tag, attrs, parent, text) {
  const e = document.createElementNS(SVGNS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (text !== undefined) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}
function buildGauge(box, opt) {
  const W = 600, X0 = 14, X1 = 586, BY = 44, BH = 14;
  const svg = el('svg', { viewBox: `0 0 ${W} 112`, role: 'img', 'aria-label': opt.aria });
  const pos = (v) => X0 + (opt.map(v) - opt.map(opt.min)) / (opt.map(opt.max) - opt.map(opt.min)) * (X1 - X0);
  opt.zones.forEach((z) => {
    el('rect', { x: pos(z.a), y: BY, width: Math.max(0, pos(z.b) - pos(z.a)), height: BH, fill: z.color, opacity: 0.9 }, svg);
  });
  opt.ticks.forEach((t) => {
    const x = pos(t.v);
    el('line', { x1: x, x2: x, y1: BY + BH, y2: BY + BH + 5, stroke: css('--muted'), 'stroke-width': 1.5 }, svg);
    el('text', { x, y: BY + BH + 20, 'text-anchor': t.anchor || 'middle', fill: css('--muted'), 'font-size': 13, 'font-family': css('--f-body') }, svg, t.label);
  });
  (opt.marks || []).forEach((m) => {
    const x = pos(m.v);
    el('line', { x1: x, x2: x, y1: 24, y2: BY, stroke: css('--fg'), 'stroke-width': 1, opacity: 0.5 }, svg);
    el('text', { x, y: 18, 'text-anchor': m.anchor || 'middle', fill: css('--fg'), 'font-size': 14, 'font-family': css('--f-body'), 'font-weight': 700 }, svg, m.label);
  });
  if (opt.band) {
    const a = pos(opt.band.a), b = pos(opt.band.b), y = BY + BH + 36;
    el('path', { d: `M${a} ${y - 6} V${y} H${b} V${y - 6}`, fill: 'none', stroke: css('--buoy'), 'stroke-width': 2 }, svg);
    el('text', { x: (a + b) / 2, y: y + 16, 'text-anchor': 'middle', fill: css('--buoy'), 'font-size': 13, 'font-family': css('--f-body'), 'font-weight': 700 }, svg, opt.band.label);
  }
  const needle = el('g', {}, svg);
  el('path', { d: `M0 ${BY - 4} L-7 ${BY - 15} H7 Z`, fill: css('--fg') }, needle);
  el('line', { x1: 0, x2: 0, y1: BY - 4, y2: BY + BH + 3, stroke: css('--fg'), 'stroke-width': 3 }, needle);
  box.innerHTML = ''; box.appendChild(svg);
  return (v) => { const x = clamp(pos(clamp(v, opt.min, opt.max)), X0, X1); needle.setAttribute('transform', `translate(${x.toFixed(1)} 0)`); };
}

// ==================================================================
//        Laboratorio modular: registro, carga de módulos, teclas, demo
// ==================================================================
// Para agregar un módulo: una línea aquí + un archivo en modulos/ que llame a Lab.register(id, { html, init }).
// init() devuelve { step(ms), draw(), readout(), every, reset(), lr, ud, key(k), demo: [fn, ...] }.
const LAB_VERSION = '2026-11-27';
const LAB_MODS = [
  { id: 're', slug: 'remolinos', name: 'Remolinos',       files: ['modulos/remolinos.js'],                  demo: 4 },
  { id: 'iw', slug: 'olas',      name: 'Olas escondidas', files: ['modulos/capas.js', 'modulos/olas.js'],   demo: 3 },
  { id: 'ri', slug: 'mezcla',    name: '¿Se mezcla?',     files: ['modulos/capas.js', 'modulos/mezcla.js'], demo: 3 },
  { id: 'ro', slug: 'rossby',    name: 'La Tierra gira',  files: ['modulos/rossby.js'],                     demo: 3 },
];
const Lab = {
  defs: {}, api: {}, active: null, paused: false,
  register(id, def) { this.defs[id] = def; },
};
function nudgeRange(id, d) { const r = $(id); r.value = clamp(+r.value + d, +r.min, +r.max); r.dispatchEvent(new Event('input')); }
function setRange(id, v) { const r = $(id); r.value = v; r.dispatchEvent(new Event('input')); }

(() => {
'use strict';
const main = $('main'), tabsEl = $('#lab-tabs');
const byId = (id) => LAB_MODS.find((m) => m.id === id);
LAB_MODS.forEach((m, i) => {
  const b = document.createElement('button');
  b.className = 'tab'; b.id = 'tab-' + m.id; b.setAttribute('role', 'tab'); b.setAttribute('aria-selected', 'false');
  b.innerHTML = `<kbd>${i + 1}</kbd>${m.name}`;
  b.addEventListener('click', () => setTab(m.id));
  tabsEl.appendChild(b);
});
$('#lab-keys').innerHTML = LAB_MODS.map((m, i) => `<kbd>${i + 1}</kbd>`).join(' ') + ' módulo';
const loadingEl = document.createElement('div'); loadingEl.className = 'loading'; loadingEl.hidden = true; loadingEl.textContent = 'Cargando módulo…';
main.appendChild(loadingEl);

const scripts = {};
function loadScript(src) {
  return scripts[src] || (scripts[src] = new Promise((ok, fail) => {
    const s = document.createElement('script'); s.src = src + '?v=' + LAB_VERSION;
    s.onload = ok; s.onerror = () => { delete scripts[src]; fail(new Error(src)); };
    document.head.appendChild(s);
  }));
}
function show() {
  LAB_MODS.forEach((m) => {
    $('#tab-' + m.id).setAttribute('aria-selected', String(m.id === Lab.active));
    const sec = $('#mod-' + m.id); if (sec) sec.hidden = m.id !== Lab.active;
  });
  loadingEl.hidden = !!Lab.api[Lab.active];
}
async function setTab(id) {
  const m = byId(id); if (!m) return;
  Lab.active = id; show();
  try { history.replaceState(null, '', '#' + m.slug); } catch (e) { /* file:// in some browsers */ }
  if (Lab.api[id] || m.loading) return;
  m.loading = true;
  try {
    for (const f of m.files) await loadScript(f);
    main.insertAdjacentHTML('beforeend', Lab.defs[id].html.trim());
    $('#mod-' + id).hidden = false;           // visible while it sets itself up (canvas sizes)
    Lab.api[id] = Lab.defs[id].init();
  } catch (e) {
    loadingEl.textContent = 'No se pudo cargar el módulo. Revisa que la carpeta "modulos" esté junto a esta página.';
    console.error(e);
  }
  m.loading = false; show();
}
Lab.setTab = setTab;

function fullscreen() {
  try {
    if (document.fullscreenElement) document.exitFullscreen().catch(() => {});
    else document.documentElement.requestFullscreen().catch(() => toast('#toast-' + Lab.active, 'Usa F11 para pantalla completa'));
  } catch (e) { /* not available */ }
}
$('#btn-fs').addEventListener('click', fullscreen);

document.addEventListener('keydown', (e) => {
  if (e.altKey || e.ctrlKey || e.metaKey) return;
  const onRange = document.activeElement && document.activeElement.type === 'range';
  const k = e.key, api = Lab.api[Lab.active];
  const n = '123456789'.indexOf(k);
  if (n >= 0 && n < LAB_MODS.length) { setTab(LAB_MODS[n].id); return; }
  if (k === 'f' || k === 'F') { fullscreen(); return; }
  if (!api) return;
  if (k === ' ') { e.preventDefault(); Lab.paused = !Lab.paused; toast('#toast-' + Lab.active, Lab.paused ? 'Pausa' : 'Continúa'); }
  else if (k === 'r' || k === 'R') api.reset();
  else if (api.key && api.key(k)) return;
  else if (onRange) return;
  else if (k === 'ArrowRight' || k === 'ArrowLeft') { e.preventDefault(); nudgeRange(api.lr, k === 'ArrowRight' ? 4 : -4); }
  else if (k === 'ArrowUp' || k === 'ArrowDown') {
    e.preventDefault();
    const up = k === 'ArrowUp';
    if (typeof api.ud === 'function') api.ud(up); else nudgeRange(api.ud, (up ? 1 : -1) * (api.udStep || 8));
  }
});

// demo (attract) mode after inactivity: 12 s per phase, module after module
const IDLE_MS = 45000, PHASE_S = 12;
const TOTAL = LAB_MODS.reduce((s, m) => s + m.demo, 0);
let lastInput = performance.now(), demo = false, demoStart = 0, demoPhase = -1;
function showDemoChip(on) { LAB_MODS.forEach((m) => { const c = $('#demo-' + m.id); if (c) c.hidden = !on || Lab.active !== m.id; }); }
function wake() { lastInput = performance.now(); if (demo) { demo = false; showDemoChip(false); } }
['pointermove', 'pointerdown', 'keydown', 'wheel'].forEach((ev) => window.addEventListener(ev, wake, { passive: true }));
function runDemo(now) {
  if (!demo && now - lastInput > IDLE_MS) { demo = true; demoStart = now; demoPhase = -1; Lab.paused = false; }
  if (!demo) return;
  const phase = Math.floor(((now - demoStart) / 1000) % (TOTAL * PHASE_S) / PHASE_S);
  if (phase === demoPhase) return;
  let p = phase, mod = LAB_MODS[0];
  for (const m of LAB_MODS) { if (p < m.demo) { mod = m; break; } p -= m.demo; }
  if (Lab.active !== mod.id) setTab(mod.id);
  const api = Lab.api[mod.id]; if (!api) return;          // still loading: try again next frame
  demoPhase = phase;
  if (api.demo && api.demo[p]) api.demo[p]();
  showDemoChip(true);
}
Lab._demo = (on) => { if (on) { lastInput = -1e9; } else wake(); };

let frameNo = 0;
function frame(now) {
  requestAnimationFrame(frame);
  runDemo(now);
  frameNo++;
  const api = Lab.api[Lab.active]; if (!api) return;
  if (!Lab.paused) api.step(9);
  api.draw();
  if (frameNo % (api.every || 3) === 0) api.readout();
}
const start = LAB_MODS.find((m) => '#' + m.slug === location.hash) || LAB_MODS[0];
setTab(start.id);
requestAnimationFrame(frame);
})();
