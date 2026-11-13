/* Race view: cumulative Rg of several runs growing together, beside their moving chains. */
const RACE_MAX = 8;
// The wiki's tile colours fail the colourblind-separation check, so lines keep this validated set.
const RACE_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'];
const DEPTH_LEVELS = 6;
const BEAD_RADIUS_NM = 0.3;
const SPINS_PER_RACE = 1.5;
const LABEL_GAP_PX = 15;
const RUNNER_IMAGE_PX = 40;
const RUNNER_IMAGE_GAP_PX = 8; // clear of the 4 px head dot and its ring
const VIEW_EXTENT_RG = 2.6; // panels fit this many Rg; a rare fully stretched frame may touch the edge
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;
const GOAL_KEY = 'raceGoalNm';
// Rg bands shaded on the race chart (nm).
// Pictures (natural size in px) sit once, in a column right of the plot, centred on their band.
const RG_REGIONS = [
  {name: 'Brain', from: 5, to: Infinity, color: '#c2408f', picture: {href: '/images/neuron.png', w: 512, h: 982}},
  {name: 'Skin', from: 4, to: 5, color: '#c27a2c', picture: {href: '/images/skin.png', w: 254, h: 336}},
  {name: 'Bone', from: -Infinity, to: 4, color: '#5a6b85', picture: {href: '/images/Good_bone.png', w: 512, h: 1149}},
];
const REGION_PICTURE_COL_PX = 52; // both pictures stand upright, so a narrow column fits them
const REGION_PICTURE_MAX_H_PX = 96;

/* One region picture, as large as fits the column and the band, centred on (cx, cy). */
function regionPicture({href, w, h, sideways}, cx, cy, bandH) {
  const [shownW, shownH] = sideways ? [h, w] : [w, h];
  const scale = Math.min(REGION_PICTURE_COL_PX / shownW, Math.min(bandH, REGION_PICTURE_MAX_H_PX) / shownH);
  // Before turning, the image box is the natural shape; rotating about its centre lays it down.
  const iw = w * scale, ih = h * scale;
  const image = svg('image', {href, x: cx - iw / 2, y: cy - ih / 2, width: iw, height: ih, preserveAspectRatio: 'xMidYMid meet'});
  if (sideways) image.setAttribute('transform', `rotate(90 ${cx} ${cy})`);
  return image;
}

/* Scoring. Each run carries the target Rg its designer entered:
     Rg points        = 900 - 100 × (nm away from that run's target)
     variation points = by how many different listed blocks the sequence uses (server: rg.VARIATION_POINTS)
     tissue bonus     = 50 if the sequence has the adhesion domain for the tissue its target lies in
   Highest total wins. Runs made before targets use the race's target box; if any run still has
   no target, the largest Rg wins instead. The constants and tissueOf live in ui.html, shared
   with the run panel; RG_REGIONS bounds must match its TISSUE_BOUNDS. */
const goalOf = e => (e.ownGoal !== null ? e.ownGoal : race.goal);
const scoring = () => race && race.entrants.every(e => goalOf(e) !== null);
const goalNm = () => (scoring() ? true : null);
const tissuePts = e => (e.adhesionTissues.includes(tissueOf(goalOf(e))) ? TISSUE_BONUS : 0);
const rgPoints = (e, v) => RG_MAX_POINTS - Math.abs(v - goalOf(e)) * RG_POINTS_PER_NM;
const totalPoints = (e, v) => rgPoints(e, v) + e.variationPts + tissuePts(e);
const fmtPts = (e, v) => `${Math.round(totalPoints(e, v))} pts`;
// Negative when entrant a (at Rg va) places ahead of b: higher total, or larger Rg without targets.
const better = (a, va, b, vb) => (goalNm() === null ? vb - va : totalPoints(b, vb) - totalPoints(a, va));
function readGoal(input) {
  const v = parseFloat(input.value);
  return input.value.trim() !== '' && Number.isFinite(v) && v > 0 ? v : null;
}

let race = null;

async function startRace(ids) {
  stopRace();
  document.body.classList.add('racing');
  const host = $('race');
  host.hidden = false;
  host.replaceChildren(el('p', {className: 'meta', textContent: 'Loading trajectories…'}));
  window.scrollTo(0, 0);
  try {
    const entrants = await Promise.all(ids.map(loadEntrant));
    if (!document.body.classList.contains('racing')) return;
    buildRace(entrants);
  } catch (err) {
    const back = el('button', {type: 'button', textContent: 'Back to runs'});
    back.addEventListener('click', stopRace);
    host.replaceChildren(el('div', {className: 'note bad', textContent: `⚠ Couldn't load this race: ${err.message}`}), back);
  }
}

function stopRace() {
  if (race) cancelAnimationFrame(race.raf);
  race = null;
  document.body.classList.remove('racing');
  $('race').hidden = true;
  $('race').replaceChildren();
  drawChart(); // in case the window resized while its panel was hidden
}

async function loadEntrant(id, slot) {
  const [d, res] = await Promise.all([api(`/api/jobs/${id}`), fetch(`/api/jobs/${id}/coords`)]);
  if (!res.ok) throw new Error(`no trajectory for ${jobTitle(d)}`);
  const coords = new Float32Array(await res.arrayBuffer());
  const beads = d.n_residues, frames = d.series.time_ns.length;
  if (coords.length !== frames * beads * 3) throw new Error(`trajectory size mismatch for ${jobTitle(d)}`);
  let maxRadius = 0;
  for (let i = 0; i < coords.length; i += 3) {
    const r = coords[i] ** 2 + coords[i + 1] ** 2 + coords[i + 2] ** 2;
    if (r > maxRadius) maxRadius = r;
  }
  return {
    id, slot, name: jobTitle(d), simulatedNs: d.simulated_ns, hasImage: d.has_image, beads, frames, coords, maxRadius: Math.sqrt(maxRadius),
    t: d.series.time_ns, cum: d.series.cumulative_nm,
    skip: d.result.frames_discarded, frameNs: d.result.frame_ns,
    final: d.result.rg_mean, sem: d.result.rg_sem,
    ownGoal: d.target_rg ?? null, types: d.blocks.types.length, variationPts: d.blocks.variation_pts,
    adhesionTissues: d.blocks.adhesion_tissues,
    meta: `${d.team ? `${d.team} · ` : ''}${d.n_residues} residues · ${d.target_rg == null ? 'no target' : `target ${d.target_rg} nm (${tissueOf(d.target_rg).toLowerCase()})`}`,
  };
}

/* Repeat runs of one sequence need telling apart: add the run length, then a number. */
function distinguishNames(entrants) {
  const repeated = e => entrants.filter(o => o.name === e.name).length > 1;
  entrants.filter(repeated).forEach(e => { e.pending = `${e.name} (${fmtNs(e.simulatedNs)} ns)`; });
  entrants.forEach(e => { if (e.pending) e.name = e.pending; });
  const seen = {};
  entrants.filter(repeated).forEach(e => { e.pending = `${e.name} #${(seen[e.name] = (seen[e.name] || 0) + 1)}`; });
  entrants.forEach(e => { if (e.pending) e.name = e.pending; });
}

function buildRace(entrants) {
  entrants.forEach(e => { e.color = RACE_COLORS[e.slot]; });
  distinguishNames(entrants);

  const play = el('button', {type: 'button', className: 'primary race-play', textContent: 'Start race'});
  const restart = el('button', {type: 'button', textContent: 'Restart'});
  const back = el('button', {type: 'button', textContent: 'Back to runs'});
  const scrub = el('input', {type: 'range', min: 0, max: 1000, value: 0, className: 'race-scrub'});
  scrub.setAttribute('aria-label', 'Race position');
  const length = el('select', {className: 'race-length'},
    ...[15, 30, 60, 120].map(s => el('option', {value: s, textContent: `${s} s race`})));
  length.value = 30;
  length.setAttribute('aria-label', 'Race length');
  // Only runs made before targets existed need this box.
  const goal = el('input', {type: 'number', step: 'any', min: 0, className: 'race-goal', placeholder: 'Target for older runs'});
  goal.setAttribute('aria-label', 'Target radius of gyration in nm for runs that were made without one');
  goal.title = 'Some runs here were made without a target Rg. Enter one to score them; leave empty to race for the largest Rg.';
  goal.hidden = entrants.every(e => e.ownGoal !== null);
  try { goal.value = localStorage.getItem(GOAL_KEY) || ''; } catch {}
  const clock = el('span', {className: 'race-clock'});
  const banner = el('div', {className: 'race-banner', hidden: true});
  banner.setAttribute('aria-live', 'polite');
  const chart = el('div', {className: 'race-chart'});
  const cards = el('div', {className: 'race-cards'});

  for (const e of entrants) {
    e.canvas = el('canvas');
    e.rankEl = el('span', {className: 'race-rank', textContent: '–'});
    e.valueEl = el('span', {className: 'race-value', textContent: 'warming up'});
    const key = el('i', {className: 'race-key'});
    key.style.background = e.color;
    e.card = el('article', {className: 'race-card'},
      el('header', {}, e.rankEl, el('span', {className: 'race-name'}, key, e.name)),
      e.canvas,
      e.valueEl,
      el('p', {className: 'race-meta', textContent: e.meta}));
    cards.append(e.card);
  }

  const notes = [];
  const lengths = entrants.map(e => e.beads), ends = entrants.map(e => e.t[e.t.length - 1]);
  if (new Set(lengths).size > 1) notes.push(`Chain lengths differ (${Math.min(...lengths)}–${Math.max(...lengths)} residues), and longer chains have a larger Rg regardless of sequence.`);
  if (new Set(ends).size > 1) notes.push('Run lengths differ. A run that ends early holds its final value.');

  $('race').replaceChildren(
    el('div', {className: 'race-controls'}, play, restart, length, goal, scrub, clock, back),
    banner,
    el('div', {className: 'race-track'},
      el('section', {className: 'card race-chart-card'},
        el('h2', {textContent: 'Cumulative Rg (nm)'}),
        chart,
        el('p', {className: 'meta', textContent: 'Each line starts after that run\'s equilibration frames.'}),
        ...notes.map(text => el('p', {className: 'meta', textContent: text}))),
      el('section', {},
        cards,
        el('p', {className: 'meta', textContent: 'Chains are the real trajectory with the view slowly turning, all at the same scale. The ring has a radius equal to the cumulative Rg.'}))));

  race = {
    entrants, chart, banner, clock, play, scrub, length, goalInput: goal, goal: readGoal(goal),
    time: 0, tEnd: Math.max(...ends), playing: false, last: null, raf: 0,
    extent: Math.min(
      Math.max(...entrants.map(e => e.maxRadius)),
      VIEW_EXTENT_RG * Math.max(...entrants.map(e => e.final))),
  };

  play.addEventListener('click', () => setPlaying(!race.playing));
  restart.addEventListener('click', () => { race.time = 0; setPlaying(true); });
  back.addEventListener('click', stopRace);
  goal.addEventListener('input', () => {
    race.goal = readGoal(goal);
    try { localStorage.setItem(GOAL_KEY, goal.value); } catch {}
    hideWinner();
    layoutRace();
    drawRace();
  });
  scrub.addEventListener('input', () => {
    setPlaying(false);
    hideWinner();
    race.time = (scrub.value / 1000) * race.tEnd;
    drawRace();
  });

  layoutRace();
  drawRace();
}

function setPlaying(on) {
  if (on) hideWinner();
  if (on && race.time >= race.tEnd) race.time = 0;
  race.playing = on;
  race.last = null;
  race.play.textContent = on ? 'Pause' : race.time > 0 ? 'Resume' : 'Start race';
  cancelAnimationFrame(race.raf);
  if (on) race.raf = requestAnimationFrame(tick);
  else drawRace();
}

function tick(now) {
  if (!race || !race.playing) return;
  if (race.last !== null) {
    // Capped so returning to a background tab doesn't skip the finish.
    const dt = Math.min(now - race.last, 100) / 1000;
    race.time = Math.min(race.time + (dt / +race.length.value) * race.tEnd, race.tEnd);
  }
  race.last = now;
  drawRace();
  if (race.time >= race.tEnd) {
    race.playing = false;
    race.play.textContent = 'Race again';
    showWinner();
    return;
  }
  race.raf = requestAnimationFrame(tick);
}

/* Value lookups: t[i] is the end of frame i, so frame position = time / frameNs - 1. */
const framePos = (e, time) => Math.max(0, Math.min(e.frames - 1, time / e.frameNs - 1));
function cumAt(e, time) {
  if (time < e.t[e.skip]) return null;
  const f = framePos(e, time), i = Math.floor(f), j = Math.min(i + 1, e.frames - 1);
  return e.cum[i] + (e.cum[j] - e.cum[i]) * (f - i);
}

/* Layout: everything that depends on size */
function layoutRace() {
  const {entrants, chart, tEnd, extent} = race;
  const W = chart.clientWidth, H = Math.max(320, Math.min(540, W * 0.68));
  // The right margin leaves room for a picture riding ahead of a line that reaches the end.
  const imageRoom = entrants.some(e => e.hasImage) ? RUNNER_IMAGE_GAP_PX + RUNNER_IMAGE_PX : 0;
  const pictureRoom = RG_REGIONS.some(r => r.picture) ? REGION_PICTURE_COL_PX + 8 : 0;
  const m = {l: 44, r: 16 + imageRoom + pictureRoom, t: 12, b: 42};
  const finite = entrants.flatMap(e => e.cum.filter(v => v !== null));
  const goals = scoring() ? entrants.map(e => ({e, v: goalOf(e)})) : [];
  finite.push(...goals.map(g => g.v));
  let lo = Math.min(...finite), hi = Math.max(...finite);
  const pad = Math.max((hi - lo) * 0.06, 0.05);
  const {ticks, decimals} = niceTicks(lo - pad, hi + pad, 7);
  lo = Math.min(lo - pad, ticks[0]); hi = Math.max(hi + pad, ticks[ticks.length - 1]);
  const x = v => m.l + (v / tEnd) * (W - m.l - m.r);
  const y = v => m.t + (1 - (v - lo) / (hi - lo)) * (H - m.t - m.b);

  const root = svg('svg', {viewBox: `0 0 ${W} ${H}`, height: H, tabindex: 0, role: 'img',
    'aria-label': 'Cumulative radius of gyration of each run over simulated time'});
  for (const v of ticks) {
    root.append(svg('line', {x1: m.l, x2: W - m.r, y1: y(v), y2: y(v), stroke: 'var(--grid)'}));
    root.append(svg('text', {x: m.l - 8, y: y(v) + 4, 'text-anchor': 'end', 'font-size': 12, fill: 'var(--text-3)'}, v.toFixed(decimals)));
  }
  for (const v of niceTicks(0, tEnd, Math.max(3, Math.floor(W / 110))).ticks) {
    root.append(svg('text', {x: x(v), y: H - m.b + 18, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--text-3)'}, v));
  }
  root.append(svg('text', {x: (m.l + W - m.r) / 2, y: H - 8, 'text-anchor': 'middle', 'font-size': 12, fill: 'var(--text-2)'}, 'Simulated time (ns)'));

  // Pictures sit under the head dots, so a dot is never covered.
  const labels = svg('g'), lines = svg('g'), images = svg('g'), heads = svg('g');
  for (const e of entrants) {
    e.points = e.t.map((v, i) => (i < e.skip ? '' : `${x(v).toFixed(1)},${y(e.cum[i]).toFixed(1)}`));
    e.line = svg('polyline', {fill: 'none', stroke: e.color, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round'});
    e.head = svg('circle', {r: 4, fill: e.color, stroke: 'var(--surface)', 'stroke-width': 2, visibility: 'hidden'});
    e.label = svg('text', {'font-size': 12, 'font-weight': 600, fill: 'var(--text)', stroke: 'var(--surface)', 'stroke-width': 3, 'paint-order': 'stroke', visibility: 'hidden'});
    lines.append(e.line); heads.append(e.head); labels.append(e.label);
    if (e.hasImage) {
      e.image = svg('image', {href: `/api/jobs/${e.id}/image`, width: RUNNER_IMAGE_PX, height: RUNNER_IMAGE_PX,
        preserveAspectRatio: 'xMidYMid meet', visibility: 'hidden'});
      images.append(e.image);
    }
  }
  const cross = svg('line', {y1: m.t, y2: H - m.b, stroke: 'var(--text-3)', visibility: 'hidden'});
  // Tissue regions by Rg, drawn faintly behind everything else.
  const regions = svg('g');
  const yTop = m.t, yBottom = H - m.b;
  for (const r of RG_REGIONS) {
    const top = Math.max(yTop, Math.min(yBottom, y(Math.min(r.to, hi))));
    const bottom = Math.max(yTop, Math.min(yBottom, y(Math.max(r.from, lo))));
    if (bottom - top < 1) continue;
    regions.append(
      svg('rect', {x: m.l, y: top, width: W - m.l - m.r, height: bottom - top, fill: r.color, 'fill-opacity': 0.14}));
    if (r.picture) regions.append(regionPicture(r.picture, W - REGION_PICTURE_COL_PX / 2 - 4, (top + bottom) / 2, bottom - top));
    regions.append(
      svg('text', {x: W - m.r - 6, y: top + 15, 'text-anchor': 'end', 'font-size': 12, 'font-weight': 600,
        fill: r.color, stroke: 'var(--surface)', 'stroke-width': 3, 'paint-order': 'stroke'}, r.name));
  }
  const target = svg('g');
  // One dashed line per distinct target; in the run's colour when runs aim at different targets.
  const shared = new Set(goals.map(g => g.v)).size === 1;
  const drawn = new Set();
  for (const {e, v} of goals) {
    if (drawn.has(v)) continue;
    drawn.add(v);
    const gy = y(v), ink = shared ? 'var(--text)' : e.color;
    target.append(
      svg('line', {x1: m.l, x2: W - m.r, y1: gy, y2: gy, stroke: ink, 'stroke-width': 1.5, 'stroke-dasharray': '6 4'}),
      svg('text', {x: m.l + 6, y: gy - 6, 'font-size': 12, 'font-weight': 600, fill: ink, stroke: 'var(--surface)',
        'stroke-width': 3, 'paint-order': 'stroke'}, `🎯 ${shared ? 'Target' : shortName(e.name)} ${v.toFixed(2)} nm`));
  }
  root.prepend(regions); // under the grid lines and tick labels
  root.append(target, lines, cross, images, heads, labels);
  const tip = el('div', {className: 'race-tip'});
  chart.replaceChildren(root, tip);
  Object.assign(race, {x, y, W, H});

  root.addEventListener('pointermove', ev => {
    // Hover can't look past the current race time, so it never spoils the finish.
    const px = ev.clientX - root.getBoundingClientRect().left;
    const time = Math.max(0, Math.min(race.time, ((px - m.l) / (W - m.l - m.r)) * tEnd));
    const cx = x(time);
    cross.setAttribute('x1', cx); cross.setAttribute('x2', cx); cross.setAttribute('visibility', 'visible');
    const rows = entrants.map(e => ({e, v: cumAt(e, Math.min(time, e.t[e.t.length - 1]))}))
      .filter(r => r.v !== null).sort((a, b) => better(a.e, a.v, b.e, b.v));
    tip.replaceChildren(el('span', {className: 'race-tip-time', textContent: `${time.toFixed(2)} ns`}),
      ...rows.map(({e, v}) => {
        const k = el('i', {className: 'race-tip-key'});
        k.style.background = e.color;
        return el('div', {}, k, el('b', {textContent: v.toFixed(2)}), goalNm() === null ? '' : ` (${fmtPts(e, v)})`, ` ${e.name}`);
      }));
    if (!rows.length) tip.append(el('div', {textContent: 'Still in equilibration'}));
    tip.style.display = 'block';
    tip.style.left = (cx + 12 + tip.offsetWidth > W ? cx - tip.offsetWidth - 12 : cx + 12) + 'px';
    tip.style.top = (m.t + 4) + 'px';
  });
  root.addEventListener('pointerleave', () => {
    cross.setAttribute('visibility', 'hidden');
    tip.style.display = 'none';
  });

  const surface = getComputedStyle(document.documentElement).getPropertyValue('--surface').trim();
  for (const e of entrants) {
    const dpr = devicePixelRatio || 1, w = e.canvas.parentElement.clientWidth - 22, h = 240;
    e.canvas.width = w * dpr; e.canvas.height = h * dpr;
    e.canvas.style.width = '100%'; e.canvas.style.height = h + 'px';
    e.ctx = e.canvas.getContext('2d');
    e.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    e.view = {w, h, scale: (Math.min(w, h) / 2 - 4) / extent};
    e.px = new Float32Array(e.beads * 3);
    e.order = Array.from({length: e.beads}, (_, i) => i);
    e.sprites = makeSprites(e.color, surface, Math.max(2, BEAD_RADIUS_NM * e.view.scale), dpr);
    e.shades = e.sprites.map(s => s.shade);
  }
}

/* Colours */
const hexRgb = h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
const mix = (a, b, f) => a.map((v, i) => Math.round(v + (b[i] - v) * f));
const rgb = c => `rgb(${c[0]},${c[1]},${c[2]})`;

/* One pre-drawn shaded ball per depth level, so hundreds of beads stay cheap to draw. */
function makeSprites(color, surface, radius, dpr) {
  const base = hexRgb(color), bg = hexRgb(surface), white = [255, 255, 255], black = [0, 0, 0];
  const sprites = [];
  for (let level = 0; level < DEPTH_LEVELS; level++) {
    const depth = level / (DEPTH_LEVELS - 1); // 0 far, 1 near
    const c = mix(base, bg, (1 - depth) * 0.6);
    const r = radius * (0.8 + 0.4 * depth);
    const canvas = document.createElement('canvas');
    canvas.width = canvas.height = Math.ceil(2 * r * dpr) + 2;
    const g = canvas.getContext('2d'), mid = canvas.width / 2, R = r * dpr;
    const grad = g.createRadialGradient(mid - R * 0.35, mid - R * 0.35, R * 0.1, mid, mid, R);
    grad.addColorStop(0, rgb(mix(c, white, 0.45)));
    grad.addColorStop(0.55, rgb(c));
    grad.addColorStop(1, rgb(mix(c, black, 0.25)));
    g.fillStyle = grad;
    g.beginPath(); g.arc(mid, mid, R, 0, 7); g.fill();
    sprites.push({canvas, size: canvas.width / dpr, shade: rgb(c)});
  }
  return sprites;
}

/* One animation frame */
function drawChain(e, time, angle) {
  const {ctx, view, px, order, coords, beads, sprites} = e;
  const {w, h, scale} = view;
  const f = framePos(e, time), a = Math.floor(f), b = Math.min(a + 1, e.frames - 1), mixT = f - a;
  const offA = a * beads * 3, offB = b * beads * 3;
  const cos = Math.cos(angle), sin = Math.sin(angle), cx = w / 2, cy = h / 2;
  const extent = race.extent;

  for (let i = 0; i < beads; i++) {
    const j = i * 3;
    const X = coords[offA + j] + (coords[offB + j] - coords[offA + j]) * mixT;
    const Y = coords[offA + j + 1] + (coords[offB + j + 1] - coords[offA + j + 1]) * mixT;
    const Z = coords[offA + j + 2] + (coords[offB + j + 2] - coords[offA + j + 2]) * mixT;
    px[j] = cx + (X * cos + Z * sin) * scale;
    px[j + 1] = cy + Y * scale;
    px[j + 2] = -X * sin + Z * cos;
  }
  order.sort((p, q) => px[p * 3 + 2] - px[q * 3 + 2]);

  ctx.clearRect(0, 0, w, h);
  const cum = cumAt(e, time);
  if (cum !== null) {
    ctx.strokeStyle = getComputedStyle(document.documentElement).getPropertyValue('--text-3');
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, cum * scale, 0, 7); ctx.stroke();
  }
  ctx.lineCap = 'round';
  ctx.lineWidth = Math.max(1, BEAD_RADIUS_NM * scale * 0.7);
  // Far to near, so nearer beads cover farther ones.
  for (const i of order) {
    const j = i * 3;
    const depth = Math.max(0, Math.min(1, (px[j + 2] / extent + 1) / 2));
    const s = sprites[Math.round(depth * (DEPTH_LEVELS - 1))];
    if (i < beads - 1) {
      ctx.strokeStyle = s.shade;
      ctx.beginPath(); ctx.moveTo(px[j], px[j + 1]); ctx.lineTo(px[j + 3], px[j + 4]); ctx.stroke();
    }
    ctx.drawImage(s.canvas, px[j] - s.size / 2, px[j + 1] - s.size / 2, s.size, s.size);
  }
}

/* "Stein: S12-T12-E12" labels the line as "Stein"; the full name stays on the card. */
const SHORT_NAME_MAX = 16;
function shortName(name) {
  const short = name.split(':')[0].trim();
  return short.length > SHORT_NAME_MAX ? short.slice(0, SHORT_NAME_MAX - 1) + '…' : short;
}

const ordinal = n => n + (['th', 'st', 'nd', 'rd'][n] || 'th');
const setText = (node, text) => { if (node.textContent !== text) node.textContent = text; };

function drawRace() {
  if (!race) return;
  const {entrants, time, tEnd, x, y} = race;
  const angle = reduceMotion ? 0 : (time / tEnd) * SPINS_PER_RACE * 2 * Math.PI;

  const standing = [];
  for (const e of entrants) {
    const last = e.t[e.t.length - 1], at = Math.min(time, last), v = cumAt(e, at);
    e.current = v;
    if (v === null) {
      e.line.setAttribute('points', '');
      e.head.setAttribute('visibility', 'hidden');
      if (e.image) e.image.setAttribute('visibility', 'hidden');
    } else {
      const k = Math.floor(framePos(e, at));
      const hx = x(at).toFixed(1), hy = y(v).toFixed(1);
      e.line.setAttribute('points', `${e.points.slice(e.skip, k + 1).join(' ')} ${hx},${hy}`);
      e.head.setAttribute('cx', hx); e.head.setAttribute('cy', hy);
      e.head.setAttribute('visibility', 'visible');
      if (e.image) {
        const top = Math.max(0, Math.min(race.H - RUNNER_IMAGE_PX, y(v) - RUNNER_IMAGE_PX / 2));
        e.image.setAttribute('x', x(at) + RUNNER_IMAGE_GAP_PX);
        e.image.setAttribute('y', top);
        e.image.setAttribute('visibility', 'visible');
      }
      standing.push(e);
    }
    drawChain(e, at, angle);
  }

  standing.sort((a, b) => better(a, a.current, b, b.current));
  // When pictures overlap, the better-placed one goes on top.
  const orderKey = standing.map(e => e.id).join();
  if (orderKey !== race.orderKey) {
    race.orderKey = orderKey;
    for (const e of [...standing].reverse()) if (e.image) e.image.parentNode.append(e.image);
  }
  // Name the line heads from the leader down, skipping any that would collide.
  const placed = [];
  standing.forEach((e, i) => {
    const hx = x(Math.min(time, e.t[e.t.length - 1])), hy = y(e.current);
    if (placed.every(p => Math.abs(p - hy) >= LABEL_GAP_PX)) {
      placed.push(hy);
      const onRight = hx < race.W * 0.6;
      const ahead = e.hasImage ? RUNNER_IMAGE_GAP_PX + RUNNER_IMAGE_PX : 0;
      e.label.setAttribute('x', onRight ? hx + 9 + ahead : hx - 9);
      e.label.setAttribute('y', onRight ? hy + 4 : hy - 9);
      e.label.setAttribute('text-anchor', onRight ? 'start' : 'end');
      setText(e.label, `${shortName(e.name)} ${e.current.toFixed(2)}`);
      e.label.setAttribute('visibility', 'visible');
    } else {
      e.label.setAttribute('visibility', 'hidden');
    }
    setText(e.rankEl, ordinal(i + 1));
    setText(e.valueEl, goalNm() === null ? `${e.current.toFixed(2)} nm cumulative`
      : `${fmtPts(e, e.current)} · ${e.current.toFixed(2)} nm · ${bonusText(e)}`);
    e.card.classList.toggle('leading', i === 0);
  });
  for (const e of entrants) {
    if (e.current === null) {
      e.label.setAttribute('visibility', 'hidden');
      setText(e.rankEl, '–');
      setText(e.valueEl, 'warming up');
      e.card.classList.remove('leading');
    }
  }

  race.scrub.value = Math.round((time / tEnd) * 1000);
  setText(race.clock, `${time.toFixed(1)} / ${tEnd.toFixed(1)} ns`);
  showResult(time >= tEnd);
}

const rgText = e => `${e.final.toFixed(2)} ± ${e.sem.toFixed(2)} nm`;
function bonusText(e) {
  let text = `+${e.variationPts} variation (${e.types} block type${e.types === 1 ? '' : 's'})`;
  if (tissuePts(e)) text += `, +${TISSUE_BONUS} ${tissueOf(goalOf(e)).toLowerCase()} adhesion`;
  return text;
}
const scoreText = e => (goalNm() === null ? rgText(e)
  : `${fmtPts(e, e.final)} (Rg ${rgText(e)}, target ${goalOf(e).toFixed(2)} nm, ${bonusText(e)})`);
const finishOrder = () => [...race.entrants].sort((a, b) => better(a, a.final, b, b.final));
// The gap is within twice the combined standard error of Rg, converted to points when scoring.
const isPhotoFinish = ([first, second]) => Math.abs(better(first, first.final, second, second.final))
  <= 2 * Math.hypot(first.sem, second.sem) * (goalNm() === null ? 1 : RG_POINTS_PER_NM);

function showResult(finished) {
  const {banner} = race;
  if (!finished) { banner.hidden = true; return; }
  const order = finishOrder(), [first, second] = order, value = scoreText;
  const targets = goalNm() === null ? [] : [...new Set(race.entrants.map(goalOf))];
  const target = targets.length === 1 ? ` Target ${targets[0].toFixed(2)} nm.` : '';
  setText(banner, (!isPhotoFinish(order)
    ? `🏆 ${first.name} finishes on top with ${value(first)}, ahead of ${second.name} at ${value(second)}.`
    : `📸 Photo finish: ${first.name} (${value(first)}) leads ${second.name} (${value(second)}), but the gap is within the noise.`) + target);
  banner.hidden = false;
}

function showWinner() {
  const order = finishOrder(), [first] = order;
  if (!first.hasImage) return;
  hideWinner();
  const overlay = el('div', {className: 'race-winner', id: 'raceWinner'},
    el('img', {src: `/api/jobs/${first.id}/image`, alt: `Picture of ${first.name}`}),
    el('div', {}, `${isPhotoFinish(order) ? '📸 Leads a photo finish' : '🏆 Winner'}: ${first.name}, ${scoreText(first)}`,
      el('small', {textContent: 'Click anywhere or press Esc to close'})));
  overlay.addEventListener('click', hideWinner);
  $('race').append(overlay);
}
function hideWinner() {
  const overlay = $('raceWinner');
  if (overlay) overlay.remove();
}
document.addEventListener('keydown', e => { if (e.key === 'Escape') hideWinner(); });

let raceResizeTimer;
window.addEventListener('resize', () => {
  clearTimeout(raceResizeTimer);
  raceResizeTimer = setTimeout(() => { if (race) { layoutRace(); drawRace(); } }, 100);
});
