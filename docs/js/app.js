/* app.js — controllers, router, init */

/* ---- browsers ---- */
const OGP = makeBrowser({
  volcanoId: 'ogp-volcano', tableId: 'ogp-table', detailId: 'ogp-detail', idKey: 'id',
  columns: [colGene, colId, colDesc, colFC, colAdjP], getData: DATA.oglcnacProtein,
});
const OGS = makeSiteBrowser();
const GA = makeOGalNAc();
const WP = makeBrowser({
  volcanoId: 'wp-volcano', tableId: 'wp-table', detailId: 'wp-detail', idKey: 'id',
  columns: [colGene, colId, colDesc, colFC, colAdjP], getData: DATA.wpProtein,
});

/* ---- whole-proteome lazy loader ---- */
let wpStarted = false;
async function loadWP() {
  if (wpStarted) return;
  wpStarted = true;
  try {
    await DATA.wpProtein();
    el('wp-loading').hidden = true;
    el('wp-split').hidden = false;
    await WP.render();
  } catch (e) {
    el('wp-loading').textContent = 'Could not load whole-proteome data.';
  }
}

/* ---- router ---- */
const PLOT_IDS = {
  'oglcnac-proteins': ['ogp-volcano', 'ogp-detail-box'],
  'oglcnac-sites': ['ogs-struct'],
  'ogalnac': ['ga-volcano', 'ga-detail-box'],
  'proteome': ['wp-volcano', 'wp-detail-box'],
};
const inited = { overview: true };

function route() {
  let h = (location.hash.replace('#', '') || 'overview');
  if (!el('view-' + h)) h = 'overview';
  document.querySelectorAll('.view').forEach(v => { v.hidden = v.id !== 'view-' + h; });
  document.querySelectorAll('.nav a').forEach(a => a.classList.toggle('active', a.dataset.view === h));
  el('sidebar').classList.add('collapsed');
  window.scrollTo(0, 0);

  if (!inited[h]) {
    inited[h] = true;
    if (h === 'oglcnac-proteins') OGP.render();
    else if (h === 'oglcnac-sites') OGS.render();
    else if (h === 'ogalnac') GA.render();
    else if (h === 'proteome') loadWP();
    else if (h === 'spectra') renderSpectra();
    else if (h === 'downloads') renderDownloads();
  } else {
    // re-fit plots that were sized while hidden
    (PLOT_IDS[h] || []).forEach(id => { const d = el(id); if (d && d.data) Plotly.Plots.resize(d); });
  }
}

/* ---- search ---- */
let indexPromise = null;
function ensureIndex() { return indexPromise || (indexPromise = buildSearchIndex()); }

function wireSearch() {
  const input = el('gene-search');
  const box = el('search-results');
  input.addEventListener('focus', ensureIndex);
  input.addEventListener('input', async () => { await ensureIndex(); runSearch(input.value); });
  input.addEventListener('keydown', e => {
    if (e.key === 'Enter') {
      const first = box.querySelector('.sr-item');
      if (first && !box.hidden) { openPassport(first.dataset.gene); box.hidden = true; input.blur(); }
    } else if (e.key === 'Escape') { box.hidden = true; input.blur(); }
  });
  document.addEventListener('click', e => {
    if (!e.target.closest('.search')) box.hidden = true;
  });
}

/* ---- passport modal ---- */
function wirePassport() {
  const close = () => { el('passport').hidden = true; };
  el('passport-close').onclick = close;
  el('passport').addEventListener('click', e => { if (e.target.id === 'passport') close(); });
  document.addEventListener('keydown', e => { if (e.key === 'Escape') close(); });
}

/* ---- draggable sidebar ---- */
const NAV_MIN = 200, NAV_MAX = 460, NAV_DEFAULT = 252;
function setNavWidth(px, persist) {
  const w = Math.min(NAV_MAX, Math.max(NAV_MIN, px));
  document.documentElement.style.setProperty('--nav-w', w + 'px');
  if (persist) { try { localStorage.setItem('og-nav-w', w + 'px'); } catch (e) {} }
  return w;
}
function wireResizer() {
  const rz = el('sidebar-resizer');
  if (!rz) return;
  let dragging = false, raf = false;
  rz.addEventListener('pointerdown', e => {
    dragging = true;
    try { rz.setPointerCapture(e.pointerId); } catch (_) {}
    document.body.classList.add('resizing');
    e.preventDefault();
  });
  rz.addEventListener('pointermove', e => {
    if (!dragging) return;
    setNavWidth(e.clientX, false);
    if (!raf) { raf = true; requestAnimationFrame(() => { raf = false; window.dispatchEvent(new Event('resize')); }); }
  });
  const end = e => {
    if (!dragging) return;
    dragging = false;
    try { rz.releasePointerCapture(e.pointerId); } catch (_) {}
    document.body.classList.remove('resizing');
    setNavWidth(parseInt(getComputedStyle(document.documentElement).getPropertyValue('--nav-w'), 10) || NAV_DEFAULT, true);
    window.dispatchEvent(new Event('resize'));
  };
  rz.addEventListener('pointerup', end);
  rz.addEventListener('pointercancel', end);
  rz.addEventListener('dblclick', () => { setNavWidth(NAV_DEFAULT, true); window.dispatchEvent(new Event('resize')); });
}

/* ---- init ---- */
async function init() {
  try {
    const meta = await DATA.meta();
    THRESH.lfc = meta.lfc_thresh;
    THRESH.adjp = meta.adjp_thresh;
  } catch (e) { /* keep defaults */ }

  // cell tabs
  makeCellTabs('ogp', c => { OGP.setCell(c); OGP.render(); });
  makeCellTabs('ogs', c => { OGS.setCell(c); OGS.render(); });
  makeCellTabs('ga', c => { GA.setCell(c); GA.render(); });
  makeCellTabs('wp', c => { WP.setCell(c); WP.render(); });

  // O-GalNAc level toggle
  const lt = el('ga-level-tabs');
  lt.querySelectorAll('button').forEach(b => b.onclick = () => {
    lt.querySelectorAll('button').forEach(x => x.classList.remove('on'));
    b.classList.add('on');
    GA.setLevel(b.dataset.level);
  });

  // nav toggle (mobile)
  el('nav-toggle').onclick = () => el('sidebar').classList.toggle('collapsed');

  wireSearch();
  wirePassport();
  wireResizer();

  await renderOverview();

  window.addEventListener('hashchange', route);
  route();
}

let _booted = false;
function boot() { if (_booted) return; _booted = true; init(); }
document.addEventListener('DOMContentLoaded', boot);
if (document.readyState !== 'loading') boot();
