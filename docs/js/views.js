/* views.js — view rendering, tables, search, gene passport */

const el = id => document.getElementById(id);

/* ---------------- cell tabs ---------------- */
function makeCellTabs(target, onChange) {
  const c = document.querySelector(`.cell-tabs[data-target="${target}"]`);
  if (!c) return;
  c.innerHTML = CELLS.map((cell, i) =>
    `<button data-cell="${cell}" class="${i === 0 ? 'on' : ''}">` +
    `<span class="cd" style="background:${CELL_COLOR[cell]}"></span>${cell}</button>`).join('');
  c.querySelectorAll('button').forEach(b => {
    b.onclick = () => {
      c.querySelectorAll('button').forEach(x => x.classList.remove('on'));
      b.classList.add('on');
      onChange(b.dataset.cell);
    };
  });
}

/* ---------------- Tabulator column presets ---------------- */
const colGene = {
  title: 'Gene', field: 'gene', width: 116, minWidth: 92,
  formatter: c => `<span class="tab-gene">${c.getValue() || '—'}</span>`,
};
const colFC = {
  title: 'log₂FC', field: 'lfc', hozAlign: 'center', width: 108, sorter: 'number',
  formatter: c => chip(c.getValue(), c.getRow().getData().sig),
};
const colAdjP = {
  title: 'adj.P', field: 'adjP', hozAlign: 'right', width: 100, sorter: 'number',
  formatter: c => `<span class="tab-mono">${fmtP(c.getValue())}</span>`,
};
const colId = {
  title: 'UniProt', field: 'id', width: 100,
  formatter: c => `<span class="tab-mono">${c.getValue() || ''}</span>`,
};
const colDesc = { title: 'Protein', field: 'desc', minWidth: 92, widthGrow: 3,
  formatter: c => c.getValue() || '' };
const colSite = {
  title: 'Site', field: 'site', width: 104, hozAlign: 'left',
  formatter: c => { const d = c.getRow().getData(); return `<span class="tab-mono">${d.res}${d.pos}</span>`; },
};

/* ---------------- generic protein / site browser ---------------- */
function makeBrowser(opts) {
  // opts: { volcanoId, tableId, detailId, getData, columns, idKey }
  let cell = CELLS[0], table = null, rows = [], selected = null, ready = false;

  async function render() {
    const all = await opts.getData();
    rows = byCell(all, cell);
    selected = null;
    drawVolcano(el(opts.volcanoId), rows, {
      cell, onClick: rec => select(rec[opts.idKey]),
    });
    buildTable();
    if (opts.detailId) el(opts.detailId).innerHTML =
      '<p class="detail-empty">Select a protein — on the plot or in the table — to see its TMT profile.</p>';
    ready = true;
  }

  function buildTable() {
    if (table) table.destroy();
    table = new Tabulator('#' + opts.tableId, {
      data: rows, height: '600px', layout: 'fitColumns', index: opts.idKey,
      columns: opts.columns, selectableRows: 1,
      initialSort: [{ column: 'adjP', dir: 'asc' }],
      placeholder: 'No entries',
    });
    // re-fit columns once the vertical scrollbar is present (avoids overflow)
    table.on('tableBuilt', () => table.redraw(true));
    table.on('rowClick', (e, row) => select(row.getData()[opts.idKey]));
  }

  function select(id) {
    selected = id;
    const rec = rows.find(r => r[opts.idKey] === id);
    if (!rec) return;
    drawVolcano(el(opts.volcanoId), rows, { cell, highlight: id, onClick: r => select(r[opts.idKey]) });
    if (table) { table.selectRow(id); table.scrollToRow(id, 'center', false).catch(() => {}); }
    showDetail(rec);
  }

  function showDetail(rec) {
    if (!opts.detailId) return;
    const box = el(opts.detailId);
    const site = rec.site ? ` ${rec.res}${rec.pos}` : '';
    box.innerHTML =
      `<div class="detail-card">
        <div class="dc-head">
          <span class="dc-gene">${rec.gene || rec.id}${site}</span>
          <span class="dc-name">${rec.name || rec.desc || ''}</span>
          <span class="dc-links">
            <a href="${LINK.uniprot(rec.id)}" target="_blank" rel="noopener">UniProt ↗</a>
            <a href="${LINK.glygen(rec.id)}" target="_blank" rel="noopener">GlyGen ↗</a>
          </span>
        </div>
        <div class="dc-sub" style="font-family:var(--f-mono);font-size:12.5px;color:var(--ink-2);margin-bottom:4px">
          ${cell} · log₂FC ${fmtFC(rec.lfc)} · adj.P ${fmtP(rec.adjP)}
        </div>
        <div class="dc-box" id="${opts.detailId}-box"></div>
      </div>`;
    const ok = drawBoxplot(el(opts.detailId + '-box'), rec, cell);
    if (!ok) el(opts.detailId + '-box').innerHTML =
      '<p class="detail-empty">No per-replicate intensities recorded for this entry.</p>';
  }

  return {
    setCell(c) { cell = c; },
    render,
    get cell() { return cell; },
    isReady: () => ready,
  };
}

/* ---------------- O-GlcNAc site browser (struct + NGL) ---------------- */
function makeSiteBrowser() {
  let cell = CELLS[0], table = null, rows = [], selected = null;

  async function render() {
    rows = byCell(await DATA.oglcnacSite(), cell);
    selected = null;
    drawStruct();
    buildTable();
    el('ngl-caption').textContent = 'Select a site — on the plot or in the table — to load its AlphaFold structure.';
  }
  function drawStruct() {
    drawStructStrip(el('ogs-struct'), rows, { selected, onClick: rec => select(rec.site) });
  }
  function buildTable() {
    if (table) table.destroy();
    table = new Tabulator('#ogs-table', {
      data: rows, height: '360px', layout: 'fitColumns', index: 'site',
      columns: [
        colGene, colSite, colFC, colAdjP,
        { title: 'pLDDT', field: 'plddt', hozAlign: 'right', width: 100, sorter: 'number',
          formatter: c => `<span class="tab-mono">${c.getValue() ?? '—'}</span>` },
        { title: 'Region', field: 'idr', hozAlign: 'center', width: 116,
          formatter: c => c.getValue() === 1
            ? '<span class="chip ns">IDR</span>'
            : '<span class="chip" style="background:#e7ecf1;color:#191f26">structured</span>' },
        { title: '2° structure', field: 'ss', minWidth: 140, widthGrow: 1, formatter: c => c.getValue() || '—' },
      ],
      initialSort: [{ column: 'adjP', dir: 'asc' }],
    });
    table.on('tableBuilt', () => table.redraw(true));
    table.on('rowClick', (e, row) => select(row.getData().site));
  }
  function select(site) {
    selected = site;
    const rec = rows.find(r => r.site === site);
    if (!rec) return;
    drawStruct();
    if (table) { table.selectRow(site); table.scrollToRow(site, 'center', false).catch(() => {}); }
    loadStructure(rec);
  }
  return { setCell(c) { cell = c; }, render };
}

/* ---------------- O-GalNAc controller (protein / site) ---------------- */
function makeOGalNAc() {
  let level = 'protein';
  const browser = makeBrowser({
    volcanoId: 'ga-volcano', tableId: 'ga-table', detailId: 'ga-detail',
    idKey: 'id', columns: [colGene, colId, colDesc, colFC, colAdjP],
    getData: () => level === 'protein' ? DATA.ogalnacProtein() : DATA.ogalnacSite(),
  });
  // site level uses a different key + columns; rebuild browser on toggle
  let current = browser;
  async function render() { await current.render(); }
  function setLevel(lv) {
    level = lv;
    current = makeBrowser({
      volcanoId: 'ga-volcano', tableId: 'ga-table', detailId: 'ga-detail',
      idKey: lv === 'protein' ? 'id' : 'site',
      columns: lv === 'protein'
        ? [colGene, colId, colDesc, colFC, colAdjP]
        : [colGene, colSite, colDesc, colFC, colAdjP],
      getData: () => lv === 'protein' ? DATA.ogalnacProtein() : DATA.ogalnacSite(),
    });
    current.setCell(cellNow);
    render();
  }
  let cellNow = CELLS[0];
  function setCell(c) { cellNow = c; current.setCell(c); }
  return { setCell, render, setLevel };
}

/* ---------------- overview ---------------- */
async function renderOverview() {
  const meta = await DATA.meta();
  const ogp = await DATA.oglcnacProtein();

  // stat strip
  const c = meta.counts;
  const total = k => CELLS.reduce((s, cell) => s + c[k][cell].n, 0);
  const stats = [
    { num: '1,109', lab: 'O-GlcNAc proteins', sub: 'union across three cell types' },
    { num: '670', lab: 'O-GlcNAc sites', sub: 'localized, with AlphaFold structural context' },
    { num: '485', lab: 'O-GalNAc proteins', sub: 'identified across three cell types' },
    { num: '3', lab: 'human cell types', sub: 'HEK293T · HepG2 · Jurkat' },
  ];
  el('stat-strip').innerHTML = stats.map(s =>
    `<div class="stat"><div class="num">${s.num}</div><div class="lab">${s.lab}</div><div class="sub">${s.sub}</div></div>`).join('');

  // hero mini-volcanoes
  CELLS.forEach(cell => {
    drawVolcano(el('hv-' + cell), byCell(ogp, cell), { mini: true, cell });
  });
  // panel click → open protein browser at that cell
  document.querySelectorAll('.volc-cell').forEach(fig => {
    fig.onclick = () => { location.hash = '#oglcnac-proteins'; };
  });
}

/* ---------------- spectra (curated / link-out) ---------------- */
function renderSpectra() {
  const g = el('spectra-gallery');
  if (g.dataset.done) return;
  g.dataset.done = '1';
  g.innerHTML =
    `<div class="ocard" style="grid-column:1/-1">
      <h3>Featured annotated spectra</h3>
      <p>Every localized site is backed by an annotated HCD / EThcD glycopeptide spectrum
      (MS1 isolation → oxonium + b/y sequence → c/z localization). A curated set is shown here;
      the full collection of raw and searched files is on
      <a href="https://www.ebi.ac.uk/pride/archive/projects/PXD073249" target="_blank" rel="noopener">PRIDE PXD073249</a>.</p>
    </div>`;
  SPECTRA.forEach(s => {
    const card = document.createElement('div');
    card.className = 'spectrum-card';
    card.innerHTML =
      `<div class="sc-head"><div class="sc-gene">${s.gene} ${s.site}</div>
        <div class="sc-meta">${s.note}</div></div>
      ${s.img ? `<a href="${s.img}" target="_blank" rel="noopener"><img src="${s.img}" alt="${s.gene} ${s.site} annotated ${s.mode} spectrum" loading="lazy"></a>`
              : `<div style="height:220px;display:flex;align-items:center;justify-content:center;color:var(--ink-3)">preview coming</div>`}
      <div class="sc-foot">${s.mode} annotation · open image for detail</div>`;
    g.appendChild(card);
  });
}

/* ---------------- downloads ---------------- */
function renderDownloads() {
  const b = el('downloads-body');
  if (b.dataset.done) return;
  b.dataset.done = '1';
  const tables = [
    ['S1', 'O-GlcNAc proteins — identification', '775 / 676 / 692 in HEK293T / HepG2 / Jurkat'],
    ['S2', 'O-GlcNAc proteins — abundance changes', 'log₂(Tuni/Ctrl), limma'],
    ['S3', 'Whole proteome — identification', ''],
    ['S4', 'Whole proteome — abundance changes', ''],
    ['S5', 'O-GlcNAc sites — identification', 'localized, site probability ≥ 0.75'],
    ['S6', 'O-GlcNAc sites — abundance changes', ''],
    ['S7', 'O-GalNAc proteins — identification', '485 proteins across three cell types'],
    ['S8', 'O-GalNAc proteins — abundance changes', ''],
    ['S9', 'O-GalNAc sites — identification', ''],
    ['S10', 'O-GalNAc sites — abundance changes', ''],
    ['S11', 'S-GlcNAc — cysteine sites', ''],
  ];
  b.innerHTML =
    `<div class="dl-group"><h3>Article &amp; raw data</h3><div class="dl-list">
      <a class="dl-item" href="https://www.ebi.ac.uk/pride/archive/projects/PXD073249" target="_blank" rel="noopener">
        <span class="dl-tag">PRIDE</span><span><span class="dl-name">PXD073249</span><br><span class="dl-desc">Raw &amp; searched MS files</span></span></a>
      <a class="dl-item" href="https://doi.org/10.1021/acs.analchem.6c00972" target="_blank" rel="noopener">
        <span class="dl-tag">DOI</span><span><span class="dl-name">Anal. Chem. 2026, 98, 15689</span><br><span class="dl-desc">Article + full Supporting Information</span></span></a>
    </div></div>
    <div class="dl-group"><h3>Supporting tables</h3><div class="dl-list">
      ${tables.map(([s, name, d]) =>
        `<a class="dl-item" href="tables/supporting_table_${s}.xlsx" download>
          <span class="dl-tag">${s}</span><span><span class="dl-name">${name}</span>${d ? `<br><span class="dl-desc">${d}</span>` : ''}</span></a>`).join('')}
    </div>
    <p class="pp-note">Excel tables as published with the article (Anal. Chem. 2026, 98, 15689–15699).</p></div>`;
}

/* ================= SEARCH + PASSPORT ================= */
let SEARCH_INDEX = [];
async function buildSearchIndex() {
  const gi = await DATA.geneIndex();
  SEARCH_INDEX = Object.keys(gi).map(gene => {
    const n = gi[gene];
    return {
      gene, node: n, id: (n.m && n.m[0]) || '',
      hay: (gene + ' ' + ((n.m && n.m[0]) || '')).toLowerCase(),
      tracks: { og: !!(n.og_p || n.og_s), ga: !!(n.ga_p || n.ga_s), wp: !!n.wp },
    };
  });
}
function layerSummary(n) {
  const bits = [];
  if (n.og_p) bits.push('O-GlcNAc');
  if (n.og_s) bits.push(n.og_s.length + ' site' + (n.og_s.length > 1 ? 's' : ''));
  if (n.ga_p || n.ga_s) bits.push('O-GalNAc');
  if (n.wp) bits.push('WP');
  return bits.join(' · ');
}
function runSearch(q) {
  q = q.trim().toLowerCase();
  const box = el('search-results');
  if (!q) { box.hidden = true; return; }
  const starts = [], has = [];
  for (const e of SEARCH_INDEX) {
    const i = e.hay.indexOf(q);
    if (i === 0 || e.gene.toLowerCase().startsWith(q)) starts.push(e);
    else if (i > 0) has.push(e);
    if (starts.length > 40) break;
  }
  const hits = starts.concat(has).slice(0, 14);
  if (!hits.length) { box.innerHTML = '<div class="sr-empty">No protein matches “' + q + '”.</div>'; box.hidden = false; return; }
  box.innerHTML = hits.map((e, i) =>
    `<div class="sr-item${i === 0 ? ' active' : ''}" data-gene="${e.gene}" role="option">
      <span class="sr-gene">${e.gene}</span>
      <span class="sr-tracks">
        ${e.tracks.og ? '<i style="background:var(--og)"></i>' : ''}
        ${e.tracks.ga ? '<i style="background:var(--ga)"></i>' : ''}
      </span>
      <span class="sr-meta">${layerSummary(e.node)}</span>
    </div>`).join('');
  box.hidden = false;
  box.querySelectorAll('.sr-item').forEach(it =>
    it.onclick = () => { openPassport(it.dataset.gene); box.hidden = true; el('gene-search').blur(); });
}

function ppCell(node, cell) {
  const x = node[cell];
  if (!x) return '<span class="chip na">—</span>';
  return chip(x.lfc, x.sig);
}
function ppSiteCell(sites, cell) {
  const inCell = sites.filter(s => s.cell === cell);
  if (!inCell.length) return '<span class="chip na">—</span>';
  const shown = inCell.slice(0, 4).map(s => {
    const sg = sigOf(s.lfc, s.adjP);
    const cls = sg === 1 ? 'up' : sg === -1 ? 'down' : 'ns';
    return `<span class="chip ${cls}" title="log₂FC ${fmtFC(s.lfc)}">${s.res}${s.pos}</span>`;
  }).join(' ');
  const more = inCell.length > 4 ? ` <span class="chip na">+${inCell.length - 4}</span>` : '';
  return `<div style="display:flex;flex-wrap:wrap;gap:3px;justify-content:center">${shown}${more}</div>`;
}
function ppRow(label, dotClass, cellsFn) {
  return `<tr><td class="rowlab">${dotClass ? `<span class="dot ${dotClass}"></span>` : ''}${label}</td>` +
    CELLS.map(c => `<td>${cellsFn(c)}</td>`).join('') + '</tr>';
}

async function openPassport(gene) {
  const gi = await DATA.geneIndex();
  const n = gi[gene];
  if (!n) return;
  const m = n.m || ['', '', ''];
  const id = m[0];
  const rows = [];
  if (n.og_p) rows.push(ppRow('O-GlcNAc protein', 'og', c => ppCell(n.og_p, c)));
  if (n.og_s) rows.push(ppRow('O-GlcNAc sites', 'og', c => ppSiteCell(n.og_s, c)));
  if (n.ga_p) rows.push(ppRow('O-GalNAc protein', 'ga', c => ppCell(n.ga_p, c)));
  if (n.ga_s) rows.push(ppRow('O-GalNAc sites', 'ga', c => ppSiteCell(n.ga_s, c)));
  if (n.wp) rows.push(ppRow('Whole proteome', '', c => ppCell(n.wp, c)));

  el('passport-body').innerHTML =
    `<div class="pp-gene">${gene}</div>
     <div class="pp-name">${m[1] || ''}</div>
     <div class="pp-id">${id ? `${id} · <a href="${LINK.uniprot(id)}" target="_blank" rel="noopener">UniProt ↗</a>
        · <a href="${LINK.alphafold(id)}" target="_blank" rel="noopener">AlphaFold ↗</a>
        · <a href="${LINK.glygen(id)}" target="_blank" rel="noopener">GlyGen ↗</a>` : ''}</div>
     <div class="pp-section-label">log₂ (Tuni / Ctrl) by layer and cell type</div>
     <table class="pp-matrix">
       <thead><tr><th class="rowlab"></th>
         ${CELLS.map(c => `<th><span class="pp-cellhead"><span class="cd" style="background:${CELL_COLOR[c]}"></span>${c}</span></th>`).join('')}
       </tr></thead>
       <tbody>${rows.join('')}</tbody>
     </table>
     <p class="pp-note"><span class="chip up">red</span> up · <span class="chip down">blue</span> down ·
       (|log₂FC| &gt; ${THRESH.lfc}, adj.P &lt; ${THRESH.adjp}). “—” = not detected in that layer/cell.</p>`;
  el('passport').hidden = false;
}
