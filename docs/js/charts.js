/* charts.js — Plotly plots + NGL structure viewer */

const nlog10 = x => -Math.log10(x);

function baseLayout(over) {
  return Object.assign({
    paper_bgcolor: '#ffffff', plot_bgcolor: '#ffffff',
    font: { family: FONT.mono, size: 14, color: '#586472' },
    margin: { l: 50, r: 16, t: 16, b: 42 },
    hovermode: 'closest',
    hoverlabel: { bgcolor: '#191f26', bordercolor: '#191f26',
                  font: { family: FONT.mono, size: 14, color: '#fff' } },
    showlegend: false,
  }, over || {});
}
const PLAIN_CFG = { displayModeBar: false, responsive: true };
const FULL_CFG = { displaylogo: false, responsive: true,
  modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d'] };

/* ---------------- volcano ---------------- */
function drawVolcano(div, records, opts) {
  opts = opts || {};
  const valid = records.filter(r => r.lfc != null && r.adjP != null && r.adjP > 0);
  const groups = { '0': [], '-1': [], '1': [] };
  valid.forEach(r => { groups[String(sigOf(r.lfc, r.adjP) || 0)].push(r); });

  const order = ['0', '-1', '1'];
  const colorOf = { '0': DIR.ns, '-1': DIR.down, '1': DIR.up };
  const pointRecs = [];
  const traces = order.map(g => {
    const rs = groups[g];
    pointRecs.push(rs);
    return {
      type: 'scattergl', mode: 'markers',
      x: rs.map(r => r.lfc), y: rs.map(r => nlog10(r.adjP)),
      marker: { size: opts.mini ? 3.5 : 6, color: colorOf[g],
                opacity: g === '0' ? 0.45 : 0.8,
                line: { width: 0 } },
      customdata: rs.map(r => [r.gene, r.adjP, r.res ? `${r.res}${r.pos}` : '']),
      hovertemplate: opts.mini ? null :
        '<b>%{customdata[0]}</b> %{customdata[2]}<br>log₂FC %{x:.2f}' +
        '<br>adj.P %{customdata[1]:.2g}<extra></extra>',
      hoverinfo: opts.mini ? 'skip' : undefined,
    };
  });

  // highlight ring
  if (opts.highlight) {
    const h = valid.find(r => (r.site && r.site === opts.highlight) || r.id === opts.highlight);
    if (h && h.adjP > 0) {
      traces.push({
        type: 'scattergl', mode: 'markers',
        x: [h.lfc], y: [nlog10(h.adjP)],
        marker: { size: 15, color: 'rgba(0,0,0,0)', line: { color: '#191f26', width: 2.5 } },
        hoverinfo: 'skip',
      });
      pointRecs.push([h]);
    }
  }

  const t = nlog10(THRESH.adjp);
  const layout = baseLayout({
    xaxis: {
      title: opts.mini ? '' : { text: 'log₂ (Tuni / Ctrl)', font: { size: 15 } },
      zeroline: false, gridcolor: '#eef1f5',
      showticklabels: !opts.mini, ticks: opts.mini ? '' : 'outside',
    },
    yaxis: {
      title: opts.mini ? '' : { text: '−log₁₀ adj.P', font: { size: 15 } },
      zeroline: false, gridcolor: '#eef1f5', rangemode: 'tozero',
      showticklabels: !opts.mini,
    },
    margin: opts.mini ? { l: 24, r: 8, t: 6, b: 20 } : { l: 52, r: 16, t: 14, b: 44 },
    shapes: [
      { type: 'line', x0: THRESH.lfc, x1: THRESH.lfc, yref: 'paper', y0: 0, y1: 1,
        line: { color: '#c9d2db', width: 1, dash: 'dot' } },
      { type: 'line', x0: -THRESH.lfc, x1: -THRESH.lfc, yref: 'paper', y0: 0, y1: 1,
        line: { color: '#c9d2db', width: 1, dash: 'dot' } },
      { type: 'line', yref: 'y', y0: t, y1: t, xref: 'paper', x0: 0, x1: 1,
        line: { color: '#c9d2db', width: 1, dash: 'dot' } },
    ],
  });

  Plotly.react(div, traces, layout, opts.mini ? PLAIN_CFG : FULL_CFG);

  if (opts.onClick && !opts.mini) {
    div.removeAllListeners && div.removeAllListeners('plotly_click');
    div.on('plotly_click', ev => {
      const p = ev.points[0];
      const rec = pointRecs[p.curveNumber] && pointRecs[p.curveNumber][p.pointNumber];
      if (rec) opts.onClick(rec);
    });
  }
}

/* ---------------- TMT boxplot for one record ---------------- */
function drawBoxplot(div, rec, cell) {
  const v = rec.int;
  if (!v || v.some(x => x == null)) { Plotly.purge(div); return false; }
  const l2 = v.map(x => (x > 0 ? Math.log2(x) : null));
  const tuni = l2.slice(0, 3), ctrl = l2.slice(3, 6);
  const mk = (name, vals, color) => ({
    type: 'box', name, y: vals, x: vals.map(() => name),
    boxpoints: 'all', jitter: 0.6, pointpos: 0, marker: { color, size: 7 },
    line: { color }, fillcolor: 'rgba(0,0,0,0)', width: 0.5, hoverinfo: 'y',
  });
  const traces = [mk('Tuni', tuni, '#586472'), mk('Ctrl', ctrl, '#aeb7c0')];
  const layout = baseLayout({
    height: 220, margin: { l: 48, r: 14, t: 10, b: 26 },
    yaxis: { title: { text: 'log₂ intensity', font: { size: 14 } }, gridcolor: '#eef1f5', zeroline: false },
    xaxis: { showgrid: false },
  });
  Plotly.react(div, traces, layout, PLAIN_CFG);
  return true;
}

/* ---------------- structured vs IDR site plot ---------------- */
function drawStructStrip(div, sites, opts) {
  opts = opts || {};
  const valid = sites.filter(s => s.lfc != null);
  const groups = { struct: valid.filter(s => s.idr === 0), idr: valid.filter(s => s.idr === 1) };
  const xpos = { struct: 1, idr: 2 };
  const jit = i => (((i * 2654435761) % 1000) / 1000 - 0.5) * 0.34; // deterministic jitter

  const pointRecs = [];
  const traces = [];

  // box distributions
  ['struct', 'idr'].forEach(k => {
    const rs = groups[k];
    traces.push({
      type: 'box', x: rs.map(() => xpos[k]), y: rs.map(r => r.lfc),
      boxpoints: false, line: { color: '#c9d2db', width: 1.4 },
      fillcolor: 'rgba(200,207,214,0.12)', width: 0.5, hoverinfo: 'skip', showlegend: false,
    });
    pointRecs.push(null);
  });
  // clickable coloured points
  ['struct', 'idr'].forEach(k => {
    const rs = groups[k];
    pointRecs.push(rs);
    traces.push({
      type: 'scatter', mode: 'markers',
      x: rs.map((r, i) => xpos[k] + jit(i)), y: rs.map(r => r.lfc),
      marker: {
        size: rs.map(r => (opts.selected && r.site === opts.selected ? 13 : 7)),
        color: rs.map(r => sigColor(r.sig)),
        opacity: k === 'idr' ? 0.7 : 0.92,
        line: rs.map(r => (opts.selected && r.site === opts.selected
          ? { color: '#191f26', width: 2 } : { color: '#fff', width: 0.6 })),
      },
      customdata: rs.map(r => [r.gene, `${r.res}${r.pos}`, r.adjP]),
      hovertemplate: '<b>%{customdata[0]}</b> %{customdata[1]}<br>log₂FC %{y:.2f}' +
                     '<br>adj.P %{customdata[2]:.2g}<extra></extra>',
    });
  });

  const layout = baseLayout({
    height: 400,
    xaxis: { tickvals: [1, 2], ticktext: ['Structured', 'Disordered (IDR)'],
             range: [0.4, 2.6], showgrid: false, tickfont: { family: FONT.body, size: 15, color: '#191f26' } },
    yaxis: { title: { text: 'log₂ (Tuni / Ctrl)', font: { size: 15 } }, gridcolor: '#eef1f5', zeroline: false },
    shapes: [{ type: 'line', yref: 'y', y0: 0, y1: 0, xref: 'paper', x0: 0, x1: 1,
               line: { color: '#c9d2db', width: 1 } }],
  });
  Plotly.react(div, traces, layout, FULL_CFG);

  if (opts.onClick) {
    div.removeAllListeners && div.removeAllListeners('plotly_click');
    div.on('plotly_click', ev => {
      const p = ev.points[0];
      const rs = pointRecs[p.curveNumber];
      if (rs && rs[p.pointNumber]) opts.onClick(rs[p.pointNumber]);
    });
  }
}

/* ---------------- NGL 3D structure ---------------- */
let _nglStage = null;
function nglStage() {
  if (!_nglStage) {
    _nglStage = new NGL.Stage('ngl-stage', { backgroundColor: '#0c1116' });
    window.addEventListener('resize', () => _nglStage.handleResize());
  }
  return _nglStage;
}

function loadStructure(rec) {
  const stage = nglStage();
  const cap = document.getElementById('ngl-caption');
  cap.innerHTML = `Loading <b>${rec.gene} ${rec.res}${rec.pos}</b> (${rec.id})…`;
  stage.removeAllComponents();
  const dirColor = rec.sig === 1 ? '#e64b35' : rec.sig === -1 ? '#3c5488' : '#f0a020';

  const tryLoad = v => stage.loadFile(LINK.afPdb(rec.id, v), { ext: 'pdb' });
  tryLoad(6).catch(() => tryLoad(4)).then(comp => {
    if (!comp) return;
    comp.addRepresentation('cartoon', {
      colorScheme: 'bfactor',
      colorScale: ['#fd7e14', '#ffd43b', '#74c0fc', '#1864ab'],
      colorDomain: [50, 90], smoothSheet: true,
    });
    const sel = `${rec.pos} and .CA`;
    const selRes = `${rec.pos}`;
    comp.addRepresentation('ball+stick', { sele: selRes, color: dirColor, aspectRatio: 2, multipleBond: true });
    comp.addRepresentation('spacefill', { sele: sel, color: dirColor, opacity: 0.35, scale: 1.6 });
    comp.autoView();
    setTimeout(() => comp.autoView(selRes, 1200), 400);
    const idr = rec.idr === 1 ? 'disordered (IDR)' : 'structured';
    cap.innerHTML =
      `<b>${rec.gene} ${rec.res}${rec.pos}</b> · ${rec.id} · ${idr}` +
      ` · pLDDT ${rec.plddt ?? '—'} · log₂FC ${fmtFC(rec.lfc)} ` +
      `<a href="${LINK.alphafold(rec.id)}" target="_blank" rel="noopener">AlphaFold ↗</a>`;
  }).catch(() => {
    cap.innerHTML = `Could not load AlphaFold model for <b>${rec.id}</b>. ` +
      `<a href="${LINK.alphafold(rec.id)}" target="_blank" rel="noopener">Open on AlphaFold DB ↗</a>`;
  });
}
