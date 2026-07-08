/* data.js — fetch + cache the JSON bundle */

const _cache = {};

async function loadJSON(name) {
  if (_cache[name]) return _cache[name];
  const res = await fetch('data/' + name);
  if (!res.ok) throw new Error(`Failed to load ${name} (${res.status})`);
  const json = await res.json();
  _cache[name] = json;
  return json;
}

const DATA = {
  meta: () => loadJSON('meta.json'),
  geneIndex: () => loadJSON('gene_index.json'),
  oglcnacProtein: () => loadJSON('oglcnac_protein.json'),
  oglcnacSite: () => loadJSON('oglcnac_site.json'),
  ogalnacProtein: () => loadJSON('ogalnac_protein.json'),
  ogalnacSite: () => loadJSON('ogalnac_site.json'),
  wpProtein: () => loadJSON('wp_protein.json'),
};

/* filter a long-format layer array to one cell type */
function byCell(rows, cell) { return rows.filter(r => r.cell === cell); }
