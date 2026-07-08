/* config.js — palette, constants, formatting helpers (shared globals) */

const CELLS = ['HEK293T', 'HepG2', 'Jurkat'];

const CELL_COLOR = { HEK293T: '#4dbbd5', HepG2: '#f39b7f', Jurkat: '#00a087' };
const GLYCAN_COLOR = { og: '#f39b7f', ga: '#4dbbd5' };
const DIR = { up: '#e64b35', down: '#3c5488', ns: '#c3c9d0' };

const THRESH = { lfc: 0.5, adjp: 0.05 };   // overwritten from meta.json

const FONT = {
  body: 'Arial, Helvetica, sans-serif',
  mono: 'Arial, Helvetica, sans-serif',
};

/* ---- number formatting ---- */
function fmtFC(x) {
  if (x === null || x === undefined) return '—';
  return (x >= 0 ? '+' : '') + x.toFixed(2);
}
function fmtP(x) {
  if (x === null || x === undefined) return '—';
  if (x === 0) return '0';
  if (x < 0.001) {
    const e = x.toExponential(1); // "2.1e-9"
    return e.replace(/e([+-])(\d+)/, 'e$1$2');
  }
  return x.toFixed(3);
}
function fmtInt(x) {
  if (x === null || x === undefined) return '—';
  return x >= 1000 ? x.toLocaleString('en-US') : String(x);
}

/* significance for an lfc/adjp pair: 1 up, -1 down, 0 ns, null missing */
function sigOf(lfc, adjp) {
  if (lfc === null || adjp === null || lfc === undefined || adjp === undefined) return null;
  if (adjp < THRESH.adjp && Math.abs(lfc) > THRESH.lfc) return lfc > 0 ? 1 : -1;
  return 0;
}
function sigColor(s) { return s === 1 ? DIR.up : s === -1 ? DIR.down : DIR.ns; }

/* chip HTML for a logFC value coloured by its significance flag */
function chip(lfc, s) {
  if (lfc === null || lfc === undefined) return '<span class="chip na">—</span>';
  const cls = s === 1 ? 'up' : s === -1 ? 'down' : 'ns';
  return `<span class="chip ${cls}">${fmtFC(lfc)}</span>`;
}

/* curated featured spectra (img filled in when previews are committed) */
const SPECTRA = [
  { gene: 'PRDX6', site: 'Y89',  mode: 'EThcD', img: 'spectra/PRDX6_Y89.png',  note: 'Tyrosine O-GlcNAc · antioxidant peroxiredoxin (redox defense)' },
  { gene: 'SON',   site: 'Y270', mode: 'EThcD', img: 'spectra/SON_Y270.png',   note: 'Tyrosine O-GlcNAc · nuclear-speckle pre-mRNA splicing factor' },
  { gene: 'HPRT1', site: 'Y105', mode: 'EThcD', img: 'spectra/HPRT1_Y105.png', note: 'Tyrosine O-GlcNAc · purine-salvage enzyme (nucleotide metabolism)' },
  { gene: 'DDX17', site: 'Y580', mode: 'EThcD', img: 'spectra/DDX17_Y580.png', note: 'Tyrosine O-GlcNAc · DEAD-box RNA helicase (RNA processing)' },
  { gene: 'CPD',   site: 'T44',  mode: 'EThcD', img: 'spectra/CPD_T44.png',     note: 'O-GalNAc · trans-Golgi carboxypeptidase (secretory cargo processing)' },
  { gene: 'PTPRC', site: 'T139', mode: 'EThcD', img: 'spectra/PTPRC_T139.png',  note: 'O-GalNAc · CD45 leukocyte phosphatase (T-cell signaling)' },
];

/* external links */
const LINK = {
  uniprot: id => `https://www.uniprot.org/uniprotkb/${id}/entry`,
  glygen: id => `https://www.glygen.org/protein/${id}`,
  alphafold: id => `https://alphafold.ebi.ac.uk/entry/${id}`,
  afPdb: (id, v) => `https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v${v || 6}.pdb`,
};
