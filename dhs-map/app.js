/* East Africa DHS Data - static version.
   Map + statistics logic. Data are read from ./data/ (made by export_static.R).
   The pure functions at the top are exported for testing in Node. */
(function (root) {
  'use strict';

  const Z975 = 1.959963984540054;
  const NO_DATA = { fill: '#BDBDBD', alpha: 0.25 };

  /* ---------------------------------------------------------------- helpers */
  const num = v => (v === null || v === undefined || v === '') ? NaN : Number(v);
  const isNum = v => typeof v === 'number' && isFinite(v);
  const esc = s => String(s).replace(/[&<>"']/g, c =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

  const fmtClass = (v, n) => isNum(v) ? 'class ' + v + ' of ' + n : 'no data';
  const formatP = p => p < 0.001 ? 'p < 0.001' : 'p = ' + p.toFixed(3);

  function classColour(lut, xv, yv) {
    if (!isNum(xv) || !isNum(yv)) return NO_DATA;
    return lut[xv + ' ' + yv] || NO_DATA;
  }

  /* ------------------------------------------------------------- statistics */
  function averageRanks(v) {
    const n = v.length, idx = v.map((_, i) => i).sort((a, b) => v[a] - v[b]);
    const r = new Array(n);
    let i = 0;
    while (i < n) {
      let j = i;
      while (j + 1 < n && v[idx[j + 1]] === v[idx[i]]) j++;
      const avg = (i + j) / 2 + 1;
      for (let k = i; k <= j; k++) r[idx[k]] = avg;
      i = j + 1;
    }
    return r;
  }

  function pearson(a, b) {
    const n = a.length;
    let ma = 0, mb = 0;
    for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
    ma /= n; mb /= n;
    let sab = 0, saa = 0, sbb = 0;
    for (let i = 0; i < n; i++) {
      const da = a[i] - ma, db = b[i] - mb;
      sab += da * db; saa += da * da; sbb += db * db;
    }
    if (saa === 0 || sbb === 0) return NaN;
    return sab / Math.sqrt(saa * sbb);
  }

  // log-gamma (Lanczos) and regularized incomplete beta (continued fraction)
  function lgamma(x) {
    const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
                 -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5, ser = 1.000000000190015;
    tmp -= (x + 0.5) * Math.log(tmp);
    for (let j = 0; j < 6; j++) { y += 1; ser += cof[j] / y; }
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  function betacf(a, b, x) {
    const MAXIT = 500, EPS = 3e-15, FPMIN = 1e-300;
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }

  function ibeta(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) +
                        a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2)
      ? bt * betacf(a, b, x) / a
      : 1 - bt * betacf(b, a, 1 - x) / b;
  }

  // Spearman rho (average ranks for ties); two-sided p from the t approximation
  // (same as R's cor.test(method="spearman", exact=FALSE)); approximate 95% CI by
  // Bonett & Wright (2000): Fisher z with SE = sqrt((1 + rho^2 / 2) / (n - 3)).
  function spearmanSummary(xv, yv) {
    const xs = [], ys = [];
    for (let i = 0; i < xv.length; i++) {
      if (isNum(xv[i]) && isNum(yv[i])) { xs.push(xv[i]); ys.push(yv[i]); }
    }
    const n = xs.length;
    if (n < 5) return null;
    const rho = pearson(averageRanks(xs), averageRanks(ys));
    if (!isFinite(rho)) return null;
    let p;
    if (Math.abs(rho) >= 1) p = 0;
    else {
      const df = n - 2, t = rho * Math.sqrt(df / (1 - rho * rho));
      p = ibeta(df / (df + t * t), df / 2, 0.5);
    }
    const se = Math.sqrt((1 + rho * rho / 2) / (n - 3)), z = Math.atanh(rho);
    return { n: n, rho: rho, p: p, ci: [Math.tanh(z - Z975 * se), Math.tanh(z + Z975 * se)] };
  }

  const api = { averageRanks, pearson, ibeta, spearmanSummary, classColour, fmtClass, formatP };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  if (typeof document === 'undefined') return;

  /* ------------------------------------------------------------------- page */
  function init() {
    const status = document.getElementById('status');
    const fail = err => {
      status.hidden = false;
      status.textContent = 'Could not load the data (' + err.message + '). If you opened ' +
        'index.html straight from disk, serve the folder over http instead (see README.md).';
    };
    const getJSON = url => fetch(url).then(r => {
      if (!r.ok) throw new Error(url + ': HTTP ' + r.status);
      return r.json();
    });
    getJSON('data/meta.json')
      .then(meta => Promise.all(meta.periods.map((_, i) => getJSON('data/period_' + (i + 1) + '.geojson')))
        .then(geos => start(meta, geos)))
      .catch(fail);
  }

  function start(meta, geos) {
    const $ = id => document.getElementById(id);
    const nP = meta.periods.length, nClass = meta.nClass;
    const state = {
      x: meta.xChoices[0].value, y: null, p: 0, pal: meta.defaultPalette,
      selected: null
    };
    state.y = meta.yChoices[state.x][0].value;

    // id -> properties, per period (for the clicked-region table)
    const index = geos.map(g => {
      const m = new Map();
      g.features.forEach(f => m.set(f.properties.rid, f.properties));
      return m;
    });

    /* ---- map ---- */
    const map = L.map('map').setView([0, 38], 6);
    L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
      maxZoom: 19,
      attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    }).addTo(map);

    const lut = () => meta.palettes[state.pal];
    const val = (props, v) => num(props[v]);

    function styleFor(feature) {
      const pr = feature.properties;
      const c = classColour(lut(), val(pr, state.x), val(pr, state.y));
      const sel = pr.rid === state.selected;
      return { fillColor: c.fill, fillOpacity: 0.6 * c.alpha,
               color: '#000', weight: sel ? 3 : 1, opacity: sel ? 0.9 : 0.3 };
    }

    const layers = [];
    let current = null;
    function getLayer(i) {
      if (!layers[i]) {
        layers[i] = L.geoJSON(geos[i], {
          style: styleFor,
          onEachFeature: (feature, layer) => {
            layer.bindTooltip(feature.properties.label, { sticky: true });
            layer.on('mouseover', e => { e.target.setStyle({ weight: 2, opacity: 0.8 }); e.target.bringToFront(); });
            layer.on('mouseout', e => e.target.setStyle(styleFor(feature)));
            layer.on('click', () => { state.selected = feature.properties.rid; update(); });
          }
        });
      }
      return layers[i];
    }

    /* ---- controls ---- */
    const fillSelect = (el, opts, selected) => {
      el.innerHTML = opts.map(o => '<option value="' + esc(o.value) + '"' +
        (o.value === selected ? ' selected' : '') + '>' + esc(o.label) + '</option>').join('');
    };
    fillSelect($('x_var'), meta.xChoices, state.x);
    fillSelect($('viridis_pal'), meta.paletteChoices, state.pal);
    const setYChoices = () => {
      const ch = meta.yChoices[state.x];
      if (!ch.some(o => o.value === state.y)) state.y = ch[0].value;
      fillSelect($('y_var'), ch, state.y);
    };
    setYChoices();

    $('period').max = nP - 1;
    $('period_labels').innerHTML = meta.periods.map(p => '<span>' + esc(p) + '</span>').join('');

    $('x_var').onchange = e => { state.x = e.target.value; setYChoices(); update(); };
    $('y_var').onchange = e => { state.y = e.target.value; update(); };
    $('viridis_pal').onchange = e => { state.pal = e.target.value; update(); };
    $('period').oninput = e => { stopPlay(); state.p = Number(e.target.value); update(); };
    $('county_close').onclick = () => { state.selected = null; update(); };

    // Play steps through the periods once (no looping)
    let timer = null;
    function stopPlay() {
      if (timer) { clearInterval(timer); timer = null; }
      $('play').textContent = '▶ Play';
    }
    $('play').onclick = () => {
      if (timer) { stopPlay(); return; }
      if (state.p >= nP - 1) state.p = 0;
      update();
      $('play').textContent = '❚❚ Pause';
      timer = setInterval(() => {
        if (state.p >= nP - 1) { stopPlay(); return; }
        state.p++; update();
      }, 2500);
    };

    /* ---- panels ---- */
    function renderLegend() {
      const L_ = lut();
      let html = '';
      for (let y = nClass; y >= 1; y--) {
        for (let x = 1; x <= nClass; x++) {
          const c = L_[x + ' ' + y] || NO_DATA;
          html += '<div class="cell" style="background:' + c.fill + ';opacity:' + c.alpha + '"></div>';
        }
      }
      const grid = $('legend_grid');
      grid.style.gridTemplateColumns = 'repeat(' + nClass + ', 1fr)';
      grid.innerHTML = html;
      $('legend_x').innerHTML = esc(meta.axisLabels[state.x]).replace(/\n/g, '<br>') + ' &rarr;';
      $('legend_y').innerHTML = esc(meta.axisLabels[state.y]).replace(/\n/g, '<br>') + ' &rarr;';
    }

    function renderAssoc() {
      const fs = geos[state.p].features;
      const a = spearmanSummary(fs.map(f => val(f.properties, state.x)),
                                fs.map(f => val(f.properties, state.y)));
      $('assoc').innerHTML = a === null
        ? '<h5>Across-region association</h5><p>Too few regions with data (or no variation) to compute a correlation.</p>'
        : '<h5>Across-region association</h5><p>Spearman &rho; = ' + a.rho.toFixed(2) +
          ' (approx. 95% CI ' + a.ci[0].toFixed(2) + ' to ' + a.ci[1].toFixed(2) + '), ' +
          formatP(a.p) + ', n = ' + a.n + ' regions</p>' +
          '<p class="note">Ecological association: it describes how region-level classes co-vary, ' +
          'and does not show that education or stigma changes outcomes for individuals. ' +
          'Regions from all countries are pooled, so differences between countries contribute to ' +
          'the correlation. Classes are coarse (ties are common), regions count equally regardless ' +
          'of survey sample size, and p-values are approximate and not adjusted for the many ' +
          'variable/period combinations you can browse.</p>';
    }

    function renderCounty() {
      const panel = $('county');
      if (state.selected === null) { panel.hidden = true; return; }
      let label = state.selected;
      for (let i = 0; i < nP && index[i]; i++) {
        const pr = index[i].get(state.selected);
        if (pr) { label = pr.label; break; }
      }
      const rows = meta.periods.map((p, i) => {
        const pr = index[i].get(state.selected) || {};
        const xv = val(pr, state.x), yv = val(pr, state.y);
        const c = classColour(lut(), xv, yv);
        return '<tr' + (i === state.p ? ' class="now"' : '') + '><td>' + esc(p) + '</td><td>' +
          fmtClass(xv, nClass) + '</td><td>' + fmtClass(yv, nClass) + '</td><td><span class="sw" ' +
          'style="background:' + c.fill + ';opacity:' + c.alpha + '"></span></td></tr>';
      }).join('');
      $('county_title').textContent = label;
      $('county_table').innerHTML = '<thead><tr><th>Period</th><th>' + esc(meta.varLabels[state.x]) +
        '</th><th>' + esc(meta.varLabels[state.y]) + '</th><th></th></tr></thead><tbody>' + rows + '</tbody>';
      panel.hidden = false;
    }

    function update() {
      const l = getLayer(state.p);
      if (current && current !== l) map.removeLayer(current);
      if (!map.hasLayer(l)) l.addTo(map);
      current = l;
      l.setStyle(styleFor);                       // recolour in place, no geometry re-sent

      $('period').value = state.p;
      Array.prototype.forEach.call($('period_labels').children, (s, i) =>
        s.classList.toggle('on', i === state.p));
      $('relationship').innerHTML = esc(meta.relLabels[state.x]) + ' &#10132; ' + esc(meta.relLabels[state.y]);
      renderLegend(); renderAssoc(); renderCounty();
    }

    update();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})(typeof window !== 'undefined' ? window : globalThis);
