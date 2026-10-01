/*
 * colors_matrix.js
 *
 * Turns the long auto-generated list of colors (one <dl> per color, each
 * containing a <span class="colored-square">) into a compact swatch matrix:
 * a compact grid where hue runs left -> right and lightness runs top -> bottom
 * (vivid | muted | gray blocks).
 * Clicking a tile opens a popup with the color name, its Fortran value
 * (e.g. rgb_color(255,0,0)), hex, copy buttons and a link to the full entry.
 *
 * Nothing in the .f90 or .rst needs to change. The original detailed entries
 * are kept (collapsed in a <details> block), so every anchor such as
 * #f/sf_colors/red keeps working, and clicking a swatch jumps to its entry.
 *
 * Usage (conf.py):   html_js_files = ['colors_matrix.js']
 * (copy this file into your _static/ folder)
 */
(function () {
  'use strict';

  /* ---------- color classification (pure functions) ---------- */

  function hexToRgb(hex) {
    var h = hex.replace('#', '');
    if (h.length === 3) h = h.split('').map(function (c) { return c + c; }).join('');
    return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)];
  }

  function rgbToHsl(r, g, b) {
    r /= 255; g /= 255; b /= 255;
    var max = Math.max(r, g, b), min = Math.min(r, g, b), d = max - min;
    var l = (max + min) / 2, s = 0, h = 0;
    if (d > 0) {
      s = d / (1 - Math.abs(2 * l - 1));
      if (max === r) h = ((g - b) / d) % 6;
      else if (max === g) h = (b - r) / d + 2;
      else h = (r - g) / d + 4;
      h *= 60; if (h < 0) h += 360;
    }
    return { h: h, s: s, l: l, v: max, chroma: d * 255 };
  }

  /* sRGB (0-255) -> OKLCH: perceptual lightness L (0-1), chroma C, hue h (deg) */
  function rgbToOklch(r, g, b) {
    function lin(v) { v /= 255; return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); }
    var R = lin(r), G = lin(g), B = lin(b);
    var l = Math.cbrt(0.4122214708 * R + 0.5363325363 * G + 0.0514459929 * B);
    var m = Math.cbrt(0.2119034982 * R + 0.6806995451 * G + 0.1073969566 * B);
    var s = Math.cbrt(0.0883024619 * R + 0.2817188376 * G + 0.6299787005 * B);
    var L = 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s;
    var a = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
    var bb = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
    var h = Math.atan2(bb, a) * 180 / Math.PI; if (h < 0) h += 360;
    return { L: L, C: Math.sqrt(a * a + bb * bb), h: h };
  }

  var SMALL_LIST = 24;    // sections with at most this many colors are shown as a plain row
  var GRAY_C = 0.012;     // OKLCH chroma below this -> neutral gray block
  var MUTED_C = 0.09;     // below this (but above gray) -> muted/dusty block
  var HUE_START = 20;     // hue (deg) where the left-to-right rainbow starts (red)

  /* merge aliases (gray/grey) and compute hsl; returns [{hex,names,ids,vals,rgb,h,l,chroma}] */
  function prepare(items, noMerge) {
    var byHex = {};
    items.forEach(function (it, i) {
      var k = it.hex.toUpperCase() + (noMerge ? '#' + i : '');
      (byHex[k] = byHex[k] || { hex: it.hex.toUpperCase(), names: [], ids: [], vals: [] });
      byHex[k].names.push(it.name);
      byHex[k].ids.push(it.id);
      byHex[k].vals.push(it.value || '');
    });
    return Object.keys(byHex).map(function (k) {
      var e = byHex[k], rgb = hexToRgb(e.hex), c = rgbToHsl(rgb[0], rgb[1], rgb[2]);
      var o = rgbToOklch(rgb[0], rgb[1], rgb[2]);
      e.rgb = rgb; e.l = o.L; e.c = o.C; e.h = (o.h - HUE_START + 360) % 360;
      return e;
    });
  }

  /* Cut a list into columns of nearly equal height, longer columns first so
     any missing cells form a tidy step at the bottom right (not scattered
     holes). `rows` is a target; up to 3 extra rows are tried to find the height
     that leaves the fewest empty cells without adding columns. Each column is
     sorted light -> dark. */
  function chunk(list, rows) {
    var n = list.length, bestK = 0, bestHoles = Infinity;
    if (!n) return [];
    for (var r = rows; r <= rows + 3; r++) {
      var k = Math.ceil(n / r), holes = k * Math.ceil(n / k) - n;
      if (holes < bestHoles) { bestHoles = holes; bestK = k; }
    }
    var k = bestK, h = Math.ceil(n / k), long_ = n - k * (h - 1);   // columns of height h
    var cols = [], pos = 0;
    for (var i = 0; i < k; i++) {
      var len = i < long_ ? h : h - 1;
      cols.push(list.slice(pos, pos + len).sort(function (a, b) {
        return (b.l - a.l) || (a.h - b.h) || (a.names[0] < b.names[0] ? -1 : 1);
      }));
      pos += len;
    }
    return cols;
  }

  /* Matrix layout: hue runs left -> right, lightness top -> bottom.
     Three blocks, each sorted by hue and cut into equally tall columns:
       vivid colors | muted/dusty colors | grays
     Keeping the muted colors apart stops their unstable hue from making them
     look out of place among the vivid ones. Returns [vividCols, mutedCols, grayCols]. */
  function split(entries) {
    var out = [[], [], []];
    entries.forEach(function (e) { out[e.c < GRAY_C ? 2 : e.c < MUTED_C ? 1 : 0].push(e); });
    return out;
  }

  function layout(entries, rows) {
    var parts = split(entries);
    return [
      chunk(parts[0].sort(function (a, b) { return (a.h - b.h) || (b.l - a.l); }), rows),
      chunk(parts[1].sort(function (a, b) { return (a.h - b.h) || (b.l - a.l); }), rows),
      chunk(parts[2].sort(function (a, b) { return b.l - a.l; }), rows)
    ];
  }

  /* smallest number of rows such that all blocks fit in `avail` tile columns */
  function pickRows(entries, avail) {
    var sizes = split(entries).map(function (p) { return p.length; });
    var nb = sizes.filter(Boolean).length;
    var rows = Math.max(1, Math.ceil(entries.length / Math.max(1, avail)));
    function cols(r) { return sizes.reduce(function (n, k) { return n + Math.ceil(k / r); }, 0) + (nb - 1); }
    while (rows < entries.length && cols(rows) > avail) rows++;
    return rows;
  }

  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { prepare: prepare, layout: layout, pickRows: pickRows };
  }
  if (typeof document === 'undefined') return;

  /* ---------- DOM part ---------- */

  var CSS =
    '.cm-wrap{margin:1em 0}' +
    '.cm-info{min-height:2.4em;margin:0 0 .5em;font-family:monospace;font-size:.9em}' +
    '.cm-grid{display:flex;gap:14px;align-items:flex-start}' +
    '.cm-block{display:grid;grid-auto-flow:column;gap:2px}' +
    '.cm-block.cm-row{display:flex;flex-wrap:wrap;gap:2px}' +
    '.cm-sw{display:block;width:26px;height:26px;border-radius:3px;' +
    'box-shadow:inset 0 0 0 1px rgba(0,0,0,.1);transition:transform .08s}' +
    '.cm-sw:hover,.cm-sw:focus{transform:scale(1.2);position:relative;z-index:1;' +
    'box-shadow: 0 3px 12px rgba(0, 0, 0, 0.50);outline:none}' +
    '.cm-sw.cm-sel{box-shadow: 0 3px 12px rgba(0, 0, 0, 0.50);position:relative;z-index:1}' +
    '.cm-pop{position:absolute;z-index:1000;min-width:230px;max-width:340px;padding:12px 14px;' +
    'background:#fff;color:#222;border:1px solid #ccc;border-radius:6px;' +
    'box-shadow:0 4px 18px rgba(0,0,0,.25);font-size:.9em}' +
    '.cm-pop-top{display:flex;gap:10px;align-items:center;margin-bottom:8px}' +
    '.cm-pop-chip{flex:0 0 44px;height:44px;border-radius:4px;box-shadow:inset 0 0 0 1px rgba(0,0,0,.15)}' +
    '.cm-pop-name{font-weight:bold;word-break:break-word}' +
    '.cm-pop-row{font-family:monospace;margin:3px 0;word-break:break-all}' +
    '.cm-pop-row span{color:#777;font-family:sans-serif;font-size:.85em;display:inline-block;width:3.6em}' +
    '.cm-pop-actions{margin-top:8px;display:flex;gap:6px;flex-wrap:wrap;align-items:center}' +
    '.cm-pop-actions button{cursor:pointer;border:1px solid #bbb;background:#f6f6f6;border-radius:3px;' +
    'padding:2px 8px;font-size:.85em}' +
    '.cm-pop-actions button:hover{background:#e8e8e8}' +
    '.cm-pop-close{position:absolute;top:2px;right:8px;cursor:pointer;font-size:1.3em;color:#888;' +
    'border:0;background:none}' +
    '.cm-details{margin:1em 0}' +
    '.cm-details>summary{cursor:pointer;font-weight:bold;margin-bottom:.5em}';

  function addStyle() {
    var s = document.createElement('style');
    s.textContent = CSS;
    document.head.appendChild(s);
  }

  function readItem(dl) {
    var sq = dl.querySelector('.colored-square');
    var dt = dl.querySelector('dt.sig');
    var nm = dl.querySelector('.sig-name');
    if (!sq || !dt || !nm) return null;
    var m = /background-color:\s*(#[0-9a-fA-F]{3,6})/.exec(sq.getAttribute('style') || '');
    if (!m) return null;
    var value = '';
    Array.prototype.forEach.call(dl.querySelectorAll('dl.field-list > dt'), function (t) {
      if (/^Default/.test(t.textContent) && t.nextElementSibling) value = t.nextElementSibling.textContent.trim();
    });
    return { name: nm.textContent.trim(), hex: m[1], id: dt.id, value: value };
  }

  /* ---------- popup ---------- */

  var pop = null, popTarget = null;

  function closePop() {
    if (pop) { pop.remove(); pop = null; }
    if (popTarget) { popTarget.classList.remove('cm-sel'); popTarget = null; }
  }

  function row(label, text) {
    var d = document.createElement('div'); d.className = 'cm-pop-row';
    var l = document.createElement('span'); l.textContent = label;
    d.appendChild(l); d.appendChild(document.createTextNode(text));
    return d;
  }

  function copyText(txt, btn) {
    var done = function () {
      var old = btn.textContent; btn.textContent = 'Copied!';
      setTimeout(function () { btn.textContent = old; }, 1200);
    };
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(txt).then(done, done);
    } else {
      var ta = document.createElement('textarea'); ta.value = txt;
      document.body.appendChild(ta); ta.select();
      try { document.execCommand('copy'); } catch (e) {}
      ta.remove(); done();
    }
  }

  function showPop(tile, e) {
    closePop();
    popTarget = tile; tile.classList.add('cm-sel');
    var value = e.vals.filter(Boolean)[0] || ('rgb_color(' + e.rgb.join(',') + ')');

    pop = document.createElement('div'); pop.className = 'cm-pop';
    pop.setAttribute('role', 'dialog');
    var x = document.createElement('button'); x.className = 'cm-pop-close';
    x.textContent = '\u00d7'; x.setAttribute('aria-label', 'Close'); x.onclick = closePop;
    var top = document.createElement('div'); top.className = 'cm-pop-top';
    var chip = document.createElement('div'); chip.className = 'cm-pop-chip';
    chip.style.backgroundColor = e.hex;
    var nm = document.createElement('div'); nm.className = 'cm-pop-name';
    nm.textContent = e.names.join(' / ');
    top.appendChild(chip); top.appendChild(nm);
    pop.appendChild(x); pop.appendChild(top);
    pop.appendChild(row('value', value));
    pop.appendChild(row('hex', e.hex));
    pop.appendChild(row('rgb', e.rgb.join(', ')));

    var act = document.createElement('div'); act.className = 'cm-pop-actions';
    var b1 = document.createElement('button'); b1.textContent = 'Copy name';
    b1.onclick = function () { copyText(e.names[0], b1); };
    var b2 = document.createElement('button'); b2.textContent = 'Copy value';
    b2.onclick = function () { copyText(value, b2); };
    var b3 = document.createElement('button'); b3.textContent = 'Go to entry';
    b3.onclick = function () {
      var id = e.ids[0], t = document.getElementById(id);
      var d = t && t.closest('details.cm-details'); if (d) d.open = true;
      closePop(); if (t) { t.scrollIntoView(); history.replaceState(null, '', '#' + id); }
    };
    act.appendChild(b1); act.appendChild(b2); act.appendChild(b3);
    pop.appendChild(act);
    document.body.appendChild(pop);

    // position next to the tile, flipping to the left / clamping if needed
    var r = tile.getBoundingClientRect(), pw = pop.offsetWidth, ph = pop.offsetHeight;
    var left = r.right + 10;
    if (left + pw > window.innerWidth - 8) left = Math.max(8, r.left - pw - 10);
    var topPx = Math.min(Math.max(8, r.top - 10), Math.max(8, window.innerHeight - ph - 8));
    pop.style.left = (left + window.pageXOffset) + 'px';
    pop.style.top = (topPx + window.pageYOffset) + 'px';
  }

  document.addEventListener('keydown', function (ev) { if (ev.key === 'Escape') closePop(); });
  document.addEventListener('click', function (ev) {
    if (pop && !pop.contains(ev.target) && !(ev.target.closest && ev.target.closest('.cm-sw'))) closePop();
  });

  function buildSection(section) {
    var dls = Array.prototype.filter.call(section.children, function (el) {
      return el.tagName === 'DL' && el.classList.contains('variable');
    });
    var items = [], used = [];
    dls.forEach(function (dl) {
      var it = readItem(dl);
      if (it) { items.push(it); used.push(dl); }
    });
    if (!items.length) return;

    // the "a, b, c" list of links Sphinx puts before the entries
    var nodes = used.slice();
    var firstP = Array.prototype.filter.call(section.children, function (el) {
      return el.tagName === 'P' && el.querySelector('.xref.f-var');
    })[0];
    if (firstP) nodes.unshift(firstP);

    // short lists (the default colors) stay as plain separate tiles, in page order:
    // no hue sorting, no merging of colors that share a value
    var small = items.length <= SMALL_LIST;
    var entries = prepare(items, small);
    var wrap = document.createElement('div'); wrap.className = 'cm-wrap';
    var info = document.createElement('div'); info.className = 'cm-info';
    info.textContent = 'Hover a swatch to see its name; click it for details.';
    var grid = document.createElement('div'); grid.className = 'cm-grid';
    wrap.appendChild(info); wrap.appendChild(grid);

    function tile(e) {
      var a = document.createElement('a');
      a.className = 'cm-sw';
      a.href = '#' + e.ids[0];
      a.style.backgroundColor = e.hex;
      var label = e.names.join(' / ') + '  ' + e.hex + '  rgb_color(' + e.rgb.join(',') + ')';
      a.title = label;
      a.setAttribute('aria-label', label);
      var show = function () { info.textContent = label; };
      a.addEventListener('mouseenter', show);
      a.addEventListener('focus', show);
      a.addEventListener('click', function (ev) {
        ev.preventDefault();
        if (popTarget === a) closePop(); else showPop(a, e);
      });
      return a;
    }

    function render() {
      closePop();
      grid.textContent = '';
      if (small) {
        var row = document.createElement('div'); row.className = 'cm-block cm-row';
        entries.forEach(function (e) { row.appendChild(tile(e)); });
        grid.appendChild(row);
        return;
      }
      var avail = Math.max(4, Math.floor((wrap.clientWidth || 700) / 28) - 1);
      var lay = layout(entries, pickRows(entries, avail));
      lay.forEach(function (cols) {
        if (!cols.length) return;
        var block = document.createElement('div'); block.className = 'cm-block';
        var nrows = 0;
        cols.forEach(function (col, ci) {
          nrows = Math.max(nrows, col.length);
          col.forEach(function (e, ri) {
            var t = tile(e);
            t.style.gridColumn = String(ci + 1); t.style.gridRow = String(ri + 1);
            block.appendChild(t);
          });
        });
        block.style.gridTemplateRows = 'repeat(' + nrows + ', auto)';
        grid.appendChild(block);
      });
    }

    var details = document.createElement('details'); details.className = 'cm-details';
    var sum = document.createElement('summary');
    sum.textContent = 'Detailed entries (' + items.length + ')';
    details.appendChild(sum);

    section.insertBefore(wrap, nodes[0]);
    section.insertBefore(details, nodes[0]);
    nodes.forEach(function (n) { details.appendChild(n); });

    render();
    var timer, lastW = wrap.clientWidth;
    window.addEventListener('resize', function () {
      clearTimeout(timer);
      timer = setTimeout(function () {
        if (wrap.clientWidth !== lastW) { lastW = wrap.clientWidth; render(); }
      }, 150);
    });
  }

  // open the collapsed list if the URL (or a click) targets an entry in it
  function openForHash() {
    if (!location.hash) return;
    var t = document.getElementById(decodeURIComponent(location.hash.slice(1)));
    var d = t && t.closest && t.closest('details.cm-details');
    if (d) { d.open = true; t.scrollIntoView(); }
  }

  function init() {
    if (!document.querySelector('dl.variable .colored-square')) return;  // other pages: no-op
    addStyle();
    var sections = [];
    Array.prototype.forEach.call(document.querySelectorAll('.colored-square'), function (sq) {
      var s = sq.closest('section');
      if (s && sections.indexOf(s) < 0) sections.push(s);
    });
    sections.forEach(buildSection);
    window.addEventListener('hashchange', openForHash);
    openForHash();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
