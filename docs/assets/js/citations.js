// SPDX-FileCopyrightText: 2026 DOPAMINE contributors
// SPDX-License-Identifier: AGPL-3.0-only
// Renders pandoc-style [@key] / [@key1; @key2] citations client-side into
// numbered links plus a References section, since GitHub Pages/Jekyll has
// no citeproc. Source markdown keeps the [@key] syntax so build.sh/Pandoc
// (which does run citeproc) still works unchanged.
// Runs after MathJax's own typeset pass finishes: both scripts walk/mutate
// the same text nodes on DOMContentLoaded, and racing MathJax's scan makes
// it silently skip a large fraction of the page's equations.
(function () {
  function whenMathJaxReady(fn) {
    if (window.MathJax && window.MathJax.startup && window.MathJax.startup.promise) {
      window.MathJax.startup.promise.then(fn);
    } else {
      fn();
    }
  }

  document.addEventListener("DOMContentLoaded", function () {
    whenMathJaxReady(renderCitations);
  });

  function renderCitations() {
    var bib = window.DOPAMINE_REFERENCES || {};
    var container =
      document.querySelector("#main-content") ||
      document.querySelector(".main-content") ||
      document.body;

    var order = [];
    var numberOf = {};
    var CITE_RE = /\[(@[A-Za-z0-9_]+(?:;\s*@[A-Za-z0-9_]+)*)\]/g;

    function citeNumber(key) {
      if (!(key in numberOf)) {
        numberOf[key] = order.length + 1;
        order.push(key);
      }
      return numberOf[key];
    }

    function isSkippable(node) {
      var el = node.parentNode;
      while (el) {
        var tag = el.tagName;
        if (
          tag === "SCRIPT" ||
          tag === "STYLE" ||
          tag === "CODE" ||
          tag === "PRE" ||
          tag === "MJX-CONTAINER"
        ) {
          return true;
        }
        el = el.parentNode;
      }
      return false;
    }

    var walker = document.createTreeWalker(container, NodeFilter.SHOW_TEXT, null, false);
    var textNodes = [];
    var node;
    while ((node = walker.nextNode())) {
      CITE_RE.lastIndex = 0;
      if (!isSkippable(node) && CITE_RE.test(node.nodeValue)) {
        textNodes.push(node);
      }
    }

    textNodes.forEach(function (node) {
      var text = node.nodeValue;
      var frag = document.createDocumentFragment();
      var lastIndex = 0;
      var match;
      CITE_RE.lastIndex = 0;
      while ((match = CITE_RE.exec(text))) {
        frag.appendChild(document.createTextNode(text.slice(lastIndex, match.index)));

        var keys = match[1].split(";").map(function (s) {
          return s.trim().replace(/^@/, "");
        });

        frag.appendChild(document.createTextNode("["));
        var sup = document.createElement("sup");
        keys.forEach(function (key, i) {
          if (i > 0) sup.appendChild(document.createTextNode(","));
          var a = document.createElement("a");
          a.href = "#ref-" + key;
          a.textContent = String(citeNumber(key));
          sup.appendChild(a);
        });
        frag.appendChild(sup);
        frag.appendChild(document.createTextNode("]"));

        lastIndex = match.index + match[0].length;
      }
      frag.appendChild(document.createTextNode(text.slice(lastIndex)));
      node.parentNode.replaceChild(frag, node);
    });

    if (order.length === 0) return;

    var section = document.createElement("section");
    section.className = "dopamine-references";
    var heading = document.createElement("h2");
    heading.textContent = "References";
    section.appendChild(heading);

    var list = document.createElement("ol");
    order.forEach(function (key) {
      var entry = bib[key];
      var li = document.createElement("li");
      li.id = "ref-" + key;
      if (!entry) {
        li.textContent = key;
      } else {
        var lead = entry.authors + ', "' + entry.title + '," ' + entry.venue + ", " + entry.year + ".";
        li.appendChild(document.createTextNode(lead + " "));
        if (entry.doi) {
          var link = document.createElement("a");
          link.href = "https://doi.org/" + entry.doi;
          link.textContent = "doi:" + entry.doi;
          li.appendChild(link);
        }
      }
      list.appendChild(li);
    });
    section.appendChild(list);
    container.appendChild(section);
  }
})();
