"use strict";
(() => {
  const data = window.ATLAS_DATA;
  const el = id => document.getElementById(id);
  const escape = value => String(value).replace(/[&<>"']/g, character => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[character]));
  const number = value => Number(value).toLocaleString("en-US");
  const bytes = value => value >= 1e9 ? `${(value / 1e9).toFixed(2)} GB` : `${(value / 1e6).toFixed(2)} MB`;
  const icons = () => window.lucide?.createIcons();
  const download = (content, filename, type = "text/tab-separated-values") => {
    const url = URL.createObjectURL(new Blob([content], {type}));
    const link = document.createElement("a"); link.href = url; link.download = filename; link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
  if (!data || data.audit_status !== "pass") { el("error").hidden = false; document.querySelectorAll(".view,.heading,.notice,.tabs").forEach(item => item.hidden = true); return; }
  const chroms = Object.keys(data.records_by_chromosome).sort((left, right) => (left === "X" ? 23 : +left) - (right === "X" ? 23 : +right));
  el("total-rows").textContent = number(data.rows);
  el("ok-rows").textContent = number(data.status_counts.ok);
  el("scored-percent").textContent = `${(100 * data.status_counts.ok / data.rows).toFixed(4)}% of allele rows`;
  el("total-chunks").textContent = number(data.chunks);
  el("total-bytes").textContent = bytes(data.bytes);
  const labels = {ok:"Successfully scored", ambiguous_context:"Ambiguous context", unsupported_allele:"Unsupported allele"};
  el("status-breakdown").innerHTML = Object.entries(data.status_counts).sort(([left], [right]) => left === "ok" ? -1 : right === "ok" ? 1 : left.localeCompare(right)).map(([status, count]) => `<div><dt>${escape(labels[status] || status)}</dt><dd>${number(count)}</dd></div>`).join("");
  const maxRecords = Math.max(...Object.values(data.records_by_chromosome));
  el("chrom-chart").innerHTML = chroms.map(chrom => `<a href="#regions" class="chrom-row" data-chrom="${escape(chrom)}" title="Chromosome ${escape(chrom)}: ${number(data.records_by_chromosome[chrom])} records"><span class="chrom-name">${escape(chrom)}</span><div class="bar-track"><div class="bar-fill" style="width:${data.records_by_chromosome[chrom] / maxRecords * 100}%"></div></div><span class="chrom-value">${number(data.records_by_chromosome[chrom])}</span></a>`).join("");
  el("region-chrom").innerHTML = '<option value="all">All chromosomes</option>' + chroms.map(chrom => `<option value="${escape(chrom)}">chr${escape(chrom)}</option>`).join("");
  el("variant-model").innerHTML = '<option value="all">All</option>' + data.models.map(model => `<option value="${escape(model.id)}">${escape(model.id)} · ${escape(model.biosample)}</option>`).join("");
  el("availability").textContent = data.download_base ? "Public shard URLs configured in manifest" : "Not yet published · metadata and examples available";
  el("example-scope").textContent = `${number(data.examples.length)} exported example alleles · ${data.example_scope}. This is not a whole-genome variant search.`;
  el("fingerprints").innerHTML = [["Frozen plan SHA256",data.plan_sha256],["Final audit SHA256",data.audit_sha256],["Audit implementation SHA256",data.audit_code_sha256],["Source release","1000 Genomes · 2019-03-12 biallelic SNV/INDEL · GRCh38"],["Biological validation","Not established"]].map(([key, value]) => `<dt>${escape(key)}</dt><dd>${escape(value)}</dd>`).join("");
  const pageSize = 25;
  let regionPage = 0, variantPage = 0;
  function models() {
    const query = el("model-search").value.trim().toLowerCase();
    const rows = data.models.filter(model => `${model.id} ${model.biosample}`.toLowerCase().includes(query));
    el("model-rows").innerHTML = rows.map(model => `<tr><td class="model-label"><strong>${escape(model.id)}</strong><small>${escape(model.biosample)}</small></td><td class="metric-cell">${model.roc_auc.toFixed(4)}<div class="metric-bar"><span style="width:${model.roc_auc * 100}%"></span></div></td><td class="metric-cell">${model.average_precision.toFixed(4)}<div class="metric-bar ap"><span style="width:${model.average_precision * 100}%"></span></div></td><td>${number(model.counts.test["1"])} / ${number(model.counts.test["0"])}</td><td><div class="encode-links"><a href="https://www.encodeproject.org/files/${escape(model.dhs)}/" target="_blank" rel="noopener">DHS · ${escape(model.dhs)}</a><a href="https://www.encodeproject.org/files/${escape(model.h3k27ac)}/" target="_blank" rel="noopener">H3K27ac · ${escape(model.h3k27ac)}</a></div></td><td>${model.audit_warning ? '<span class="warn">ENCODE ERROR audits retained</span>' : model.small_sample ? '<span class="warn">Small positive set</span>' : '<span class="muted">—</span>'}</td></tr>`).join("") || '<tr><td colspan="6" class="empty">No matching models</td></tr>';
    el("model-count").textContent = `${rows.length} of ${data.models.length} models`;
  }
  function interval() {
    const query = el("region-query").value.trim().replace(/,/g, "");
    if (!query) { el("region-error").textContent = ""; return null; }
    const match = /^(?:chr)?(\d{1,2}|X):(\d+)(?:-(\d+))?$/i.exec(query);
    if (!match || !chroms.includes(match[1].toUpperCase()) || +match[2] < 1 || +(match[3] || match[2]) < +match[2] || !Number.isSafeInteger(+(match[3] || match[2]))) {
      el("region-error").textContent = "Enter a valid GRCh38 interval, for example 22:20000001-25000000."; return false;
    }
    el("region-error").textContent = "";
    return {chrom:match[1].toUpperCase(),start:+match[2],end:+(match[3] || match[2])};
  }
  function pager(prefix, total, page) {
    const pages = Math.max(1, Math.ceil(total / pageSize));
    el(`${prefix}-page`).textContent = `${page + 1} / ${pages}`;
    el(`${prefix}-prev`).disabled = page === 0;
    el(`${prefix}-next`).disabled = page + 1 >= pages;
    el(`${prefix}-count`).textContent = total ? `${number(page * pageSize + 1)}–${number(Math.min(total, (page + 1) * pageSize))} of ${number(total)}` : "0 results";
  }
  function regions() {
    const range = interval();
    const chrom = el("region-chrom").value;
    const rows = range === false ? [] : data.regions.filter(region => (chrom === "all" || region.chrom === chrom) && (!range || region.chrom === range.chrom && region.start <= range.end && region.end >= range.start));
    regionPage = Math.min(regionPage, Math.max(0, Math.ceil(rows.length / pageSize) - 1));
    el("region-rows").innerHTML = rows.slice(regionPage * pageSize, (regionPage + 1) * pageSize).map(region => `<tr><td>chr${escape(region.chrom)}:${number(region.start)}–${number(region.end)}</td><td>${number(region.rows)}</td><td>${number(region.status_counts.ok || 0)}</td><td>${number(region.rows - (region.status_counts.ok || 0))}</td><td>${region.shards}</td><td>${bytes(region.bytes)}</td><td><button class="icon-button region-export" data-index="${region.index}" title="Download region manifest" aria-label="Download manifest for chr${escape(region.chrom)}:${region.start}-${region.end}"><i data-lucide="download" aria-hidden="true"></i></button></td></tr>`).join("") || '<tr><td colspan="7" class="empty">No matching regions</td></tr>';
    pager("region", rows.length, regionPage); icons();
  }
  const kind = row => row.ref.length === 1 && row.alt.length === 1 ? "snv" : row.alt.length > row.ref.length ? "insertion" : row.alt.length < row.ref.length ? "deletion" : "substitution";
  const variantId = row => row.id && row.id !== "." ? row.id : `chr${row.chrom}:${row.pos}:${row.ref}:${row.alt}`;
  function filteredVariants() {
    const query = el("variant-search").value.trim().toLowerCase().replace(/,/g, "");
    const selectedKind = el("variant-kind").value;
    const selectedModel = el("variant-model").value;
    const contexts = data.models.filter(model => selectedModel === "all" || model.id === selectedModel);
    const rows = data.examples.filter(row => `${variantId(row)} ${row.pos} ${row.ref}>${row.alt}`.toLowerCase().includes(query) && (selectedKind === "all" || kind(row) === selectedKind)).flatMap(row => contexts.map(model => ({...row, variant:variantId(row), context:model.id, biosample:model.biosample, ref_score:row[`${model.id}_ref_score`], alt_score:row[`${model.id}_alt_score`], delta:row[`${model.id}_delta`]})));
    const delta = row => row.status !== "ok" || row.delta === "" ? null : Number(row.delta);
    rows.sort((left, right) => {
      const first = delta(left), second = delta(right);
      if (el("variant-sort").value === "position") return +left.pos - +right.pos;
      if (first === null || second === null) return first === null ? second === null ? 0 : 1 : -1;
      return el("variant-sort").value === "magnitude" ? Math.abs(second) - Math.abs(first) : el("variant-sort").value === "increase" ? second - first : first - second;
    });
    return rows;
  }
  function variants() {
    const rows = filteredVariants();
    variantPage = Math.min(variantPage, Math.max(0, Math.ceil(rows.length / pageSize) - 1));
    el("variant-rows").innerHTML = rows.slice(variantPage * pageSize, (variantPage + 1) * pageSize).map(row => {
      const delta = Number(row.delta), valid = row.status === "ok";
      const score = suffix => valid ? Number(row[`${suffix}_score`]).toFixed(6) : "—";
      return `<tr><td>${escape(row.variant)}</td><td class="model-label"><strong>${escape(row.context)}</strong><small>${escape(row.biosample)}</small></td><td>chr${escape(row.chrom)}:${number(row.pos)}</td><td class="allele">${escape(row.ref)} → ${escape(row.alt)}</td><td>${score("ref")}</td><td>${score("alt")}</td><td><div class="effect"><span class="effect-value ${delta >= 0 ? "positive" : "negative"}">${valid ? (delta >= 0 ? "+" : "") + delta.toFixed(6) : "—"}</span><span class="effect-track"><span style="${delta >= 0 ? "left:50%" : "right:50%"};width:${valid ? Math.min(50, Math.abs(delta) * 50) : 0}%;background:${delta >= 0 ? "var(--teal)" : "var(--red)"}"></span></span></div></td><td>${escape(row.status)}</td></tr>`;
    }).join("") || '<tr><td colspan="8" class="empty">No matching example variants</td></tr>';
    pager("variant", rows.length, variantPage);
  }
  function navigate() {
    const requested = location.hash.slice(1) || "overview";
    const tab = ["overview","models","regions","variants","provenance"].includes(requested) ? requested : "overview";
    document.querySelectorAll(".view").forEach(view => view.hidden = view.id !== tab);
    document.querySelectorAll("[data-tab]").forEach(link => { link.classList.toggle("active", link.dataset.tab === tab); if (link.dataset.tab === tab) link.setAttribute("aria-current","page"); else link.removeAttribute("aria-current"); });
    icons();
  }
  el("model-search").addEventListener("input", models);
  for (const id of ["region-chrom","region-query"]) el(id).addEventListener(id === "region-query" ? "input" : "change", () => { regionPage = 0; regions(); });
  el("reset-regions").addEventListener("click", () => { el("region-chrom").value = "all"; el("region-query").value = ""; regionPage = 0; regions(); });
  document.querySelectorAll("[data-chrom]").forEach(link => link.addEventListener("click", () => { el("region-chrom").value = link.dataset.chrom; el("region-query").value = ""; regionPage = 0; regions(); }));
  el("region-rows").addEventListener("click", event => { const button = event.target.closest(".region-export"); if (button) { const region = data.regions.find(item => item.index === +button.dataset.index); download(JSON.stringify({region, catalog:"shards.tsv", audit_sha256:data.audit_sha256, download_base:data.download_base},null,2),`atlas-region-${region.index}.json`,"application/json"); } });
  for (const id of ["variant-model","variant-search","variant-kind","variant-sort"]) el(id).addEventListener(id === "variant-search" ? "input" : "change", () => { variantPage = 0; variants(); });
  el("region-prev").addEventListener("click", () => { regionPage--; regions(); });
  el("region-next").addEventListener("click", () => { regionPage++; regions(); });
  el("variant-prev").addEventListener("click", () => { variantPage--; variants(); });
  el("variant-next").addEventListener("click", () => { variantPage++; variants(); });
  el("export-variants").addEventListener("click", () => { const rows = filteredVariants(); const headers = ["chrom","pos","id","variant","ref","alt","context","biosample","ref_score","alt_score","delta","status"]; download([headers.join("\t"),...rows.map(row => headers.map(key => row[key]).join("\t"))].join("\n") + "\n", "atlas-example-variants.tsv"); });
  addEventListener("hashchange", navigate);
  models(); regions(); variants(); navigate(); icons();
})();