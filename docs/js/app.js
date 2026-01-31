/* Dockline – Static Demo Frontend */
"use strict";

// Base path for static data files
const DATA_BASE = "data";

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let viewer = null;
let currentComplex = null;
let currentRank = 1;
let surfaceVisible = false;
let labelsVisible = false;
let confidenceChart = null;
let cachedProteinText = null;
let currentMetrics = null;

// Cache for loaded data
const dataCache = {
  complexes: null,
  proteins: {},
  poses: {},
  metrics: {},
  ensembles: {},
  comparison: null,
};

// Current complex info (includes color)
let currentComplexInfo = null;

// ---------------------------------------------------------------------------
// Initialisation
// ---------------------------------------------------------------------------
document.addEventListener("DOMContentLoaded", () => {
  initViewer();
  loadComplexes();

  document.getElementById("complexSelect").addEventListener("change", onComplexChange);
  document.getElementById("rankSlider").addEventListener("input", onRankChange);
  document.getElementById("btnSurface").addEventListener("click", toggleSurface);
  document.getElementById("btnLabels").addEventListener("click", toggleLabels);
  document.getElementById("btnReset").addEventListener("click", resetView);
  document.getElementById("comparisonModal").addEventListener("show.bs.modal", loadComparison);
});

function initViewer() {
  const el = document.getElementById("viewer3d");
  viewer = $3Dmol.createViewer(el, {
    backgroundColor: "white",
    antialias: true,
  });
}

// ---------------------------------------------------------------------------
// Welcome overlay
// ---------------------------------------------------------------------------
function hideWelcome() {
  document.getElementById("welcomeOverlay").classList.add("hidden");
  document.getElementById("viewerControls").style.display = "";
}

function showWelcome() {
  document.getElementById("welcomeOverlay").classList.remove("hidden");
  document.getElementById("viewerControls").style.display = "none";
}

// ---------------------------------------------------------------------------
// Complex list loading (from static JSON)
// ---------------------------------------------------------------------------
async function loadComplexes() {
  const sel = document.getElementById("complexSelect");
  try {
    if (!dataCache.complexes) {
      const res = await fetch(`${DATA_BASE}/complexes.json`);
      dataCache.complexes = await res.json();
    }
    const data = dataCache.complexes;

    sel.innerHTML = '<option value="">-- Select a complex --</option>';
    data.forEach(c => {
      const opt = document.createElement("option");
      opt.value = c.name;
      opt.textContent = c.label;
      opt.dataset.numPoses = c.num_poses;
      opt.dataset.ligand = c.ligand || "";
      opt.dataset.color = c.color || "#6c757d";
      sel.appendChild(opt);
    });
  } catch (e) {
    sel.innerHTML = '<option value="">Error loading complexes</option>';
    console.error(e);
  }
}

// ---------------------------------------------------------------------------
// Complex & rank selection
// ---------------------------------------------------------------------------
async function onComplexChange() {
  const sel = document.getElementById("complexSelect");
  const name = sel.value;
  const badgeDiv = document.getElementById("currentLigandBadge");

  if (!name) {
    showWelcome();
    document.getElementById("metricsCard").style.display = "none";
    document.getElementById("chartCard").style.display = "none";
    if (badgeDiv) badgeDiv.innerHTML = "";
    currentComplex = null;
    cachedProteinText = null;
    return;
  }

  currentComplex = name;
  cachedProteinText = null;

  // Update rank slider max and get ligand info
  const opt = sel.selectedOptions[0];
  const numPoses = parseInt(opt.dataset.numPoses) || 10;
  const ligand = opt.dataset.ligand || "";
  const color = opt.dataset.color || "#6c757d";

  // Store current complex info
  currentComplexInfo = { name, ligand, color };

  // Display ligand badge
  if (badgeDiv && ligand) {
    badgeDiv.innerHTML = `<span class="badge" style="background-color: ${color};">${ligand}</span>`;
  }

  const slider = document.getElementById("rankSlider");
  slider.max = numPoses;
  slider.value = 1;
  currentRank = 1;
  document.getElementById("rankLabel").textContent = "1";

  hideWelcome();
  document.getElementById("metricsCard").style.display = "";
  document.getElementById("chartCard").style.display = "";

  await loadScene(name, 1);
  loadMetrics(name, 1);
  loadEnsemble(name);
}

async function onRankChange() {
  const rank = parseInt(document.getElementById("rankSlider").value);
  document.getElementById("rankLabel").textContent = rank;
  currentRank = rank;
  if (!currentComplex) return;

  await loadScene(currentComplex, rank);
  loadMetrics(currentComplex, rank);
}

// ---------------------------------------------------------------------------
// Data fetching (static files)
// ---------------------------------------------------------------------------
async function fetchProtein(name) {
  const cacheKey = name;
  if (dataCache.proteins[cacheKey]) {
    return dataCache.proteins[cacheKey];
  }
  const res = await fetch(`${DATA_BASE}/${name}/protein.pdb`);
  const text = await res.text();
  dataCache.proteins[cacheKey] = text;
  return text;
}

async function fetchPose(name, rank) {
  const cacheKey = `${name}_${rank}`;
  if (dataCache.poses[cacheKey]) {
    return dataCache.poses[cacheKey];
  }
  const res = await fetch(`${DATA_BASE}/${name}/pose_${rank}.sdf`);
  const text = await res.text();
  dataCache.poses[cacheKey] = text;
  return text;
}

async function fetchMetrics(name, rank) {
  const cacheKey = `${name}_${rank}`;
  if (dataCache.metrics[cacheKey]) {
    return dataCache.metrics[cacheKey];
  }
  const res = await fetch(`${DATA_BASE}/${name}/metrics_${rank}.json`);
  const data = await res.json();
  dataCache.metrics[cacheKey] = data;
  return data;
}

async function fetchEnsemble(name) {
  if (dataCache.ensembles[name]) {
    return dataCache.ensembles[name];
  }
  const res = await fetch(`${DATA_BASE}/${name}/ensemble.json`);
  const data = await res.json();
  dataCache.ensembles[name] = data;
  return data;
}

async function fetchComparison() {
  if (dataCache.comparison) {
    return dataCache.comparison;
  }
  const res = await fetch(`${DATA_BASE}/comparison.json`);
  const data = await res.json();
  dataCache.comparison = data;
  return data;
}

// ---------------------------------------------------------------------------
// 3Dmol.js scene management
// ---------------------------------------------------------------------------
async function loadScene(name, rank) {
  viewer.clear();
  viewer.removeAllLabels();

  const [protText, poseText] = await Promise.all([
    fetchProtein(name),
    fetchPose(name, rank),
  ]);

  cachedProteinText = protText;

  // Model 0: Protein
  viewer.addModel(protText, "pdb");
  viewer.setStyle({ model: 0 }, {
    cartoon: { color: "spectrum", opacity: 0.85 }
  });

  // Model 1: Ligand
  viewer.addModel(poseText, "sdf");
  viewer.setStyle({ model: 1 }, {
    stick: { radius: 0.2, colorscheme: "default" },
    sphere: { scale: 0.25, colorscheme: "default" }
  });

  if (surfaceVisible) {
    addSurface();
  }

  if (labelsVisible) {
    addResidueLabels();
  }

  viewer.zoomTo({ model: 1 });
  viewer.zoom(0.8);
  viewer.render();
}

function addSurface() {
  viewer.addSurface($3Dmol.SurfaceType.VDW, {
    opacity: 0.15, color: "lightblue"
  }, { model: 0, within: { distance: 8, sel: { model: 1 } } });
}

function addResidueLabels() {
  try {
    viewer.setStyle(
      { model: 0, within: { distance: 4.5, sel: { model: 1 } } },
      {
        cartoon: { color: "spectrum", opacity: 0.85 },
        stick: { radius: 0.1, colorscheme: "default" }
      }
    );
    viewer.addResLabels(
      { model: 0, atom: "CA", within: { distance: 5.0, sel: { model: 1 } } },
      {
        font: "Arial",
        fontSize: 11,
        fontColor: "white",
        backgroundColor: "rgba(40,40,40,0.75)",
        backgroundOpacity: 0.75,
        borderThickness: 0,
        padding: 2,
      }
    );
  } catch (e) {
    console.warn("addResLabels not supported:", e);
  }
}

function toggleSurface() {
  surfaceVisible = !surfaceVisible;
  const btn = document.getElementById("btnSurface");
  btn.classList.toggle("active", surfaceVisible);
  if (currentComplex) loadScene(currentComplex, currentRank);
}

function toggleLabels() {
  labelsVisible = !labelsVisible;
  const btn = document.getElementById("btnLabels");
  btn.classList.toggle("active", labelsVisible);
  if (currentComplex) loadScene(currentComplex, currentRank);
}

function resetView() {
  if (viewer) {
    viewer.zoomTo({ model: 1 });
    viewer.zoom(0.8);
    viewer.render();
  }
}

// ---------------------------------------------------------------------------
// Metrics panel + viewer overlay
// ---------------------------------------------------------------------------
async function loadMetrics(name, rank) {
  const panel = document.getElementById("metricsPanel");
  panel.innerHTML = '<p class="text-muted small">Loading metrics...</p>';

  try {
    const m = await fetchMetrics(name, rank);
    currentMetrics = m;
    if (m.error && !m.confidence) {
      panel.innerHTML = `<p class="text-danger small">${m.error}</p>`;
      return;
    }
    renderMetrics(m);
    updateViewerOverlay(m, rank);
  } catch (e) {
    panel.innerHTML = '<p class="text-danger small">Failed to load metrics.</p>';
  }
}

function renderMetrics(m) {
  const panel = document.getElementById("metricsPanel");
  let html = "";

  if (m.confidence !== null && m.confidence !== undefined) {
    html += metricRow("Confidence", m.confidence.toFixed(2));
  }

  if (m.rmsd_to_crystal !== null && m.rmsd_to_crystal !== undefined) {
    const cls = m.rmsd_to_crystal < 2 ? "rmsd-green" :
                m.rmsd_to_crystal < 4 ? "rmsd-yellow" : "rmsd-red";
    html += metricRow("RMSD to Crystal",
      `<span class="${cls}">${m.rmsd_to_crystal.toFixed(2)} &Aring;</span>`);
  }

  html += metricRow("H-bonds", m.hbonds ? m.hbonds.length : 0);
  if (m.hbonds && m.hbonds.length) {
    html += '<div class="contact-list">';
    m.hbonds.forEach(h => {
      html += `<span>${h.protein_residue}:${h.protein_atom} (${h.distance}&Aring;)</span> `;
    });
    html += "</div>";
  }

  html += metricRow("Salt Bridges", m.salt_bridges ? m.salt_bridges.length : 0);
  if (m.salt_bridges && m.salt_bridges.length) {
    html += '<div class="contact-list">';
    m.salt_bridges.forEach(s => {
      html += `<span>${s.residue} (${s.distance}&Aring;)</span> `;
    });
    html += "</div>";
  }

  html += metricRow("Hydrophobic Contacts", m.hydrophobic_contacts || 0);
  html += metricRow("Buried SA", (m.buried_surface_area || 0) + " &Aring;&sup2;");

  if (m.close_contacts && m.close_contacts.length) {
    html += metricRow("Close Contact Residues", m.close_contacts.length);
    html += '<div class="contact-list">';
    m.close_contacts.forEach(c => { html += `<span>${c}</span> `; });
    html += "</div>";
  }

  panel.innerHTML = html;
}

function metricRow(label, value) {
  return `<div class="metric-row"><span class="metric-label">${label}</span><span class="metric-value">${value}</span></div>`;
}

function updateViewerOverlay(m, rank) {
  let overlay = document.getElementById("viewerInfoOverlay");
  if (!overlay) {
    overlay = document.createElement("div");
    overlay.id = "viewerInfoOverlay";
    overlay.style.cssText = "position:absolute;top:8px;left:8px;z-index:20;pointer-events:none;";
    document.getElementById("viewer3d").parentElement.appendChild(overlay);
  }

  let badges = "";
  badges += `<span class="badge bg-dark me-1">Rank ${rank}</span>`;
  if (m.confidence !== null && m.confidence !== undefined) {
    badges += `<span class="badge bg-primary me-1">Conf: ${m.confidence.toFixed(2)}</span>`;
  }
  if (m.rmsd_to_crystal !== null && m.rmsd_to_crystal !== undefined) {
    const cls = m.rmsd_to_crystal < 2 ? "bg-success" :
                m.rmsd_to_crystal < 4 ? "bg-warning text-dark" : "bg-danger";
    badges += `<span class="badge ${cls} me-1">RMSD: ${m.rmsd_to_crystal.toFixed(2)} \u00C5</span>`;
  }
  if (m.hbonds) {
    badges += `<span class="badge bg-info text-dark me-1">H-bonds: ${m.hbonds.length}</span>`;
  }
  overlay.innerHTML = badges;
}

// ---------------------------------------------------------------------------
// Ensemble / confidence chart
// ---------------------------------------------------------------------------
async function loadEnsemble(name) {
  const divSummary = document.getElementById("diversitySummary");
  divSummary.textContent = "";

  try {
    const data = await fetchEnsemble(name);
    renderConfidenceChart(data);

    if (data.pairwise_rmsd) {
      const p = data.pairwise_rmsd;
      divSummary.innerHTML =
        `Pose diversity (pairwise RMSD): mean ${p.mean} &Aring;, ` +
        `min ${p.min} &Aring;, max ${p.max} &Aring;`;
    }
  } catch (e) {
    console.error("Ensemble load error:", e);
  }
}

function renderConfidenceChart(data) {
  const canvas = document.getElementById("confidenceChart");
  const ctx = canvas.getContext("2d");
  const scores = (data.confidence_scores || []);

  const labels = scores.map(s => `R${s.rank}`);
  const values = scores.map(s => s.confidence != null ? -s.confidence : 0);
  const colors = scores.map((_, i) => i === 0 ? "#198754" : "#6c757d");

  if (confidenceChart) {
    confidenceChart.destroy();
    confidenceChart = null;
  }

  confidenceChart = new Chart(ctx, {
    type: "bar",
    data: {
      labels,
      datasets: [{
        label: "Confidence (-score)",
        data: values,
        backgroundColor: colors,
        borderRadius: 3,
      }],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      animation: { duration: 300 },
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: item => {
              const raw = scores[item.dataIndex].confidence;
              return `Confidence: ${raw != null ? raw.toFixed(2) : "N/A"}`;
            }
          }
        }
      },
      scales: {
        y: {
          beginAtZero: true,
          title: { display: true, text: "Score (negated)" },
        },
        x: {
          ticks: { font: { size: 10 } },
        },
      },
      onClick: (_evt, elements) => {
        if (elements.length > 0) {
          const idx = elements[0].index;
          const rank = scores[idx].rank;
          const slider = document.getElementById("rankSlider");
          slider.value = rank;
          slider.dispatchEvent(new Event("input"));
        }
      }
    }
  });
}

// ---------------------------------------------------------------------------
// Comparison modal
// ---------------------------------------------------------------------------
async function loadComparison() {
  const wrap = document.getElementById("comparisonTableWrap");
  wrap.innerHTML = '<p class="text-muted">Loading...</p>';

  try {
    const rows = await fetchComparison();
    if (!rows.length) {
      wrap.innerHTML = '<p class="text-muted">No complexes to compare.</p>';
      return;
    }
    renderComparisonTable(rows);
  } catch (e) {
    wrap.innerHTML = '<p class="text-danger">Failed to load comparison.</p>';
  }
}

function renderComparisonTable(rows) {
  const wrap = document.getElementById("comparisonTableWrap");
  let html = '<div class="table-responsive"><table class="table table-sm table-striped align-middle">';
  html += "<thead><tr>";
  html += "<th>Ligand</th><th>Complex</th><th>Confidence</th><th>RMSD</th>";
  html += "<th>H-bonds</th><th>Salt Bridges</th><th>Hydro.</th>";
  html += "<th>Buried SA</th>";
  html += "</tr></thead><tbody>";

  rows.forEach(m => {
    const conf = m.confidence !== null && m.confidence !== undefined
      ? m.confidence.toFixed(2) : "N/A";
    const rmsd = m.rmsd_to_crystal !== null && m.rmsd_to_crystal !== undefined
      ? m.rmsd_to_crystal.toFixed(2) + " \u00C5" : "\u2014";
    const ligand = m.ligand || "";
    const color = m.color || "#6c757d";
    html += "<tr>";
    html += `<td><span class="badge" style="background-color: ${color};">${ligand}</span></td>`;
    html += `<td class="fw-bold">${m.label}</td>`;
    html += `<td>${conf}</td>`;
    html += `<td>${rmsd}</td>`;
    html += `<td>${m.hbonds ? m.hbonds.length : 0}</td>`;
    html += `<td>${m.salt_bridges ? m.salt_bridges.length : 0}</td>`;
    html += `<td>${m.hydrophobic_contacts || 0}</td>`;
    html += `<td>${m.buried_surface_area || 0} \u00C5\u00B2</td>`;
    html += "</tr>";
  });

  html += "</tbody></table></div>";
  wrap.innerHTML = html;
}
