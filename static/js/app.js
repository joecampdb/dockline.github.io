/* Dockline – frontend logic */
"use strict";

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let viewer = null;
let currentComplex = null;
let currentRank = 1;
let surfaceVisible = false;
let labelsVisible = false;
let confidenceChart = null;
let cachedProteinText = null;   // avoid re-fetching protein on rank change
let currentMetrics = null;      // latest metrics for overlay

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
  document.getElementById("samplesSlider").addEventListener("input", e => {
    document.getElementById("samplesLabel").textContent = e.target.value;
  });
  document.getElementById("dockingForm").addEventListener("submit", onDockSubmit);
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
// Complex list loading (no auto-select)
// ---------------------------------------------------------------------------
async function loadComplexes() {
  const sel = document.getElementById("complexSelect");
  try {
    const res = await fetch("/api/complexes");
    const data = await res.json();
    // Keep the placeholder option, add complexes after it
    sel.innerHTML = '<option value="">-- Select a complex --</option>';
    data.forEach(c => {
      const opt = document.createElement("option");
      opt.value = c.name;
      opt.textContent = c.label;
      opt.dataset.numPoses = c.num_poses;
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
  if (!name) {
    showWelcome();
    document.getElementById("metricsCard").style.display = "none";
    document.getElementById("chartCard").style.display = "none";
    currentComplex = null;
    cachedProteinText = null;
    return;
  }

  currentComplex = name;
  cachedProteinText = null;

  // Update rank slider max
  const opt = sel.selectedOptions[0];
  const numPoses = parseInt(opt.dataset.numPoses) || 10;
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

  // Reload full scene with cached protein (avoids 3Dmol model-index bugs)
  await loadScene(currentComplex, rank);
  loadMetrics(currentComplex, rank);
}

// ---------------------------------------------------------------------------
// 3Dmol.js scene management
// ---------------------------------------------------------------------------
async function loadScene(name, rank) {
  viewer.clear();
  viewer.removeAllLabels();

  // Fetch protein (use cache if same complex) + pose
  let protPromise;
  if (cachedProteinText && currentComplex === name) {
    protPromise = Promise.resolve(cachedProteinText);
  } else {
    protPromise = fetch(`/api/complex/${name}/protein`).then(r => r.text());
  }

  const [protText, poseText] = await Promise.all([
    protPromise,
    fetch(`/api/complex/${name}/pose/${rank}`).then(r => r.text()),
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

  // Surface (if toggled on)
  if (surfaceVisible) {
    addSurface();
  }

  // Residue labels (if toggled on)
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
  // Label nearby protein residues (CA atoms within 5A of ligand)
  try {
    // Show binding-site residues as sticks
    viewer.setStyle(
      { model: 0, within: { distance: 4.5, sel: { model: 1 } } },
      {
        cartoon: { color: "spectrum", opacity: 0.85 },
        stick: { radius: 0.1, colorscheme: "default" }
      }
    );
    // Add residue-name labels
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
  // Reload scene to apply (3Dmol surface add/remove is fragile)
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
    const res = await fetch(`/api/complex/${name}/metrics/${rank}`);
    const m = await res.json();
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

  // Confidence
  if (m.confidence !== null && m.confidence !== undefined) {
    html += metricRow("Confidence", m.confidence.toFixed(2));
  }

  // RMSD to crystal
  if (m.rmsd_to_crystal !== null && m.rmsd_to_crystal !== undefined) {
    const cls = m.rmsd_to_crystal < 2 ? "rmsd-green" :
                m.rmsd_to_crystal < 4 ? "rmsd-yellow" : "rmsd-red";
    html += metricRow("RMSD to Crystal",
      `<span class="${cls}">${m.rmsd_to_crystal.toFixed(2)} &Aring;</span>`);
  }

  // H-bonds
  html += metricRow("H-bonds", m.hbonds ? m.hbonds.length : 0);
  if (m.hbonds && m.hbonds.length) {
    html += '<div class="contact-list">';
    m.hbonds.forEach(h => {
      html += `<span>${h.protein_residue}:${h.protein_atom} (${h.distance}&Aring;)</span> `;
    });
    html += "</div>";
  }

  // Salt bridges
  html += metricRow("Salt Bridges", m.salt_bridges ? m.salt_bridges.length : 0);
  if (m.salt_bridges && m.salt_bridges.length) {
    html += '<div class="contact-list">';
    m.salt_bridges.forEach(s => {
      html += `<span>${s.residue} (${s.distance}&Aring;)</span> `;
    });
    html += "</div>";
  }

  // Hydrophobic contacts
  html += metricRow("Hydrophobic Contacts", m.hydrophobic_contacts || 0);

  // Buried SA
  html += metricRow("Buried SA", (m.buried_surface_area || 0) + " &Aring;&sup2;");

  // Close contacts
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
  // Remove existing overlay if any
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
// Ensemble / confidence chart (negated scores → bars go UP)
// ---------------------------------------------------------------------------
async function loadEnsemble(name) {
  const divSummary = document.getElementById("diversitySummary");
  divSummary.textContent = "";

  try {
    const res = await fetch(`/api/complex/${name}/ensemble`);
    const data = await res.json();
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
  // Negate so bars grow upward (DiffDock scores are negative; higher = better)
  const values = scores.map(s => s.confidence != null ? -s.confidence : 0);

  // Gradient: rank 1 is green, others fade to grey
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
    const res = await fetch("/api/comparison");
    const rows = await res.json();
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
  html += "<th>Complex</th><th>Confidence</th><th>RMSD</th>";
  html += "<th>H-bonds</th><th>Salt Bridges</th><th>Hydro. Contacts</th>";
  html += "<th>Buried SA</th><th>Close Residues</th>";
  html += "</tr></thead><tbody>";

  rows.forEach(m => {
    const conf = m.confidence !== null && m.confidence !== undefined
      ? m.confidence.toFixed(2) : "N/A";
    const rmsd = m.rmsd_to_crystal !== null && m.rmsd_to_crystal !== undefined
      ? m.rmsd_to_crystal.toFixed(2) + " \u00C5" : "\u2014";
    html += "<tr>";
    html += `<td class="fw-bold">${m.label}</td>`;
    html += `<td>${conf}</td>`;
    html += `<td>${rmsd}</td>`;
    html += `<td>${m.hbonds ? m.hbonds.length : 0}</td>`;
    html += `<td>${m.salt_bridges ? m.salt_bridges.length : 0}</td>`;
    html += `<td>${m.hydrophobic_contacts || 0}</td>`;
    html += `<td>${m.buried_surface_area || 0} \u00C5\u00B2</td>`;
    html += `<td>${m.close_contacts ? m.close_contacts.length : 0}</td>`;
    html += "</tr>";
  });

  html += "</tbody></table></div>";
  wrap.innerHTML = html;
}

// ---------------------------------------------------------------------------
// Docking job submission & polling
// ---------------------------------------------------------------------------
async function onDockSubmit(e) {
  e.preventDefault();
  const form = e.target;
  const formData = new FormData(form);

  try {
    const res = await fetch("/api/dock", { method: "POST", body: formData });
    const data = await res.json();
    if (data.error) {
      alert("Error: " + data.error);
      return;
    }
    form.reset();
    document.getElementById("samplesLabel").textContent = "10";

    // Show progress toast
    showJobToast(data.job_id, formData.get("complex_name"));
    pollJob(data.job_id);
  } catch (err) {
    alert("Submission failed: " + err.message);
  }
}

function showJobToast(jobId, complexName) {
  const container = document.getElementById("jobToastContainer");
  const div = document.createElement("div");
  div.className = "job-toast";
  div.id = `job-${jobId}`;
  div.innerHTML = `
    <div class="d-flex justify-content-between align-items-center mb-1">
      <span class="job-title">${complexName}</span>
      <span class="badge bg-secondary" id="job-badge-${jobId}">queued</span>
    </div>
    <div class="progress">
      <div class="progress-bar progress-bar-striped progress-bar-animated"
           id="job-bar-${jobId}" style="width:5%"></div>
    </div>
    <small class="text-muted" id="job-step-${jobId}">Waiting...</small>
  `;
  container.appendChild(div);
}

const STEP_PROGRESS = {
  queued: 5,
  preparing: 15,
  starting: 20,
  embeddings: 35,
  loading_model: 50,
  diffusion: 70,
  scoring: 85,
  saving: 95,
  complete: 100,
};

async function pollJob(jobId) {
  const interval = setInterval(async () => {
    try {
      const res = await fetch(`/api/dock/${jobId}/status`);
      const data = await res.json();

      const badge = document.getElementById(`job-badge-${jobId}`);
      const bar = document.getElementById(`job-bar-${jobId}`);
      const step = document.getElementById(`job-step-${jobId}`);

      if (!badge) { clearInterval(interval); return; }

      badge.textContent = data.status;
      step.textContent = data.progress_step || data.status;

      const pct = STEP_PROGRESS[data.progress_step] || STEP_PROGRESS[data.status] || 10;
      bar.style.width = pct + "%";

      if (data.status === "complete") {
        clearInterval(interval);
        badge.className = "badge bg-success";
        bar.classList.remove("progress-bar-animated");
        bar.style.width = "100%";
        step.textContent = "Complete!";
        // Reload complex list and auto-select the new result
        await loadComplexes();
        const sel = document.getElementById("complexSelect");
        sel.value = data.complex_name;
        sel.dispatchEvent(new Event("change"));
      } else if (data.status === "failed") {
        clearInterval(interval);
        badge.className = "badge bg-danger";
        bar.classList.remove("progress-bar-animated");
        bar.classList.add("bg-danger");
        step.textContent = data.error || "Failed";
      }
    } catch (e) {
      console.error("Poll error:", e);
    }
  }, 3000);
}
