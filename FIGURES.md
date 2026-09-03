# FIGURES.md: Figure Generation Doctrine

Operating standard for generating reproducible, publication-grade figures. The rules below establish a consistent visual hierarchy, color palette, and layout across all generated plots.

---

## 0. The stack

```txt
# requirements (pin known-good versions in the real repo)
matplotlib
ultraplot            # proplot successor; styling + multi-panel. NOT proplot (dead since 2023)
seaborn              # use the seaborn.objects grammar interface, OR plotnine. Pick ONE.
adjustText           # ggrepel-equivalent label repeller (fallback)
textalloc            # faster + more deterministic label placement (default for CI)
cmcrameri            # perceptually-uniform continuous colormaps (batlow/vik/roma)
palettable           # named categorical palettes (Tol, Okabe-Ito, ColorBrewer)
PyComplexHeatmap     # annotated/clustered heatmaps with metadata tracks
statannotations      # significance brackets on seaborn plots (pin a version, maintenance is lumpy)
matplotlib-label-lines  # direct line labeling, kills legends
```

**Engine ownership rule:** Decide per figure type which layout engine owns the Figure object. Do not mix plotnine and ultraplot in the same figure. Prefer matplotlib's native `subplot_mosaic` for layout because it is zero-dependency, deterministic, and diffable in pull requests.

---

## 1. Import-time preamble (paste at top of every figure module)

```python
import matplotlib
matplotlib.use("Agg")                       # headless: no backend coin-flip in CI

import matplotlib.pyplot as plt
from matplotlib import rcParams

# --- editable-text export (the pro tell) ---
rcParams["pdf.fonttype"] = 42               # TrueType, text stays editable in Illustrator
rcParams["ps.fonttype"]  = 42               # default Type 3 = per-letter shapes = uneditable
rcParams["svg.fonttype"] = "none"           # keep text as text, not outlined paths

# --- determinism (byte-stable output, no git churn) ---
rcParams["svg.hashsalt"] = "figures"        # any fixed string; else SVG element IDs randomize each run
# seed EVERY rng that touches a figure (jitter, bootstrap, strip/swarm plots)
# strip timestamps on save:  fig.savefig(p, metadata={"Date": None})  (PDF/SVG)

# --- font: register a bundled font, never trust system fallback ---
# from matplotlib import font_manager
# font_manager.fontManager.addfont("assets/FiraSans-Regular.ttf")
# rcParams["font.family"] = "Fira Sans"
```

---

## 2. House style (`house.mplstyle`, applied at import)

```ini
# house.mplstyle  ->  plt.style.use("house.mplstyle")

# size-first: set figsize per-figure to the journal column width, never resize after.
# single col ~3.5in (88mm), double col ~7.2in (180mm). points then mean what they say.
figure.figsize      : 3.5, 2.6
figure.dpi          : 150
savefig.dpi         : 600
savefig.bbox        : tight
savefig.facecolor   : white          # avoid surprise-transparent backgrounds

# typography: one family, one compressed size hierarchy, no stray bold
font.family         : sans-serif
font.size           : 8
axes.titlesize      : 8
axes.titleweight    : regular
axes.labelsize      : 8
xtick.labelsize     : 7
ytick.labelsize     : 7
legend.fontsize     : 7

# data-ink: drop top/right spines, thin lines, inward minor ticks
axes.spines.top     : False
axes.spines.right   : False
axes.linewidth      : 0.8
xtick.direction     : in
ytick.direction     : in
xtick.minor.visible : True
ytick.minor.visible : True
xtick.major.size    : 3.0
xtick.major.width   : 0.8
xtick.minor.size    : 1.8
ytick.major.size    : 3.0
ytick.major.width   : 0.8

# grid: off by default; if on, faint and behind data
axes.grid           : False
grid.linewidth      : 0.5
grid.alpha          : 0.3

# legend: frameless, doesn't sit on data by default
legend.frameon      : False

# layout solver (replaces tight_layout; handles colorbars + outside legends)
figure.constrained_layout.use : True

# colorblind-safe categorical cycle (Okabe-Ito)
axes.prop_cycle : cycler('color', ['0072B2','E69F00','009E73','CC79A7','56B4E9','D55E00','F0E442','000000'])
```

Things the style can't express, do in the wrapper: trim spines to the data range
(`seaborn.despine(trim=True)` or set spine bounds), and kill the corner offset/sci-notation
label by folding the scale into the axis label (e.g. "Intensity (x10^6)").

---

## 3. Text/chart overlap: it's FIVE problems, not one

Agents throw `adjustText` at all of them; it only solves #2. Match tool to category.

| Symptom | Real cause | Fix |
|---|---|---|
| Subplots collide, labels clip at edge, titles overrun | layout engine | `constrained_layout=True` at figure creation. Half of all overlap dies here. |
| Point labels overlap each other/points (gene names on volcano) | label placement | `textalloc` (deterministic, fast) → `adjustText` fallback. **Label top-N only.** |
| Axis tick labels crammed | too many ticks | Horizontal bar chart (best), or `MaxNLocator`, or wrap/truncate. NOT 45deg rotation. |
| Legend sits on the data | legend placement | Move out with `bbox_to_anchor`, or direct-label with `labellines` and delete the legend. |
| Bar value labels collide | manual text placement | `ax.bar_label(bars, padding=3)`. Never hand-place with `ax.text` + offsets. |
| Overplotted point cloud (ink is the problem) | density, not labels | `hexbin` / `mpl-scatter-density` / `datashader`, label only the few that matter. |

### Label placement: the rule that actually fixes volcanoes

`adjustText` is a force-directed solver: each label is a body with repulsive forces from
other labels and points, plus a spring to its anchor, iterated until settled. It **degrades
badly past ~30-50 labels** and is iterative (slow, not perfectly deterministic: fights CI
reproducibility). So the lever is NOT the solver, it's N. A volcano with 12 labeled hits looks
designed; the same plot with 200 labels is static no matter what solver runs.

`textalloc` uses a non-iterative placement strategy: faster, more deterministic. Make it the
default for CI; keep adjustText as the fallback for cases textalloc can't place.

```python
# pattern: filter to top-N by effect+significance, THEN place
def label_top(ax, x, y, names, n=12, by=None):
    import numpy as np
    order = np.argsort(by if by is not None else (np.abs(x) * y))[::-1][:n]
    xs, ys, ls = x[order], y[order], [names[i] for i in order]
    try:
        import textalloc as ta
        ta.allocate(ax, xs, ys, ls,
                    x_scatter=x, y_scatter=y,
                    textsize=7, linecolor="0.5", linewidth=0.5)
    except Exception:
        from adjustText import adjust_text
        texts = [ax.text(xi, yi, l, fontsize=7) for xi, yi, l in zip(xs, ys, ls)]
        adjust_text(texts, x=x, y=y, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5),
                    expand=(1.2, 1.4))
```

---

## 4. Color = hierarchy, not decoration

- Continuous: perceptually-uniform only (`cmcrameri` batlow/vik, or viridis). Never jet: it
  invents structure that isn't in the data.
- Categorical: colorblind-safe, max ~6 hues (Okabe-Ito in the style above, or Tol via palettable).
- **Highest-impact single move:** grey out everything, color only the subject. Non-significant
  volcano points in light grey, hits in color. Saturated-everywhere reads as noise.
- Encode redundantly (color + marker/linestyle) so it survives B&W and colorblindness.

---

## 5. Export

- Vector (PDF/SVG) for all line art. 300+ dpi only matters for raster (photos).
- **Big point clouds:** rasterize ONLY the data layer (`artist.set_rasterized(True)`) at high
  dpi, keep axes/text/annotations vector. Small file, crisp text, no 50k-vector-object monster.
- `figure.dpi` (screen) is separate from `savefig.dpi` (what ships). Set both.

---

## 6. Modular Figure Specification Layer

Avoid writing unconstrained plotting code across individual scripts. Use a thin, validated figure specification layer-such as `volcano(df, ...)`, `ranked_bar(...)`, and `clustered_heatmap(...)`-that emits constrained figures. Calling scripts supply the data specification while the wrapper enforces typography, spacing, and palettes. This ensures repository-wide consistency: modifying the house style updates every figure across the pipeline.

Megatensor wrappers live in `src/megatensor/viz/specs.py`.

---

## 7. Pre-ship checklist (the designed-vs-generated tells)

- [ ] figsize set to final journal width; not resized afterward
- [ ] one font family, one compressed size hierarchy, no stray bold
- [ ] `pdf.fonttype=42`, `svg.fonttype='none'` → text editable in vector editor
- [ ] top/right spines gone; remaining spines trimmed to data range
- [ ] no floating corner offset/sci-notation label; scale folded into axis label
- [ ] context greyed, color only on the subject; ≤6 categorical hues; redundant encoding
- [ ] bars start at zero; consistent rounding/sig figs; units in every axis label
- [ ] legend off the data (moved out or direct-labeled)
- [ ] label_top_N enforced so the placer never sees an intractable count
- [ ] title states the finding, not the variable name
- [ ] deterministic: Agg backend, fixed hashsalt, seeded RNG, stripped metadata
- [ ] dense clouds rasterized as a layer; everything else vector
