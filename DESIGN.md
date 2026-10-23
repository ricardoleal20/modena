# MODENA design system

Identity and interface guidance · 2026-10-23

MODENA means **Molecular Dynamics Engine for Numerical Atomistics**. Its visual
identity should feel like a reliable scientific instrument: precise geometry,
legible numbers, quiet surfaces, and explicit evidence. The engine is an early
Rust library and CLI with Python and MCP interfaces. The workbench in
[`docs/ui-prototype/index.html`](docs/ui-prototype/index.html) is a standalone
static prototype, not a shipped simulation frontend.

This document introduces a **proposed visual system**. Logo development remains
local and deferred; no logo artwork is included in this change. It records the
prototype's existing choices separately. It does not claim that
the proposed tokens, accessibility behaviors, or interface states are already
implemented. Product and validation statements remain governed by
[`README.md`](README.md) and [`docs/VALIDATED_SCOPE.md`](docs/VALIDATED_SCOPE.md).

## 1. Visual direction

- **Precise, open geometry.** The angular M and periodic-cell corners express
  numerical structure. The central particle supplies a single controlled accent.
- **Data before decoration.** Give configuration, units, diagnostics, provenance,
  and errors the strongest hierarchy. Keep brand illustration outside plots.
- **Restrained color.** Deep ink and paper form the structure; teal identifies
  intentional actions and selection. Reserve semantic colors for actual states.
- **Scientific restraint.** A completed run, low energy drift, or attractive plot
  is not proof of physical validity. Show the scope and method behind a result.

The symbol is an abstract identity, not a chemical structure, trajectory,
periodic-boundary algorithm diagram, or promise that all box geometries are
supported. Do not present its connecting strokes as simulated molecular bonds.

## 2. Logo development: deferred

An angular M, open periodic-cell corners, and a single particle are an explored
identity concept. The artwork and its local previews are outside this change.
Use ordinary system-font text reading “MODENA” until reviewed artwork is added
in a separate contribution. Geometry, file inventory, clear space, minimum
reproduction sizes, and final usage rules remain to be established there.

A future logo contribution should include light, inverse, monochrome, and small
variants with a reproducible preview. Verify each at its intended display size,
retain accessible titles and descriptions, and document authorship and license.
Avoid font dependencies, embedded raster images, scripts, and external resources
in SVG artwork. Inverse artwork needs an explicit dark host surface.

For future HTML logo images, use `alt="MODENA"`; when adjacent text already names
MODENA, use `alt=""` to avoid duplicate speech. Give inline SVG title/description
IDs unique prefixes when embedding multiple copies. Decorative copies should
use `aria-hidden="true"` and omit an announced image role.

## 3. Existing prototype: observed baseline

These values are present in the prototype. They are **observations, not new
requirements or an accessibility certification**.

| Existing token or pattern | Current value |
|---|---|
| Page / panel / hover | `--bg: #FAFAFA`, `--panel: #FFFFFF`, `--panel-hover: #F4F4F5` |
| Text | `--fg: #18181B`, `--fg-strong: #0A0A0A`, `--fg-secondary: #71717A`, `--fg-muted: #A1A1AA` |
| Rules | `--border: #E5E5E5`, `--border-strong: #D4D4D8` |
| States | `--green: #16A34A`, `--amber: #D97706`, `--red: #DC2626` |
| State tints | `#DCFCE7`, `#FEF3C7`, `#FEE2E2` |
| Radius | Small 6 px; default 10 px; large 12 px |
| Typography | 14 px base, 1.5 line height; system sans with named fallbacks; system monospace for measurements |
| Layout | 52 px top bar; 24 px bar/panel padding; 32 px tab content padding; 24 px main gaps; 1200 px content maximum |
| Charts | Mostly grayscale strokes; generated sample series and particles |
| Interaction | Tab changes, local play/pause appearance, segmented selections, expandable reasoning examples |

The prototype includes synthetic measurements, a sample campaign, a confidence
percentage, and illustrative agent reasoning. Treat these as **demonstration
content**. They are not run evidence, calibrated confidence, or backend
capabilities. Some small muted labels and color-only cues need review before
production. The prototype removes default outlines, so a future implementation
must add a visible focus treatment before following that reset.

## 4. Proposed color tokens

Use semantic names so application code can change theme without changing meaning.
This proposal uses light surfaces as the default. The inverse colors cover
controlled dark contexts; they do not constitute a complete dark application
theme.

| Proposed token | Hex | Role |
|---|---|---|
| `--color-ink` | `#10262B` | Primary text, strong icons, brand structure |
| `--color-accent` | `#087F78` | Primary action, link, selection indicator, logo particle on light |
| `--color-accent-inverse` | `#5ED8C1` | Logo particle and restrained accents on ink only |
| `--color-canvas` | `#F7FAF8` | Page background |
| `--color-surface` | `#FFFFFF` | Data and configuration surfaces |
| `--color-text-secondary` | `#4D666B` | Secondary readable text |
| `--color-control-border` | `#78908E` | Boundary when required to identify a control |
| `--color-rule` | `#D8E3DF` | Decorative separators; never a required control boundary |
| `--color-warning` | `#9D4E09` | Warning text/icon paired with a label |
| `--color-error` | `#B63D47` | Failed or invalid state paired with a label |
| `--color-information` | `#315F9B` | Informational state or second plot series |

Do not use teal as an automatic success signal. “Run complete” describes job
execution; “validation passed” requires a named check and its threshold. Use
text and a distinct icon alongside any status color. Links remain underlined
in prose. White text may be used on solid teal buttons. Mint is unsuitable for
text, required icons, or fine marks on white.

### Contrast verification

Ratios below use the WCAG relative luminance formula for **opaque sRGB colors**:
linearize each normalized channel (`c/12.92` when `c ≤ 0.04045`, otherwise
`((c + 0.055)/1.055)^2.4`), then compute `L = 0.2126R + 0.7152G + 0.0722B` and
`(L_lighter + 0.05) / (L_darker + 0.05)`. Values are rounded to two decimals for
display; pass decisions use unrounded results.

| Foreground / background | Ratio | Intended use |
|---|---:|---|
| Ink / white | 15.74:1 | Body text |
| Ink / paper | 14.98:1 | Body text |
| Secondary / white | 6.13:1 | Secondary text |
| Secondary / paper | 5.83:1 | Secondary text |
| Teal / white | 4.86:1 | Text, links, meaningful graphics |
| Teal / paper | 4.63:1 | Text, links, meaningful graphics |
| White / teal | 4.86:1 | Button text |
| Paper / ink | 14.98:1 | Inverse text and logo |
| Mint / ink | 9.05:1 | Inverse logo accent |
| Control border / white | 3.40:1 | Required control boundary, not body text |
| Control border / paper | 3.23:1 | Required control boundary, not body text |
| Warning / white | 5.95:1 | Labeled warning text |
| Error / white | 5.60:1 | Labeled error text |
| Information / white | 6.48:1 | Labeled information text |
| Rule / white | 1.31:1 | Decoration only |
| Mint / white | 1.74:1 | Do not use for meaningful foreground |

Normal text targets 4.5:1; large text 3:1. Meaningful control cues and graphics
target 3:1 against adjacent surfaces. These follow the W3C explanations of
[text contrast](https://www.w3.org/WAI/WCAG22/Understanding/contrast-minimum.html)
and [non-text contrast](https://www.w3.org/WAI/WCAG22/Understanding/non-text-contrast.html).
Color arithmetic alone does not establish interface accessibility. Opacity,
thin strokes, overlays, hover states, focus rings, and actual backgrounds need
verification in the implemented view.

## 5. Typography and numeric notation

Use locally available system fonts; do not fetch fonts for the interface.

```css
--font-sans: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
--font-mono: ui-monospace, "SFMono-Regular", Menlo, Consolas, monospace;
```

The font stacks reference host fonts without redistributing font files; their
licenses remain with their respective owners. This proposal adds no external
font dependency. If a future contribution bundles a font,
include the exact font source, version, redistribution license, and required
license files before adopting it. Do not assume a named prototype fallback is
bundled or installed.

| Role | Suggested size / line height | Weight and treatment |
|---|---|---|
| Page title | 24 / 32 px | 600; sentence case |
| Section title | 18 / 26 px | 600 |
| Body and form input | 14–16 / 21–24 px | 400; use 16 px for narrow/touch layouts |
| Field label and table header | 13 / 20 px | 500–600 |
| Supporting annotation | 12 / 18 px minimum | Secondary color; never hide critical state here |
| Measured value and code | 13–14 / 20–21 px | Monospace; tabular numerals |

Use `font-variant-numeric: tabular-nums` on comparison tables; right-align numeric
columns and align precision within each column. Retain full precision in exported
data while choosing a meaningful display precision. Keep signs and exponents
visible. Label unavailable values as “Not recorded” or explain the em dash; never
render missing values as zero. Avoid all-caps prose and decorative letterspacing
in controls. The plain MODENA name may remain uppercase.

## 6. Spacing, surfaces, and layout

Use a 4 px base scale: **4, 8, 12, 16, 24, 32, 48, 64**. Place related label/value
pairs 4–8 px apart, adjacent controls 8–12 px apart, surface contents 16–24 px
from their edges, and major sections 24–32 px apart. Use 4 px radii for compact
controls and 8 px for surfaces in future branded UI. These are proposed values;
they do not change the prototype's 6/10/12 px radii.

Prefer a strong heading and quiet rule to repeated card nesting. White panels
sit on paper with a 1 px decorative rule. Use elevation only for overlays that
actually float above content. Avoid glossy particles, gradients, glow, and
decorative blur around scientific data. Reserve pure ink fields for identity
proofs or controlled inverse contexts.

Keep configuration before execution, execution state next to run identity, and
diagnostics with their provenance. At wide sizes, plots and configuration may
sit side by side. At narrow sizes, stack by task order; prevent page-level
horizontal scrolling. Tables may scroll inside a labeled region with sticky
headers, while run name, status, and essential controls remain available. Do not
hide validation warnings to make a layout fit. Let explanatory prose reflow to
roughly 60–80 characters per line.

## 7. Scientific plots, units, and provenance

- Put quantity **and unit convention** in axes, headings, and exports. The engine
  computes internally in reduced units; real-unit conversion occurs at the I/O
  boundary. See [`units.rs`](crates/modena-core/src/units.rs). Label reduced
  quantities explicitly, for example `T* (kBT/ε)` and `t* (t/τ)`; distinguish
  temperature from thermostat target. Real outputs may use Å, kcal/mol, amu, K,
  ps, and atm as appropriate. Never silently reinterpret a reduced value as K.
- Distinguish totals from per-particle values, instantaneous measurements from
  averages, and energy drift from statistical uncertainty. Show the actual
  normalization used for drift and the calculation window. Never add an
  invented confidence score to an agent observation.
- Keep trace colors stable for the same quantity within a run comparison. Start
  with teal `#087F78`, blue `#315F9B`, and ink `#10262B` on white. Pair each series
  with a direct label, distinct dash or marker, and a text/table alternative.
  These colors are not guaranteed distinguishable from each other or under every
  color-vision condition; do not rely on color alone. More series require a
  reviewed palette and additional encodings.
- Use at least 1.5–2 CSS px meaningful plot strokes. Separate decorative grids
  from readable ticks. Do not use the faint separator color for axes or labels
  required to read the chart. Check plot contrast against actual plot surfaces.
- State axis ranges, scale type, sampling stride, and any smoothing or decimation.
  Never smooth away instability. Retain gaps for missing/invalid measurements.
  Give uncertainty bands their method and sample window; do not draw them when
  uncertainty has not been estimated.
- Preserve run ID, engine version/commit, input/configuration, unit system, seed
  when applicable, potential parameters, ensemble, timestep, step count, and
  output paths with a result. Label demo and derived content visibly. Keep an
  agent's suggestion distinguishable from a measured observation.
- A particle viewer needs a box/length scale and explicit coloring legend. Mark
  bonds only when the loaded topology provides them. Show camera orientation and
  periodic wrapping behavior. Playback frame/time is separate from wall time or
  engine progress. Export a still with the same legends and provenance.

## 8. States, accessibility, and motion

Future interfaces must provide complete idle, validating, running, paused,
completed, failed, cancelled, empty, and unavailable states wherever supported
by the real backend. Label unsupported operations honestly. Confirmation of
an action must describe the action actually taken, not a cosmetic UI toggle.

Keep native labels, semantic buttons, keyboard interaction, and logical focus
order. Show an inline error beside the responsible field with a useful repair
instruction. Restore focus after dialogs. Give plots a concise text summary
and a route to inspect/export the underlying values. Announce major run state
changes without announcing every simulation step. Avoid automatic focus jumps.

Use a visible **2 px teal focus ring with a 2 px surface-colored gap** on light
surfaces; use mint on ink. The gap prevents the ring from blending into a teal
control. Check focus visibility across all host backgrounds and states. Support
200% text zoom and reflow without losing controls or warnings. Prefer 40–44 px
control hit areas; meet at least the applicable
[WCAG 2.2 target-size requirement](https://www.w3.org/WAI/WCAG22/Understanding/target-size-minimum.html).
Do not copy the prototype's outline reset without a replacement.

Motion should explain a state change. Use 120–180 ms opacity or small transform
transitions; avoid animated layout shifts, springy panels, and pulsing numbers.
Honor `prefers-reduced-motion: reduce` by removing decorative transitions and
autoplay. A trajectory viewer starts paused unless playback is explicitly chosen,
always offers pause and frame stepping, and keeps reduced-motion mode usable
with still frames. Polling or data refresh must not animate every new value.
Any future logo should remain static.

## 9. Delivery and license

This contribution delivers design guidance only. Logo artwork, previews,
favicons, and app icons remain deferred. A future artwork contribution must
provide its own inventory, source files, usage rules, and reproduction checks.

This design document is contributed under the repository's
[MIT license](LICENSE), copyright © 2026 MODENA contributors. It does not bundle
third-party images, font files, or icon packages.
