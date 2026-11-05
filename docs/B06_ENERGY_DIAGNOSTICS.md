# B06: energy error and sampling diagnostics

B06 implements the diagnostic work proposed in the [core research backlog](research/2026-11-05-core-improvements/phase3_analysis/ACCEPTANCE-BACKLOG.md). It separates integration error, expected energy exchange and observable-specific sampling estimates. It does not change forces, integration schemes or controllers. These diagnostics do not establish material accuracy, ergodicity or parity with another engine.

## Integration statistics

`Simulation::run_streaming` owns a constant-memory `EnergyMonitor`. It observes the initial state and **every integration step**, after thermal and pressure coupling, regardless of `log_every` or trajectory cadence. Statistics cover the current invocation, including a checkpoint continuation; they are not accumulated across checkpoint files. Its time origin is zero at invocation start, while ordinary output keeps the global step/time labels.

All diagnostic energies and times use **reduced units**, including for a real-unit input. Ordinary thermodynamic CSV rows and `RunEvent::Thermo.energies` still use the selected display units. For samples `E_i` at reduced times `t_i`, including `i=0`, the report contains:

| Field | Definition |
| --- | --- |
| `samples`, `duration` | Initial state plus accepted observations; elapsed reduced time |
| `delta_energy` | `E_last - E_0` |
| `delta_per_atom` | `(E_last - E_0) / N` |
| `max_abs_deviation` | `max_i abs(E_i - E_0)` |
| `energy_range` | `max_i E_i - min_i E_i` |
| `rms_fluctuation` | `sqrt(mean((E_i - mean(E))^2))`, population RMS around the sample mean |
| `slope` | OLS energy versus physical reduced time, in reduced energy / reduced time; null before a nonzero time interval |
| `reference_energy` | Initial kinetic energy `KE_0`, null when zero |
| `normalized_change` | `delta_energy / KE_0`, or null |
| `normalized_max_abs_deviation` | `max_abs_deviation / KE_0`, or null |

Centered online moments avoid cancellation from summing large energy offsets. Adding a constant to potential energy preserves the new error metrics up to floating-point representation. No algorithm can recover increments already lost when forming the input energies. Initial kinetic energy is a declared scale, not a universal tolerance; cold initial states produce unavailable normalized values. Nonfinite or overflowing observations abort the run rather than issuing a favorable verdict.

A returning endpoint can hide an excursion. A large excursion can also be bounded rather than a secular trend. The slope is a descriptive estimate, without a significance test or an uncertainty interval; compare excursions, RMS and slopes across timesteps at a fixed physical duration.

## Advisory policy and ensembles

The runner uses these serialized verdicts:

| Verdict | Meaning |
| --- | --- |
| `uncalibrated` | NVE has a normalized metric, but no tolerances were selected; this is the default |
| `healthy`, `warning`, `drifting` | NVE maximum excursion / initial kinetic energy compared with caller-selected thresholds |
| `not_applicable` | Thermal or pressure coupling exchanges physical energy; no NVE energy verdict applies |
| `unavailable` | Normalization or a finite normalized metric is unavailable |
| `insufficient_data` | NVE has no observed integration interval yet |

`drifting` is a legacy label for an error above a threshold; it does not prove a trend. Select and justify tolerances for the fixture and precision before interpreting the classification. For example, this is the **syntax**, not a recommended tolerance for every material:

```toml
[run.energy_error_thresholds]
warning = 0.01
drifting = 0.05
```

Both values must be finite and satisfy `0 < warning < drifting`. The selected policy is echoed in `diagnostics.advisory_thresholds`. It changes diagnostics only, so it is excluded from checkpoint physical identity. NVT/NPT remain `not_applicable` even with this section. Rust users can set `RunSpec.energy_error_thresholds` or call `EnergyMonitor::with_thresholds` directly.

The old number `(E_last - E_0) / max(abs(E_0), 1e-12)` remains in `ConservationReport.relative_drift` and in the CLI's NVE `conservation drift` line for existing parsers. It depends on the arbitrary energy zero and does **not** select the runner's verdict. Legacy `ConservationReport::new`, `pending`, `observe` and `ConservationVerdict::from_drift` retain their endpoint classifier and historical `1e-4/1e-3` thresholds; their applicability is `unspecified` and their diagnostics are null. Use the monitor or runner report for contextual diagnostics.

## CLI, Python and MCP

```bash
modena run --diagnostics-json -i examples/argon_nve.toml
```

The optional flag appends one `energy diagnostics = {JSON}` line containing the whole report. The human summary shows the legacy endpoint number separately from the normalized excursion and time-series statistics.

B01 additionally emits a separate `run metadata = {JSON}` line under the same
flag; the energy-diagnostics line remains unchanged. Python appends optional
`RunResult.metadata` (None with older CLIs). See [B01 migration](B01_EXTERNAL_PROVIDER.md#checkpoints-and-public-api-migration).

Python's `RunResult.energy_diagnostics` contains that report. `RunResult.conservation` exposes the scoped verdict in the existing title-case style, including `Uncalibrated`, `NotApplicable`, `Unavailable` and `InsufficientData`. Configure optional tolerances with `Run(energy_error_thresholds={"warning": 0.01, "drifting": 0.05})`. With an older CLI, the wrapper retries without the new flag only after its exact argument-parser rejection, before a simulation starts; the old verdict remains available and `energy_diagnostics` is null. Physics failures are not retried.

MCP live status consumes the runner's full-rate report, including at sparse log cadence; it does not reconstruct energy errors from displayed thermodynamic samples. `detect_anomaly` emits the legacy `energy_drift` finding only for NVE above explicit tolerances, with the actual metric, policy and diagnostics in its evidence. A failed run generates a `run_failure` finding. A monitor failure reports the offending step and error even when the last logged sample was finite. No findings does not certify sampling or accuracy. MCP's ordinary accumulated sample means remain uncorrected for burn-in/autocorrelation.

Rust API changes: `RunEvent::Thermo` includes `conservation`, `RunEvent::Failed` includes `step`, `RunSpec` includes optional `energy_error_thresholds`, and `ConservationReport` includes applicability/diagnostics. Downstream exhaustive patterns, struct literals and enum matches must accommodate these additions. CSV columns and existing endpoint values are unchanged.

## Observable sampling analysis

`bench/analyze.py` is a Python 3.11+ standard-library postprocessor. The user supplies a uniformly sampled observable and the time between stored rows, in the CSV's time units. Burn-in counts **rows**, not integration steps. For the argon example, rows are spaced by `log_every * dt = 0.05`:

```bash
python3 bench/analyze.py examples/argon_nve_energy.csv \
  --observable potential --discard 20 --sample-interval 0.05 \
  --output target/argon-potential-analysis.json
```

Choose burn-in for the studied protocol. `--discard auto` optionally searches at most 21 starts in the first half of the series, maximizing a coarse batch-means estimate of retained effective samples. This is a bounded observable-specific heuristic, **not** PyMBAR's algorithm or proof of global equilibration. Selection uncertainty is not included in the subsequent mean band.

The retained series reports:

- Mean, unbiased sample variance and statistical inefficiency `g = max(1, 1 + 2 sum rho_k)`, retaining positive adjacent autocorrelation pairs until the first nonpositive pair or the lag cap. The default cap is `min(1000, floor(N/5))`, at least one and at most `N-1`; `--max-lag` can override it.
- `N_eff = N/g`, bounded to `[1,N]`; standard error `sqrt(variance/N_eff)` and an approximate normal 95% band, conditional on stationarity and sufficient effective samples.
- Integrated correlation time with the explicit convention `tau_int = g * sample_interval / 2`, including the zero-lag half contribution.
- Batch-mean standard error with block size `ceil(2g)`, block count and unused tail; first-half versus second-half mean difference.
- Warnings for a capped positive autocorrelation tail or fewer than 50 effective samples. Constant series have `N_eff=1` and unavailable uncertainty/correlation estimates rather than fictitious independent samples.

Finite-series autocorrelation estimates can miss slow modes or transitions; a long apparently stationary window need not explore all relevant states. Block size is a diagnostic choice, not a guarantee of independent blocks. The approximate mean band does not include force-field, discretization, reference-model or burn-in-selection errors. Invalid/nonfinite samples are rejected rather than silently removed. Input and analyzer hashes preserve provenance.

## Scientific harness

```bash
cargo build --locked --release -p modena-core --example scientific_validation
python3 bench/validate.py --output target/b06-scientific-smoke
# Optional overrides: --discard auto, --discard ROWS, --max-lag LAGS, --extended
```

The harness keeps its existing scalar expectations and tolerances. Schema version 2 additionally records kinetic-energy series for CSVR and constrained rigid rotation, separately for seeds 7, 42 and 91, plus full-rate two-body NVE diagnostics at three timesteps with the same duration. To bound artifact size, time series are uniformly thinned to at most about 20,000 rows; their stride, first step, time interval and fixed 1,000-step equilibration discard are recorded. Scalar regression moments still use every production step. Effective sample estimates refer to the stored series.

`validation.stdout` retains raw series and observations. `validation.json` adds per-series analysis and the observed min/max of replica means, mean of replica means and between-replica standard error. The min/max band is not a confidence interval. Runs remain separate: correlated frames are not pooled as independent replicas. Distinct seeds identify runs; statistical independence and stationarity still require justification. Old scalar-only harness payloads remain accepted. The `passed` field checks scalar tolerances and successful finite analysis, not ensemble validity; inspect the warnings separately.

Implementation and executed checks are recorded in [B06 verification](B06_VERIFICATION.md). Earlier research and scientific results retain their original revisions and scope.
