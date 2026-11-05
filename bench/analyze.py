#!/usr/bin/env python3
"""Observable-specific sampling diagnostics; Python >=3.11, standard library only.

These are finite-series estimates, not ensemble or ergodicity certificates.
The optional coarse burn-in search uses batch means; it is not PyMBAR's method.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import statistics


def _finite_values(values):
    if not isinstance(values, (list, tuple)) or len(values) < 2:
        raise ValueError('a series needs at least two finite numeric samples')
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
           for x in values):
        raise ValueError('series samples must be finite numbers; none are silently dropped')
    return [float(x) for x in values]


def _center(values):
    # Subtract a nearby origin before summing to retain precision for energies
    # with large arbitrary offsets. No estimator can recover bits absent in input.
    origin = values[0]
    shifted = [x - origin for x in values]
    mean_shift = math.fsum(shifted) / len(values)
    centered = [x - mean_shift for x in shifted]
    ss = math.fsum(x * x for x in centered)
    if not math.isfinite(ss):
        raise ValueError('series variance must be finite')
    return origin + mean_shift, centered, ss


def _batch_effective(values):
    _, centered, ss = _center(values)
    if ss == 0.0:
        return 1.0
    size = max(2, math.isqrt(len(values)))
    blocks = [math.fsum(centered[i:i + size]) / size
              for i in range(0, len(values) - size + 1, size)]
    g = max(1.0, size * statistics.variance(blocks) / (ss / len(values)))
    return max(1.0, len(values) / g)


def _select_discard(values):
    # Search at most 21 starts through the first half. The bounded coarse search
    # is a heuristic for this observable, not a claim of global equilibration.
    candidates = sorted({i * (len(values) // 2) // 20 for i in range(21)
                         if len(values) - i * (len(values) // 2) // 20 >= 16})
    if not candidates:
        return 0, [0]
    best = max(candidates, key=lambda start: _batch_effective(values[start:]))
    return best, candidates


def analyze_series(values, *, discard=0, sample_interval=1.0, max_lag=None):
    """Report a mean, autocorrelation-adjusted uncertainty and window diagnostics.

    Statistical inefficiency g = 1+2*sum(rho_k), retaining positive adjacent
    pairs until the first nonpositive pair or max_lag. g is conservatively
    floored at 1, so negative correlation never yields N_eff > N. The reported
    tau_int = g*sample_interval/2 includes the zero-lag half contribution.
    Normal 95% bands assume a stationary series and enough effective samples.
    """
    values = _finite_values(values)
    if isinstance(sample_interval, bool) or not isinstance(sample_interval, (int, float)) \
            or not math.isfinite(sample_interval) or sample_interval <= 0:
        raise ValueError('sample_interval must be finite and positive')
    policy, search = 'explicit', None
    if discard == 'auto':
        discard, candidates = _select_discard(values)
        policy = 'auto_batch_means_heuristic'
        search = {'candidates': candidates, 'maximum_discard_fraction': 0.5,
                  'criterion': 'coarse batch-means estimate of retained effective samples',
                  'global_equilibration_guaranteed': False}
    if isinstance(discard, bool) or not isinstance(discard, int) or discard < 0 \
            or discard > len(values) - 2:
        raise ValueError('discard must leave at least two samples')
    retained = values[discard:]
    n = len(retained)
    if max_lag is not None and (isinstance(max_lag, bool) or not isinstance(max_lag, int) or max_lag < 1):
        raise ValueError('max_lag must be a positive integer')
    cap = min(n - 1, max_lag if max_lag is not None else max(1, min(1000, n // 5)))
    mean, centered, ss = _center(retained)
    variance = ss / (n - 1)
    warnings = []
    g, neff, sem, interval, tau = None, 1.0, None, None, None
    last_lag, truncated = 0, False
    if ss == 0.0:
        warnings.append('constant_series')
    else:
        g = 1.0
        for lag in range(1, cap + 1, 2):
            # Biased finite-sample autocovariance divided by gamma_0 = ss/n.
            pair = math.fsum(centered[i] * centered[i + lag] for i in range(n - lag)) / ss
            end = min(lag + 1, cap)
            if end > lag:
                pair += math.fsum(centered[i] * centered[i + end] for i in range(n - end)) / ss
            last_lag = end
            if pair <= 0.0:
                break
            g += 2.0 * pair
            if end == cap:
                truncated = True
        g = max(1.0, g)
        neff = max(1.0, min(float(n), n / g))
        sem = math.sqrt(variance / neff)
        interval = [mean - 1.959963984540054 * sem, mean + 1.959963984540054 * sem]
        tau = 0.5 * g * sample_interval
        if truncated:
            warnings.append('autocorrelation_window_truncated')
        if neff < 50:
            warnings.append('few_effective_samples')
    block_size = max(1, math.ceil(2.0 * (g or n)))
    blocks = [math.fsum(centered[i:i + block_size]) / block_size
              for i in range(0, n - block_size + 1, block_size)]
    block_sem = (statistics.stdev(blocks) / math.sqrt(len(blocks))) if len(blocks) >= 2 else None
    half = n // 2
    half_difference = (math.fsum(centered[half:]) / (n - half)
                       - math.fsum(centered[:half]) / half)
    result = dict(input_samples=len(values), discarded=discard, samples=n, discard_policy=policy,
                  equilibration_search=search, sample_interval=sample_interval,
                  mean=mean, variance=variance, statistical_inefficiency=g, effective_samples=neff,
                  standard_error=sem, mean_interval_95=interval,
                  interval_method='normal approximation using estimated g; conditional on stationarity',
                  autocorrelation_time=tau, autocorrelation_time_convention='tau_int=g*sample_interval/2',
                  autocorrelation_last_lag=last_lag, max_lag=cap,
                  block_size=block_size, block_count=len(blocks),
                  block_discarded_tail=n - len(blocks) * block_size,
                  block_mean_standard_error=block_sem, half_mean_difference=half_difference,
                  warnings=warnings, validates_ensemble=False)
    # Never export Infinity/NaN as a successful analysis artifact.
    try:
        json.dumps(result, allow_nan=False)
    except (ValueError, OverflowError) as error:
        raise ValueError('sampling statistics must be finite') from error
    return result


def analyze_records(records, *, discard=None, max_lag=None):
    """Analyze recorded replicas separately; never pool correlated frames as replicas."""
    if not isinstance(records, list):
        raise ValueError('timeseries must be a list')
    series, groups = [], {}
    for record in records:
        name, seed, ensemble = record['name'], record['seed'], record['ensemble']
        if not isinstance(name, str) or not name or not isinstance(ensemble, str):
            raise ValueError('each series needs a name and ensemble')
        if isinstance(seed, bool) or not isinstance(seed, int):
            raise ValueError('each replica needs an integer seed')
        interval = record['sample_interval']
        summary = analyze_series(record['values'], discard=record.get('discard', 0) if discard is None else discard,
                                 sample_interval=interval, max_lag=max_lag)
        group = groups.setdefault((name, ensemble), [])
        if any(row['seed'] == seed for row in group):
            raise ValueError(f'duplicate replica seed for {name}')
        if group and group[0]['sample_interval'] != interval:
            raise ValueError(f'replicas of {name} must share a sample interval')
        row = dict(name=name, seed=seed, ensemble=ensemble, sample_interval=interval, summary=summary)
        series.append(row)
        group.append(row)
    replicas = []
    for (name, ensemble), rows in groups.items():
        means = [row['summary']['mean'] for row in rows]
        mean, _, ss = _center(means) if len(means) > 1 else (means[0], [0.0], 0.0)
        replicas.append(dict(name=name, ensemble=ensemble, replicas=len(rows),
                             seeds=[row['seed'] for row in rows],
                             mean_of_replica_means=mean,
                             replica_mean_range=[min(means), max(means)],
                             between_replica_standard_error=math.sqrt(ss / (len(rows) - 1) / len(rows))
                             if len(rows) > 1 else None,
                             band_kind='observed min/max of replica means, not a confidence interval',
                             independence='distinct seeds identify runs; independence and stationarity are assumptions'))
    return {'series': series, 'replicas': replicas}


def _discard_argument(value):
    if value == 'auto':
        return value
    try:
        count = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError('discard must be a nonnegative integer or auto') from error
    if count < 0:
        raise argparse.ArgumentTypeError('discard must be nonnegative')
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='uniformly sampled thermodynamic CSV')
    parser.add_argument('--observable', required=True, help='numeric CSV column')
    parser.add_argument('--discard', type=_discard_argument, default=0, help='sample count or auto (heuristic)')
    parser.add_argument('--sample-interval', type=float, required=True, help='physical time between CSV rows')
    parser.add_argument('--max-lag', type=int)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        if args.output and args.output.resolve() == args.input.resolve():
            raise ValueError('analysis output must differ from the input CSV')
        raw = args.input.read_bytes()
        rows = list(csv.DictReader(raw.decode('utf-8').splitlines()))
        values = [float(row[args.observable]) for row in rows]
        for column in ('step', 'time'):
            if rows and column in rows[0]:
                coordinates = [float(row[column]) for row in rows]
                differences = [b - a for a, b in zip(coordinates, coordinates[1:])]
                if differences and (not all(math.isfinite(x) and x > 0 for x in differences)
                                    or not all(math.isclose(x, differences[0], rel_tol=1e-8, abs_tol=1e-12)
                                               for x in differences)):
                    raise ValueError(f'CSV needs uniform increasing {column} spacing')
                if column == 'time' and differences and not math.isclose(
                        differences[0], args.sample_interval, rel_tol=1e-8, abs_tol=1e-12):
                    raise ValueError('sample_interval disagrees with CSV time spacing')
        summary = analyze_series(values, discard=args.discard, sample_interval=args.sample_interval,
                                 max_lag=args.max_lag)
        result = dict(observable=args.observable, summary=summary, provenance=dict(
            created_utc=datetime.now(timezone.utc).isoformat(), input=str(args.input.resolve()),
            input_sha256=hashlib.sha256(raw).hexdigest(),
            analyzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            time_units='same units as declared sample_interval / CSV time'))
        rendered = json.dumps(result, indent=2, allow_nan=False) + '\n'
        if args.output:
            args.output.write_text(rendered)
        else:
            print(rendered, end='')
    except (ValueError, TypeError, KeyError, OSError, OverflowError) as error:
        parser.error(str(error))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
