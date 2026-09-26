# Acta — Comité de rendimiento vs LAMMPS (consenso final)

Fecha: 2026-09-26. Convocado por el jefe de investigación (orquestador) bajo
`/goal`: iterar rondas de presentación → discusión cruzada → auditoría hasta
consenso grupal explícito **sin dudas remanentes**.

## Votación final (4/4 YEA)

| Expositor | Dominio | Voto | Correcciones/registros al acta |
|---|---|---|---|
| a7c3abfd (kernel) | kernel/micro | **YEA** | proyección single-core corregida a ~1.4–1.7× (1t: K×B×C, no A); iónico ~8–9× (no 2×); retracta superlinearidad de cache; post-B flag A/B en G2 |
| afb96abe (ewald/scale) | electrostatics | **YEA** | conteo real 16.9M sincos/step (no 33.8M) → ~10.1 ns/sincos (~30 ciclos, latency-bound); Ewald-b → ~8.6×; tablas trig por-eje (recomendadas) → ~2×; columna `rebuild_count` en G4 |
| ae5cc3f3 (gpu) | heterogéneo/Apple-GPU | **YEA** | flag ×1000 resuelto por matemática: 341.99 s/2000 = **171 ms/step** (no µs) = coincide al 3.4% con el estimado 177 ms de R1 — no hay discrepancia; fila GPU de Stage 0 = provisional; "GPU small-N" (typo "GUI"); hardening OpenCL `neigh no` opcional no-bloqueante; carta del nodo Stage 5 incluye rama cell-list (3.5–5.5 d) |
| ad1a9b44 (rig) | compilador/medición | **YEA** | condicional: G1 verifica byte-identidad post-cambio (ya en el pipeline); flags "0%" solo arm64, "medido pre-monomorfización" o A/B post-B |

## Hechos medidos (auditoría, main e9ef83a, 2026-09-26)

1. **Gap single-core reconciliado y plano**: 2.06× (864) / 2.10× (2048) a 1
   hilo; 2.55× / 4.64× con default-thread. La superlinearidad de cache no
   existe a 1t+skin — era artefacto del bench de kernel (rebuild 52–68%).
2. **Default-threading se autodaña**: rayon-11t es 31% peor que 1t a 864.
3. **Gap iónico medido ~34×**: LAMMPS 9.93 s vs MODENA 341.99 s (mismo deck,
   216 átomos, 2000 steps) = 5.0 vs 171 ms/step.
4. **Bit-determinismo**: byte-idéntico a 1t; no determinista a default threads.
5. **Flags**: 0% en arm64 (reales ±3%); lto=thin + codegen-units=1 ya activos.
6. Ewald: 16.9M sincos/step × ~10.1 ns ≈ 171 ms/step (todo el step es
   transcendentales).

## Decisiones convergidas (sin dudas)

- **Rechazados con matemática**: f32-pair en CPU (FD-error ~600× sobre
  tolerancia del contrato); tabulación erfc (C⁰-lin = fuga de deriva); overlap
  CPU↔GPU (cero beneficio en LJ); f32-CPU mixed precision.
- **Diferidos post-paper con precondiciones registradas**: SPME (virial
  strain-form en grilla + protocolo matched-κ/precisión; O(N log N));
  cell-list GPU (spike de gather-throughput en G3); cluster-pair GROMACS
  (internal); locality sort H; ESP.
- **El claim GPU** (único "win" opcional) queda ligado a G3: drift f32 vs
  clasificador Healthy, timing fence-sync, Kahan-bajo-naga; wording exacto con
  caveat `neigh no` de lammps#3531 + disclosure f32 en la misma tabla.

## Pipeline final (gates con criterio numérico)

| Stage | Contenido | Esfuerzo | Gate (criterio) |
|---|---|---|---|
| 0 | rig `bench.sh` + `profile.sh`; regenerar TODOS los números commiteados (incl. fila iónica; fila GPU provisional marcada); stamps (SHA/rustc/flags/threads/lmp 22Jul2025-update6) | 0.5–1 d | G0: CSV fresco firmado; tabla corregida (16.9M/10.1 ns) |
| 1 | thread plumbing: serial bajo umbral + `with_min_len` + scratch reuse | 1–2 d | G1: byte-ids 1t; ≥1.25× (864) y ≥2× (2048) default; drift Healthy; embedding MCP/GUI |
| 2 | kernel PR: B monomorfizar + C fusionar (gated) + K stencil-half (dims≤3) | 1–2 d | G2: bench prebuilt (resuelve 2048-anomaly); rebuild≤15%; byte-compare B; post-B flag A/B |
| 3 | GPU-e (Kahan + buffers persistentes) + ladder GPU-d 3 motores + drift | 1–1.5 d | G3: drift f32 vs Healthy; timing fence; Kahan-naga; spike gather; buffer edge-cases; overhead floor |
| 4 | Ewald-b (±k + phase reuse + kvec cache) como correctness/higiene | 1–2 d | G4: ≥1.8× iónico; gate 1e-10 con **agrupamiento A/B fijado**; `rebuild_count`; retener el camino original como referencia |
| 5 | freeze + nodo de decisión del paper | — | win (si G3 gana) vs gap (1t 1.4–1.7× tras K×B×C) — la rama cell-list pre-paper es admisible |

## Erratas aritméticas que no deben sobrevivir en el paper

- "33.8M sincos/step" → **16.9M**; "~5 ns/sincos" → **~10.1 ns**.
- Proyecciones: single-core "1.2–1.5×" → **~1.4–1.7× (1t, K×B×C)**; iónico
  "~2×" → **~8.6× (Ewald-b)** o **~2× solo con tablas trig** (recomendado,
  +0.5–1 d) / **O(N log N) con SPME (post-paper)**.
- CSV 9.70 s stale → regenerar; fila GPU provisional etiquetada.

## Registro de dudas final

Todas las dudas declaradas en R2 (4×~60, clusters) tienen estado: **resueltas
por medición**, **resueltas por scoping**, o **gates con criterio numérico**
(G1–G5). Cero dudas de investigación abiertas. La condición del `/goal`
(consenso tras auditoría, sin dudas) se declara cumplida a las 4/4 YEA.