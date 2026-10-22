# Verificación y reproducción

Fecha: 3 de octubre de 2026. Repo: `/Users/ricardoleal20/Documents/Code/Misc/modena`. HEAD: `b68abd112339aa6970261f470a572da6633bc617`. No se modificaron fuentes de producción; los artefactos de auditoría se añadieron bajo este directorio. Los archivos no rastreados que ya existían no se alteraron.

## Entorno

- Apple M3 Pro, 11 CPUs lógicas; macOS 26.6.2 arm64.
- `rustc 1.97.1 (8bab26f4f 2026-07-14)`; toolchain estable configurada por el repo.
- CLI MODENA compilado en release, perfil del workspace: opt-level=3, LTO thin, codegen-units=1.
- `lmp_serial`: 22 Jul 2025, Update 6; MPI STUBS, sin OpenMP habilitado, f64.
- La documentación externa consultada corresponde a las páginas oficiales actuales; se distingue de la versión local del ejecutable.
- Benchmark LJ con carga compartida de escritorio (load average cercano a 4), tres repeticiones por fila. No hubo builds/tests de esta auditoría durante su medición. No se garantiza máquina ociosa ni aislamiento de procesos del usuario.

## Pruebas ejecutadas

| Verificación | Resultado y alcance |
|---|---|
| `rtk cargo test --workspace` | 152 passed, 1 ignored, 19 suites, 176.83 s. La salida resumida la produjo RTK; no se conservó aquí el stdout completo de esta ejecución inicial. |
| `rtk cargo build --release --bin modena` | Compilación exitosa, 24.59 s. |
| Probes de API pública | Bussi, BAOAB con control y NPT/cache con control; se ejecutan correctamente y reproducen defectos físicos. Véanse PROBE-RESULTS.txt y probes.rs. |
| MPI a dos ranks | Un test único pasó en ambos procesos, 11.66 s. No valida migración multinodo. Véase MPI-RESULTS.txt. |
| GPU | Seis tests pasaron, 1.43 s, compilación previa 24.75 s. No aparecen avisos de ausencia de adapter en stdout/stderr. Véase GPU-RESULTS.txt. Cubre el backend LJ y sus tolerancias, no toda física molecular. |
| Cliente Python | Smoke real de check/run contra CLI release: 32 átomos, 100 pasos, 11 filas thermo y verdict Warning; véase PYTHON-SMOKE.txt. Comprueba transporte/parser, no una validación física adicional. La suite pytest del cliente no se ejecutó: pytest no estaba instalado en el Python disponible. |
| Benchmark LJ nuevo | Tres tamaños, 1000 pasos, 200 warmup descartados, tres repeticiones. CLI release contra LAMMPS serial; JSON, raw logs y script conservados. |

Pasar la suite no invalida los contraejemplos: los probes agregan escenarios independientes que no estaban cubiertos por los asserts existentes. No se repararon los defectos en esta auditoría.

## Repetir suite y probes

Desde la raíz del repo:

```bash
rtk cargo test --workspace
rtk proxy python3 docs/audits/2026-10-03-modena-vs-lammps/run_probes.py
```

El launcher compila el core con features por defecto, identifica la biblioteca exacta desde los mensajes estructurados de Cargo, compila probes.rs contra ella y conserva la procedencia en PROBE-BUILD.txt. Los probes utilizan el API público; no insertan ni cambian tests dentro del motor. Las series tienen correlación temporal y un seed fijo; se usan para exhibir errores y contrastarlos con predicciones/control, no para una estimación formal de ESS o errores estándar.

## Repetir benchmark

Ejecutar después de terminar cualquier compilación o test:

```bash
rtk cargo build --release --bin modena
rtk proxy python3 docs/audits/2026-10-03-modena-vs-lammps/benchmark_lj.py
```

El script especifica pasos, misma densidad/modelo/unidades y DOF inicial, ejecuta un proceso separado de 200 pasos para calentamiento del entorno, alterna el orden de engines, comprueba exit codes y número de pasos, y conserva stdout/stderr y drift entre endpoints. El calentamiento no equilibra las trayectorias medidas, que reinician sus decks. Mide wall time de proceso completo. Las velocidades de los dos motores usan RNG diferentes. `modena_default` dispone de varios hilos, mientras que `lammps_serial` usa un core; la igualdad de recursos sólo se aproxima en la comparación `modena_1t`/serial. No se midieron mejores configuraciones aceleradas de LAMMPS ni propiedades de equilibrio.

## Repetir pruebas opcionales

```bash
rtk cargo test -p modena-core --features gpu --test gpu_pair -- --nocapture
rtk proxy cargo test -p modena-core --features mpi --test mpi_decomp --no-run --message-format=json > /tmp/modena-audit-mpi-build.jsonl
```

Para MPI, identificar `executable` del mensaje `compiler-artifact` cuyo `target.name` sea `mpi_decomp` y lanzar ese binario:

```bash
rtk proxy mpiexec --oversubscribe -n 2 /ruta/al/binario/mpi_decomp-hash --nocapture
```

MPI-RESULTS.txt conserva la ruta exacta usada en esta ejecución. No ejecutar dos procesos Cargo como sustituto del binario MPI: deben participar los dos ranks en el mismo conjunto de llamadas colectivas.

## Trazabilidad del análisis

- SCOPE.md: protocolo, alcance y revisión crítica inicial.
- SOURCES.md: 24 entradas primarias verificadas, criterio de inclusión y limitaciones.
- PHYSICS-NOTES.md / PERFORMANCE-NOTES.md: evidencia y líneas del código.
- SYNTHESIS.md / CHECKPOINT-2.md: síntesis y contraste crítico previo al borrador.
- REPORT.md: informe integrado.
- REVIEW.md / CHECKPOINT-3.md: revisión editorial/integridad y crítica final.

Análisis asistido por IA con revisión de agentes. Los resultados ejecutados se distinguen de las inferencias y de campañas históricas.
