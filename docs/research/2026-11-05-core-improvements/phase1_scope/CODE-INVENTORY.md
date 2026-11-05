# Inventario inicial para delimitar la investigación

Estado: observaciones locales de Fase 1; no constituye la síntesis bibliográfica ni una campaña de rendimiento. Revisión inspeccionada: `f166399070b2d6ebfd91f51e5934c132f04d46c0`. El árbol ya contenía archivos sin seguimiento en `assets/brand/archive/`, `docs/app-design/` y un plan de aplicación; no se modificaron.

## Base que ya existe

La revisión parte de la [cobertura vigente](../../../../docs/VALIDATED_SCOPE.md) y los [resultados del endurecimiento científico](../../../../docs/scientific-hardening/RESULTS.md). La auditoría archivada de octubre describe defectos que fueron corregidos después. No se trasladarán al nuevo estudio como defectos actuales sin reproducción.

| Observación comprobada por inspección | Evidencia local | Pregunta a investigar después de confirmar el alcance |
| --- | --- | --- |
| Verlet con skin, lista de medio par, rebuild por desplazamiento y cambios de caja ya existen. Las celdas triclínicas usan enumeración de todos los pares al reconstruir. | `crates/modena-core/src/neighbor.rs:239`, `:300`, `:330` | ¿Qué estructura de vecinos/localidad mejora el coste manteniendo exhaustividad geométrica? |
| El camino de pares usa especialización LJ, umbral serial y reducción Rayon con buffers `3N` por partición; aún asigna buffers durante la evaluación paralela. | `crates/modena-core/src/force_model.rs:113`, `:236`, `:302` | ¿Conviene reutilizar scratch, cambiar partición, vectorizar por bloques o emplear listas completas por átomo? Medir antes de elegir. |
| Las coordenadas son un buffer contiguo `x0,y0,z0,x1,y1,z1,...`, aunque los comentarios lo llaman SoA. | `crates/modena-core/src/system.rs:268`, `:319` | ¿Ofrece una vista SoA o AoSoA beneficios suficientes para compensar permutaciones y compatibilidad? |
| Las exclusiones están indexadas; SPME reutiliza planes FFT, grids y scratch. | `crates/modena-core/src/neighbor.rs:67`, `crates/modena-core/src/spme.rs:140`, `:280`, `:362` | Estas optimizaciones ya están; buscar el siguiente cuello con perfil de fase, especialmente real-space, spread, gather y FFT. |
| `Virial` contiene solo `w`; no hay tensor de esfuerzo en el contrato del proveedor. | `crates/modena-core/src/force_model.rs:33` | ¿Qué contrato de derivadas de celda y contribuciones de restricciones permite presión física y futuras celdas variables? |
| La presión con restricciones omite impulsos; el test de rotor deja el límite explícito y el validador rechaza reescalado restringido. | `crates/modena-core/tests/sampling_hardening.rs:198`, `crates/modena-core/src/config/validation.rs:256` | ¿Cómo incorporar el esfuerzo de restricciones antes de habilitar NPT molecular? |
| SPME tiene malla de orden cuatro y selección basada en `kspace_cutoff`; rechaza celdas triclínicas. | `crates/modena-core/src/spme.rs:82`, `:234` | ¿Cómo exponer tolerancia de error y validar convergencia de fuerza, virial y traslación? |
| GPU incluye all-pairs y lista por átomo, pero sube posiciones y hace readback bloqueante de resultados en cada cómputo. | `crates/modena-core/src/gpu_comm.rs:261`, `:780`, `:809`, `:834` | ¿Qué límite de backend permite mantener el paso completo en device y reportar observables por lotes? |
| El driver MPI envía todos los átomos propios a cada vecino; reconstruye sistema local, lista y adyacencia en cada evaluación. Empaqueta IDs enteros mediante `f64`. | `crates/modena-core/src/mpi_comm.rs:104`, `:135`, `:170`, `:203` | ¿Cómo pasar a estado local persistente, halo geométrico, migración y payload tipado con validación multirank? |
| `ForceModel` es extensible, pero `ModelSpec` y el runner construyen potenciales de pares más términos moleculares/electrostáticos. No se encontró un proveedor de inferencia ML o un adaptador ASE/QE. | `crates/modena-core/src/force_model.rs:39`, `crates/modena-core/src/io/input.rs:356`, `crates/modena-core/src/runner.rs:1006` | ¿Cuál es la ruta mínima comprobable para integrar modelos ML y referencias DFT sin confundir sus dominios? |
| `System.species` es índice de modelo; `elements` se documenta cosmético. No existe campo explícito de números atómicos. | `crates/modena-core/src/system.rs:264` | ¿Cómo distinguir identidad química, tipos de force field y metadatos electrónicos requeridos por modelos? |
| El informe de conservación compara energía inicial/final con umbrales fijos calibrados en LJ. | `crates/modena-core/src/conservation.rs:19`, `:72` | ¿Qué diagnóstico incluye error acotado, pendiente temporal, escala de energía y el ensemble pertinente? |

Cada fila es una observación o pregunta de diseño. No establece el porcentaje de runtime, la magnitud de una mejora futura ni equivalencia frente a otro motor. Las fuentes oficiales y los trabajos primarios se incorporarán en Fase 2, tras la confirmación requerida por la habilidad.

Elaborado con asistencia de IA mediante inspección de código y documentación local.
