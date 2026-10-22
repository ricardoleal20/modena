# Verificación de la entrega científica de MODENA

La implementación completa se verificó sobre el árbol limpio
`e51ff1f47bc72a08296f984604add1f2649df46e`, rama
`ricardoleal20/modena-scientific-hardening`, partiendo de
`b68abd112339aa6970261f470a572da6633bc617`. El commit de cierre posterior sólo
añade documentación y evidencia; no modifica el código verificado. Las cinco
tareas del plan están completas y las revisiones independientes quedaron sin
hallazgos críticos o importantes abiertos.

Cada uno de los 26 comandos finales terminó con código 0 y registra HEAD,
estado limpio, argumentos, duración y stdout/stderr sin filtrar en
[summary.json](evidence/root-verification/summary.json). Son comprobaciones
locales ejecutadas, no una afirmación de que CI remoto haya corrido.

## Comprobaciones ejecutadas

| Comprobación | Resultado | Evidencia |
|---|---|---|
| `cargo test --workspace` | 266 aprobadas, 0 fallidas; 1 doctest previamente ignorada | [stdout](evidence/root-verification/workspace.stdout), [stderr](evidence/root-verification/workspace.stderr) |
| `cargo fmt --all -- --check` | Aprobado | [registro](evidence/root-verification/fmt.json) |
| `cargo clippy --workspace --all-targets -- -D warnings` | Aprobado | [registro](evidence/root-verification/clippy-default.json) |
| Clippy anterior con `--features modena-core/gpu` y con `--features modena-core/mpi` | Ambos aprobados | [GPU](evidence/root-verification/clippy-gpu.json), [MPI](evidence/root-verification/clippy-mpi.json) |
| `cargo build --workspace --release` | CLI, MCP y core compilados | [registro](evidence/root-verification/release.json) |
| `cargo build --release -p modena-core --example scientific_validation` | Compilado | [registro](evidence/root-verification/scientific-build.json) |
| `cargo bench -p modena-bench --no-run` | Los tres ejecutables de benchmark compilados | [stderr](evidence/root-verification/bench-compile.stderr) |
| Python: `-m unittest discover -s python/modena/tests -v`, con `MODENA_BINARY` apuntando al CLI release verificado | 10 aprobadas, 0 omitidas | [registro](evidence/root-verification/python.json), [salida](evidence/root-verification/python.stderr) |
| Python: `-m unittest discover -s bench/tests -v` | 10 aprobadas, 0 omitidas; incluye los fallos de lanzamiento GPU | [salida](evidence/root-verification/tooling.stderr) |
| `bench/validate.py --output <directorio>` | 67/67 observaciones aprobadas, límites y bases independientes registrados | [JSON](evidence/root-verification/scientific/validation.json) |
| `cargo test -p modena-core --features gpu --test gpu_pair -- --nocapture --test-threads=1` | 7 aprobadas con adaptador disponible; 0 omitidas | [salida](evidence/root-verification/gpu-adapter.stdout) |
| `cargo test -p modena-core --features mpi --test mpi_decomp --no-run --message-format=json`; `mpiexec --oversubscribe -n 2 <ejecutable exacto> --nocapture` | Un test distinto aprobado en dos procesos reales; 0 omitidos | [artefacto exacto y hash](evidence/root-verification/mpi-executable.json), [salida](evidence/root-verification/mpi-2rank.stdout) |
| Probes de la auditoría original compilados contra el core release final y ejecutados con `RAYON_NUM_THREADS=1` | 8/8 comprobaciones analíticas aprobadas | [observaciones](evidence/root-verification/original-probe-observations.json), [salida](evidence/root-verification/original-probes.stdout) |
| Los cinco [ejemplos químicos](../../examples/chemical/README.md), cada uno mediante `target/release/modena check -i <config>` y `run -i <config>` | 10/10 acciones aprobadas: molecular NVT, continuación, iones SPME, MC NPT y agua a volumen fijo | [registros de comandos](evidence/root-verification/summary.json) |

MPI seleccionó el único `compiler-artifact` de Cargo cuyo target era
`mpi_decomp`, tipo `test` y feature `mpi`; no eligió el primer archivo de un glob.
Las salidas de ambos procesos representan el mismo test distribuido. Los tests
GPU verifican el backend opcional disponible; no demuestran aceleración general
del CLI ni una campaña de rendimiento GPU.

La compilación de los probes utilizó `rustc --edition=2021 -O -C lto=thin
-C codegen-units=1`, `--extern modena_core=target/release/libmodena_core.rlib` y
`-L dependency=target/release/deps`. Sus ocho observaciones incluyen los momentos
CSVR, Maxwell con/sin restricciones ajenas y coherencia de energía/fuerzas tras
el barostato. La precisión impresa de los probes y las muestras correlacionadas
no constituyen intervalos de incertidumbre formales.

## Entorno y procedencia

- Apple M3 Pro, macOS 26.6.2 arm64, Rust 1.97.1; Open MPI 5.0.10.
- Python 3.12.14 del runtime existente, sin instalar dependencias:
  `/Users/ricardoleal20/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`.
  El Python 3.9 del sistema no cumple el mínimo 3.11 del cliente.
- Los hashes de CLI, MCP, core y exportador numérico están en el resumen; el
  ejecutable MPI tiene su registro separado. No se archivan binarios compilados.
- Los originales del coordinador `verify-command.py`, `run-mpi-verification.py`
  y `collect-verification.py` se guardan bajo `evidence/root-verification/tools`.
  Son scripts históricos de la sesión: esperaban estar junto a un directorio
  `root-verification` y ejecutarse desde la raíz del worktree. Los JSON conservan
  los comandos efectivos para repetir las comprobaciones en otra ruta.
- [manifest.json](evidence/root-verification/manifest.json) identifica por SHA256
  los archivos finales archivados. El snapshot anterior a la corrección final
  queda separado en `evidence/root-verification-pre-final-fix`; sus 252 tests no
  se suman a los 266 finales.

La campaña de rendimiento retenida ejecutó 24 procesos sobre el commit limpio
`11d7683841786699e77d02d7c55c2d8af9cdac40`. Sus binarios, recursos, inputs y
ventana medida son los documentados en [RESULTS.md](RESULTS.md). No se volvió
a medir ni se atribuyen sus hashes al código final: la corrección posterior
afecta DOF/inicialización SPME cargada, aleaciones y catálogo MCP; el workload LJ
medido conserva su política COM3 y sus kernels. En esa campaña de un hilo,
MODENA tardó entre 1.945 y 2.655 veces lo que LAMMPS serial. Para N=20000,
RSS mediano fue 427212800 frente a 22003712 bytes.

## Revisión independiente y cierre

La [revisión completa](evidence/implementation/final-review.md) detectó la
política COM incorrecta para SPME cargado, el pánico de una aleación pequeña
válida y el catálogo MCP incompleto. La
[revisión de las correcciones](evidence/implementation/final-rereview.md)
declaró los tres hallazgos atendidos y aprobó el código, sin nuevos hallazgos
críticos o importantes. Se preservan las reproducciones RED y las regresiones
GREEN en el [reporte del fix](evidence/implementation/final-fix-report.md).

El único hallazgo menor de esa relectura, M2, era un fixture omitido del archivo
Git. Se cerró añadiendo por fuerza el `charges.xyz` original, sin alterar bytes,
con SHA256 `f69364f16fc3a62970cc8782efe69ce9fa1eaa3a5c4eb2351d7e4df5351b14a9`;
los doce miembros del manifiesto del probe coinciden. La
[nota de reproducción](evidence/final-review-probe/README.md) explica la ruta
histórica que conserva su fuente. El [ledger final](evidence/implementation/progress.md)
registra las decisiones y el cierre. El directorio temporal SDD de este plan se
retiró después de archivar y comprobar la evidencia; rama y worktree se conservan.

## Compatibilidad y límites científicos

La conversión correcta es `q*=q_e*sqrt(C/(sigma*epsilon))`. Resultados
electrostáticos anteriores en unidades reales con `sigma!=1` necesitan
recalcularse; la continuación normal rechaza identidades físicas incompatibles.
Unidades reducidas y `sigma=1` no cambian por esta corrección.

Los checkpoints SPME cargados que declaraban COM3 ahora rechazan continuación
normal. `state_only` importa mecánica, reevalúa los grados de libertad y crea
controladores nuevos; no promete la misma trayectoria estocástica. Checkpoints
COM0 elegibles y de pares/Ewald conservan sus reglas de identidad. La nueva
repartición determinista de aleaciones puede cambiar conteos y asignaciones
sembradas respecto al antiguo redondeo independiente.

La entrega respalda el nicho clásico delimitado: pares/iones y moléculas con
modelos existentes a volumen fijo, más MC NPT isotrópico sin restricciones. La
presión con restricciones sigue siendo diagnóstica; NPT rígido/anisótropo,
ensembles líquidos convergidos, validación química amplia, física reactiva o
cuántica y producción HPC general permanecen fuera del alcance aceptado.
Las pruebas y mediciones no establecen completitud ni paridad global con LAMMPS.
