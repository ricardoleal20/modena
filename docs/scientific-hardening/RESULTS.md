# MODENA: resultados del endurecimiento científico

La entrega es una base de dinámica molecular clásica para pares, iones y modelos
moleculares existentes a volumen fijo. No demuestra paridad general con LAMMPS,
completitud química, física cuántica/reactiva/muchos cuerpos ni producción HPC.
Las cinco tareas y la corrección de la revisión final están terminadas. La
[verificación final](VERIFY.md) registra 266 pruebas Rust aprobadas, una doctest
ya ignorada, 20 pruebas Python/herramientas, 67 observaciones numéricas, GPU con
adaptador y MPI con dos procesos reales; ninguna comprobación falló.

## Inventario implementado

- CSVR con momentos canónicos; BAOAB restringido con ruido tangente correcto;
  propagación rígida que conserva rotación libre; grados de libertad físicos
  explícitos, COM conservado y restricciones; errores al fallar la proyección.
- Fuerzas actualizadas tras cambios de caja/intercambio; tablas Hermite conjuntas
  con `dU/dr=-F`; viriales Ewald/SPME por derivada de deformación; exclusiones
  indexadas y espacios de trabajo SPME reutilizables. Pares/Ewald directos admiten
  celdas triclínicas; SPME sigue limitado a cajas ortorrómbicas.
- NPT canónico por Monte Carlo isotrópico sin restricciones. Berendsen sirve
  para equilibrar; el pistón heredado es experimental. Acoplamiento anisótropo
  completo y NPT rígido/restringido se rechazan. La presión molecular omite el
  esfuerzo de los impulsos de restricción y sigue siendo diagnóstica.
- Checkpoints con estado estocástico/extendido e identidad física; importación
  molecular tipada y políticas separadas LJ/Coulomb según el grafo; validación
  estricta; cliente Python y MCP con unidades/metadatos/promedios reales.
- Herramientas de benchmark que fallan ante errores y conservan entradas,
  stdout/stderr, código de salida, pasos, precisión, hilos, semilla, hardware,
  hashes y tiempos/RSS. Criterion separa listas cacheadas y reconstrucción e
  incluye escaleras LJ, iónicas y moleculares SPME. CI añade humo numérico acotado.

## Corrección dimensional y compatibilidad

La resolución del ledger amplió Task4 para corregir `q*=q*sqrt(C/(sigma*epsilon))`.
El factor anterior `sqrt(C*sigma/epsilon)` multiplicaba la energía Coulomb física
por sigma²: no era una convención. Los casos reducidos o sigma1 no cambian.
Trayectorias electrostáticas reales anteriores con sigma!=1 tienen otro
Hamiltoniano; la continuación normal verifica factor/cargas y rechaza identidades
incompatibles. No se promete migración histórica. Las energías/timings reales
iónicos antiguos no validan el modelo corregido. [Referencia de unidades](UNITS.md).
La afirmación antigua de WATER_MODEL de que la discrepancia no podía ser un
error de unidades queda retirada; relajación estructural tampoco explica por sí
misma deriva NVE. Los originales quedan archivados sin modificar.

## Corrección final de COM y mezclas pequeñas

La malla finita SPME tiene energía dependiente de una traslación continua común
respecto a la cuadrícula. Sus fuerzas analíticas conservan `F=-grad U`; la suma
puede ser distinta de cero, incluso en una celda neutra. Con cargas no nulas,
VV/CSVR/NHC y MC usan COM0 y `3N-restricciones`; las velocidades iniciales incluyen
los modos COM. Langevin mantiene COM0. Pares/Ewald directo y SPME sin cargas usan
COM3 sólo cuando el protocolo conserva el momento cero inicial/importado.
Agua rígida cargada con SPME usa `6Nw` DOF, sin resta COM3. La presión de
restricciones sigue siendo diagnóstica.

La continuación v2 rechaza metadatos COM3 para SPME cargado. `state_only` conserva
la mecánica importada y reevalúa COM según modelo, controlador, masas y momento,
también después de proyectar restricciones. No cambió el formato de checkpoint
ni la firma `modena-physical-v2`: la comparación existente de metadatos físicos
ahora se hace contra una política elegible. Checkpoints de pares/Ewald y de
SPME/Langevin COM0 siguen siendo elegibles con la misma identidad. No se promete
una continuación idéntica de los antiguos estados SPME COM3.

Las poblaciones `random_alloy` usan cuotas piso y mayores residuos (empates por
orden de especies) antes del shuffle sembrado. Se admiten mezclas pequeñas y
fracciones cero; el total es exactamente N. Algunos conteos/asignaciones difieren
de la regla antigua de redondeo independiente; la semilla sigue siendo repetible
con la regla actual. Entradas y capacidades inválidas regresan `Result::Err` en
el constructor de aleaciones. `list_potentials` deriva las cinco familias y sus
campos del mismo esquema MCP que usan las solicitudes.

La [evidencia de la corrección final](evidence/implementation/final-fix-report.md)
conserva las reproducciones RED del revisor, sus hashes y los tests RED/GREEN.
Los 67 resultados y 24 procesos LJ históricos permanecen intactos, con sus
commits/hashes originales. El force-energy SPME no cambió: las filas estáticas
SPME siguen describiendo ese modelo. La inicialización/DOF de trayectorias
SPME cargadas sí cambió; no se presenta la campaña histórica como una nueva
validación de esas trayectorias. Los workloads LJ con COM3 no cambian por esta
política, ni por la repartición de especies de aleaciones que no usan.

## Evidencia numérica observada

El [JSON de la campaña de Task5](evidence/scientific-release/validation.json) contiene
67 observaciones con valor esperado, tolerancia y base independiente; todas
pasan. Semillas 7/42/91, 40000 muestras por semilla, descarte 1000, f64, un hilo.
Son regresiones acotadas y momentos correlacionados, no incertidumbres formales
ni observables de líquido convergidos.
La [repetición final sobre e51ff1f](evidence/root-verification/scientific/validation.json)
también pasó las 67 observaciones con los mismos límites; el cuadro siguiente
conserva las mediciones de Task5 y su procedencia.

| Prueba | Medición | Expectativa independiente y límite |
|---|---|---|
| CSVR K, f12 T1 | medias 6.00088/5.98484/5.98738; varianzas 5.99023/6.05303/5.98428 | Gamma: media=varianza=6; tolerancias .085/.25 |
| BAOAB átomo libre, con/sin dímero restringido | 3 semillas, 3 ejes, media y varianza registradas | Maxwell media0, varianza T/m=.5; límites .025/.035 |
| Rotación rígida, gamma2 dt.01 | medias 1.53338/1.50309/1.51553; varianzas 1.53344/1.48569/1.50677 | Gamma(3/2,1); bandas OU conservadoras .22652/.65956 |
| Restricciones rígidas | residual máximo <1e-10 en las 3 semillas | distancias analíticas sqrt(3), límite1e-10 |
| NVE par LJ, tiempo físico2 | bandas4.53572e-5/1.13386e-5/2.83474e-6 en dt.004/.002/.001 | error O(dt²), ratios4.00024/3.99987 dentro4±.3 |
| SPME malla15³/20³, caja9, kappa.8, rc4.2 | energía -.39908485/-.39806559; error máximo de fuerza .00233876/.000594538 | Ewald directo -.39773102; error energía/fuerza <.005 |
| Virial por deformación Ewald/SPME | valores y FD de caja+posiciones en JSON | diferencia absoluta<1e-7, h1e-6 |
| Seam estocástico JSON BAOAB13+17 | error máximo0 | trayectoria continua idéntica, tolerancia0 |

La malla inicial de kcut3.2 **falló** el límite .005: error de energía .00856 y
fuerza .01242. [Salida fallida retenida](evidence/scientific-smoke/validation.json).
La aceptación utiliza kcut4.8/6.4 (mallas15³/20³); esto evidencia control de
resolución y no se oculta como defecto resuelto de producción.

Los [probes originales ejecutados durante Task5](evidence/raw/original-probes-current.json)
miden T CSVR=.99933785 y varianza=.16712408 (expectativas1 y1/6), varianza Maxwell
libre1.00073504 con/sin restricciones y diferencias de caché tras barostato0 a
la precisión impresa. La compilación directa requirió `-C lto=thin` para leer el
bitcode Rust1.97 con el linker Apple; se conserva el intento fallido y el éxito.
La [ejecución final de los mismos probes](evidence/root-verification/original-probe-observations.json)
pasó las ocho comprobaciones analíticas sobre el core de e51ff1f. Las etiquetas
antiguas de defecto en el stdout original se conservan como texto histórico.

## Reproducción y comprobaciones de Task5

```sh
cargo build --release --bin modena
cargo build --release -p modena-core --example scientific_validation
python3 bench/validate.py --output target/scientific-smoke
python3 bench/validate.py --extended --output target/scientific-extended
python3 -m unittest discover -s bench/tests -v
cargo bench -p modena-bench --no-run
cargo bench -p modena-bench --bench force_kernel -- --test
cargo bench -p modena-bench --bench integrator_step -- --test
```

Se usó Python3.12.14 del runtime Codex sin instalar dependencias; la ruta exacta
está en la procedencia. La campaña inicial verificó 5 pruebas ejecutables RED/GREEN
para salida fallida, metadatos/pasos/render numérico, ausencia de ficheros de
trayectoria/CSV, resultados numéricos no finitos/fuera de tolerancia y binario
faltante. Criterion ejecutó18 casos de fuerza y9 de integración con éxito.
Los logs Criterion incluyen el aviso no fatal «Gnuplot not found, using plotters
backend»: se utilizó Plotters para gráficos; no es un fallo de física o compilación.
Los bytes originales se conservan. El benchmark de integración incluye la
destrucción/deallocación de buffers al finalizar el closure, además del paso;
no hay un claim numérico Criterion publicado que dependa de esta ventana.
La ronda de revisión añadió 5 fixtures GPU ejecutables RED/GREEN: binario ausente,
proceso fallido, errores del sistema al lanzar binario/timer y éxito sin timer.
Los 5 pasan; la evidencia de fallos conserva hash nulo cuando no está disponible,
código de salida, logs y JSON. Sin timer, CPU/RSS son null. No se repitieron
campañas CPU/GPU ni se cambiaron los kernels medidos.
Compilar benches, fmt y clippy de los objetivos modificados pasaron. El ejemplo
GPU compiló con feature gpu y ahora declara COM3 antes de inicializar velocidades.
Task5 no duplicó campañas generales: las comprobaciones frescas de workspace,
Python completo, GPU con adaptador y MPI2 corresponden al coordinador después
de la revisión independiente. Todas pasaron sobre e51ff1f y están archivadas en
[VERIFY.md](VERIFY.md) y [el resumen de comandos](evidence/root-verification/summary.json).

## Timings actuales

La [campaña actual](evidence/matched-lj/benchmark.json) midió el commit limpio
`11d7683841786699e77d02d7c55c2d8af9cdac40`: Apple M3 Pro,11 CPUs lógicas,
Rust1.97.1, Python3.12.14, LAMMPS22Jul2025 Update6 serial, un hilo/rank por motor.
24/24 procesos terminaron con código0 (6 warmup+18 medidos). Se pausaron builds,
tests y trabajo pesado del equipo durante las mediciones; el escritorio compartido
no se convierte por eso en hardware dedicado. Whole-process: inicialización más
integración, sin archivos de trayectoria/termo; stdout capturado. Warmup200 pasos
en proceso separado,1000 pasos medidos,3 repeticiones, semilla87287. LJ reducido
shifted-force sigma=epsilon=m1,rc2.5,skin.3,a1.6796,dt.005,T1,DOF3N-3; los RNGs
de los motores difieren. Hashes de binarios/entradas, comando y logs completos se
conservan. Las razones actuales están entre condiciones comparables; no implican
aceleración frente al baseline histórico.

| N | MODENA1t ms/paso (mediana) | LMP serial ms/paso (mediana) | Razón MODENA/LMP | RSS MODENA/LMP bytes (mediana) |
|---:|---:|---:|---:|---:|
|864|.601879|.309412|1.945|3883008 /12140544|
|2048|1.738660|.736351|2.361|5111808 /12500992|
|20000|17.326401|6.525943|2.655|427212800 /22003712|

MODENA sigue más lento en esta campaña. El RSS N20000 también es una brecha
medida importante (~427MB frente a22MB); no se atribuye a una causa sin perfilar.
Los valores CPU/tiempos brutos están en el JSON; memoria no disponible sería null,
no un cero fabricado. Deriva relativa MODENA/LMP: N864: 3.171e-5/2.80594e-5,
N2048: 1.919e-5/1.73939e-5, N20000: 1.356e-5/1.45231e-5; no se atribuye una diferencia
en microestados a superioridad del algoritmo. [Resumen numérico](evidence/matched-lj/summary.json).
El ejemplo GPU usa otra ventana de tiempo y f32: se compiló, pero Task5 no midió
una campaña GPU ni presenta una aceleración GPU/CLI con scopes diferentes.

La [base histórica](evidence/original-audit/LJ-RESULTS.json) preserva medianas
LMP/MODENA1t en ms/paso: N864 .2586267/.3877270, N2048 .6465104/1.2711487,
N20000 6.8488693/13.3516402. Las mediciones default multicore son condiciones
diferentes y no victorias con iguales recursos. La inicialización histórica
usaba otras DOF; no se deriva una aceleración antes/después de microestados,
DOF o cargas de escritorio distintos. [Metodología](../../bench/BENCHMARKS.md).

## Brechas competitivas restantes

Faltan ensembles líquidos convergidos, validación física molecular de presión,
NPT rígido/anisótropo completo, matrices amplias de parámetros químicos y
estática real-unit iónica/agua comparada con LAMMPS tras la corrección dimensional.
No hay nuevas funciones reactivas, cuánticas, muchos cuerpos o producción HPC.
GPU sigue siendo offload LJ f32 con readback al host; domain MPI rechaza políticas
especiales y cambios de caja. No son aceleración general del CLI/MCP. Los tiempos
no establecen superioridad numérica: ambos kernels shifted-force son consistentes
y Velocity-Verlet tiene error acotado, no conservación exacta.

## Archivo durable

[evidence/original-audit](evidence/original-audit/REPORT.md) contiene fuentes,
probes, baseline y logs de la auditoría; [implementation](evidence/implementation/progress.md)
conserva reports, revisiones, ruling y logs de tareas1–4. Las narrativas antiguas
de benchmark/validación/agua están en legacy-docs. La evidencia describe su
commit/hash/dirty state: outputs previos al
commit de herramientas no se presentan como campañas de un árbol limpio.
También se conservan las revisiones de Task5 y final, su relectura de los fixes,
el ledger de cierre y todos los logs de la verificación final. El fixture XYZ
omitido inicialmente por una regla de Git se añadió con su hash original; la
[nota del probe](evidence/final-review-probe/README.md) permite preparar su ruta
histórica sin modificar el código ni las salidas archivadas.
