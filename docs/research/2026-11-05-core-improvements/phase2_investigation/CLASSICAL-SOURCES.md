# Fuentes anotadas: MD clásica, algoritmos y validación física

**Estado:** Fase 2, bibliografía dirigida; no contiene síntesis ni recomendaciones de Fase 3. **Consulta:** 2026-11-05 UTC. Se leyeron contenidos abiertos, no solo resultados del buscador. Los casos C1–C3 y las puertas G1–G4 proceden de [CASES.md](CASES.md) y [METHODOLOGY.md](../phase1_scope/METHODOLOGY.md). La correspondencia con casos expresa pertinencia, no validación de MODENA ni rendimiento medido.

Hay **12 registros bibliográficos**: ocho publicaciones de métodos/software con publicación en revista identificada y cuatro fuentes oficiales de implementación. CL01 agrupa tres secciones del mismo manual. Los manuscritos de autores se distinguen del texto editorial publicado; no se presupone identidad completa entre versiones. No se trasladan aceleraciones externas al checkout local. Los metadatos no comprobados se señalan o se omiten. El [ledger de búsqueda](CLASSICAL-SEARCH-LOG.md) conserva consultas, acceso y exclusiones.

**Deduplicación con el corpus paralelo:** CL12 es la misma obra que V1 en [VALIDATION-SOURCES.md](VALIDATION-SOURCES.md), no evidencia independiente adicional. CL11 se solapa temáticamente con V5/V6, con una derivación de otro manual; V7 cubre por separado la aceptación del barostato MC molecular. Este corpus no establece que todo algoritmo NPT necesite virial instantáneo.

## CL01 — LAMMPS: vecinos, geometría y distintos cutoffs

LAMMPS developers. (s. f.). *LAMMPS manual: Neighbor lists; neighbor command; atom_modify command* (stable, git info `30Sep2026 update1`). Recuperado el 5 de noviembre de 2026, de [Neighbor lists](https://docs.lammps.org/stable/Developer_par_neigh.html), [neighbor](https://docs.lammps.org/stable/neighbor.html) y [atom_modify](https://docs.lammps.org/stable/atom_modify.html).

- **Tipo y acceso:** documentación oficial; HTML completo abierto; versión declarada en las tres páginas. DOI: no aplica.
- **Localizadores verificados:** Neighbor lists, descripción inicial y stencils triclínicos; neighbor, Description y Note; atom_modify, `sort` y notas de redondeo. Casos: C1 y controles de vecinos de C2/C3.
- **Anotación:** Documenta listas Verlet paginadas, migración, ordenamiento espacial y stencils distintos para cajas triclínicas. `multi` emplea bins por colecciones cuando los cutoffs difieren; añade coste de preparación y no cubre automáticamente subestilos superpuestos sobre los mismos tipos. El ordenamiento puede cambiar redondeo y trayectorias, y su beneficio depende del sistema. Respalda mecanismos y excepciones, incluidos casos pequeños donde enumerar pares puede convenir. No demuestra que cambiar skin, sorting o `multi` mejore MODENA; tampoco autoriza asumir equivalencia bit a bit tras permutar átomos.

## CL02 — Representación empaquetada de coordenadas en GROMACS

GROMACS developers. (2026). *atomdata.h: Functionality for per-atom data in the nbnxm module* (tag `v2026.3`) [Documentación de interfaz en el repositorio oficial]. [Fuente oficial](https://raw.githubusercontent.com/gromacs/gromacs/v2026.3/src/gromacs/nbnxm/atomdata.h).

- **Tipo y acceso:** fuente oficial de implementación con comentarios de interfaz; archivo completo abierto, sin ejecución ni compilación. DOI: no aplica.
- **Localizadores verificados:** `nbatXYZ`, `nbatXYZQ`, `nbatX4`, `nbatX8`; `c_packX4`, `c_packX8`, `atom_to_x_index`, `getCoordinate`; comentarios de `nbnxn_atomdata_copy_x_to_nbat_x`. Caso: C1.
- **Anotación:** Declara formatos xyz/xyzq y paquetes de cuatro u ocho componentes; la función de índice y la lectura de coordenadas muestran bloques de x, y y z. Describir estos últimos como AoSoA es una interpretación de la disposición observable, no terminología atribuida al artículo de 2013. La interfaz documenta copias desde coordenadas ordinarias al orden de corto alcance, también en GPU. Respalda que existen representaciones especializadas y conversión explícita. No muestra qué formato selecciona cada ejecución, ni mide coste de conversión o beneficio para el buffer actual de MODENA.

## CL03 — Pares de clusters y SIMD

Páll, S., & Hess, B. (2013). A flexible algorithm for calculating pair interactions on SIMD architectures. *Computer Physics Communications, 184*, 2641–2650. [DOI](https://doi.org/10.1016/j.cpc.2013.06.003). [Manuscrito de autores, arXiv:1306.1737v1](https://arxiv.org/pdf/1306.1737).

- **Tipo y acceso:** artículo publicado en revista; manuscrito aceptado de autores consultado íntegramente como PDF de 17 páginas. La ficha de arXiv confirma aceptación y DOI; la web editorial devolvió 403.
- **Localizadores verificados:** §§2–2.5, agrupamiento espacial, máscaras de exclusiones y coste de construcción; implementación descrita en GROMACS 4.6. Casos: C1 y fuerza de C2/C3.
- **Anotación:** Formula clusters espaciales de tamaño fijo y calcula sus interacciones cruzadas para reutilizar datos y adaptar el trabajo al ancho SIMD. La lista por clusters incluye interacciones adicionales que resultan nulas, de modo que menor eficiencia algorítmica puede acompañar mejor utilización del hardware. Describe también exclusiones mediante máscaras y compatibilidad con buffer Verlet. Respalda el mecanismo de bloques y sus compromisos. No demuestra que autovectorizar un bucle existente produzca el mismo resultado, ni que cualquier tamaño de bloque o potencial conserve su ventaja en otra arquitectura.

## CL04 — GROMACS: paso residente y ejecución heterogénea

Páll, S., Zhmurov, A., Bauer, P., Abraham, M., Lundborg, M., Gray, A., Hess, B., & Lindahl, E. (2020). Heterogeneous parallelization and acceleration of molecular dynamics simulations in GROMACS. *The Journal of Chemical Physics, 153*, 134110. [DOI](https://doi.org/10.1063/5.0018516). [Manuscrito arXiv:2006.09167v2](https://arxiv.org/pdf/2006.09167).

- **Tipo y acceso:** artículo publicado en revista; PDF completo del manuscrito de autores, 15 páginas. Identidad editorial contrastada con [PubMed](https://pubmed.ncbi.nlm.nih.gov/33032406/).
- **Localizadores verificados:** §III, reducción de fuerzas por hilo; §IV y Fig. 2, offload de fuerzas frente a iteración completa; §V, clusters y pruning. Casos: C1–C3.
- **Anotación:** Distingue transferir coordenadas y fuerzas cada paso de trasladar también integración, restricciones y propiedad del estado a GPU. El CPU conserva búsqueda y descomposición infrecuentes y puede ejecutar fuerzas especiales. Describe reducción dispersa por hilo, clusters y poda de listas con dependencias explícitas. Respalda residencia como propiedad del flujo completo. También presenta rutas híbridas de PME: FFT pequeñas o hardware legado pueden favorecer trabajo en CPU. No respalda que eliminar toda participación del host sea siempre óptimo ni estima el beneficio para el backend de MODENA.

## CL05 — HOOMD-blue: cobertura GPU y límites por método

Anderson, J. A., Glaser, J., & Glotzer, S. C. (2020). HOOMD-blue: A Python package for high-performance molecular dynamics and hard particle Monte Carlo simulations. *Computational Materials Science, 173*, 109363. [DOI](https://doi.org/10.1016/j.commatsci.2019.109363). [Manuscrito arXiv:1308.5587v2](https://arxiv.org/pdf/1308.5587).

- **Tipo y acceso:** artículo publicado en revista; manuscrito completo abierto, 7 páginas, revisión de 2019 que describe HOOMD-blue v2.6. Año/DOI/volumen corroborados en [documentación oficial de citas 7.2.0](https://hoomd-blue.readthedocs.io/en/v7.2.0/citing.html); la editorial devolvió 403.
- **Localizadores verificados:** §§4–5, backend, operaciones GPU, excepciones y BVH. Caso: C1.
- **Anotación:** Describe un backend C++/CUDA en el que estructuras de datos, fuerzas, integradores y comunicación MPI pueden ejecutarse en GPU. Explicita excepciones de v2.6 y dependencia de densidad, partículas y force field; una interfaz Python no implica transportar el estado en cada operación. Compara listas por celdas y BVH en regímenes de tamaños distintos. Respalda que cobertura de operaciones y selección de vecinos pertenecen al diseño GPU. No prueba cobertura actual de todos los métodos, ni convierte sus umbrales de tamaño o benchmarks históricos en criterios válidos para MODENA.

## CL06 — OpenMM: extensibilidad y precisión

Eastman, P., Swails, J., Chodera, J. D., McGibbon, R. T., Zhao, Y., Beauchamp, K. A., Wang, L.-P., Simmonett, A. C., Harrigan, M. P., Stern, C. D., Wiewiora, R. P., Brooks, B. R., & Pande, V. S. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Computational Biology, 13*(7), e1005659. [Artículo y DOI](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005659).

- **Tipo y acceso:** artículo primario revisado por pares; HTML editorial completo abierto. La copia PMC presentó comprobación de navegador y no se utilizó como evidencia abierta.
- **Localizadores verificados:** Custom forces, Custom integrators, Precision modes y Results/Energy conservation. Casos: C1–C3.
- **Anotación:** Presenta fuerzas derivadas de expresiones energéticas, integradores compilados desde operaciones y una arquitectura extensible por plugins. Contrasta precisión simple, mixta y doble mediante simulaciones NVE replicadas, y advierte que el termostato puede ocultar errores. Respalda separar exactitud de integración, acumulación de energía y cálculo de fuerzas; la comparación de fuerzas/energías entre interfaces también está descrita. Los resultados corresponden a OpenMM 7, sistemas y hardware concretos. No demuestran una precisión universalmente suficiente ni documentan por sí solos todos los detalles de residencia del OpenCL actual.

## CL07 — OpenMM: estado de dispositivo e integración

OpenMM developers. (s. f.). *The OpenCL Platform*. En *OpenMM Developer Guide* (versión mostrada: 8.6). Recuperado el 5 de noviembre de 2026, de [guía oficial](https://docs.openmm.org/latest/developerguide/06_opencl_platform.html).

- **Tipo y acceso:** documentación oficial; HTML completo abierto. URL `latest`, versión registrada al consultar; fecha individual de publicación desconocida. DOI: no aplica.
- **Localizadores verificados:** §§6.1–6.2, OpenCLContext/OpenCLArray; §6.6, Reordering of Particles; §6.7, OpenCLIntegrationUtilities. Casos: C1/C3.
- **Anotación:** El contexto guarda posiciones, velocidades, fuerzas, cola y scratch; las matrices abstraen memoria de dispositivo y sus transferencias. Documenta integración y restricciones sobre desplazamientos, así como listas de pares de bloques de 32 partículas. El reordenamiento interno conserva una relación con los índices del System y exige respetar su definición. Respalda un estado persistente con utilidades comunes y una identidad lógica distinta del orden físico. No afirma ausencia absoluta de transferencias, no cubre aquí paralelismo multidispositivo y no permite atribuir automáticamente su representación a CUDA, Metal o wgpu.

## CL08 — SPME: energía, gradientes y tensor de virial

Essmann, U., Perera, L., Berkowitz, M. L., Darden, T., Lee, H., & Pedersen, L. G. (1995). A smooth particle mesh Ewald method. *The Journal of Chemical Physics, 103*(19), 8577–8593. [DOI](https://doi.org/10.1063/1.470117). [Texto primario consultado](https://biomolmd.org/mw/images/e/e0/Spme.pdf).

- **Tipo y acceso:** artículo primario publicado; copia del PDF editorial completo, 17 páginas. El DOI redirigió a un recurso editorial no accesible; contenido verificado en el PDF, con cabecera de revista, fecha y autores.
- **Localizadores verificados:** §§II–IV, energía/celda, B-splines y fuerzas analíticas; §IV, discusión de energía frente a momento. Caso: C2.
- **Anotación:** Reformula PME mediante B-splines y deriva gradientes de la energía aproximada y contribuciones al tensor de virial. La diferenciación de la asignación de carga es parte de la definición de la fuerza, no un detalle intercambiable sin consecuencias. El propio artículo declara que esta formulación no conserva momento total y describe una corrección de fuerza media en sus pruebas. Respalda consistencia energía–fuerza dentro del funcional discretizado. No garantiza invariancia traslacional exacta respecto de la malla, ni error prefijado para cualquier elección de orden, grid o parámetro de separación.

## CL09 — Comparar asignación, diferenciación y conservación

Deserno, M., & Holm, C. (1998). How to mesh up Ewald sums. I. A theoretical and numerical comparison of various particle mesh routines. *The Journal of Chemical Physics, 109*, 7678. [DOI](https://doi.org/10.1063/1.477414). [Manuscrito arXiv:cond-mat/9807099v1](https://arxiv.org/pdf/cond-mat/9807099).

- **Tipo y acceso:** artículo publicado en revista; manuscrito completo de 18 páginas y ficha de arXiv con referencia editorial/DOI abiertos. Se conserva solo página inicial comprobada; intervalo editorial no verificado.
- **Localizadores verificados:** discusión de diferenciación y conservación, pp. 7–8 del PDF; Investigating the accuracy y comparaciones numéricas. Caso: C2.
- **Anotación:** Compara PME, SPME y P3M descomponiendo asignación de carga, función de influencia y diferenciación. Distingue gradiente analítico de la energía, diferenciación en Fourier y operadores de diferencia; expone condiciones suficientes de conservación de momento para interpolación simétrica y derivadas centradas. También señala que corregir globalmente la fuerza media apenas corrige fluctuaciones locales del error. Respalda comparar coste a precisión equivalente y especificar qué conservación se exige. No establece superioridad universal de una etiqueta de método ni equivalencia de fuerzas cuando se combinan componentes distintos.

## CL10 — Error RMS de P3M y límites de sus estimadores

Deserno, M., & Holm, C. (1998). How to mesh up Ewald sums. II. An accurate error estimate for the particle-particle-particle-mesh algorithm. *The Journal of Chemical Physics, 109*, 7694. [DOI](https://doi.org/10.1063/1.477415). [Manuscrito arXiv:cond-mat/9807100v1](https://arxiv.org/pdf/cond-mat/9807100).

- **Tipo y acceso:** artículo publicado en revista; manuscrito completo de 9 páginas y ficha de arXiv con referencia editorial/DOI abiertos. Intervalo de páginas editorial no comprobado.
- **Localizadores verificados:** Eq. 22 y aproximación analítica posterior; discusión final y Fig. 7 sobre polielectrolitos. Caso: C2.
- **Anotación:** Deriva una estimación del error RMS de fuerzas P3M y una aproximación útil para seleccionar el parámetro de separación junto con malla y orden. La aproximación descrita se restringe a celda cúbica, igual número de nodos por dirección y diferenciación ik. La Fig. 7 conserva evidencia adversa: inhomogeneidad fuerte de cargas produce desviaciones sistemáticas respecto al error real. Respalda control explícito de error acompañado por comprobación del régimen. No constituye una cota universal ni un estimador directamente válido para SPME analítico, celdas triclínicas o cargas netas.

## CL11 — Virial de restricciones y periodicidad

GROMACS developers. (2026). *Some implementation details: Virial from SHAKE*. En *GROMACS 2026.4 Reference Manual*. [Sección oficial](https://manual.gromacs.org/2026.4/reference-manual/details.html#virial-from-shake).

- **Tipo y acceso:** documentación oficial, HTML completo abierto; versión 2026.4. DOI: no aplica.
- **Localizadores verificados:** Single Sum Virial in GROMACS, The intra-molecular shift, Virial from Covalent Bonds y Virial from SHAKE, ecuaciones 484–487. Caso: C3.
- **Anotación:** Explicita que las fuerzas de restricción contribuyen al virial y describe cómo recuperarlas desde la corrección de posiciones para leapfrog. La recuperación posterior tiene un límite geométrico: sin manejo de periodicidad no sirve para el caso general y el virial debe acumularse dentro de SHAKE. El capítulo conserva shifts moleculares y términos de contorno en otras contribuciones. Respalda que satisfacer distancias no basta para presión completa. No proporciona una fórmula portable sin adaptar la convención de signo, integrador, tensor y tratamiento periódico del motor que la implemente.

## CL12 — Validación física de integrador y ensemble

Merz, P. T., & Shirts, M. R. (2018). Testing for physical validity in molecular simulations. *PLOS ONE, 13*(9), e0202764. [Artículo y DOI](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202764).

- **Tipo y acceso:** artículo primario revisado por pares; texto editorial HTML completo abierto; scripts y datos enlazados por los autores, no ejecutados en esta investigación.
- **Localizadores verificados:** Integrator validation; Kinetic energy distribution y separation; Ensemble validation; Discussion. Casos: C1–C3.
- **Anotación:** Presenta controles de fluctuaciones energéticas frente al timestep, distribución cinética, equipartición por subconjuntos y relaciones entre ensembles. Muestra que temperatura media correcta e incluso distribución cinética total adecuada pueden ocultar redistribución errónea entre grados de libertad. El análisis estadístico exige considerar correlación de muestras y el estimador cinético a timestep finito. Respalda pruebas físicas complementarias a comparación puntual de fuerzas. Superar estos controles no prueba ergodicidad, convergencia de observables materiales ni ausencia de errores en todas las combinaciones de parámetros y controladores.

## Límites de este corpus

Es una selección dirigida de métodos y fuentes oficiales en inglés, con trabajos fundacionales deliberadamente antiguos. No es una revisión exhaustiva ni evaluación clínica; no se asignan niveles clínicos de evidencia. No se descargaron datasets ni se ejecutaron benchmarks externos. La comprobación aquí es lectura de contenido y metadatos visibles; no se realizó auditoría de todos los resultados, revisión de cada ecuación ni triangulación automática de índices bibliográficos. Las cuestiones de rendimiento y aceptación del core local siguen sin medición nueva. Documento preparado con asistencia de IA.
