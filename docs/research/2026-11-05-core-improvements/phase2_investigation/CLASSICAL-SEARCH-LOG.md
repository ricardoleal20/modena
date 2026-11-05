# Ledger de búsqueda: MD clásica y algoritmos

**Fase:** 2, investigación bibliográfica. **Fecha:** 2026-11-05 UTC. **Última comprobación de reloj durante la sesión:** 05:30:53 UTC. La búsqueda empezó tras la confirmación de Fase 1 comunicada por el orquestador. Se consultaron [METHODOLOGY.md](../phase1_scope/METHODOLOGY.md), [CODE-INVENTORY.md](../phase1_scope/CODE-INVENTORY.md) y [CASES.md](CASES.md) antes de seleccionar evidencia.

## Parámetros y criterios aplicados

- **Medio:** buscador web y apertura directa de fuentes primarias; sitios oficiales LAMMPS, GROMACS, OpenMM y HOOMD-blue; artículos editoriales y manuscritos de autores depositados en arXiv. Sin bases bibliográficas de suscripción ni búsqueda exhaustiva de Scopus/Web of Science.
- **Temas:** vecinos/localidad/múltiples cutoffs; clusters y empaquetado; residencia GPU; SPME y estimación de error; virial de restricciones y validación de ensemble.
- **Ventana:** métodos fundacionales sin corte inferior y documentación disponible al 2026-11-05. Idioma de los materiales leídos: inglés. No se seleccionó por magnitud de aceleración.
- **Inclusión:** texto primario o documentación oficial con contenido abierto que sostenga un mecanismo, condición matemática, capacidad o límite relacionado con C1–C3.
- **Exclusión:** blogs, foros, Wikipedia y revisiones secundarias como evidencia; papers tangenciales; snippets como prueba final; comparaciones sin identidad de modelo/precisión como aceleraciones transferibles.
- **Cribado:** pertinencia por título/resumen, después apertura de contenido y lectura de las secciones anotadas. Duplicación resuelta por título, autores, DOI y relación explícita entre manuscrito y revista. Sin `literature_corpus[]` de usuario recibido ni consultas de deduplicación a Semantic Scholar/OpenAlex/Crossref; no se atribuye una comprobación negativa a esos índices.

No se reporta flujo PRISMA. El buscador mostró selecciones limitadas de hits y no un universo exportable ni totales reproducibles por consulta. La lista siguiente registra hits **observados**, no su número total, y no pretende incluir cada resultado irrelevante mostrado. Los lotes de múltiples consultas pueden mezclar sus resultados: la correspondencia se registra por tema y no como conteo individual.

## Consultas y resultados observados

| ID | Consulta enviada | Hits pertinentes observados y decisión |
| --- | --- | --- |
| Q01 | `site.docs.lammps.org Developer_par_neigh atom_modify sort neighbor multi` | Manual oficial de Neighbor lists, atom_modify, neighbor y neigh_modify. Abiertas las tres primeras, incorporadas como secciones de un único manual en CL01; no se tomó el snippet como prueba. |
| Q02 | `Pall Hess 2013 flexible algorithm calculating pair interactions SIMD full paper` | Editorial CPC y arXiv:1306.1737; PDF y ficha arXiv abiertos, CL03. Editorial devolvió 403. |
| Q03 | `Pall 2020 heterogeneous parallelization acceleration molecular dynamics GROMACS GPU resident` | PubMed 33032406 y arXiv:2006.09167; PDF y ficha PubMed abiertos, CL04. Copias de USC y suplementario Zenodo descubiertos; no usados como evidencia adicional. |
| Q04 | `HOOMD blue Anderson Glaser Glotzer 2020 general purpose particle simulation arxiv` | arXiv:1308.5587 y 1412.3387. Elegido el manuscrito revisado de 1308.5587; el segundo no se anotó separadamente porque el foco de este registro era cobertura del backend, no escalado MPI. |
| Q05 | `site:docs.openmm.org developers guide CUDA OpenCL platforms GPU data CPU GPU transfers` | El lote no aportó un hit oficial útil atribuido a esta cadena; se abrió después la URL oficial de OpenCL directamente. No se inventa un total cero de resultados. |
| Q06 | `Smooth particle mesh Ewald method Essmann 1995 full text interpolation derivatives` | PDF del artículo en biomolmd.org y copias secundarias ResearchGate/Scribd. Abierto y leído el PDF editorial primario; CL08. Las páginas secundarias no se usaron como sustento. |
| Q07 | `Deserno Holm How to mesh up Ewald sums II accurate error estimates particle mesh methods full text` | arXiv:cond-mat/9807100, enlace al trabajo I y PDF alojado en matsci.org. Abiertos manuscritos y fichas arXiv de I/II; CL09/CL10. No se necesitó la copia del foro. |
| Q08 | `Merz Shirts testing for physical validity molecular simulation 2018 ensemble constraints virial` | PLOS ONE 0202764, PDF editorial y PMC6126824. Abierto el HTML editorial completo; CL12. |
| Q09 | `OpenMM 7 rapid development high performance algorithms molecular dynamics 2017 PMC` | PMC5549999, PLOS Computational Biology 1005659, PubMed y eScholarship. PMC presentó comprobación de navegador al abrir; se abrió PLOS y se leyó el artículo completo, CL06. |
| Q10 | `site:glotzerlab.engin.umich.edu HOOMD-blue 2020 109363 Anderson Glaser` | En el lote se observó un PDF antiguo v1; no se confundió con el manuscrito v2 consultado. Se hizo Q13 para confirmar identidad editorial. |
| Q11 | `site:manual.gromacs.org current reference manual pressure virial constraint` | Manuales oficiales PDF; se refinó con Q14 y se abrió el HTML versionado de implementación, CL11. |
| Q12 | `site:manual.gromacs.org AoSoA cluster layout` | No apareció un hit primario útil en el lote observado; se refinó con Q15–Q18. No se atribuyó AoSoA al texto de Páll sin comprobar la implementación. |
| Q13 | `"HOOMD-blue" "109363"` | Editorial, documentación oficial de citas y registro ORNL. Abierta la documentación oficial para año, volumen y DOI; PDF arXiv completo para contenido, CL05. Editorial devolvió 403. |
| Q14 | `site:manual.gromacs.org "constraint" "virial" "Pressure"`; `site:manual.gromacs.org "Constraint virial"` | Some implementation details de varias versiones, Molecular Dynamics 2021 y notas 2021.6. Se abrió details 2026.4; CL11. La nota de bug histórico solo fue un hit de búsqueda y no se usa como defecto actual ni como evidencia abierta. |
| Q15 | `site:manual.gromacs.org "structure-of-arrays"` | El lote no devolvió evidencia oficial útil para la afirmación buscada. Continuó lectura de fuente oficial. |
| Q16 | `site:github.com/gromacs/gromacs "AoSoA"` | El lote mostró resultados ajenos al proyecto y demostraciones genéricas de AoSoA. Excluidos por falta de pertinencia/primariedad para GROMACS. |
| Q17 | `site:github.com/gromacs/gromacs "X4" "cluster" "layout"` | Sin evidencia oficial útil en la selección observada. No se usa como ausencia de implementación. |
| Q18 | `GROMACS nbatXYZQ nbatX4 nbnxn_atomdata_t coordinates xstride` | Espejos Fossies revelaron el nombre exacto de la interfaz; solo descubrimiento. Se abrió directamente el archivo del repositorio oficial en tag v2026.3; CL02. |

Además se usó búsqueda dentro de textos abiertos para `layout`, `stored`, `structure`, `AoS`, `resident`, `memory`, `entire`, `x0`, `CustomIntegrator`, `CPU and GPU`, `momentum`, `cubic`, `constraint`, `constraints`, `conservation` y `Drift`. Hubo búsquedas sin coincidencia y una cadena de prueba irrelevante (`xxxxx`); ninguna ausencia se convirtió en afirmación científica. Las aperturas con posición dentro del documento sirvieron para revisar el contexto de cada coincidencia.

## Registro de acceso y versiones

| Fuente final | Evidencia abierta y versión | Acceso incompleto, discrepancia o alcance |
| --- | --- | --- |
| CL01 | HTML de tres secciones del manual LAMMPS, `stable`, cabecera `30Sep2026 update1`. | `stable` es móvil; se conserva el identificador observado. Vecinos/multi/sort son secciones del mismo manual, no estudios independientes. |
| CL02 | Archivo oficial `src/gromacs/nbnxm/atomdata.h`, tag v2026.3, texto completo. | Intentos de v2026.4 en raw y GitHub devolvieron cache miss; no se afirma ausencia del tag. La palabra AoSoA es interpretación explícita de índices y paquetes. |
| CL03 | arXiv PDF 1306.1737v1 y ficha con aceptación en CPC y DOI. | Editorial 403; el footer del PDF lleva fecha de generación 2024, mientras ficha y cabecera identifican 2013. No se trató como publicación nueva de 2024. |
| CL04 | arXiv PDF 2006.09167v2, 2020; PubMed para identidad editorial. | Manuscrito de autores, no comparación ejecutada del GROMACS actual. |
| CL05 | arXiv PDF 1308.5587v2, revisado 2019; documentación oficial de citas 7.2.0. | El paper describe v2.6 y se publicó en 2020; no se confunden esos números con el software vigente. Editorial 403. |
| CL06 | HTML editorial PLOS, DOI 10.1371/journal.pcbi.1005659, 2017. | PMC devolvió comprobación de navegador; fuente sustituida por HTML editorial completo. |
| CL07 | OpenMM Dev Guide 8.6, capítulo 6 OpenCL, URL `latest`. | Intento del supuesto capítulo CUDA `04_cuda_platform.html` no accesible; no se dedujeron detalles de CUDA de ese fallo. Fecha individual desconocida. |
| CL08 | PDF primario editorial de 17 páginas alojado en biomolmd.org, revista/volumen/año/autores visibles. | DOI 10.1063/1.470117 redirigió a recurso editorial inaccesible; se conservó el acceso alternativo. No se asumió revisión de versión posterior. |
| CL09 | PDF arXiv cond-mat/9807099v1 y ficha con referencia JCP 109 (1998), 7678, DOI. | No se comprobó intervalo de páginas editorial; no se inventó. |
| CL10 | PDF arXiv cond-mat/9807100v1 y ficha con JCP 109 (1998), 7694, DOI. | El PDF muestra una fecha interna de generación 2018; la ficha de depósito/referencia identifica 1998. No se usó esa fecha interna como año de publicación. |
| CL11 | HTML versionado GROMACS 2026.4, Virial from SHAKE y contexto de periodicidad. | Algunas ecuaciones HTML tienen notación difícil de extraer; no se transcribió la fórmula como implementación validada. |
| CL12 | HTML editorial completo PLOS ONE, 2018. | Datos/scripts enlazados no ejecutados ni auditados en esta revisión. |

## Cobertura y sesgos de selección observables

**DISTRIBUTIONAL_SKEW_ADVISORY — método:** ocho de ocho artículos de revista identificados son trabajos de algoritmos/software y evaluación computacional (8/8 conocidos); las otras cuatro entradas son implementación oficial. La concentración es coherente con la pregunta técnica solicitada; no representa diversidad de validación experimental material. No se infirió geografía de autores o instituciones como geografía de evidencia.

**DISTRIBUTIONAL_SKEW_ADVISORY — antigüedad:** los ocho artículos de revista identificados son anteriores a 2021 (8/8 con año conocido). Se retienen deliberadamente por fundamentos matemáticos y mecanismos; las cuatro entradas oficiales aportan versiones consultadas de 2026/8.6, con año de publicación individual desconocido cuando corresponde. Esto no demuestra ausencia de trabajos recientes.

Los artículos de arquitectura proceden de desarrolladores de los motores estudiados; aportan conocimiento directo pero sus benchmarks no son replicación independiente. Hay dependencia de familias de autores: Páll/Hess y Deserno/Holm tienen más de una entrada. La evidencia adversa se registró dentro de las anotaciones: trabajo redundante por clusters, sobrecoste de preparación, dependencia de hardware, precisión, falta de conservación de momento de cierta diferenciación, inhomogeneidad de cargas y límites periódicos del virial. No se calcularon efectos combinados, aceleraciones esperadas ni prioridades para MODENA.

**Cierre del cribado:** 12 registros finales anotados, con ocho publicaciones de revista identificadas y cuatro entradas oficiales; no equivalen a 12 experimentos independientes. No se conocen totales de descubrimiento o exclusión. Corpus preparado con asistencia de IA y lectura directa de fuentes.

**Alias entre entregables de Fase 2:** CL12 = V1 de [VALIDATION-SOURCES.md](VALIDATION-SOURCES.md). Descontar esta duplicación al contar obras del corpus conjunto. CL11 y V5/V6 abordan restricciones/presión desde manuales diferentes y no son replicaciones experimentales. Se incorporó por referencia V7, ya investigado por el orquestador, sin repetir búsquedas de barostatos MC.
