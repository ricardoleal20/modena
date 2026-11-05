# Registro de búsqueda: DFT y ASE

Fase 2, consulta **2026-11-05 UTC**, iniciada antes de las 05:30 UTC. Revisión técnica dirigida conforme a [METHODOLOGY.md](../phase1_scope/METHODOLOGY.md) y los casos [C4/C5](CASES.md). El [catálogo anotado](DFT-SOURCES.md) contiene ocho registros bibliográficos; uno agrupa tres capítulos del mismo manual. No se usaron cifras estimadas de resultados del buscador ni un flujo PRISMA.

## Parámetros y criterios previos

**Preguntas:** alcance electrónico de PWscf; BOMD frente a CPMD; disponibilidad y convergencia de E/F/estrés; pseudopotenciales, malla k y unidades; contrato ASE, identidad atómica y PBC; socket i-PI y acceso a proveedores ML.

**Fuentes objetivo:** documentación y código oficiales de Quantum ESPRESSO y ASE. **Idioma:** inglés, con registro en español. **Ventana:** sin corte inicial; se prefiere la documentación oficial servida en la fecha de consulta y se registra la versión visible. Los años de resultados o de rastreo del buscador no son metadatos bibliográficos.

**Inclusión:** acceso abierto al contenido original; vínculo directo con las preguntas; versión/título localizables; afirmación respaldada por el pasaje o código inspeccionado. **Exclusión:** resúmenes de buscador como prueba, tutoriales ajenos cuando existe documentación oficial, foros, duplicados de un mismo manual y referencias no comprobadas. No se persiguió una cuota de artículos revisados por pares: este subencargo pregunta por contratos y capacidades de software y solicita fuentes primarias oficiales.

**Corpus previo:** no se recibió `literature_corpus[]` ni Material Passport. Se usaron el método y los casos locales como contexto, no como fuentes externas que acrediten capacidades de QE/ASE. La deduplicación se hizo por URL canónica, título y versión. No se consultó Semantic Scholar: los registros seleccionados son manuales/código, sin DOI asignado ni supuesto identificador de artículo. No se declaró que una búsqueda no ejecutada hubiera fallado.

## Búsquedas y navegación efectivamente ejecutadas

Todas las consultas se ejecutaron en el buscador web disponible, en el orden indicado. Las cadenas se conservan literalmente. Se inspeccionó cada grupo de resultados para escoger originales; los hits de otras procedencias no se contaron como corpus incluido.

| ID | Consulta literal o acción | Resultado y decisión trazable |
| --- | --- | --- |
| Q1 | `site.quantum-espresso.org Doc pw_user_guide convergence forces stress molecular dynamics` | Descubrió secciones oficiales de capacidades, dinámica y troubleshooting; se recuperó después el PDF completo 7.5 para DFT-01. |
| Q2 | `site.quantum-espresso.org Doc INPUT_PW calculation md scf tstress tprnfor` | Recuperó la especificación PW 7.5, incluida como DFT-02; se abrieron los parámetros relevantes. |
| Q3 | `site.ase-lib.org espresso calculator implemented_properties stress units` | No permitió resolver la ruta oficial en este primer grupo; se amplió la consulta en Q5. No se atribuye un conteo de resultados. |
| Q4 | `site.ase-lib.org socketio i-PI calculator` | Aparecieron espejos; no se incluyeron como autoridad. Se resolvió la página oficial mediante navegación posterior. |
| N1 | Aperturas de `https://ase-lib.org/ase/calculators/espresso.html`, `.../calculators.html`, `https://ase-lib.org/ase/atoms.html`, `.../socketio/socketio.html` | Las cuatro rutas devolvieron 404. Se corrigió el host a `docs.ase-lib.org` tras Q5–Q7. |
| N2 | Apertura de `https://www.quantum-espresso.org/Doc/pw_user_guide/pw_user_guide.html` | No accesible mediante la herramienta. Se utilizó el PDF oficial y las secciones HTML individuales descubiertas. |
| Q5 | `ASE atomic simulation environment official documentation calculators espresso implemented_properties stress` | Identificó el dominio oficial actual `docs.ase-lib.org`; las páginas del manual fueron abiertas directamente. |
| Q6 | `ASE documentation atoms get_stress Voigt order eV Angstrom units` | Descubrió código/API de Atoms; se inspeccionaron el capítulo de Atoms y Units oficiales para DFT-05. |
| Q7 | `ASE documentation socket I/O Quantum Espresso i-PI limitations` | Descubrió la página socket bajo la ruta histórica DTU; se recuperó su equivalente oficial actual, DFT-08. |
| Q8 | `Quantum Espresso pw user guide 3.3 Electronic structure calculations convergence k points cutoffs pseudopotentials` | Confirmó sección electrónica y documentación fuente; se eligió el manual oficial, sin usar tutoriales o resultados de foros como respaldo. |
| N3 | Apertura directa de `https://docs.ase-lib.org/ase/calculators/espresso.html`, `.../calculators.html`, `https://docs.ase-lib.org/ase/atoms.html`, `.../socketio/socketio.html` | Contenido completo accesible y títulos comprobados; registros DFT-05, DFT-06 y DFT-08. |
| N4 | Apertura directa de `https://docs.ase-lib.org/_modules/ase/io/espresso.html` y `https://docs.ase-lib.org/_modules/ase/calculators/calculator.html` | Código documentado accesible, incluido en DFT-07 y DFT-04 respectivamente. |
| N5 | Apertura directa de `https://www.quantum-espresso.org/Doc/INPUT_CP.html` y `https://www.quantum-espresso.org/Doc/user_guide_PDF/pw_user_guide.pdf` | Versiones 7.5 visibles; DFT-03 y DFT-01. |
| N6 | Apertura directa de `https://docs.ase-lib.org/ase/units.html` | Unidades base y tiempo interno comprobados; capítulo agrupado en DFT-05. |

## Verificación de contenido y decisiones

Se leyó el texto primario de las secciones seleccionadas, no solo los fragmentos del buscador. Para navegar documentos largos se buscaron literalmente los siguientes símbolos o expresiones y se inspeccionó su contexto:

| Objeto | Localizadores efectivamente consultados | Registro |
| --- | --- | --- |
| PWscf User’s Guide | `convergence`, `k-points`, `pseudopotentials`; §§1.1, 3.3, 3.4 y 5 | DFT-01 |
| INPUT_PW | `ecutrho REAL`, `conv_thr REAL`, `assume_isolated CHARACTER`, `input_dft CHARACTER`, `Monkhorst`, `pseudopotential file`; tarjetas de estructura y especies | DFT-02 |
| INPUT_CP | `emass`, `time step for molecular dynamics`; contexto de `electron_dynamics` | DFT-03 |
| calculator.py documentado | `implemented_properties`, `all_changes`; métodos de obtención y comparación de estado | DFT-04 |
| Manual ASE | `MACE`, `pbc`, `Voigt`, `eV`; capítulos Atoms, Calculators y Units | DFT-05 |
| Espresso | `units`, `Last updated`; parámetros y listas de implementación | DFT-06 |
| io.espresso | `stress *=`, `create_units`, `pbc=`, `isinstance(kpts, np.ndarray)`; parseo de energía y comentarios térmicos | DFT-07 |
| Sockets ASE | `numbers`, `stress`; introducción, compatibilidad, ejemplo cliente y clase del servidor | DFT-08 |

Las búsquedas de `Limitations`, `Limit` y `Warning` en la página socket no devolvieron una sección con esos términos. Esto no se interpretó como ausencia de límites: se leyó el ejemplo cliente y el contrato del servidor. `implemented_properties` no apareció en el capítulo general de calculators, por lo que se verificó en su código oficial. La consulta de `pseudo_file` no encontró ese símbolo en INPUT_PW; el campo real se localizó como `PseudoPot_X`.

**Discrepancia conservada:** dentro de Espresso, la lista de funciones pendientes y la API detallada no concuerdan sobre k-points. Se contrastó con el código escritor y se registró el caso en DFT-06/07. La discrepancia no se resolvió inventando una versión de ASE ni suponiendo que toda la página estuviera actualizada homogéneamente.

**Criterio de versión:** para trasladar esa observación a una instalación se exige identificar release/commit de ASE, conservar código o hash de la función escritora y comprobar su tarjeta k con pesos sobre una entrada fijada. Sin esa evidencia se atribuye la observación exclusivamente a las páginas live consultadas. No se ejecutó tal comprobación en esta investigación.

## Selección y exclusiones

| Material encontrado | Decisión | Motivo |
| --- | --- | --- |
| Manual PWscf PDF 7.5 y secciones HTML oficiales | Incluir el manual; HTML como localizadores del mismo trabajo | Evita contar cada sección como una fuente independiente. |
| PDF PW alojado en `/wp-content/uploads/2022/03/` | No incluir como registro adicional | Copia anterior; se eligió la ruta actual del manual 7.5. |
| Entradas oficiales PW y CP 7.5 | Incluir separadamente | Contratos de ejecutables y algoritmos diferentes. |
| Manual ASE y código oficial documentado | Incluir según catálogo | Evidencia primaria del contrato y de las conversiones reales. |
| `yuzie007ase.readthedocs.io` y ruta histórica `wiki.fysik.dtu.dk` | No incluir como registros principales | Se recuperaron las páginas del host oficial actual. |
| Wikipedia, Reddit, tutoriales de curso y páginas de terceros devueltos por búsqueda | Excluir del respaldo factual | Hay originales oficiales disponibles para las preguntas concretas. |
| Artículos QE de 2009/2017 y otros papers aparecidos en Q8 | No incluir en este subconjunto | No se inspeccionó su texto completo; los contratos actuales se documentaron con fuentes oficiales. No se les atribuyen hallazgos no leídos. |

Se seleccionaron **ocho registros bibliográficos** tras lectura y comparación. El buscador no proporcionó un universo estable de candidatos; por ello no se informan totales de hits, tasas de exclusión ni diagrama PRISMA. La agrupación de DFT-05 es por manual, con tres URL explícitas y localizadores separados; no implica independencia entre capítulos.

## Cobertura y límites de la búsqueda

**DISTRIBUTIONAL_SKEW_ADVISORY:** documentación/código de software = 8/8 registros conocidos (100%). Es concentración deliberada del tipo de fuente, coherente con el subencargo de capacidades e interfaces, no un defecto estadístico. Respuesta de búsqueda: sin expansión a artículos de validación material, que pertenecen a otros ejes del estudio. Los registros provienen de dos proyectos; los múltiples capítulos y el código ASE no constituyen corroboraciones independientes de precisión científica. No se asignó distribución geográfica ni año editorial a documentos sin esos metadatos.

Limitaciones: búsqueda dirigida, sin API de bases académicas ni archivado local de cada original; documentación web mutable de ASE sin commit fijado; no consulta de releases instaladas ni pruebas de ejecución. La fecha visible de actualización del sitio se distingue de fecha de publicación. La evidencia permite describir capacidades y contratos, pero no medir errores, rendimiento, fiabilidad de un MLIP, contaminación de datasets ni convergencia de un caso material. Investigación realizada con asistencia de IA.
