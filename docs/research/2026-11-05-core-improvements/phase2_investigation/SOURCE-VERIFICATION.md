# Verificación de fuentes — Fase 2

**Fecha:** 2026-11-05 UTC. **Objeto:** calidad, identidad y respaldo de afirmaciones del corpus; no síntesis, prioridades ni validación de MODENA. Se revisaron los **38 registros** de [CLASSICAL-SOURCES.md](CLASSICAL-SOURCES.md), [DFT-SOURCES.md](DFT-SOURCES.md), [ML-SOURCES.md](ML-SOURCES.md) y [VALIDATION-SOURCES.md](VALIDATION-SOURCES.md), con sus tres `SEARCH-LOG.md` y el registro integrado de VALIDATION. **CL12 y V1 son la misma obra: quedan 37 unidades bibliográficas, no 37 estudios independientes.** Ninguna fuente se rechaza por inexistencia demostrada; varias inferencias quedan excluidas o condicionadas abajo.

## Método y alcance de la comprobación

- **R:** reapertura del original en esta invocación y examen de pasajes pertinentes; no implica auditoría completa de todas las ecuaciones/resultados.
- **M:** reapertura de metadatos en ficha primaria, editorial o registro bibliográfico; el contenido detallado se conserva de la lectura y los localizadores del bibliógrafo.
- **L:** contraste del registro, sus localizadores y ledger; no se reabrió ese original en esta invocación. No equivale a segunda lectura independiente.

Las filas marcan la comprobación más fuerte efectuada, con los límites correspondientes. El agente de esta verificación también redactó el corpus DFT en una invocación anterior: una reapertura R de DFT no se presenta como revisión por otra persona o por un investigador independiente. No hay atestación de lectura humana, adquisición local de todos los originales, ejecución de benchmarks, inspección de pesos ni prueba de inferencia.

**S2-API-UNAVAILABLE:** el orquestador comunicó el intento fallido para Merz; ML-SEARCH-LOG registra otros intentos fallidos de Graph API. No se repitieron ni inventaron IDs. Se utilizó deduplicación por DOI, título, autores y relación manuscrito/publicación.

Se intentó abrir `doi.org` para los 17 registros distintos de artículos/preprints con DOI del corpus: CL03–CL06, CL08–CL10, CL12/V1, ML03 y ML06–ML11, V2–V3. Ninguna apertura dio acceso íntegro a una landing verificable en esta herramienta: V2/V3 devolvieron 403; otras devolvieron recurso no accesible o redirecciones fallidas, incluidas rutas de PLOS que sí abrieron por URL editorial directa. **Esto no prueba DOI inexistente ni autoriza la etiqueta `DOI_VERIFIED`.** La identidad se contrasta por los originales alternativos de las tablas. No se convierte un abstract de PubMed en evidencia del método.

## Criterio de calidad

La calificación es **relativa a la afirmación**: **A**, fuente primaria apta para una definición, algoritmo o contrato acotado; **B**, evidencia primaria útil con límite de versión, evaluación o acceso; **C**, utilización que requiere una condición todavía no resuelta. La identidad bibliográfica y el respaldo del contenido son ejes separados.

La jerarquía I–VII de la habilidad se adapta de forma explícita: **III-T**, comparación computacional controlada o análisis matemático con pruebas numéricas; **VI-T**, descripción primaria de implementación; **VII-T**, documentación/especificación técnica. Estos rótulos no convierten simulaciones en ensayos clínicos ni penalizan un manual por no ser un RCT. Se valoran adecuación, revisión formal, método, datos/alcance, vigencia y conflictos. Los manuales pueden ser A para su API y no demostrar exactitud material; los resultados computacionales pueden ser B para transferencia a otro código. No se asignan puntuaciones clínicas o precisión estadística inexistentes.

**Intereses, columna I:** `D` = documentación/código mantenidos por sus desarrolladores; `P` = autores proponen o evalúan su propio método/software, con posible interés intelectual/institucional; `F` = declaración financiera concreta consultada; `N` = declaración explícita de ausencia consultada. D/P no significan conflicto financiero oculto. Donde no se leyó una declaración, no se certifica su ausencia.

## Matriz por fuente: MD clásica

| ID | Nivel; grado para el uso descrito | Comprobación de identidad/contenido | I y límite |
| --- | --- | --- | --- |
| CL01 | VII-T; A para manual, B para transferencia | **L**. Tres capítulos oficiales, LAMMPS stable `30Sep2026 update1`; identidad y secciones trazables en ledger. | D. `stable` mutable; no mide mejora en MODENA. |
| CL02 | VI-T; A para disposición observada | **L**. Archivo oficial `atomdata.h`, tag `v2026.3`; los símbolos y localizadores están documentados. | D. AoSoA es interpretación explícita del índice, no denominación atribuida al paper de 2013. |
| CL03 | III-T; B | **M**. [arXiv 1306.1737](https://arxiv.org/abs/1306.1737) confirma Páll/Hess, título, 2013, aceptación CPC y DOI relacionado. Contenido §§2–2.5: ledger del bibliógrafo. | P. Manuscrito distinto del acceso editorial; no se transfieren tamaños óptimos ni speedups. |
| CL04 | III-T; B | **M**. [arXiv 2006.09167](https://arxiv.org/abs/2006.09167) confirma ocho autores, título y v2 de 2020; DOI relacionado. Referencia editorial completa contrastada por el bibliógrafo en PubMed. | P. Arquitectura y mediciones históricas, no ejecución del GROMACS instalado. |
| CL05 | VI-T/III-T; B | **M**. [arXiv 1308.5587](https://arxiv.org/abs/1308.5587) confirma título actual, tres autores y revisión 2019. Año 2020, volumen y artículo proceden del registro oficial de citas leído por el bibliógrafo. | P. v2.6 descrita; no se identifica con el software actual ni con arXiv v1 de 2013. |
| CL06 | III-T; A para método/contraste, B para generalización | **R**. [PLOS OpenMM 7](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005659), texto editorial y condición peer-reviewed visibles. | **F/P**. La declaración identifica vínculos de V. S. Pande y J. D. Chodera con Schrödinger, de Pande con Andreessen Horowitz y empleo de K. A. Beauchamp en Counsyl. No invalida el método; los benchmarks no son replicación independiente. |
| CL07 | VII-T; A para OpenCL documentado | **L**. OpenMM Developer Guide 8.6, capítulo OpenCL y secciones precisas disponibles en ledger. | D. No autoriza atribuir estructuras a todos los backends ni ausencia total de transferencias. |
| CL08 | III-T; A para definición matemática | **R**. [PDF SPME original](https://biomolmd.org/mw/images/e/e0/Spme.pdf), cabecera con seis autores, JCP 103(19), 8577–8593, 1995; §IV y p. 8581 examinadas. | P. Copia editorial en alojamiento tercero; DOI inaccesible. La diferenciación analítica es conservativa, sin promesa de energía numérica exactamente constante a paso finito. |
| CL09 | III-T; A para distinciones, B para transferencia | **R**. [Ficha](https://arxiv.org/abs/cond-mat/9807099) y [PDF primario](https://arxiv.org/pdf/cond-mat/9807099), pp. 7–8: asignación, diferenciación, condiciones suficientes de conservación. Autores, año, JCP 109, 7678 y DOI coinciden. | P. La corrección global de fuerza media no elimina errores locales; no se inventa intervalo editorial. |
| CL10 | III-T; B | **M**. [Ficha primaria](https://arxiv.org/abs/cond-mat/9807100) confirma autores, título, JCP 109, 7694, 1998 y DOI; límites geométricos extraídos por el bibliógrafo. | P. Estimador P3M/ik en régimen específico; no cota universal SPME o triclínica. |
| CL11 | VII-T; A para algoritmo documentado | **L**. Manual GROMACS 2026.4, «Virial from SHAKE» y contexto periódico con ecuaciones/localizadores registrados. | D. Procedimiento para integrador concreto; no receta universal portable. |
| CL12 = V1 | III-T; A para pruebas, B para cobertura global | **R**. [PLOS Merz/Shirts](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202764): autores, 2018, 13(9), e0202764 y métodos visibles. | **N/P**. Declaran no tener competing interests; desarrollan la herramienta evaluada. No certifica ergodicidad ni observables materiales. |

## Matriz por fuente: DFT y ASE

| ID | Nivel; grado | Comprobación de identidad/contenido | I y límite |
| --- | --- | --- | --- |
| DFT-01 | VII-T; A | **R**. [Manual original PWscf](https://www.quantum-espresso.org/Doc/user_guide_PDF/pw_user_guide.pdf), versión 7.5, capacidades y distinción BOMD/CPMD. | D. Capacidad electrónica documentada, sin validación de una configuración C5. |
| DFT-02 | VII-T; A para especificación | **L**. INPUT_PW 7.5, campos y unidades comprobados en la invocación bibliográfica anterior; ledger completo. | D. SCF y discretización no son el mismo control; referencia convergida todavía no ejecutada. |
| DFT-03 | VII-T; A | **L**. INPUT_CP 7.5: dinámica electrónica, masa efectiva y unidades temporales documentadas. | D. `cp.x` admite más de un tratamiento electrónico; no todos sus modos equivalen a una única CPMD. |
| DFT-04 | VI-T; A | **L**. Código oficial `ase.calculators.calculator`, métodos y errores localizados. | D. Capacidad declarada no garantiza una propiedad físicamente válida para un modelo. |
| DFT-05 | VII-T; A para contrato | **L**. Tres capítulos del mismo manual, con URL y secciones separadas; no tres estudios. | D. Elementos/PBC/unidades son metadatos del contrato, no validación de cobertura ML. |
| DFT-06 | VII-T; B; C para extrapolar a instalación sin pin | **L**, contrastado con código **R** de DFT-07. Discrepancia k-points explícita y criterio de pin registrado. | D. La lista antigua de pendientes no prevalece sobre una rama implementada sin inspeccionar versión. |
| DFT-07 | VI-T; A para código publicado | **R**. [Fuente oficial ASE](https://docs.ase-lib.org/_modules/ase/io/espresso.html), `_get_stress` y lectura de estructuras: signo/conversión y periodicidad. | D. Documentación live de 2026-11-02; no corresponde a un paquete instalado fijado. |
| DFT-08 | VII-T; A para protocolo documentado | **L**. Socket I/O oficial ASE y ejemplo cliente, leídos/localizados por el bibliógrafo. | D. Persistencia del proceso no prueba reducción de coste del solver ni transmisión completa de identidad. |

## Matriz por fuente: proveedores ML

| ID | Nivel; grado | Comprobación de identidad/contenido | I y límite |
| --- | --- | --- | --- |
| ML01 | VI/VII-T; A para lo observado, C para receta conjunta sin pin | **R**. [Predictor `main`](https://raw.githubusercontent.com/facebookresearch/fairchem/main/docs/core/common_tasks/ase_calculator.md), [tasks](https://raw.githubusercontent.com/facebookresearch/fairchem/main/docs/core/uma.md), [calculator](https://raw.githubusercontent.com/facebookresearch/fairchem/main/src/fairchem/core/calculate/ase_calculator.py) y [guía MD](https://facebookresearch.github.io/fairchem/molecular-dynamics/). | D. La guía renderizada usa `turbo`; `main` lo declara retirado y usa Nitro. No se fusionan en instrucciones de una release no identificada. |
| ML02 | VII-T; A para términos visibles | **R**. [Tarjeta UMA](https://huggingface.co/facebook/UMA) y [LICENSE.md](https://raw.githubusercontent.com/facebookresearch/fairchem/main/LICENSE.md). Nombre/hash anunciado 1.2.1 comprobado en tabla; no archivo descargado. | D. Código MIT y materiales/pesos FAIR Chemistry License v1 son superficies distintas; acceso gated visible. No dictamen jurídico ni autorización aceptada. |
| ML03 | III-T; B | **R**. [Ficha UMA](https://arxiv.org/abs/2506.23971) y [v2](https://arxiv.org/html/2506.23971v2): dieciocho autores, depósito 2025, revisión 2026-03-04, entrenamiento y límites. | P/institucional Meta. Se mantiene como preprint; no se verificó publicación editorial posterior. Descripción de entrenamiento no demuestra cada modo/checkpoint actual. |
| ML04 | VII-T; B | **L**. MACE docs 0.3.13, interfaces original y ML-IAP con preparaciones y límites diferenciados. | D. Rutas distintas, beta y multigpu no se confunden; no ONNX universal. |
| ML05 | VII-T; B | **L**. NequIP latest: compilación, tipos, conversiones y grafos; URL y localizadores en ledger. | D. Dispositivo/formato/versión condicionan despliegue; no todos los pesos son intercambiables. |
| ML06 | III-T; A para mecanismo, B para extrapolación | **R**. [PMLR ICML 2025](https://proceedings.mlr.press/v267/fu25h.html) confirma autores, volumen 267 y pp. 17875–17893; [v2](https://arxiv.org/html/2502.12147v2) para diseño/conservación. | P; autores relacionados con FAIR-Chem. Localizador NVE correcto en v2: **A.1**, ya corregido en ML-SOURCES; C es eficiencia de inferencia. |
| ML07 | III-T; B | **M**. [arXiv](https://arxiv.org/abs/2210.07237) confirma siete autores, v2 y referencia TMLR 2023; benchmark detallado: lectura del bibliógrafo. | P. OpenReview presentó challenge al bibliógrafo; no se atribuye lectura editorial íntegra aquí. No evalúa UMA actual. |
| ML08 | III-T; A para contraejemplos y esquema, B para transferencia | **R**. [PMLR](https://proceedings.mlr.press/v267/bigi25a.html) confirma tres autores, ICML 2025, pp. 4384–4414; [arXiv v6](https://arxiv.org/html/2412.11569v6), §4.8 y apéndice I, explicita MTS y límites. | P. La revisión posterior a conferencia se distingue; el caso favorable a MTS no justifica conmutación arbitraria. |
| ML09 | III/VI-T; B | **M**. [arXiv DP-GEN](https://arxiv.org/abs/1910.12690) confirma siete autores, CPC 253 (2020) 107206 y DOI; ciclo con motores separados explícito. | P. Umbrales y cobre del ejemplo no son garantías para otra química. |
| ML10 | III-T; B | **M**. [PubMed 37724552](https://pubmed.ncbi.nlm.nih.gov/37724552/) confirma autores, DOI, 2023, 25(37), 25828–25837. Método: manuscrito arXiv leído por el bibliógrafo, no PubMed. | P; declaración financiera no inspeccionada. RSC inaccesible; no acceso editorial íntegro supuesto. |
| ML11 | III-T; B | **M**. [arXiv](https://arxiv.org/abs/2108.05748) confirma autores, referencia PRE 105 (2022) 015311 y DOI; título de preprint con «On the», distinguido del publicado. | P; declaración financiera no inspeccionada. Cuatro casos no demuestran que todo ensemble falle. |

## Matriz por fuente: validación y presión

| ID | Nivel; grado | Comprobación de identidad/contenido | I y límite |
| --- | --- | --- | --- |
| V1 | — | Mismo registro que CL12, no se suma. | Véase CL12. |
| V2 | III-T; A para método, B para transferencia | **M**. [arXiv Shirts](https://arxiv.org/abs/1208.0910) confirma autor, título, JCTC 9(2), 909–926, 2013 y DOI. Lectura del método: bibliógrafo. | P. DOI 403; requiere muestras/solapamiento adecuados, no diagnóstico causal de todo fallo. |
| V3 | III-T; A para método, B para cobertura global | **R**. [Manuscrito en PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4945107/), autor, título y referencia final 2016, 12(4), 1799–1805 visibles. | P. Descarte/autocorrelación de la serie observada; no garantía de transiciones lentas no vistas. |
| V4 | VII-T; A para API, B para heurística | **L**. PyMBAR rotulado `4.2.1.dev0+gcf12100f8.d20251016`; métodos/localizadores documentados. | D. No versión instalada ni equilibración global certificada. |
| V5 | VII-T; A | **R**. [LAMMPS SHAKE/RATTLE](https://docs.lammps.org/fix_shake.html), release `30Sep2026`, sección virial y restricciones. | D. Contribución relevante al estimador; no condición universal para aceptar todo MC NPT. |
| V6 | VII-T; A | **L**. LAMMPS compute pressure, stable `30Sep2026 update1`; términos y convención registrados. | D. Definición de un estimador, no prueba de implementación MODENA. |
| V7 | VII-T; A para algoritmo | **R**. [OpenMM 8.6 §19.12](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#montecarlobarostat): propuesta por centros moleculares y peso de aceptación. | D. MC molecular y estimador de presión son cuestiones diferentes; no se atribuye equivalencia al MC atómico local. |

## Veredictos de afirmaciones críticas

`VERIFIED` significa respaldo en el alcance y la versión identificados; no certificación universal. `CONDITIONAL` conserva una condición no resuelta. `EXCLUDED` rechaza **la inferencia indicada**, no la existencia o utilidad de las fuentes.

| Afirmación examinada | Veredicto | Evidencia y condición |
| --- | --- | --- |
| La diferenciación analítica de SPME deriva fuerzas de su energía aproximada, pero no impone momento total conservado. | **VERIFIED** | CL08 y CL09, ambos reabiertos, autores independientes; §IV/p. 8581 y pp. 7–8 respectivamente. No se confunde fuerza conservativa con energía exactamente constante a timestep finito. |
| Una suma de fuerzas SPME no nula demuestra por sí sola un bug del kernel analítico. | **EXCLUDED** | Contradice la propiedad declarada por CL08 y discutida por CL09; depende del método de diferenciación y del error de malla. |
| Las fuerzas de restricción contribuyen al virial/estimador de presión. | **VERIFIED** | V5 reabierto, con CL11 y V6 como referencias de otras implementaciones/localizadores. El detalle depende de integrador y periodicidad. |
| Un tensor de restricciones completo es requisito universal antes de que cualquier algoritmo NPT pueda aceptar movimientos. | **EXCLUDED** | V7 ofrece un contraejemplo explícito de MC molecular por energía/volumen. La presión reportada y un barostato dinámico conservan sus propios requisitos. |
| OpenMM MC molecular usa cambio energético, presión externa, volumen y número de moléculas en la aceptación. | **VERIFIED** | V7, §19.12. Escala centros moleculares; no requiere virial instantáneo para esa aceptación. No valida la medida/propuesta del MC atómico de MODENA. |
| Drift pequeño, temperatura media o MAE bajo bastan para validar ensemble y observables. | **EXCLUDED** | CL12/V1 y ML07; ML06/ML08 aportan controles y contraejemplos complementarios. Son fuentes de autores/implementaciones diferentes, no duplicados del mismo benchmark. |
| Suavidad del grafo/cutoff y conservatividad importan además de error estático. | **VERIFIED**, en sistemas ensayados | ML06 y ML08; no convierte NVE en prueba de ergodicidad ni garantiza toda propiedad material. |
| Dispersión de ensembles es automáticamente un margen calibrado de error fuera de distribución. | **EXCLUDED** | ML10 y ML11, estudios distintos con contextos definidos. El respaldo detallado procede de las lecturas primarias de sus bibliógrafos; aquí se reabrieron metadatos. |
| Las tareas UMA representan niveles de teoría distintos y requieren metadatos físicos. | **VERIFIED** | ML01, tabla oficial y código. OMol y OMat no constituyen una referencia energética común por compartir modelo. El encabezado histórico de cinco tasks no se usa como censo inmutable de las filas actuales. |
| Una receta única puede mezclar `main` Nitro y los ejemplos MD `turbo` sin fijar versión. | **EXCLUDED** | ML01, discrepancia reabierta. El contrato observable es utilizable; la reproducción queda **CONDITIONAL** a release/commit, checkpoint/task, ajustes de inferencia y dependencias concordantes. |
| Copiar energy a free_energy demuestra una PMF o energía libre termodinámica. | **EXCLUDED** | ML01, código calculator: alias de compatibilidad de API; no estimación de potencial de fuerza media. |
| Código MIT implica que los pesos UMA también tienen licencia MIT. | **EXCLUDED** | ML02, dos documentos de distribución distintos reabiertos. Condiciones visibles verificadas; no se aceptaron términos ni descargaron pesos. |
| E/F/estrés de un MLIP equivalen al conjunto de propiedades electrónicas de QE. | **EXCLUDED** | DFT-01 y ML01: dominios de salida diferentes; la documentación UMA remite a DFT para work functions. Integración para etiquetas no acredita reemplazo electrónico. |
| Una comparación QE–ML puede ignorar funcional, pseudopotenciales, smearing, estados electrónicos y convergencia de la referencia. | **EXCLUDED** | DFT-02 y ML01. Una comparación concreta permanece **CONDITIONAL** a identidad física y convergencia declaradas; no se ejecutó C5. |
| ASE conserva sin transformación el estrés y toda la intención de frontera del output QE. | **EXCLUDED** | DFT-07: conversión/signo y construcción de periodicidad explícitos. Los metadatos originales no se infieren de un Atoms parseado solamente. |
| ASE no admite en ninguna versión k-points explícitos. | **EXCLUDED** | DFT-06/07: discrepancia de texto frente a rama N×4; la afirmación sobre una instalación requiere el pin y comprobación definidos en el corpus. |
| Una interfaz MACE/NequIP/ASE demuestra exportación ONNX universal o igual semántica de todos los pesos. | **EXCLUDED** | ML04/05 documentan rutas positivas específicas; ausencia de un hit ONNX tampoco prueba imposibilidad. |
| Aplicar correcciones conservativas mediante MTS valida cualquier fallback instantáneo basado en incertidumbre. | **EXCLUDED** | ML08 especifica un algoritmo y examina límites; ML09/10 no proporcionan esa equivalencia dinámica. |
| Los mecanismos de vecinos, empaquetado y residencia garantizan determinada aceleración en MODENA. | **CONDITIONAL**, aceleración no demostrada | CL01–CL07 apoyan mecanismos, condiciones y alternativas; ninguna es medición del checkout local. |

## Duplicación, contradicciones y calidad editorial

- **Duplicados:** CL12=V1 por DOI/título/autores. Las fichas de arXiv, el manuscrito y la publicación de una obra son accesos relacionados, no réplicas. Se agrupan capítulos de manuales sin tratarlos como estudios independientes. V2/V3/V4 y CL12 comparten antecedentes metodológicos; no se suman como demostraciones estadísticas independientes de una implementación.
- **Contradicciones documentales:** FAIR-Chem `main`/guía MD y ASE k-points son discrepancias de versión/superficie, conservadas como tales. La discrepancia de localizador eSEN quedó corregida a A.1. Las unidades, outputs y bibliografía restantes no presentan contradicción crítica demostrada en lo inspeccionado; esto no equivale a auditoría exhaustiva.
- **Procedencia editorial:** se observaron artículos primarios en PLOS, JCP, CPC, Computational Materials Science, JCTC, PCCP, PRE y PMLR/ICML, y una referencia TMLR. UMA se mantiene como preprint. No se comprobó Scopus/WoS, COPE, DOAJ, Cabell, ISSN o identidad ORCID para todo el corpus; por tanto no hay sello de «revistas libres de prácticas predatorias». Los manuales y licencias no se evalúan como revistas. No apareció una falsedad bibliográfica demostrada que justifique retirar una fuente.
- **Conflictos:** CL06 contiene declaración financiera verificada; CL12/V1 declaración explícita de ausencia. En los demás registros, el vínculo autoral con método/producto se señala sin inventar intereses financieros ni dar por leída una declaración que no se abrió. En ML01–ML03 y los artículos de motores, la evaluación propia no constituye replicación independiente.
- **Retracciones y correcciones editoriales:** no se consultó una base de retracciones ni se revisó sistemáticamente Crossmark/avisos de cada obra. Se registra **no comprobado**, no «sin retractaciones». Un DOI inaccesible no se interpreta como retirada.

## Límites del gate

**Observación local adicional:** se leyó [DIAGNOSTIC-OBSERVATION.md](DIAGNOSTIC-OBSERVATION.md), fuera del conteo bibliográfico externo. Su ejemplo de sumar una constante al cero energético es algebraicamente compatible con la fórmula de drift relativo transcrita: cambia el denominador y puede cambiar la etiqueta sin cambiar ΔE ni las fuerzas. **VERIFIED** para esa dependencia matemática; la procedencia `conservation.rs`/runner queda como inspección local del orquestador. No es una trayectoria ejecutada por este agente ni prueba de error del integrador o del ensemble.

Esta es una verificación dirigida con reaperturas críticas, metadatos alternativos y evidencia trazada de los bibliógrafos. No es una nueva búsqueda exhaustiva ni revisión humana. Las triangulaciones explícitas separan obras de autores diferentes cuando existen; varios contratos solo tienen como autoridad al proyecto que los implementa. No se calcula un porcentaje de corroboración independiente usando capítulos o espejos como si fueran fuentes independientes.

La aceptación documental de una fuente no afirma que MODENA cumpla sus contratos, que una referencia DFT esté convergida, que un modelo permanezca en dominio o que un resultado publicado se transfiera a otra arquitectura. Las afirmaciones admisibles son las acotadas en las tablas; las condiciones pendientes siguen visibles. Elaborado con asistencia de IA.
