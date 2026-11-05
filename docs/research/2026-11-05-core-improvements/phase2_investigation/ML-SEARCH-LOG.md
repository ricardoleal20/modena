# Registro de búsqueda: potenciales aprendidos

Fase 2, `bibliography_agent`; consulta efectuada el 5 de noviembre de 2026 UTC. Salida asociada: [11 fuentes anotadas](ML-SOURCES.md). Alcance dirigido a C4–C5 y al contrato de proveedor ML; no revisión sistemática PRISMA, censo exhaustivo ni experimento. Se usó búsqueda web para descubrir y navegación a textos primarios/documentación oficial para extraer evidencia. No se inspeccionaron fases posteriores ni se modificó el core.

## Parámetros definidos antes de seleccionar

| Parámetro | Criterio |
| --- | --- |
| Pregunta | Interfaces oficiales actuales FAIR-Chem/UMA, MACE y NequIP; periodicidad, unidades y metadatos; conservación/suavidad; relación entre error estático y MD; active learning y calibración. |
| Fuentes | Repositorios y documentación de desarrolladores, model card oficial, arXiv de autores, PMLR/ICML, OpenReview/TMLR y editoriales de investigación primaria. PubMed se usa para metadatos, no como evidencia científica independiente. |
| Tiempo | Sin límite inferior para métodos pertinentes; software actual a fecha de consulta. La carpeta de investigación no actúa como corte bibliográfico. |
| Idioma | Búsqueda técnica en inglés; anotaciones en español. |
| Inclusión | Texto primario abierto u otra versión abierta identificada; pertinencia directa; existencia y metadatos comprobables. Documentación sirve para capacidad y formato, artículos para resultados del método. |
| Exclusión | Blogs, agregadores y redes como respaldo final; reviews generales como evidencia primaria; resúmenes sin texto pertinente accesible; promesas de marketing; resultados de otro modelo tratados como prueba de UMA. |
| Evidencia contraria buscada | Fallos dinámicos con MAE bajo, discontinuidades por grafos/cutoff, sesgos no conservativos, sobreconfianza de ensembles y limitaciones de exportación. |

No se recibió un corpus bibliográfico externo estructurado: el contexto preexistente fueron los documentos locales de Fase 1 y `CASES.md`. La selección conserva una unidad por artículo y agrupa páginas oficiales de una misma familia documental; no cuenta preprint y publicación como dos estudios.

## Consultas ejecutadas y observaciones

Las cadenas siguientes son las consultas efectivamente enviadas; el buscador no entregó totales completos ni estables. Se conserva el resultado cualitativo observado y la decisión posterior. No se inventan cantidades de registros identificados, cribados o excluidos.

| Lote | Consultas exactas | Resultado observado y acción |
| --- | --- | --- |
| A | `site.fair-chem.github.io UMA calculator task_name charge spin pbc direct forces conservative units`; `site.huggingface.co facebook UMA license model card`; `machine learning force fields energy conservation smoothness molecular dynamics Forces are not enough benchmark` | Aparecieron repositorio FAIR-Chem, documentación calculator, tarjeta UMA, Fu y Bigi. Se abrieron las fuentes oficiales y los manuscritos. ML01–ML03, ML07–ML08. |
| B | `site.mace-docs.readthedocs.io ASE calculator units LAMMPS export compile periodic`; `site.nequip.readthedocs.io ASE calculator nequip-compile package periodic AtomicData`; `"machine learning" "potential" "uncertainty" "calibration" ensemble molecular dynamics`; `"energy conservation" "smoothness" "machine learning" "potentials"` | Se recuperaron guías oficiales NequIP y MACE, estudios de calibración y material de conservación. Se siguieron enlaces ASE/LAMMPS y se descartaron reviews generales como soporte primario. ML04–ML05, ML10. |
| C | `"Learning smooth and expressive interatomic potentials"`; `"DP-GEN" "Computer Physics Communications" 2020`; `"uncertainty" "interatomic potentials" "overconfident" ensemble`; `"uncertainty" "potentials" "adversarial" "Schwalbe"` | Se localizaron eSEN, DP-GEN, Kahle/Zipoli, Lu et al., Bilbrey et al. y Schwalbe-Koda et al. Se retuvieron ML06, ML09 y ML11. Los otros trabajos quedaron como candidatos descritos abajo. |
| D | `site.facebookresearch.github.io/fairchem "eV/" "stress"`; `site.facebookresearch.github.io/fairchem "direct_forces"`; `site.nequip.readthedocs.io "edge_cell_shift" "cell"`; `site.mace-docs.readthedocs.io "ML-IAP"` | Se recuperó tabla oficial de unidades MD, campos del grafo NequIP y guía MACE ML-IAP. También surgió un issue de autograd FAIR-Chem; no se aceptó como demostración de un fallo vigente sin reproducción y revisión fijada. |
| E | `"Learning Smooth and Expressive Interatomic Potentials" "proceedings"`; `"Graph Neural Network Interatomic Potential Ensembles" "2023" "26100"`; `"Uncertainty quantification for neural network potential foundation models" "Bilbrey"`; `"Forces are not Enough" "openreview.net/forum"` | Se confirmó ICML para eSEN mediante PMLR; el número de página incluido como término de búsqueda no se tomó como metadato. TMLR quedó corroborado por arXiv y PDF editorial indexado. Bilbrey se mantuvo como candidato con acceso editorial incompleto. |
| F | `"10.1039/D3CP02143B" "258"`; `"10.1039/D3CP02143B" "259"`; `"Graph neural network interatomic potential ensembles" "Phys. Chem. Chem. Phys."`; `site.facebookresearch.github.io/fairchem "ONNX"` | PubMed permitió corroborar páginas 25828–25837 y DOI de Busk. La búsqueda ONNX devolvió material irrelevante, archivos de configuración genéricos y documentación general de ONNX; no se encontró con ella una exportación universal de MLIP respaldada. Esta ausencia de recuperación no prueba imposibilidad técnica. |

## Navegación y extracción registradas

| Fuente retenida | Superficie primaria consultada | Elementos comprobados |
| --- | --- | --- |
| ML01 | `facebookresearch/fairchem` → `docs/core/quickstart.md`, `uma.md`, `common_tasks/ase_calculator.md`; `src/fairchem/core/calculate/ase_calculator.py`; guía MD renderizada | Predictor/calculator, task, propiedades implementadas, carga/multiplicidad, validación de celda, unidades/Voigt y límites de los tasks. |
| ML02 | `huggingface.co/facebook/UMA`; `fairchem/main/LICENSE.md` | Identidad del propietario Meta, licencia visible, acceso gated, tabla de checkpoints y MD5 publicado; código MIT separado. |
| ML03 | arXiv `2506.23971v2`, HTML y ficha | Autores, DOI de depósito, revisiones 2025/2026; §2.4 y apéndice G. |
| ML04 | MACE `latest`, encabezado 0.3.13; páginas `ase.html`, `lammps.html`, `lammps_mliap.html` | Formatos y comandos de conversión, grafo periódico original frente a nueva ruta ML-IAP, límites de heads y GPU. |
| ML05 | NequIP `latest`: ASE, workflow y data fields | Compilación AOTInductor, tipos/especies, factores de unidades y desplazamientos periódicos del grafo. |
| ML06 | PMLR `v267/fu25h.html`; arXiv `2502.12147v2` | Autores, conferencia/páginas y evidencia de ablaciones de top-K/cutoff/envelopes. |
| ML07 | arXiv `2210.07237v2`, HTML/ficha; PDF TMLR indexado | Autores/año, TMLR, tabla de sistemas y métricas dinámicas. El foro presentó challenge. |
| ML08 | PMLR `v267/bigi25a.html`; arXiv `2412.11569v6` | Autores, conferencia/páginas; fuerzas directas, sesgos dinámicos y MTS con correcciones conservativas. |
| ML09 | arXiv `1910.12690v1`, HTML/ficha | Autores, DOI de revista y año 2020; ciclo exploración/etiquetado/entrenamiento. |
| ML10 | arXiv `2305.16325v2`, HTML/PDF; PubMed `37724552` | Autores y método en texto primario; año/DOI/volumen/páginas en registro bibliográfico. |
| ML11 | APS DOI `10.1103/PhysRevE.105.015311`; PDF arXiv `2108.05748v2` | Autores y fecha editorial; cuatro casos y alcance de sobreconfianza en manuscrito abierto. |

## Incidencias y decisiones de cribado

- **Desalineación documental FAIR-Chem:** la guía MD renderizada observada aún usa `turbo`; `ase_calculator.md` de `main` afirma que ese modo fue retirado y presenta Nitro como predeterminado. Se conservaron ambas procedencias, sin convertirlas en una receta de versión única. El lector debe fijar release/commit antes de reproducir. [Guía MD](https://facebookresearch.github.io/fairchem/molecular-dynamics/), [guía predictor `main`](https://raw.githubusercontent.com/facebookresearch/fairchem/main/docs/core/common_tasks/ase_calculator.md).
- **Tasks por checkpoint:** `uma.md` conserva una explicación de cinco datasets y añade `oc22`/`oc25` en su tabla con la marca “Only available in UMA-1.2”. No se convirtió esa tabla en una lista runtime de todos los checkpoints. No se cargó el manifiesto `dataset_to_tasks` del archivo de pesos; verificarlo queda pendiente.
- **Acceso alternativo:** `.../fairchem/main/LICENSE` devolvió 404; `LICENSE.md` resolvió. El archivo raw de la guía MD devolvió cache miss; la versión renderizada fue accesible. RSC devolvió 403 para HTML, landing y PDF de Busk; se usó arXiv para texto y PubMed para metadatos. APS expuso metadatos pero pidió acceso para cuerpo editorial; se usó manuscrito arXiv. OpenReview presentó challenge en los foros consultados; se usaron PMLR y arXiv según correspondía.
- **Candidatos no retenidos:** [Schwalbe-Koda et al., 2021](https://www.nature.com/articles/s41467-021-25342-8.pdf), active learning por muestreo adversarial, era pertinente pero redundante para el corpus acotado; no se extrajeron resultados adicionales. [Lu et al., arXiv:2309.00195](https://arxiv.org/abs/2309.00195), ENN y sobreconfianza, quedó como ampliación pendiente; no se comprobó su texto completo en esta búsqueda. [Bilbrey et al., 2025](https://www.nature.com/articles/s41524-025-01572-y), UQ para foundation models, apareció con DOI y resultados indexados, pero abrir artículo/PDF llevó a autenticación/500; no se incorporó como evidencia de texto completo consultado. Reviews de MLFF y agregadores no se emplearon para respaldar afirmaciones.
- **ONNX:** no se presentó una inferencia de imposibilidad. Las rutas positivas documentadas en ML04–ML05 son específicas de PyTorch/TorchScript/AOTInductor o ML-IAP; no se extrapolaron a una exportación universal. Un `.gitattributes` que enumera extensiones `.onnx` no demuestra que exista un checkpoint exportado.
- **Fallback:** las menciones de recuperación mediante métodos de mayor fidelidad en ML10 y de MTS en ML08 se registraron por separado. Ninguna fuente consultada establece que conmutar fuerzas instantáneamente por un umbral de desacuerdo conserve la dinámica objetivo. Esa afirmación permanece sin respaldo en este corpus.

## Deduplicación y cobertura

Se intentaron consultas a Semantic Scholar Graph API para `ARXIV:2506.23971`, `ARXIV:2502.12147`, `ARXIV:2210.07237`, `ARXIV:2412.11569` y los DOI `10.1016/j.cpc.2020.107206`, `10.1039/D3CP02143B`, `10.1103/PhysRevE.105.015311` y `10.1038/s41524-025-01572-y`, solicitando `paperId,title,year,externalIds`. La herramienta no pudo acceder a esas URLs. Se aplica degradación explícita: `s2_unresolved`; deduplicación por título, autores, arXiv y DOI. No se inventan IDs.

El corpus retenido contiene cuatro unidades documentales oficiales y siete artículos/manuscritos primarios, de los cuales seis tienen publicación revisada por pares comprobada y uno se conserva como preprint. Las siete fuentes de investigación corresponden a años 2020, 2022, 2023 (dos) y 2025 (tres); UMA incluye revisión 2026. No hay un año que alcance 70% de las fuentes con año conocido. No se infirió distribución geográfica a partir de afiliaciones. La concentración computacional/técnica es intencional por la pregunta; no representa evidencia experimental de materiales independiente ni permite un meta-análisis.

**DISTRIBUTIONAL_SKEW_ADVISORY:** método computacional/benchmark de software = 7/7 artículos primarios con método conocido (100%). Es señal de cobertura, no defecto ni validación. No se amplió a estudios experimentales porque esta tarea documenta interfaces y criterios de validación, con la transferencia material pendiente.

La incertidumbre restante incluye correspondencia exacta entre documentación mutable y release, rendimiento/paridad en el transporte de MODENA y comportamiento de un checkpoint fijado sobre C4. No se midió latencia, memoria, fuerzas, estrés, estabilidad ni calibración. Elaborado con asistencia de IA; ninguna entrada implica lectura humana o archivo original adquirido localmente.
