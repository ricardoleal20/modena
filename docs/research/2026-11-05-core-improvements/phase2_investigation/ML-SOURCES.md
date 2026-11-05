# Fuentes anotadas: proveedores ML y validación dinámica

Fase 2, rol `bibliography_agent`. Consulta: 2026-11-05 UTC. Corpus dirigido de 11 fuentes; las páginas relacionadas de una misma documentación se agrupan como una unidad documental. Contexto: [casos C4–C5](CASES.md), [subpregunta 3](../phase1_scope/RESEARCH-QUESTION.md) y puerta G5 del [método](../phase1_scope/METHODOLOGY.md). Este archivo registra evidencia por fuente; no constituye síntesis, priorización ni validación de MODENA.

Las anotaciones se basan en texto primario abierto y documentación oficial consultados mediante navegador de investigación. No se descargaron pesos, no se aceptaron condiciones de acceso y no se ejecutó inferencia. `main`/`latest` son referencias mutables: se registra lo observado, sin atribuirlo a una versión instalada. Los DOI de arXiv identifican depósitos, no revisión por pares. Véase [búsqueda y limitaciones de acceso](ML-SEARCH-LOG.md).

## ML01 — Contrato oficial de FAIR-Chem/UMA

FAIR Chemistry & Collaborators. (s. f.). *Inference using ASE and predictor interface; Hello world with UMA; UMA models; Molecular dynamics with ASE* [Documentación y código de referencia]. Meta FAIR. Recuperado el 5 de noviembre de 2026.

- **Tipo:** documentación oficial y código primario; capacidad de interfaz, no evidencia independiente de exactitud.
- **URLs:** [predictor/calculator](https://raw.githubusercontent.com/facebookresearch/fairchem/main/docs/core/common_tasks/ase_calculator.md), [quickstart](https://raw.githubusercontent.com/facebookresearch/fairchem/main/docs/core/quickstart.md), [tasks y límites](https://raw.githubusercontent.com/facebookresearch/fairchem/main/docs/core/uma.md), [unidades MD](https://facebookresearch.github.io/fairchem/molecular-dynamics/), [código calculator](https://raw.githubusercontent.com/facebookresearch/fairchem/main/src/fairchem/core/calculate/ase_calculator.py).
- **Versión observada:** `fairchem` rama `main`; quickstart usa `uma-s-1p2p1`; SHA y versión de paquete no fijados.

**Anotación.** `FAIRChemCalculator` envuelve `MLIPPredictUnit`; `task_name` determina dominio y nivel de teoría. OMol recibe carga total y multiplicidad mediante `atoms.info`; OMat documenta energía, fuerzas y estrés. La guía MD especifica posiciones/celda en Å, energía en eV, fuerzas en eV/Å y estrés en eV/Å³, con orden Voigt ASE. La periodicidad exige celda válida. Respalda el contrato de C4; limita su generalización porque estrés calculable sin entrenamiento no equivale a estrés validado, OMol se entrenó con datos aperiódicos y las energías de tasks diferentes no son una referencia común.

**Localizadores:** quickstart, flujo de cuatro pasos y advertencia final; `uma.md`, tabla de tasks; guía MD, requisitos NPT y tabla de salidas; calculator, `implemented_properties`, `calculate` y copia de `energy` a `free_energy`. Esta última es compatibilidad ASE: no certifica energía libre electrónica ni PMF. La documentación de tasks señala expresamente que UMA no proporciona work functions y remite a DFT para obtenerlas.

**Frontera por checkpoint:** la tabla de `uma.md` marca `oc22` y `oc25` como disponibles solamente en UMA-1.2, mientras describe los cinco dominios originales por separado. Esto es una afirmación de la documentación, no una inspección del manifiesto del archivo `uma-s-1p2p1.pt`. El calculator consulta las tasks y propiedades del predictor cargado; su lista efectiva permanece pendiente de comprobar sobre pesos fijados.

## ML02 — Pesos UMA, acceso y licencia del código

Meta Platforms, Inc. (2025). *UMA: A family of universal models for atoms* [Model card y FAIR Chemistry License v1, 14 de mayo de 2025]. Hugging Face. Meta Platforms, Inc. (s. f.). *FAIR-Chem MIT license* [Licencia de código].

- **Tipo:** artefactos oficiales de distribución; identificación y términos, sin interpretación jurídica.
- **URLs:** [model card y condiciones visibles](https://huggingface.co/facebook/UMA), [licencia MIT del repositorio](https://raw.githubusercontent.com/facebookresearch/fairchem/main/LICENSE.md).
- **Versión observada:** `uma-s-1.2.1` / `uma-s-1p2p1.pt`, MD5 publicado `3497615fd30a24c5b35cd3b41a682e6e`; no comprobado contra un archivo descargado.

**Anotación.** La tarjeta oficial enumera checkpoints, hashes y variantes, y conserva `uma-s-1` como archivado con un problema conocido de extensividad. El acceso a los archivos requiere aceptar condiciones del repositorio gated. Los materiales UMA se distribuyen bajo FAIR Chemistry License v1 y política de uso; el código del repositorio FAIR-Chem declara MIT. Respalda separar procedencia de pesos y licencia de integración. La disponibilidad pública del código no demuestra permiso equivalente sobre cada checkpoint, y el hash mostrado únicamente identifica lo anunciado por la tarjeta consultada.

## ML03 — Arquitectura, entrenamiento y límites de UMA

Wood, B. M., Dzamba, M., Fu, X., Gao, M., Shuaibi, M., Barroso-Luque, L., Abdelmaqsoud, K., Gharakhanyan, V., Kitchin, J. R., Levine, D. S., Michel, K., Sriram, A., Cohen, T., Das, A., Rizvi, A., Sahoo, S. J., Ulissi, Z. W., & Zitnick, C. L. (2025). *UMA: A family of universal models for atoms* [Preprint, versión 2 revisada el 4 de marzo de 2026]. arXiv. [DOI 10.48550/arXiv.2506.23971](https://doi.org/10.48550/arXiv.2506.23971).

- **Tipo:** investigación primaria; preprint en la versión consultada, no se verificó publicación editorial posterior.
- **Texto abierto:** [HTML v2](https://arxiv.org/html/2506.23971v2); [metadatos y versiones](https://arxiv.org/abs/2506.23971).

**Anotación.** Presenta una GNN con mixture of linear experts y entradas de composición, task, carga y spin. El entrenamiento tiene una etapa rápida de fuerzas directas seguida de eliminación de esa cabeza y ajuste con derivadas automáticas conservativas de energía para fuerzas y estrés. Respalda distinguir arquitectura final de estrategia de preentrenamiento. Sus benchmarks son propios de dominios y datasets definidos. El apéndice de diatómicas muestra mejoras entre versiones y límites de disociación relacionados con carga/spin globales; ni tamaño del entrenamiento ni conservación demuestran cobertura química universal.

**Localizadores:** §2.2.1, §2.4, §6 y apéndices A.1 y G. No se traslada la descripción de UMA a todos los modelos históricos de FAIR-Chem ni a todos sus modos de evaluación.

## ML04 — Interfaces oficiales de MACE y grafos periódicos

MACE Developers. (s. f.). *ASE calculator; MACE in LAMMPS; MACE in LAMMPS with ML-IAP* [Documentación, encabezado 0.3.13]. Recuperado el 5 de noviembre de 2026.

- **Tipo:** documentación oficial de interfaces; rutas documentadas, sin ejecución local.
- **URLs:** [ASE](https://mace-docs.readthedocs.io/en/latest/guide/ase.html), [interfaz original LAMMPS](https://mace-docs.readthedocs.io/en/latest/guide/lammps.html), [ML-IAP](https://mace-docs.readthedocs.io/en/latest/guide/lammps_mliap.html).

**Anotación.** MACE expone `MACECalculator` para ASE. La interfaz original de LAMMPS prepara un modelo TorchScript con `create_lammps_model.py`; `no_domain_decomposition` construye un grafo periódico y evita tratar ghosts como nodos independientes. La ruta ML-IAP requiere conversión `--format=mliap`, invoca Python y documenta aceleración GPU y descomposición para varias GPU. Respalda formatos y semántica periódica concretos. La guía marca ML-IAP como beta y sin múltiples heads; las restricciones de una ruta no deben atribuirse a la otra. Estas páginas no establecen exportación ONNX universal.

**Localizadores:** preparación del modelo y `pair_style` de cada guía; la interfaz original muestra `units metal` y orden explícito de elementos en `pair_coeff`. La recomendación de una sola GPU pertenece a esa interfaz original.

## ML05 — NequIP: tipos, unidades, empaquetado y compilación

NequIP Developers. (s. f.). *ASE; The NequIP workflow; Data fields* [Documentación oficial, `latest`]. Recuperado el 5 de noviembre de 2026.

- **Tipo:** documentación oficial de API y despliegue; versión de paquete y SHA no fijados.
- **URLs:** [ASE](https://nequip.readthedocs.io/en/latest/integrations/ase.html), [workflow](https://nequip.readthedocs.io/en/latest/guide/getting-started/workflow.html), [campos del grafo](https://nequip.readthedocs.io/en/latest/api/data_fields.html).

**Anotación.** La vía actual documenta checkpoints `.ckpt`, paquetes `.nequip.zip` y `nequip-compile --target ase` con AOTInductor para crear `.nequip.pt2`, cargado por `NequIPCalculator.from_compiled_model`. El dispositivo de compilación debe corresponder al de evaluación. El mapeo de especies químicas a tipos del modelo es explícito; las unidades dependen del dataset y admiten factores hacia eV y Å. Los campos distinguen `cell`, `pbc` y `edge_cell_shift`. Respalda contratos específicos de C4. No prueba que todos los checkpoints tengan igual semántica química, formatos intercambiables o exportación ONNX.

## ML06 — Suavidad, vecinos y conservación: eSEN

Fu, X., Wood, B. M., Barroso-Luque, L., Levine, D. S., Gao, M., Dzamba, M., & Zitnick, C. L. (2025). Learning smooth and expressive interatomic potentials for physical property prediction. En *Proceedings of the 42nd International Conference on Machine Learning* (PMLR, Vol. 267, pp. 17875–17893). [Registro editorial](https://proceedings.mlr.press/v267/fu25h.html). DOI del depósito: [10.48550/arXiv.2502.12147](https://doi.org/10.48550/arXiv.2502.12147).

- **Tipo:** artículo primario revisado por pares en ICML 2025; [texto abierto v2](https://arxiv.org/html/2502.12147v2).

**Anotación.** Las ablaciones relacionan diseño del modelo y grafo con conservación y propiedades físicas. Un límite top-K de vecinos puede producir discontinuidades cuando cambia el conjunto vecino; omitir envelopes al cutoff también perjudica conservación, pese a errores estáticos similares. Ser gradiente de una energía no basta para conservarla numéricamente con paso finito: importan suavidad y derivadas acotadas. Respalda examinar la construcción del grafo y dinámica además de MAE. Los resultados son de los modelos y regímenes ensayados; pasar su prueba NVE no certifica por sí solo un ensemble ni propiedades convergidas.

**Localizadores:** §§3 y 5.3; tabla 1; apéndice A.1 para protocolo de conservación.

## ML07 — MAE y calidad de trayectoria: Forces are not Enough

Fu, X., Wu, Z., Wang, W., Xie, T., Keten, S., Gomez-Bombarelli, R., & Jaakkola, T. (2023). Forces are not enough: Benchmark and critical evaluation for machine learning force fields with molecular simulations. *Transactions on Machine Learning Research*. [Publicación](https://openreview.net/forum?id=A8pqQipwkt). DOI del depósito: [10.48550/arXiv.2210.07237](https://doi.org/10.48550/arXiv.2210.07237).

- **Tipo:** benchmark primario revisado por pares; [texto abierto v2](https://arxiv.org/html/2210.07237v2). Registro de arXiv confirma TMLR 2023; el foro OpenReview presentó verificación de navegador.

**Anotación.** Evalúa moléculas, agua, un péptido y un material mediante trayectorias y observables, además de error de fuerzas. Los protocolos incluyen estabilidad, distribuciones de distancias, RDF, difusión y reconstrucción de superficie de energía libre en un caso con metadinámica. Encuentra desacoplamiento entre precisión estática y utilidad dinámica, con fallos de modelos seleccionados. Respalda la distinción entre prueba de snapshots y validación de observables. El benchmark no evalúa los checkpoints UMA actuales; su caso de energía libre tiene condiciones particulares y no autoriza inferir PMF de cualquier MAE.

**Localizadores:** tabla 1 y §§4–7; las referencias de cada dataset y protocolo deben mantenerse al reutilizar un caso.

## ML08 — Fuerzas no conservativas y corrección estructurada

Bigi, F., Langer, M. F., & Ceriotti, M. (2025). The dark side of the forces: Assessing non-conservative force models for atomistic machine learning. En *Proceedings of the 42nd International Conference on Machine Learning* (PMLR, Vol. 267, pp. 4384–4414). [Registro editorial](https://proceedings.mlr.press/v267/bigi25a.html). DOI del depósito: [10.48550/arXiv.2412.11569](https://doi.org/10.48550/arXiv.2412.11569).

- **Tipo:** estudio primario teórico y computacional revisado por pares; [texto abierto v6](https://arxiv.org/html/2412.11569v6), revisión de octubre de 2025, posterior a la conferencia.

**Anotación.** Analiza asimetría del Jacobiano, trabajo no conservativo y sus efectos en optimización, NVE, NVT y propiedades temporales. Los ejemplos muestran que termostatizar no elimina todos los sesgos de fuerzas directas. Examina preentrenamiento seguido de ajuste conservativo y un esquema multiple time stepping con correcciones conservativas periódicas. Respalda distinguir fuerza directa de gradiente y reconocer una alternativa híbrida formalmente especificada. El estudio no valida alternar arbitrariamente proveedores por umbral de incertidumbre; sus resultados de aceleración y calidad pertenecen al esquema y sistemas evaluados.

**Localizadores:** §§3, 4.4–4.8 y apéndice I. La versión de arXiv consultada y la publicación de conferencia se identifican separadamente.

## ML09 — Active learning con etiquetas electrónicas: DP-GEN

Zhang, Y., Wang, H., Chen, W., Zeng, J., Zhang, L., Wang, H., & E, W. (2020). DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models. *Computer Physics Communications, 253*, 107206. [DOI 10.1016/j.cpc.2020.107206](https://doi.org/10.1016/j.cpc.2020.107206).

- **Tipo:** artículo primario de método/software; [texto abierto del manuscrito v1](https://arxiv.org/html/1910.12690v1), [metadatos con referencia publicada](https://arxiv.org/abs/1910.12690).

**Anotación.** DP-GEN itera exploración, selección de configuraciones, etiquetado con métodos electrónicos y entrenamiento. Usa desviación entre fuerzas de modelos para seleccionar candidatos entre límites de confianza definidos por el usuario, y demuestra el flujo en cobre. Documenta motores separados para exploración y etiquetas, incluidos LAMMPS y Quantum ESPRESSO. Respalda un ciclo de adquisición y trazabilidad en C4–C5. Los umbrales del ejemplo no son universales, la desviación no es una garantía calibrada de error, y generar etiquetas DFT continúa requiriendo el solver y sus parámetros físicos.

**Localizadores:** §II y ejemplo de §IV. La fecha de publicación procede del registro editorial enlazado por arXiv; el HTML muestra una fecha de compilación que no se usa como año bibliográfico.

## ML10 — Calibración explícita de incertidumbre en energía y fuerzas

Busk, J., Schmidt, M. N., Winther, O., Vegge, T., & Jørgensen, P. B. (2023). Graph neural network interatomic potential ensembles with calibrated aleatoric and epistemic uncertainty on energy and forces. *Physical Chemistry Chemical Physics, 25*(37), 25828–25837. [DOI 10.1039/D3CP02143B](https://doi.org/10.1039/D3CP02143B).

- **Tipo:** estudio primario de calibración; [texto abierto v2](https://arxiv.org/html/2305.16325v2), [registro bibliográfico PubMed](https://pubmed.ncbi.nlm.nih.gov/37724552/). El sitio editorial RSC devolvió 403; los resultados científicos se extrajeron del manuscrito abierto, no de PubMed.

**Anotación.** Propone ensembles GNN con componentes epistémica y aleatoria y recalibración posterior mediante una función no lineal. Evalúa energías y fuerzas en ANI-1x y Transition1x, con geometrías alejadas del equilibrio. Respalda medir correspondencia entre incertidumbre y error observado, en vez de interpretar directamente dispersión entre modelos como margen fiable. Su calibración se evalúa sobre particiones concretas y no ofrece una garantía sobre cualquier trayectoria futura o química nueva. Menciona detener o recurrir a mayor fidelidad; no demuestra que un cambio instantáneo de fuerzas preserve Hamiltoniano o muestreo.

**Localizadores:** §§2.5–2.6, §3 y discusión; DOI, autores, año, volumen y páginas corroborados mediante el registro bibliográfico enlazado.

## ML11 — Evidencia contraria a confiar en ensembles sin calibrar

Kahle, L., & Zipoli, F. (2022). Quality of uncertainty estimates from neural network potential ensembles. *Physical Review E, 105*(1), 015311. [DOI 10.1103/PhysRevE.105.015311](https://doi.org/10.1103/PhysRevE.105.015311).

- **Tipo:** estudio primario revisado por pares; [metadatos editoriales](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.105.015311), [manuscrito abierto v2](https://arxiv.org/pdf/2108.05748v2). El cuerpo editorial requiere acceso; se consultó el manuscrito de autores.

**Anotación.** Compara incertidumbre estimada y errores reales fuera de distribución mediante cuatro estudios de ensembles de potenciales neuronales. Observa sobreconfianza y subestimación del error, y presenta evidencia favorable a determinadas aproximaciones bayesianas por muestreo de parámetros. Respalda evaluar calibración por sistema y arquitectura antes de usarla para decisiones de active learning o recuperación. No afirma que todo ensemble falle ni que un método alternativo sea universalmente fiable. Sus casos no son una evaluación directa de UMA; conservar esa frontera evita convertir una advertencia metodológica en un ranking no medido.

**Localizadores:** resumen, casos computacionales y conclusiones del manuscrito; título publicado sin el prefijo “On the” que conserva arXiv.

Elaborado con asistencia de IA. No existe atestación de lectura humana ni adquisición local del corpus; la comprobación aquí descrita es consulta de fuentes web, sujeta a los límites de acceso registrados.
