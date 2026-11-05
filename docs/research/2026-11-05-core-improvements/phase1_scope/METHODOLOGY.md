# Metodología propuesta para la investigación del core de MODENA

**Estado:** Fase 1, diseño metodológico; pendiente de confirmación del alcance por el usuario antes de investigar fuentes externas. Esta entrega define cómo evaluar mejoras; no presenta resultados nuevos de rendimiento ni una revisión bibliográfica concluida.

**Objeto:** identificar cambios concretos del algoritmo y del core de MODENA que permitan sustituir herramientas existentes en casos de uso delimitados. La pregunta definitiva y sus límites quedan en [RESEARCH-QUESTION.md](RESEARCH-QUESTION.md). La unidad de comparación será un flujo científico reproducible, con Hamiltoniano, precisión, recursos y observables especificados, en lugar del nombre completo de una plataforma.

**Base de inspección:** checkout `f166399070b2d6ebfd91f51e5934c132f04d46c0`. La procedencia de cada medición futura debe registrar también el estado del árbol y los cambios no confirmados. Las mediciones archivadas conservarán los commits y condiciones que les corresponden; no se atribuirán al checkout actual.

## Paradigma y diseño

**Paradigma pragmatista:** la decisión que debe apoyar el estudio es qué mejorar primero para resolver un trabajo científico concreto. Se combinará análisis cualitativo de arquitectura y evidencia técnica con especificaciones cuantitativas de experimentos de aceptación. La revisión identifica candidatos y sus mecanismos; las campañas que no se ejecuten se presentarán como experimentos propuestos.

**Diseño:** estudio comparativo técnico con auditoría de implementación, revisión dirigida de fuentes primarias y análisis de brechas. No se registrará como revisión sistemática exhaustiva ni metaanálisis. No hay poblaciones clínicas ni tamaños de efecto combinables; PRISMA o GRADE no deben utilizarse como sellos de calidad fuera de su propósito. El estándar de reporte será una memoria técnica trazable: preguntas, búsquedas, selección, evidencias, contradicciones, límites y protocolo de reproducción.

La distinción entre implementación, pruebas y validación material seguirá [el alcance vigente](../../../VALIDATED_SCOPE.md), que advierte que una prueba o una simulación terminada no certifican propiedades convergidas. El [README](../../../../README.md) da contexto de interfaces; [la arquitectura](../../../../ARCHITECTURE.md) y [el roadmap](../../../../ROADMAP.md) se tratarán como documentos de diseño cuya vigencia debe comprobarse en código.

## Datos y trazabilidad

| Clase | Material previsto | Uso permitido |
|---|---|---|
| Evidencia local actual | Código, tests, configuraciones, runners y contratos públicos del checkout | Demostrar qué ruta existe, qué acepta y cómo calcula; citar archivo y línea |
| Evidencia local histórica | JSON, logs y reportes con revisión y hardware identificados | Establecer antecedentes; conservar fecha, versión y alcance original |
| Documentación oficial externa | Manuales y repositorios oficiales de los motores y modelos evaluados | Confirmar interfaces, capacidades, restricciones y versiones consultadas |
| Literatura primaria | Artículos originales, publicaciones de métodos y preprints verificables | Sustentar algoritmos, condiciones de validez y mecanismos de mejora |
| Experimento nuevo, si se ejecuta | Artefactos crudos de una prueba o campaña reproducible | Responder una hipótesis acotada con límites, procedencia y resultado negativo incluido |

La investigación externa comenzará tras la confirmación de Fase 1. Antes de citar, se verificará título, autores, año, enlace o DOI, versión y que el texto respalde la afirmación asociada. Una referencia no verificable se excluye de la evidencia, incluso si aparece citada en notas anteriores. Se registrará acceso incompleto. La documentación oficial sirve para capacidad de software; un preprint se identificará como tal y no recibirá el estado de artículo revisado por pares.

Los reportes [científicos](../../../scientific-hardening/RESULTS.md), [de verificación](../../../scientific-hardening/VERIFY.md) y [de benchmarks](../../../../bench/BENCHMARKS.md) guían el diseño. Las propuestas previas de [rendimiento](../../../PERF_ROADMAP.md), [comparación de motores](../../../MATCH_ENGINES.md) y [comité](../../../COMMITTEE_ACTA.md) serán antecedentes sujetos a revalidación: no se sumarán sus proyecciones de aceleración ni se tratarán como mediciones nuevas.

## Estrategia de búsqueda propuesta

Se cubrirán trabajos fundamentales cuando expliquen un método y literatura reciente para implementaciones o modelos actuales. La ventana se cierra en la fecha de consulta; se registrarán búsquedas y versiones en un ledger. Se dará prioridad al texto primario, al manual oficial vigente y al código de la versión documentada. Los resúmenes de buscador sirven para descubrir fuentes, no para respaldar conclusiones.

| Eje | Cadenas iniciales orientativas | Fuentes primarias a consultar |
|---|---|---|
| MD clásica y rendimiento | `molecular dynamics neighbor list memory cluster pair SIMD`; `LAMMPS neighbor list atom sorting performance`; `GPU resident molecular dynamics neighbor list` | Manual y código LAMMPS; artículos originales de algoritmos; GROMACS/OpenMM/HOOMD-blue cuando aporten un mecanismo comparable |
| Física y muestreo | `molecular dynamics validation ensemble sampling timestep convergence`; `constraint virial rigid water NPT`; `smooth particle mesh Ewald force energy consistency stress` | Publicaciones originales de integradores, electrostática y restricciones; documentación de implementaciones de referencia |
| DFT e interoperabilidad | `Quantum ESPRESSO PWscf forces stress molecular dynamics convergence`; `ASE Quantum ESPRESSO calculator units stress` | Manuales y código oficiales de Quantum ESPRESSO y ASE; publicaciones del método aplicable |
| Potenciales aprendidos | `FAIR-Chem calculator molecular dynamics energy conservation forces stress`; `universal machine learning interatomic potential out of distribution uncertainty molecular dynamics`; `MACE NequIP equivariant potential conservative forces` | Documentación, repositorios, model cards y artículos originales de FAIR-Chem y de los modelos pertinentes |
| Evaluación ML | `machine learning interatomic potentials simulation stability validation`; `active learning interatomic potentials uncertainty calibration` | Estudios primarios de estabilidad, cobertura química y validación dinámica; datasets y protocolos públicos |

Se buscará también evidencia contraria: regresiones por precisión mixta, sobrecoste de listas y reconstrucciones, fallos de extrapolación química, conservación insuficiente para muestreo correcto y límites de incertidumbre de ensembles. Cada candidato debe registrar por qué se incluyó o excluyó. Los nombres y familias de modelos no se asumirán vigentes hasta comprobarlos en sus fuentes oficiales.

## Marco de análisis

1. **Inventariar rutas alcanzables.** Separar core, configuración TOML, CLI, Python, MCP y ejemplos. Una feature compilable no prueba que el runner de usuario la seleccione. Registrar requisitos de caja, topología, precisión, hardware y combinación de controladores.
2. **Definir casos de sustitución.** Para cada flujo, documentar entrada, física exigida, salida y nivel de confianza. Contrastar MD clásica, cálculo electrónico y evaluación de potenciales ML como necesidades distintas. La ausencia de una capacidad se confirma por inspección y documentación; no se infiere de un único archivo.
3. **Mapear mecanismos a código.** Cada mejora deberá señalar punto de integración, cambio del algoritmo o estructura de datos, invariantes que puede alterar, dependencias y evidencia que falta. Distinguir cuello de botella medido de hipótesis de perfilado.
4. **Triangular.** Contrastar contrato matemático, implementación y prueba independiente. El mismo algoritmo implementado dos veces no constituye necesariamente un oráculo independiente. Derivadas, identidades analíticas o referencias convergidas deben acompañar la comparación entre motores.
5. **Priorizar por dependencias.** Ordenar primero corrección y confianza científica, luego cobertura de casos de uso y rendimiento medido. Examinar la sensibilidad del orden a supuestos de esfuerzo y hardware. No multiplicar aceleraciones hipotéticas de optimizaciones relacionadas.
6. **Revisar adversarialmente.** Solicitar objeciones a alcance, síntesis y reporte final según los tres checkpoints de la habilidad. Resolver problemas que invaliden conclusiones antes de continuar; conservar límites y desacuerdos restantes.

## Matriz de evidencia y confianza

Cada hallazgo o recomendación tendrá una fila con:

| Campo | Registro exigido |
|---|---|
| Identidad | ID estable, pregunta asociada y caso de uso |
| Observación | Qué se vio, con archivo/línea o fuente y sección |
| Procedencia | Commit, versión externa, fecha de acceso y artefacto cuando exista |
| Categoría | Hecho observado; resultado histórico; inferencia técnica; propuesta por comprobar |
| Alcance | Potencial, sistema, caja, ensemble, precisión, hardware y ruta de ejecución |
| Evidencia contraria | Excepciones, resultados negativos o fuentes discrepantes |
| Confianza | Alta, media o baja, acompañada de razón explícita |
| Acción | Cambio propuesto, dependencia, riesgo y prueba que podría refutar su utilidad |

**Confianza alta:** observación directa reproducible para el alcance citado, o capacidad documentada y corroborada en la versión inspeccionada. **Media:** soporte primario pertinente con transferencia todavía no probada a MODENA. **Baja:** hipótesis arquitectónica o de rendimiento pendiente de medición. La confianza de una recomendación no aumenta por repetirla en varios documentos que dependen de la misma fuente.

## Priorización propuesta

Se evaluarán impacto, esfuerzo, riesgo y dependencia por separado; cualquier puntuación es juicio de ingeniería, no probabilidad de éxito ni resultado experimental.

| Dimensión | Criterios |
|---|---|
| Impacto científico | Corrige un sesgo; habilita un observable fiable; amplía cobertura química o un ensemble requerido |
| Impacto en uso | Completa un flujo reproducible; permite interoperar; reduce trabajo manual con errores |
| Impacto en recursos | Reduce coste o memoria bajo condiciones emparejadas; mantiene precisión y estabilidad |
| Esfuerzo | Archivos afectados, cambios de contrato, dependencias nativas, migración y trabajo de validación; expresar intervalos o tamaños relativos |
| Riesgo | Alteración del Hamiltoniano, ensemble, identidad de átomos, determinismo, restart, portabilidad o compatibilidad |
| Dependencias | Qué condición debe cumplirse antes y cuál es la alternativa si el experimento fracasa |

La tabla final deberá incluir prioridad sugerida, objetivo, cambio concreto, evidencia, aceptación y alternativa. No se asignarán plazos exactos sin conocer recursos humanos ni se prometerá paridad global por completar una lista de características.

Guardas añadidas tras el [Checkpoint 1](CHECKPOINT-1.md): antes de evaluar resultados se congelará una ficha por caso con entrada, Hamiltoniano/modelo, observable, requisito operativo, tolerancia y motivo de selección. Se conservará al menos un caso hoy bloqueado y se registrarán cambios de selección. Cada recomendación se comparará, cuando sea pertinente, con una alternativa de coste mínimo que satisfaga la necesidad mediante integración, validación o mejora del flujo sin modificar el algoritmo. Las prioridades seguirán siendo provisionales mientras no se confirme un caso material específico del usuario.

## Experimentos emparejados y puertas de aceptación

Estas puertas son especificaciones propuestas. Sus tolerancias numéricas se fijarán por caso y precisión antes de medir; no se copiará un umbral universal para sistemas con escalas distintas.

| Puerta | Diseño de verificación previsto | Evidencia de aceptación |
|---|---|---|
| G0 — Identidad y reproducibilidad | Congelar entradas y condiciones; registrar revisión, dirty state, compilador, flags, binario, hardware, hilos/ranks, semillas, precisión y unidades | Manifiesto completo; logs y fallos conservados; métricas ausentes como `null` |
| G1 — Fuerzas y celda | Coordenadas idénticas; energías por término, fuerzas y virial/tensor cuando corresponda; diferencias finitas con barrido de paso y pruebas de simetría | Errores absolutos y normalizados dentro de tolerancias predefinidas; estabilidad frente al paso de diferencia finita y a convergencia de referencia |
| G2 — Integración y muestreo | Tiempo físico fijo para barrido de timestep; distribuciones de control analíticas; energía o invariante correcto para el ensemble; restricciones y grados de libertad coherentes | Orden de convergencia compatible; momentos/distribuciones adecuados con incertidumbre y autocorrelación consideradas |
| G3 — Flujos materiales | LJ, sistema iónico y agua como familias iniciales; después química o material objetivo; coincidencia de Hamiltoniano, correcciones, condiciones y precisión | Observables relevantes —por ejemplo estructura, densidad, presión o difusión— con réplicas y análisis de convergencia; no solo una trayectoria estable |
| G4 — Rendimiento | Baseline conservado y candidato sobre mismos snapshots, pasos y recursos; microbench, integración y proceso completo reportados por separado | Distribución de tiempos, memoria y reparto de costes; mejora reproducible con física de G1–G3 preservada; casos desfavorables incluidos |
| G5 — ML y DFT | Verificar especies, unidades, periodicidad, energía, fuerza, estrés, cutoff, versión/pesos del modelo; referencia DFT convergida y datasets fuera de entrenamiento cuando sean verificables | Precisión estática, estabilidad dinámica y observables bajo dominio declarado; trazabilidad de datos y pesos; incertidumbre calibrada si se utiliza para decidir |
| G6 — Integración de producto | Ejecución desde interfaz de usuario, cancelación, fallos, restart y recuperación; misma identidad física y controladores | Pruebas del flujo completo y mensajes verificables; feature experimental claramente diferenciada de ruta validada |

Para benchmarks de MD se reutilizarán los runners documentados en [bench/BENCHMARKS.md](../../../../bench/BENCHMARKS.md) cuando su ventana responda a la pregunta. Se ampliarán solo si hace falta aislar un coste no observable. Inicializaciones con semillas iguales pero RNG diferentes no garantizan microestados iguales. Un benchmark de fuerza exigirá el mismo snapshot; un benchmark de throughput con diferentes microestados declarará esa limitación.

Una comparación de electrostática fijará precisión objetivo y parámetros convergidos, además de tiempo; usar nombres de métodos parecidos no demuestra error equivalente. Una comparación GPU registrará sincronización, transferencias y qué estado reside en dispositivo. Una campaña de escala distinguirá escalado fuerte y débil, comunicación y trabajo útil. Las ejecuciones en máquinas diferentes no demostrarán una aceleración atribuible exclusivamente al algoritmo.

## Validez, límites y ética

La validez interna depende de mantener la misma física al comparar rendimiento. La validez externa depende de cubrir regímenes variados: tamaño, densidad, química, caja, restricciones y recursos. La fiabilidad exige artefactos crudos, versiones y reproducción; la validez de inferencia exige separar resultados ejecutados de hipótesis y diseños futuros.

La revisión dirigida puede omitir alternativas y no cuantifica toda la literatura. Las capacidades externas pueden cambiar; se indicará fecha y versión. Los resultados archivados pertenecen a su entorno. Las estimaciones de esfuerzo requieren revisión del mantenedor. Un conjunto reducido de benchmarks no certifica una herramienta para toda MD, toda DFT o toda química. Un potencial ML no hereda fiabilidad en un dominio por su nombre o tamaño; la cobertura se evaluará por caso con fuentes verificadas.

El estudio utiliza código y documentación técnica, sin participantes humanos ni acceso previsto a datos personales. Se declarará el uso de IA en la investigación, se respetarán atribución y licencias, y no se fabricarán referencias, mediciones ni autorización de uso de modelos/datasets. No se requiere un plan de administración de sujetos humanos para este diseño.

**Prerregistro:** no requerido para esta revisión exploratoria. Se recomienda congelar hipótesis y criterios antes de cualquier campaña confirmatoria derivada. Declaración de artefacto completado: `not_provided`; companion: ninguno. Este documento no crea un sidecar de prerregistro ni calcula su procedencia.

**Siguiente paso:** tras confirmar pregunta, alcance y metodología, ejecutar búsqueda y verificación de fuentes de Fase 2; después producir síntesis y recomendaciones con sus revisiones. No se ha realizado investigación externa en esta entrega.
