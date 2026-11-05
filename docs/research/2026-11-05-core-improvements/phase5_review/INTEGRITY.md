# Revisión de integridad de la investigación

**Veredicto de integridad: CLEARED.** Fecha: 2026-11-05 UTC. Rol: `ethics_review_agent`, Fase 5. Se revisaron [DRAFT.md](../phase4_composition/DRAFT.md), [ACCEPTANCE-BACKLOG.md](../phase3_analysis/ACCEPTANCE-BACKLOG.md), [SOURCE-VERIFICATION.md](../phase2_investigation/SOURCE-VERIFICATION.md) y los cuatro corpus bibliográficos enlazados por el borrador. Este dictamen cubre atribución, transparencia y representación de evidencia; no valida la física, el rendimiento, una licencia de uso concreta ni la aceptación de las recomendaciones.

No se detectó fabricación de referencias, ocultación del uso de IA, apropiación de resultados como mediciones propias ni tergiversación sistemática en los pasajes comprobados. No hay correcciones obligatorias de integridad para entregar el informe con sus límites actuales. Los detalles editoriales o de localización corregibles no se clasifican como faltas críticas.

## Dimensiones examinadas

| Dimensión | Estado | Observación y alcance |
| --- | --- | --- |
| Divulgación de IA | **pass** | El método se declara asistido por IA; la declaración final enumera búsqueda, lectura asistida, comprobación documental, síntesis, redacción y crítica. Los límites identifican a los revisores como agentes de la misma sesión. No se atribuye revisión experta humana ni se confunde confirmación del alcance con revisión de fuentes. |
| Integridad de atribución | **pass**, con límites de cobertura | Se contrastaron identidad, versión y uso de 21/37 unidades bibliográficas como se detalla abajo. Las publicaciones y anexos dan crédito a autores/proyectos; los mirrors y manuscritos no se presentan como estudios adicionales. No se detectaron citas textuales extensas ni apropiación literal en los pasajes comparados. No se ejecutó un detector exhaustivo de plagio. |
| Uso dual | **pass** | Riesgo **bajo** para el contenido concreto: métodos generales de simulación, interoperabilidad y validación de software. No contiene un procedimiento operativo de daño que justifique una restricción o declaración especial. |
| Representación justa | **pass** | LAMMPS se trata como referencia de MD; QE conserva cálculo electrónico; FAIR-Chem aporta modelos y herramientas sujetos a task/checkpoint. Se mantienen alternativas de integración y uso externo. Las limitaciones de MODENA se describen sin convertirlas en prueba de fuerzas incorrectas o inferioridad universal. |
| Ética de datos y materiales | **pass** | El objeto son código, manuales y publicaciones públicas. El texto no afirma adquisición de datasets/pesos ni aceptación de términos gated. El acceso web al corpus no se transforma en permiso para redistribuir pesos. No hay participantes, datos de sujetos ni actividad con datos personales en este estudio. |
| Conflictos e intereses | **pass**, comprobación acotada | El beneficiario práctico es el desarrollo de MODENA. El gate distingue documentación de sus mantenedores y resultados propios de sus autores de replicaciones independientes. No se inventó una declaración financiera para todo el corpus. No se proporcionó una declaración de financiación del informe; no se infiere ni patrocinio ni ausencia de él. |
| Investigación con sujetos humanos | **no aplica al objeto revisado** | Es una investigación de software/literatura técnica. No se activa revisión de participantes, consentimiento, acceso a datos identificables ni selección de vía institucional. Este dictamen de integridad no otorga autorización para actividades distintas. |

## Cobertura y modalidad del contraste bibliográfico

El universo más conservador es el corpus completo al que remite el borrador: **38 registros anotados y 37 unidades distintas**, después de tratar CL12 y V1 como la misma obra. Se contrastó el uso de **21/37 unidades (56.8%)** con contenido primario pertinente, sin multiplicar capítulos, copias o versiones. Dentro de las 14 referencias académicas enumeradas o identificadas expresamente al final del borrador —13 viñetas más DP-GEN—, se cubrieron **12/14**; Páll–Hess 2013 y Páll et al. 2020 no reciben aquí un nuevo contraste de contenido original.

Hay dos modalidades y no son revisión independiente entre personas:

- **P — lectura primaria previa reutilizada:** las 11 unidades ML01–ML11 fueron consultadas por este mismo agente en su invocación bibliográfica anterior de Fase 2. Ahora se contrastaron las afirmaciones del borrador/backlog con esas lecturas y sus registros; no se reabrieron íntegramente las 11 fuentes en esta invocación. Los enlaces originales, versiones, acceso y limitaciones están en [ML-SOURCES.md](../phase2_investigation/ML-SOURCES.md) y [ML-SEARCH-LOG.md](../phase2_investigation/ML-SEARCH-LOG.md).
- **R — reapertura puntual actual:** se abrieron 10 unidades adicionales y se examinaron los pasajes indicados abajo. R significa lectura del contenido necesario para la afirmación revisada, no reproducción de resultados ni auditoría de todas las ecuaciones. No se llevó a cabo una nueva búsqueda de DOI o retractaciones.

| Unidad | Modalidad | Qué se contrastó y resultado |
| --- | --- | --- |
| ML01 | P | Predictor/calculator, task y propiedades condicionadas; coherente con el contrato propuesto. Se conserva la divergencia documental `turbo`/Nitro. |
| ML02 | P | MIT del código frente a FAIR Chemistry License de materiales/pesos y acceso gated; el borrador no equipara las licencias. |
| ML03 | P | UMA se conserva como preprint con depósito 2025 y revisión 2026; entrenamiento no se extrapola a todos los checkpoints. |
| ML04 | P | Las rutas MACE ASE/LAMMPS/ML-IAP se usan como contratos específicos, sin ONNX universal. |
| ML05 | P | Tipos, factores de unidades y campos periódicos NequIP; no se atribuye igual semántica a todo peso compilado. |
| ML06 | P | Suavidad/top-K y conservación; el uso queda acotado a mecanismo y casos ensayados. El localizador A.1 ya corregido no plantea una incidencia de integridad. |
| ML07 | P | El benchmark distingue error estático de calidad dinámica; no se presenta como evaluación de UMA actual. |
| ML08 | P | El MTS del estudio es un algoritmo especificado; el borrador no lo convierte en permiso para switching arbitrario. |
| ML09 | P | DP-GEN separa exploración, etiquetas y entrenamiento; QE sigue siendo solver para las etiquetas electrónicas. |
| ML10 | P | La calibración explícita se usa con su dominio y particiones; el fallback mencionado no se presenta como garantía de muestreo. |
| ML11 | P | Sobreconfianza en estudios definidos; no se afirma que falle todo ensemble ni se asigna ranking a UMA. |
| CL01 | R | [Manual LAMMPS, vecinos](https://docs.lammps.org/stable/Developer_par_neigh.html), versión visible `30Sep2026 update1`: lista Verlet, migración y ordenamiento. Respalda el mecanismo; no una aceleración local garantizada. Solo se reabrió este capítulo de la unidad agrupada. |
| CL08 | R | [PDF de Essmann et al.](https://biomolmd.org/mw/images/e/e0/Spme.pdf), §IV/p. 8581: fuerza como derivada y distinción entre energía/momento. El borrador evita imponer momento conservado al SPME analítico. |
| CL09 | R | [Manuscrito I](https://arxiv.org/pdf/cond-mat/9807099), pp. 7–8: diferenciación, corrección de fuerza media y comparación a precisión equivalente. No se extrapola una etiqueta de método a superioridad universal. |
| CL10 | R | [Manuscrito II](https://arxiv.org/pdf/cond-mat/9807100), aproximación tras Eq. 22 y Fig. 7: restricciones de celda/ik y desviaciones con cargas inhomogéneas. El informe conserva el límite del estimador P3M. |
| CL12 = V1 | R | [Merz–Shirts, PLOS](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202764), secciones de validación de integrador y cinética: pruebas complementarias, no certificado universal. Identidad y declaración visible de intereses también comprobadas. |
| V2 | R | [Manuscrito de Shirts](https://arxiv.org/pdf/1208.0910), resumen y discusión de solapamiento/ensembles: contraste de distribuciones con muestras suficientes. El borrador lo distingue del drift energético. |
| V3 | R | [Chodera, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4945107/), metadatos y «The essential idea», Eq. 14: maximización de muestras efectivas de una serie. La cautela del borrador sobre estados lentos no vistos no contradice ese alcance. |
| V7 | R | [OpenMM 8.6 §19.12](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#montecarlobarostat): centros moleculares y aceptación por energía/volumen/número de moléculas. El informe distingue esta propuesta de completar presión atómica. |
| DFT-01 | R | [PWscf 7.5](https://www.quantum-espresso.org/Doc/user_guide_PDF/pw_user_guide.pdf), §§1.1 y 3.3: propiedades electrónicas, fuerzas y bandas/DOS. La frontera frente al MLIP queda representada correctamente. |
| DFT-02 | R | [INPUT_PW 7.5](https://www.quantum-espresso.org/Doc/INPUT_PW.html), entradas de estado, base y SCF; `conv_thr`, `ecutwfc`/`ecutrho`. El manifiesto propuesto se presenta como requisito de referencia, sin tolerancia universal atribuida a QE. |

Las otras 16 unidades se examinaron mediante sus fichas y el gate, sin contarlas como lectura primaria de contenido de esta revisión. Las observaciones del checkout y el resultado de pruebas locales se trataron como evidencia trazada del orquestador: no se repitieron esas pruebas ni se presenta este dictamen como una segunda auditoría del código.

## Puntos de integridad examinados expresamente

**Propuestas frente a resultados.** El resumen, el método, las puertas de salida y el backlog indican que no hay campañas nuevas de comparación LAMMPS, inferencia ML o DFT. La tabla del backlog afirma que sus cambios no están implementados ni aceptados. Tamaños S/M/L/XL no se convierten en días o probabilidades. Los objetivos electrostáticos y las pruebas de rendimiento se identifican como propuestos. La demostración aritmética del diagnóstico energético se distingue de una trayectoria. No se detectó promoción de estas propuestas a mejoras medidas.

**Atribución y versiones.** Los registros de Fu/eSEN y Bigi remiten a ICML/PMLR; Forces are not Enough se identifica como TMLR 2023; UMA conserva estado de preprint. Los manuscritos y versiones posteriores no se atribuyen silenciosamente al PDF editorial. Los dos Deserno–Holm de 1998 usan las letras a/b de manera coherente en el borrador. El cuerpo se apoya en anexos bibliográficos explícitos; no hay una cita fantasma demostrada en la muestra.

**Retracciones y DOI.** No se recibió una fila canónica v1.1 `bibliographic_integrity_signals[].retraction_status` ni un resultado del finalizador de citas. Por ello esta revisión **no asigna estado de retractación** a ninguna obra, no lo deriva de campos heredados ni afirma que no existan avisos. El gate declara que no hizo una búsqueda sistemática de retractaciones, y el borrador conserva esa limitación. Las aperturas fallidas de DOI no se presentan como inexistencia ni como resolución `DOI_VERIFIED`. Este es un límite declarado, no una acusación bibliográfica ni una autorización para ampliar lo que se comprobó.

**Lectura humana y originalidad.** No hay atestación de lectura humana en los materiales examinados. La declaración de IA es fiel al flujo observado y no inventa supervisión humana. El solapamiento entre autores de métodos y sus herramientas se señala en el gate; no se cuenta como corroboración independiente. No se puede calcular una tasa de autocitas personales del informe sin una autoría humana declarada: no se presume una tasa cero ni se confunde recurrencia de un desarrollador citado con autocita del autor del informe.

## Incidencias y acciones

- **Críticas:** ninguna localizada.
- **Condicionales obligatorias:** ninguna.
- **Advisory para la entrega:** conservar la declaración de IA, la ausencia de validación material y los límites de versiones/acceso al trasladar el borrador al informe final. Si se añade posteriormente una afirmación de ausencia de retractaciones, permiso de uso de pesos, lectura humana o mejora medida, necesitará su propia evidencia; este dictamen no la proporciona.

No hay decisiones del usuario sobre ítems condicionales/bloqueados que registrar. El veredicto `CLEARED` permite continuar la entrega en las dimensiones de integridad revisadas, conservando los límites descritos. Informe preparado con asistencia de IA; no es revisión humana independiente ni validación física de MODENA.
