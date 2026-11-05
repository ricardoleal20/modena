# Pregunta de investigación: mejoras al core de MODENA

Estado: Fase 1, propuesta para confirmar antes de la búsqueda bibliográfica. Documento de alcance, no resultado de una investigación concluida.

## Tema y decisión que debe apoyar

Identificar mejoras concretas al algoritmo y a las interfaces del core de MODENA que permitan sustituir flujos de trabajo atomísticos seleccionados actualmente resueltos con LAMMPS, Quantum ESPRESSO o FAIR-Chem. La unidad de comparación será un caso de uso reproducible con requisitos explícitos de exactitud, muestreo, rendimiento y operación, no la cantidad de funciones anunciadas por cada proyecto.

El punto de partida documentado es MODENA v0.1.1 como motor de dinámica molecular clásica en desarrollo: pares, términos moleculares, restricciones, Ewald/SPME y determinados ensembles. El alcance vigente distingue cobertura de regresión de validación de propiedades materiales; MPI/GPU son experimentales y el runner habitual usa CPU/Rayon. Estas observaciones proceden de [README](../../../../README.md) y [VALIDATED_SCOPE](../../../../docs/VALIDATED_SCOPE.md), no de una nueva ejecución experimental.

## Pregunta principal

¿Qué mejoras del core de MODENA ofrecen la mayor reducción verificable de las brechas que impiden sustituir flujos de simulación atomística seleccionados, según una auditoría del código actual contrastada con fuentes primarias y criterios reproducibles de aceptación?

## Evaluación FINER

Las puntuaciones son juicios de diseño del alcance, no mediciones científicas.

| Criterio | Puntuación | Justificación |
| --- | --- | --- |
| Feasible | 4/5 | Hay código, tests y documentación local para contrastar con fuentes públicas; medir todos los backends y sistemas de referencia requeriría recursos adicionales. |
| Interesting | 5/5 | La respuesta decide qué trabajo del core puede convertir implementaciones experimentales en flujos utilizables. |
| Novel | 3/5 | La aportación prevista es una síntesis específica para este repositorio; no se presupone un algoritmo físicamente novedoso. |
| Ethical | 5/5 | Investigación de software y literatura pública, sin participantes humanos; exige atribución y evitar afirmaciones de exactitud sin evidencia. |
| Relevant | 5/5 | Produce prioridades, ubicaciones de cambio y criterios para decidir si MODENA cubre un trabajo real. |
| **Promedio** | **4.4/5** | Supera el umbral de alcance de la habilidad: promedio ≥3 y ninguna dimensión <2. |

## Preguntas candidatas consideradas

| Candidata | F | I | N | E | R | Promedio | Decisión |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Priorizar mejoras del core por brecha verificable en casos de uso seleccionados | 4 | 5 | 3 | 5 | 5 | 4.4 | Seleccionada: combina utilidad inmediata y una comparación falsable. |
| ¿Qué optimizaciones reducen el tiempo por paso del kernel de pares actual a precisión fija? | 5 | 4 | 2 | 5 | 4 | 4.0 | Viable y concreta, pero demasiado estrecha para el objetivo de herramienta. |
| ¿Cómo puede MODENA reemplazar completamente las tres herramientas? | 1 | 4 | 2 | 5 | 4 | 3.2 | Descartada: “completamente” carece de fronteras operativas; factibilidad por debajo del mínimo. |

La segunda candidata permite experimentos locales claros, pero podría optimizar un kernel sin resolver las necesidades del usuario. La tercera presupone que productos con funciones distintas pueden reducirse a un único objetivo de paridad y no fija un experimento que permita declarar éxito.

## Límites del alcance

**Dentro del alcance**

- **Población/objeto:** revisión actual del repositorio MODENA, especialmente `crates/modena-core`, más sus interfaces CLI/Python/MCP cuando condicionen el uso del core.
- **Dominio:** dinámica molecular clásica; corrección de fuerzas, energía, virial, restricciones y ensembles; búsqueda de vecinos, electrostática, memoria y paralelismo; contratos para proveedores de fuerzas y potenciales ML.
- **Referencias de producto:** investigar por separado MD clásica en LAMMPS, estructura electrónica/DFT en Quantum ESPRESSO y modelos/interoperabilidad ML en FAIR-Chem. Esta separación guía la búsqueda y deberá sustentarse con documentación oficial en Fase 2.
- **Casos iniciales para evaluar:** pares tipo LJ; sistemas iónicos con electrostática; moléculas/agua con restricciones; un flujo de fuerzas ML conectado a un motor de integración. La selección final dependerá de evidencia y coste de reproducción.
- **Horizonte documental:** fuentes primarias fundacionales y recientes, disponibles al realizar la búsqueda, con fecha de consulta y versión cuando sean recuperables. El nombre de la carpeta es un identificador de sesión, no un corte bibliográfico.
- **Entorno:** evaluación de software sin restricción geográfica; mediciones limitadas al hardware y dependencias efectivamente disponibles.
- **Salida prevista:** recomendaciones priorizadas, referencia al código, fundamento, riesgos, dependencia entre cambios y un protocolo de validación por recomendación.

**Fuera del alcance de esta investigación**

- Implementar los cambios propuestos o prometer equivalencia de producción antes de medirla.
- Reescribir un solver DFT completo como requisito implícito del core de MD. Se evaluará la frontera entre integración, orquestación y sustitución de cálculos.
- Afirmar superioridad global a partir de drift NVE pequeño, tests aprobados, tiempos históricos o benchmarks que no igualen el problema físico.
- Meta-análisis estadístico de resultados heterogéneos; revisión sistemática PRISMA; validación exhaustiva de todas las químicas o todos los modelos ML.
- Entrenamiento de un modelo fundacional, adquisición de nuevos datasets masivos o campañas extensivas en hardware no disponible.

**Supuestos a comprobar**

1. El objetivo de “reemplazo” puede expresarse mediante casos de uso y criterios de salida; no se presupone paridad funcional universal.
2. Una arquitectura de proveedores de fuerzas podría permitir integrar modelos o motores externos. Su viabilidad, costes y límites se comprobarán en código y fuentes; no se da por implementada.
3. Las afirmaciones de los documentos históricos no se considerarán vigentes sin contrastarlas con la revisión actual y sus evidencias.
4. Los resultados de tests, mediciones nuevas, inferencias del auditor e hipótesis de mejora se etiquetarán por separado.

## Subpreguntas

1. ¿Qué brechas de corrección física y validación del core bloquean primero los casos de uso seleccionados?
2. ¿Qué cambios de algoritmos, estructuras de datos y ejecución ofrecen mejoras medibles a precisión y problema físico comparables?
3. ¿Qué contratos de proveedores de fuerzas e interoperabilidad permiten cubrir flujos de MD/ML y conectar referencias DFT con la menor ampliación justificada del core?

## Vinculaciones de subpreguntas (`sub_question_bindings`, Schema 1)

| Subpregunta | Ejes heredados de la pregunta principal | Desviaciones |
| --- | --- | --- |
| 1 | Objeto=MODENA actual y superficies dependientes; dominio=MD clásica y frontera MD/ML/DFT; casos=LJ, iones, moléculas/agua y flujo ML; horizonte=fuentes disponibles en la búsqueda; geografía=no aplica; entorno=recursos disponibles. | Ninguna. |
| 2 | Objeto=MODENA actual y superficies dependientes; dominio=MD clásica y frontera MD/ML/DFT; casos=LJ, iones, moléculas/agua y flujo ML; horizonte=fuentes disponibles en la búsqueda; geografía=no aplica; entorno=recursos disponibles. | Ninguna. |
| 3 | Objeto=MODENA actual y superficies dependientes; dominio=MD clásica y frontera MD/ML/DFT; casos=LJ, iones, moléculas/agua y flujo ML; horizonte=fuentes disponibles en la búsqueda; geografía=no aplica; entorno=recursos disponibles. | Ninguna. |

## Condición para avanzar

La habilidad `deep-research` exige confirmar el alcance después de Fase 1. Este documento no autoriza por sí mismo la búsqueda ni presupone aceptación del usuario. La continuación corresponde al orquestador, después del Methodology Blueprint y del Checkpoint 1.

Documento elaborado con asistencia de IA a partir de los archivos locales enlazados. No incluye referencias externas aún ni afirma haber ejecutado benchmarks.
