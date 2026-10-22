# Revisión editorial e integridad de la auditoría MODENA–LAMMPS

Fecha: 3 de octubre de 2026. Documento revisado: [REPORT.md](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/REPORT.md), sobre el código `b68abd112339aa6970261f470a572da6633bc617`. Se aplicaron los criterios de Phase 5 de deep-research, calibrados para una auditoría técnica privada, no para un artículo sometido a una revista. Un mismo revisor aplicó las dos rúbricas; esto no representa dos revisiones humanas independientes. No se modificó el informe ni código de producción.

## Revisión editorial

**Dictamen: Accept. Puntuación ponderada: 4.1/5.** No se localizaron problemas críticos ni mayores que invaliden las conclusiones. Las precisiones menores identificadas durante la revisión fueron incorporadas antes de este dictamen.

| Dimensión | Peso | Nota | Justificación |
|---|---:|---:|---|
| Contribución y utilidad | 20% | 4.5 | Une inventario del código, fallos físicos reproducibles y nueva medición; convierte la pregunta amplia en decisiones verificables por dominio. |
| Rigor metodológico | 25% | 4.0 | Los probes tienen controles y predicción analítica; el benchmark conserva entradas y salidas. Reconoce series correlacionadas, carga compartida y sólo tres repeticiones. |
| Suficiencia de evidencia | 25% | 4.0 | Los tres defectos centrales tienen código y contraejemplos. Se separan implementaciones de integración operativa y mediciones de inferencias. La validación completa de ensembles y el rendimiento en cluster siguen abiertos, como declara el texto. |
| Coherencia argumental | 15% | 4.0 | El dictamen se deriva de hallazgos específicos; las prioridades exigen criterios de salida. No exige replicar todo LAMMPS para sostener investigación en un nicho. |
| Claridad y precisión | 15% | 4.0 | Las tablas, enlaces y límites facilitan revisar la evidencia. Conserva algunos términos de ingeniería en inglés, apropiados para el destinatario técnico. No se impone APA a un informe de auditoría. |

Fortalezas concretas:

1. La revisión de Bussi reconcilia la errata del autor con A7 y una implementación independiente; evita convertir un párrafo incorrecto del paper en autoridad final.
2. Los probes distinguen un contraejemplo de corrección de una campaña estadística de validación: no presentan estados correlacionados como muestras independientes ni extrapolan la magnitud del barostato a cualquier NPT.
3. La comparación LJ distingue un hilo frente a serial y varios hilos frente a serial. La mediana y los cocientes publicados se contrastaron con los tiempos individuales del JSON y coinciden. No se afirma una victoria con igualdad de recursos.
4. El informe trata agua y escalabilidad con límites explícitos: concordancia de temperatura no implica causalidad única; la deriva histórica no certifica un líquido; el prototipo MPI no establece rendimiento multinodo.

Revisiones menores resueltas:

| Pasaje | Observación | Resultado de la revisión |
|---|---|---|
| SOURCES, Q01 | La descripción inicial restringía xTB a cajas ortogonales. La documentación admite cajas ortogonales y triclínicas fijas. | Corregido; el informe no había basado su dictamen en esa restricción errónea. |
| SOURCES, método | Las URLs release/stable son mutables y pueden mostrar distintas generaciones; no debía asignarse una versión uniforme al corpus. | Corregido; se distingue además el ejecutable local 22Jul2025 Update 6. |
| REPORT, metodología LJ | Los 200 pasos descartados ocurren en un proceso separado. No equilibran las trayectorias posteriormente medidas. | Corregido y explícito en el informe. |

Sugerencia opcional: en la frase MPI de [REPORT.md:131](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/REPORT.md:131), reemplazar «evita migración» por «mantiene ownership inicial fijo y no prueba transporte de migración». El test declara la intención de evitar migración; no contiene una demostración general de ausencia de cruces. La distinción no altera la conclusión sobre integración distribuida pendiente.

El informe es entregable con estos límites. Las campañas futuras propuestas no son requisitos para entregar esta auditoría: son requisitos para las afirmaciones más fuertes de robustez y competitividad que hoy no se realizan.

## Revisión de integridad de investigación asistida por IA

**Dictamen: CLEARED.** Este dictamen cubre atribución, transparencia y representación de la evidencia; no constituye certificación científica o institucional.

| Dimensión | Estado | Evaluación |
|---|---|---|
| Divulgación de IA | Pass | REPORT declara auditoría asistida por IA y revisión de agentes; SOURCES precisa búsqueda, extracción y síntesis. No se presenta como revisión humana independiente. |
| Integridad de atribución | Pass | Fuentes primarias identificables; paráfrasis pertinentes; correcciones menores de corpus resueltas. No se encontraron referencias fabricadas ni plagio en las fuentes comprobadas. |
| Uso dual | Pass | Riesgo bajo: auditoría de métodos y software MD. No contiene instrucciones operativas que habiliten daño concreto. |
| Representación justa | Pass | Reconoce fortalezas de MODENA, paquetes/restricciones de LAMMPS y comparaciones incompletas. No atribuye superioridad universal a ningún motor. |
| Ética de datos | Pass | Código del workspace autorizado, documentación pública y datos sintéticos de simulación. No se identificó uso de datos personales. |
| Conflicto de interés | Pass | La finalidad es evaluar la plataforma del solicitante; el trabajo no se presenta como evaluación institucional externa ni como ensayo comercial imparcial certificado. |
| Sujetos humanos | No aplica al contenido auditado | No se analizaron participantes ni datos humanos. No se emite determinación de una vía institucional. |

Verificación de referencias: se contrastaron 18 de las 24 entradas del corpus (75%): L01–L13, L16, V01, V02, V04 y Q01. Se comprobó identidad y contenido pertinente en fuentes primarias; las páginas combinadas cuentan como una entrada del corpus. La matemática CSVR se contrastó además con la reconciliación de SOURCES y los resultados reproducibles. No se recibió un resultado canónico de retraction-status del citation gate; esta revisión no afirma certificar ese estado para cada artículo.

Las dos discrepancias menores de fuentes se resolvieron editorialmente, sin decisiones u overrides del usuario. No hay condiciones pendientes de integridad ni motivo para abrir un flujo de aprobación adicional. La revisión no valida todas las combinaciones de modelos, las propiedades experimentales de sustancias ni predicciones de aceleración futura.
