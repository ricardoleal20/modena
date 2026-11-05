# Devil’s Advocate — Checkpoint 3

**Veredicto: PASS.** El [borrador completo](../phase4_composition/DRAFT.md) incorpora los ajustes de Checkpoint 2 y mantiene sus conclusiones dentro de la evidencia. No hay problemas Critical o Major que bloqueen la revisión final.

**Base de revisión:** borrador, [backlog](../phase3_analysis/ACCEPTANCE-BACKLOG.md), [Checkpoint 2](../phase3_analysis/CHECKPOINT-2.md) y [verificación de fuentes](../phase2_investigation/SOURCE-VERIFICATION.md). No se efectuó búsqueda nueva, ejecución de código ni segunda verificación íntegra de originales.

## Hallazgo concreto

**DA3-1 — La posición de “nueva química” en el orden puede confundirse con una decisión tardía del dominio. Severidad: Minor.**

La tabla “Orden de trabajo y puertas de salida” coloca “Nueva química y adquisición ML validada” en el punto 5, después de GPU/MPI. El contexto explica que el orden depende del caso, pero un lector podría interpretar que debe optimizar primero y seleccionar después el trabajo científico. Sin material objetivo no se sabe si throughput LJ, presión de agua o un modelo ML resuelve su necesidad. La tabla no invalida la conclusión porque el texto reconoce la falta de dominio/presupuesto y ofrece caminos alternativos.

**Ajuste accionable de redacción:** añadir antes de la tabla una decisión inicial breve: elegir el flujo/material, propiedades requeridas y recursos disponibles; esa selección determina qué filas aplican. Distinguir la elección temprana del dominio de la implementación posterior de nueva física o adquisición. No requiere una nueva pregunta al usuario ni una nueva investigación para cerrar este informe: puede quedar registrada como condición de la futura campaña.

## Verificación de vulnerabilidades principales

| Punto | Evaluación |
| --- | --- |
| Vertical mínimo antes de arquitectura general | Incorporado: snapshot, identidad/unidades, proveedor oficial persistente, E/F y ausencias. Puede ser solución suficiente; no obliga a portar pesos o generalizar todos los backends. |
| C4/C5 pendientes | Declarado: material/modelo/referencia no están fijados y paridad de software no equivale a eficacia química. |
| Speedups y plazos | No se prometen. Los mecanismos de CPU/GPU/MPI se presentan como candidatos y se exigen perfiles del paso completo. |
| DFT frente a ML | Se preserva el solver para propiedades electrónicas y etiquetas convergidas. Un proveedor de fuerzas no se presenta como reemplazo general de QE. |
| Presión y MC molecular | Se distinguen el estimador con impulsos de restricciones y la aceptación MC por energía/volumen/medida. No se recomienda habilitar el MC atómico actual quitando una validación. |
| SPME: energía y momento | Se conserva la diferencia entre funcional analítico, error de malla y momento total. No hay corrección automática de fuerza media ni trasplante de una cota P3M/ik. |
| Conservatividad y ensemble | Drift, gradientes, muestreo y observables tienen controles distintos. La termostatización no se usa como garantía para fuerzas no conservativas. |
| Correcciones históricas | No se presentan como defectos actuales: scratch FFT, skin, medio par y endurecimiento previo se reconocen como existentes. |
| Independencia y calidad de fuentes | El texto deduplica obras y declara revisores IA, fuentes mutables, acceso incompleto y ausencia de réplica de papers. |
| Originalidad y relevancia | La contribución es una aplicación técnica al checkout, con puntos de cambio y aceptación; no se anuncia un nuevo algoritmo científico. El backlog permite refutar propuestas y conservar alternativas de menor intervención. |

## Contrargumento más fuerte

> El informe aún no demuestra qué trabajo científico real puede migrar el usuario a MODENA. Las recomendaciones pueden producir un motor mejor diseñado sin alcanzar una sustitución útil si el coste dominante es parametrización, cobertura química o validación material. Usar directamente un calculator y motor existentes podría satisfacer antes la necesidad.

El borrador admite esa alternativa y limita su dictamen a una ruta de investigación e implementación revisable. El valor demostrado aquí es reducir incertidumbre sobre qué probar y dónde cambiar el código; el valor de adopción y la reducción de coste científico siguen por medir. Esa distinción debe conservarse al presentar el resultado.

## Prueba final de inferencia

- Sin las fuentes de arquitectura GPU, permanecen los límites locales y la recomendación de medir antes de reescribir; no cae una conclusión de aceleración porque no se formula.
- La conclusión opuesta —posponer cambios de core y usar integración externa— es compatible con el informe y con varios criterios de salida.
- No se infiere validez para toda química, todo hardware ni una campaña no ejecutada.
- La declaración de IA y las limitaciones de validación son explícitas; no se detectó una preocupación Critical de integridad que requiera confirmación adicional.

**Cierre del checkpoint:** satisfechos los tres checkpoints de Devil’s Advocate. Incorporar DA3-1 durante la revisión final, junto con las observaciones editoriales y de integridad. Este PASS revisa la argumentación, no certifica corrección numérica o aptitud de producción. Elaborado con asistencia de IA.
