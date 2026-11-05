# Devil’s Advocate — Checkpoint 1

**Veredicto: PASS.** Pregunta y método son adecuados para una investigación exploratoria que priorice cambios verificables del core. No prueban que MODENA ya pueda sustituir ningún flujo concreto. No se identificaron problemas Critical ni Major que impidan presentar el alcance para confirmación.

**Material revisado:** [pregunta](RESEARCH-QUESTION.md), [metodología](METHODOLOGY.md) e [inventario local](CODE-INVENTORY.md). Revisión del repositorio declarada por esos documentos: `f166399070b2d6ebfd91f51e5934c132f04d46c0`. Esta revisión no ejecutó tests, benchmarks ni búsquedas externas.

## Fortalezas del diseño

- La pregunta evita asumir paridad total: define la comparación por flujos científicos y separa las necesidades de MD clásica, evaluación de modelos ML y referencia DFT. La metodología permite concluir que conviene integrar una herramienta externa en vez de recrear su función.
- El inventario distingue avances ya presentes de hipótesis nuevas: skin y medio par, especialización LJ, planes y scratch FFT, y exclusiones indexadas ya constan como implementados. No se reconstruye una lista de pendientes a partir de documentos antiguos.
- El método exige igualar física y precisión, mantiene procedencia de resultados históricos y diferencia corrección, muestreo y rendimiento. Las puertas G0–G6 permiten refutar una propuesta sin redefinir el éxito después.

## Problemas concretos

### DA1 — La selección de casos aún permite favorecer al motor

**Severidad: Minor. Tipo: selección de casos / moving goalposts.**

**Ubicación:** “Casos iniciales para evaluar” de la pregunta y “Definir casos de sustitución” de la metodología.

Las familias LJ, iones, agua y ML son razonables para explorar, pero no fijan todavía sistema, observables ni necesidad de usuario. Si los casos se eligen después de ver qué soporta MODENA, una tabla favorable podría demostrar cobertura de ejemplos cómodos y presentarse erróneamente como avance hacia sustitución práctica. No invalida el diseño porque las conclusiones previstas son exploratorias y el método ya exige tolerancias anteriores a las mediciones.

**Ajuste accionable:** al abrir Fase 2, congelar una ficha por caso con entrada, Hamiltoniano/modelo, observable, condición operativa mínima, criterio de aceptación y motivo de selección. Conservar al menos un caso actualmente bloqueado. Registrar cualquier cambio de selección y declarar las prioridades provisionales mientras falte un caso de uso confirmado por el usuario.

### DA2 — La priorización necesita una alternativa de coste mínimo

**Severidad: Minor. Tipo: sesgo de intervención.**

**Ubicación:** “Priorización propuesta” y “Impacto en uso” de la metodología.

La tabla final contempla alternativas cuando una optimización fracasa, pero no exige comparar cada propuesta con mantener el core actual y cubrir la necesidad mediante integración o mejora del flujo. Una investigación enfocada en algoritmos puede sobrevalorar cambios profundos aunque la barrera de adopción sea entrada química, parametrización, interoperabilidad o evidencia de validación.

**Ajuste accionable:** añadir a cada recomendación la alternativa “sin modificar el algoritmo” —cuando tenga sentido— y comparar qué necesidad satisface, coste de validación y dependencia externa. Si el cambio de core no supera esa alternativa para un caso definido, reducir su prioridad. No convertir una preferencia por Rust, GPU o potenciales ML en evidencia de utilidad.

## Contrargumento más fuerte

> Aunque MODENA obtuviera kernels más rápidos y fuerzas correctas, podría seguir sin sustituir el trabajo científico del usuario: cobertura química, parámetros fiables, observables, herramientas de preparación y análisis e integración con referencias pueden determinar la utilidad. Mejorar el core es una hipótesis sobre la vía de adopción, no una demostración de que esa vía sea la de mayor impacto.

El diseño resiste este contrargumento porque incluye las interfaces cuando condicionan al core, exige casos de uso y no promete equivalencia global. Debe conservarse la posibilidad de recomendar integración o validación antes que nuevos algoritmos.

## Pruebas de resistencia del alcance

| Prueba | Resultado |
| --- | --- |
| ¿La pregunta sigue siendo contestable si no se ejecuta una nueva campaña extensa? | Sí: puede producir prioridades justificadas y protocolos. Las ganancias seguirían siendo hipótesis, no resultados medidos. |
| ¿Es creíble la respuesta opuesta: no reescribir un subsistema? | Sí: integración, reutilización y mejora de validación pueden ser la recomendación principal. |
| ¿Permite concluir que MODENA no sustituye una necesidad DFT o ML concreta? | Sí. No obliga a recomendar un solver DFT ni entrenamiento de un modelo fundacional. |
| ¿Los resultados serían generalizables a toda MD o toda química? | No. El alcance y la metodología lo excluyen expresamente. |
| ¿Se reabren como actuales defectos históricos ya corregidos? | No en los documentos revisados: el inventario exige reproducción antes de tratarlos como defectos vigentes. |
| ¿Se justifica la relevancia práctica? | Sí, si la síntesis enlaza cada mejora con un flujo bloqueado o un coste real; DA1 exige hacer explícito ese vínculo. |

## Límites y siguiente paso

No hay corpus bibliográfico todavía; no es posible evaluar su cobertura o sesgo en este checkpoint. La confirmación del alcance no constituye aceptación de una conclusión ni autorización para afirmar rendimiento, exactitud material o sustitución de producción. Tras confirmar, Fase 2 debe aplicar DA1 y DA2 como guardas de análisis, verificar fuentes y preparar los Checkpoints 2 y 3.

Revisión elaborada con asistencia de IA. No se presentan fuentes externas, mediciones nuevas ni afirmaciones de autoridad científica independiente.
