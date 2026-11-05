# Devil’s Advocate — Checkpoint 2

**Veredicto: PASS.** La síntesis respeta el alcance exploratorio y las condiciones de las fuentes verificadas. No se identificaron problemas Critical ni Major que obliguen a detener la composición. Las prioridades son propuestas de ingeniería; no son resultados de sustitución ni aceleraciones demostradas.

**Material revisado:** [SYNTHESIS.md](SYNTHESIS.md), [SOURCE-VERIFICATION.md](../phase2_investigation/SOURCE-VERIFICATION.md), [CASES.md](../phase2_investigation/CASES.md), [CODE-INVENTORY.md](../phase1_scope/CODE-INVENTORY.md), [RESEARCH-QUESTION.md](../phase1_scope/RESEARCH-QUESTION.md) y el Checkpoint 1 de alcance. Se aplicó la definición de `devils_advocate_agent`. No se realizó una búsqueda nueva ni una verificación independiente de todos los originales.

## Problemas concretos

### DA2-1 — “P1 compartida” puede inducir una generalización prematura del contrato

**Severidad: Minor. Tipo: inferencia arquitectónica / coste de intervención.**

**Ubicación:** interpretación “Un núcleo común útil necesita declarar capacidades…” y tabla “Prioridades y sensibilidad”, fila P1 compartida.

Las fuentes y el inventario sustentan que faltan identidad química y contratos explícitos de propiedades; no demuestran que un rediseño general de evaluación deba preceder a cualquier adaptador útil. La síntesis llama inferencia a esta elección y admite el proceso externo como alternativa, por lo que no incurre en una afirmación de rendimiento falsa. Sin embargo, el rótulo “compartida” puede leerse como una dependencia extensa e inevitable incluso para el caso LJ, que no necesita un proveedor ML.

**Ajuste accionable para el reporte:** definir como primer entregable una sección vertical mínima: un snapshot con identidad química/unidades explícitas, un proveedor externo persistente, resultados E/F y propiedades ausentes declaradas, y paridad frente a la interfaz oficial. Separar esa prueba del diseño general de backends o del portado de inferencia a Rust. Mantener P1 condicionada al flujo objetivo y permitir concluir que el adaptador externo basta. No se exige modificar la síntesis antes de redactar, si el reporte conserva esa distinción.

### DA2-2 — C4/C5 siguen siendo familias de evaluación, no un benchmark químico fijado

**Severidad: Minor. Tipo: selección / validez externa.**

**Ubicación:** C4 y C5 de CASES; “La vía ML/DFT más corta…” y “Prioridades y sensibilidad”.

Los casos congelan qué propiedades y metadatos hacen falta, pero el material, checkpoint y task de C4 y la referencia concreta de C5 aún no están elegidos. Esta apertura está declarada y no invalida el estudio. Sí limita la fuerza de una recomendación sobre qué modelo incorporar primero: elegir el que mejor resulte después de probar varios podría introducir el sesgo que la congelación de casos busca evitar.

**Ajuste accionable para el reporte:** distinguir dos decisiones: primero demostrar el contrato de interoperabilidad sin afirmar eficacia química; después fijar material, dominio, checkpoint/task, referencia DFT, snapshots de evaluación y criterios antes de comparar candidatos por calidad o coste. Registrar descartes y mantener al menos un resultado negativo. Con hardware y presupuesto desconocidos, no convertir el orden tentativo en un plan con retorno o plazos prometidos.

## Aspectos comprobados del razonamiento

| Riesgo solicitado | Evaluación del contenido revisado |
| --- | --- |
| Selección favorable y ausencia de controles | Hay cinco familias, una deliberadamente bloqueada, mecanismos de refutación y alternativas de menor intervención. Falta materializar C4/C5; se conserva en DA2-2. |
| Falsa equivalencia LAMMPS/QE/ML | La síntesis mantiene MD, evaluación ML y propiedades electrónicas como necesidades distintas. Fuerzas ML no se presentan como sustitución de bandas, orbitales o densidad. |
| Presupuesto o hardware desconocidos | Se reconoce expresamente que impiden ranking de retorno universal y plazos. GPU y MPI están condicionados al caso y perfil. |
| Aceleraciones sin medir | Ninguna cifra de speedup se transfiere a MODENA. Scratch, AoSoA, stencil y residencia son hipótesis; la evidencia local demuestra estructuras y transferencias, no porcentaje de runtime. |
| Independencia aparente de fuentes | La verificación deduplica CL12/V1, no cuenta accesos como réplicas y declara dependencia autoral/proyecto. Alta confianza en un contrato no equivale a replicación externa del rendimiento. |
| Presión restringida y MC | C3 y la síntesis distinguen virial para estimador/control dinámico de aceptación MC molecular por energía, volumen y medida correcta. No recomiendan eliminar sin más el rechazo del MC atómico actual. |
| Conservatividad, SPME y ensemble | Se separan F=-∇U, momento, discretización temporal y muestreo. No se propone corregir automáticamente fuerza neta ni validar un ensemble por drift pequeño. |
| Switching ML/DFT | Un esquema MTS publicado no se usa para justificar fallback arbitrario basado en incertidumbre. El dominio y la calibración siguen pendientes. |
| Repetir defectos históricos | La síntesis reconoce las correcciones de unidades, muestreo, tablas y restart; no las presenta como bugs actuales. El cero energético se plantea como límite de diagnóstico, no como prueba de integrador incorrecto. |

## Contrargumento más fuerte

> La investigación identifica mecanismos razonables, pero todavía no demuestra que reescribir el core acerque más a un reemplazo útil que conectar MODENA a los motores y modelos existentes. Sin material objetivo, presupuesto ni campaña física emparejada, la prioridad más defendible puede ser un flujo estrecho de interoperabilidad y validación, incluso si produce pocos cambios algorítmicos.

La síntesis soporta esta lectura: ofrece alternativas externas, condiciona las prioridades al caso y reconoce mediciones ausentes. El reporte debe preservar esa posibilidad como conclusión legítima, sin presentar integración como una etapa inferior obligatoria antes de reimplementar.

## Pruebas de resistencia

| Prueba | Resultado |
| --- | --- |
| Retirar el paper de arquitectura más favorable: ¿se mantiene el diagnóstico local? | Sí para la estructura observada y los vacíos de interfaz; disminuye el respaldo de una estrategia concreta, y el speedup seguiría sin medir. |
| ¿La alternativa “no reescribir GPU/MPI ahora” es creíble? | Sí; la tabla la admite y exige caso/hardware/perfil. |
| ¿La alternativa “conservar un motor externo para DFT o NPT” es aceptada? | Sí; no contradice los límites de sustitución propuestos. |
| ¿Un modelo conservativo elimina la necesidad de estudiar observables? | No; la síntesis lo rechaza y evita atribuir a benchmarks previos una evaluación de pesos UMA no ensayados. |
| ¿Las conclusiones se trasladan a toda química o todo hardware? | No; están limitadas por casos y evidencia no ejecutada. |
| ¿El cambio del cero de energía invalida la trayectoria? | No; afecta al diagnóstico relativo. La síntesis conserva esta separación. |

## Límites del checkpoint

Este PASS evalúa coherencia entre síntesis, fuentes documentadas y alcance. No certifica DOI, reproducibilidad de artículos, calidad de pesos ni corrección numérica de MODENA. No hubo nueva búsqueda contraria en esta invocación; se revisaron las alternativas y exclusiones registradas por los bibliógrafos y la verificación de fuentes. Los Checks 1 y 2 están satisfechos; falta el Checkpoint 3 sobre el reporte completo.

Elaborado con asistencia de IA. No se presenta como revisión humana independiente.
