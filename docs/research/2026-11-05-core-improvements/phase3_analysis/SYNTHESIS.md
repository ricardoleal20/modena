# Síntesis: mejoras del core y sustitución por casos de uso

Fase 3. Base local `f166399070b2d6ebfd91f51e5934c132f04d46c0`. Alcance confirmado. Entradas: [pregunta](../phase1_scope/RESEARCH-QUESTION.md), [casos congelados](../phase2_investigation/CASES.md), los cuatro corpus anotados y [verificación de fuentes](../phase2_investigation/SOURCE-VERIFICATION.md). Esta síntesis formula inferencias y decisiones de ingeniería. No aporta nuevas trayectorias, propiedades materiales ni aceleraciones medidas.

## Matriz de evidencia

Los IDs remiten a [MD clásica](../phase2_investigation/CLASSICAL-SOURCES.md), [ML](../phase2_investigation/ML-SOURCES.md), [DFT/ASE](../phase2_investigation/DFT-SOURCES.md) y [validación](../phase2_investigation/VALIDATION-SOURCES.md). Alta = soporte directo para un mecanismo o contrato. La utilidad o aceleración concreta en MODENA sigue sin medir.

| Tema | Convergencia entre fuentes | Evidencia local | Confianza y vacío |
| --- | --- | --- | --- |
| Exactitud, muestreo y observables son puertas diferentes | CL12/V1, V2, ML06–08 | Suite registrada de 284 aprobadas/1 ignorada; alcance vigente limita propiedades de líquido | Alta para separar pruebas; falta campaña material emparejada |
| Vecinos y memoria deben diseñarse para el backend | CL01–04, CL07, ML04–06 | Medio par/skin existentes, buffers xyz intercalados, scratch 3N por partición, rebuild triclínico all-pairs | Alta sobre estructuras, baja sobre speedup local |
| Energía, fuerzas y esfuerzo requieren semántica explícita | CL08–11, V5–7, DFT-04/07, ML01/03/08 | Virial escalar, presión restringida diagnóstica, proveedor actual de acumulación | Alta sobre contratos, media sobre arquitectura propuesta |
| La residencia GPU afecta al paso completo | CL04–07 | Upload de posiciones/readback por cómputo e integración host; selector TOML ausente | Alta sobre barrera arquitectónica, baja sobre retorno de reescritura |
| ML amplía cobertura dentro de un dominio y referencia definidos | ML01–11, DFT-01/02/05 | Sin modelo ML en ModelSpec, elementos cosméticos y tipos FF numéricos | Alta sobre brecha de interfaz; eficacia química no evaluada |
| QE mantiene una función electrónica y de etiquetado | DFT-01–03/06–08, ML01/09 | No solver electrónico ni adaptador externo integrado | Alta sobre diferencia funcional; coste/beneficio del conector por medir |

Las observaciones locales y sus archivos se encuentran en [CODE-INVENTORY](../phase1_scope/CODE-INVENTORY.md), [alcance vigente](../../../VALIDATED_SCOPE.md) y [verificación](../phase1_scope/VERIFICATION.md). Las fuentes secundarias y las proyecciones de planes antiguos no se usan como pruebas de aceleración actual.

## Interpretación integrada

**Un núcleo común útil necesita declarar capacidades físicas y computacionales.** La interfaz ASE distingue propiedades soportadas, estado físico y resultados ausentes (DFT-04/05); los modelos ML necesitan además tipos químicos, metadatos electrónicos y grafos periódicos (ML01/04/05). Al unir esa evidencia con `ForceModel`, `ModelSpec` y `System` locales, la inferencia es ampliar un contrato general de evaluación antes de añadir una colección de backends específicos. La composición debe declarar si un término es potencial completo o corrección. Sumar LJ/Coulomb a un modelo total sin definir un residual puede contar energía dos veces. Esta última advertencia es una inferencia matemática de composición, no un fallo reproducido de MODENA.

**La corrección existente debe preservarse y ampliarse.** Los defectos de muestreo, unidades, tablas y reinicio de la auditoría de octubre tienen correcciones documentadas en el checkout. El siguiente trabajo no es anunciarlos otra vez, sino extender controles a los casos nuevos. V1/V2 y ML06–08 justifican pasar de igualdad estática a convergencia temporal, distribución y observable. La observación [del cero energético](../phase2_investigation/DIAGNOSTIC-OBSERVATION.md) añade una mejora concreta: el indicador relativo actual cambia al sumar una constante a U. Conviene informar ΔE, pendiente/error acotado y una normalización física declarada independiente de ese cero. El ensemble debe determinar qué invariante o estadística tiene sentido. Propuesta, no implementación.

**CPU: minimizar tráfico y asignación antes de elegir instrucciones.** `pair_pass` ya especializa LJ y evita paralelismo en listas pequeñas. El siguiente candidato local es reutilizar scratch y medir la reducción de buffers de tamaño 3N por partición. Después se puede evaluar orden espacial, bloques AoSoA y cluster-pairs (CL01–04). El beneficio SIMD puede exigir hacer más pares aritméticos para reducir cargas/ramas (CL03): contar solo pares por segundo no resuelve la comparación. Un reordenamiento exige remapear topología, excepciones, IDs y restart. Para celdas triclínicas, reemplazar all-pairs requiere un stencil geométricamente demostrado, no trasladar sin prueba el ortorrómbico. Todo es hipótesis local de rendimiento pendiente de perfil.

**SPME: optimizar una curva de error/coste.** CL08–10 muestran que asignación, influencia y diferenciación definen el método. La malla analítica local preserva F=-∇U discretizado y puede no conservar momento. No se recomienda restar fuerza media automáticamente ni copiar un estimador P3M/ik como si certificara SPME. El candidato es exponer un objetivo de error de fuerza definido, calibrarlo con snapshots y Ewald convergido, y medir real-space, spread, FFT y gather. Los planes/scratch FFT ya son reutilizables. Cambiar el operador de diferenciación constituye una decisión de modelo que requiere nueva validación, aunque conserve el nombre PME.

**Agua restringida tiene dos caminos distintos.** CL11/V5/V6 apoyan completar el esfuerzo de restricciones para reportar presión atómica y futuro control dinámico. V7 aporta una alternativa más acotada para NPT: MC molecular que preserva geometría interna y usa la medida/Jacobiano correctos. El segundo no necesita virial instantáneo para aceptar cambios, pero tampoco arregla el primer estimador. La enmienda explícita de C3 preserva esta distinción. No basta con quitar el rechazo del validador actual.

**GPU y MPI exigen estado persistente y pruebas de trayectoria.** CL04/05/07 apoyan mover más operaciones al dispositivo, manteniendo rutas híbridas cuando convengan. El driver MPI local manda el conjunto completo propio a vecinos y reconstruye vista/lista por evaluación. La inferencia arquitectónica es introducir estado local persistente, halo geométrico y migración real. La prueba multirank existente conserva sistema global y no demuestra migración distribuida. Nada de ello prueba que una reescritura supere LAMMPS o merezca prioridad antes de un caso químico útil.

**La vía ML/DFT más corta es interoperar y comprobar, después especializar.** DFT-04–08 y ML01/04/05 permiten diseñar un evaluador MODENA para ASE y un proveedor externo persistente hacia calculators oficiales. Son direcciones distintas. Una conexión i-PI necesita un handshake de identidad fuera del transporte de coordenadas; el protocolo no transmite números atómicos (DFT-08). QE puede producir etiquetas convergidas para entrenar/validar. Sustituir algunas fuerzas DFT por ML dentro de un dominio no sustituye bandas, orbitales o densidad electrónica. Para adquisición, ML09–11 apoyan un ciclo de exploración/etiquetado/reentrenamiento con calibración, y ML08 impide justificar un switching arbitrario a partir de un esquema MTS particular.

## Tensiones y resoluciones propuestas

T1. **Conservación de energía frente a momento:** CL08/09 no se contradicen al describir operadores distintos. El contrato debe identificar funcional y diferenciación. No corregir silenciosamente un invariante destruyendo otro.

T2. **Residencia GPU frente a trabajo híbrido:** CL04/05 muestran arquitecturas y regímenes diferentes. Minimizar la dependencia host por paso es un mecanismo razonable; suprimir todo trabajo CPU no es un objetivo universal.

T3. **Calibración de ensembles frente a sobreconfianza:** ML10 calibra predicciones en particiones determinadas; ML11 documenta sobreconfianza fuera de distribución en otros casos. Es compatible que un ensemble sea útil tras validación y que no ofrezca una garantía general. Medir cobertura por dominio.

T4. **Modelo conservativo frente a calidad dinámica:** ML03 describe un modelo final conservativo, mientras ML07 exige observables y estabilidad. No hay contradicción: ese benchmark anterior no evalúa los checkpoints UMA inspeccionados. Una propiedad de arquitectura no implica validación del caso C4.

Las divergencias documentales FAIR-Chem `main` Nitro/guía `turbo` y ASE k-points permanecen sin resolver para una instalación concreta. La acción es fijar y comprobar release/commit/pesos/opciones, no inventar compatibilidad. La presión atómica incompleta frente a MC molecular por energía es una diferencia entre estimador y algoritmo, no un desacuerdo entre motores. [Gate de fuentes](../phase2_investigation/SOURCE-VERIFICATION.md).

```yaml
cross_paper_tensions:
  - pair_id: CP-001
    paper_a: CL08
    paper_b: CL09
    candidate_basis: shared construct/outcome/measure
    overlap_topic: energia y momento en electrostatica de malla
    a_finding: fuerzas analiticas derivadas de energia sin momento total necesariamente conservado
    a_evidence_pointer: CLASSICAL-SOURCES CL08, seccion IV
    b_finding: asignacion y diferenciacion cambian las propiedades de conservacion
    b_evidence_pointer: CLASSICAL-SOURCES CL09, paginas 7-8 del manuscrito
    pair_assessment: conditional_difference
    resolution_status: resolved_in_synthesis
    resolution_pointer: SYNTHESIS > Tensiones y resoluciones propuestas > T1
    scholar_confirmation: pending
  - pair_id: CP-002
    paper_a: CL04
    paper_b: CL05
    candidate_basis: shared RQ subtopic
    overlap_topic: residencia GPU y reparto CPU
    a_finding: paso residente con rutas PME hibridas segun regimen
    a_evidence_pointer: CLASSICAL-SOURCES CL04, secciones IV-V
    b_finding: cobertura GPU con excepciones y dependencias del metodo
    b_evidence_pointer: CLASSICAL-SOURCES CL05, secciones 4-5
    pair_assessment: conditional_difference
    resolution_status: resolved_in_synthesis
    resolution_pointer: SYNTHESIS > Tensiones y resoluciones propuestas > T2
    scholar_confirmation: pending
  - pair_id: CP-003
    paper_a: ML10
    paper_b: ML11
    candidate_basis: shared construct/outcome/measure
    overlap_topic: fiabilidad de incertidumbre de ensembles
    a_finding: calibracion evaluada en datasets y particiones concretas
    a_evidence_pointer: ML-SOURCES ML10, secciones 2.5-2.6 y 3
    b_finding: sobreconfianza y error subestimado en cuatro casos
    b_evidence_pointer: ML-SOURCES ML11, casos y conclusiones
    pair_assessment: conditional_difference
    resolution_status: resolved_in_synthesis
    resolution_pointer: SYNTHESIS > Tensiones y resoluciones propuestas > T3
    scholar_confirmation: pending
  - pair_id: CP-004
    paper_a: ML03
    paper_b: ML07
    candidate_basis: shared RQ subtopic
    overlap_topic: potenciales ML y utilidad dinamica
    a_finding: ajuste final conservativo y benchmarks de dominios definidos
    a_evidence_pointer: ML-SOURCES ML03, seccion 2.4 y limites de evaluacion
    b_finding: MAE estatico y calidad dinamica pueden desacoplarse
    b_evidence_pointer: ML-SOURCES ML07, secciones 4-7
    pair_assessment: no_material_conflict
    resolution_status: not_applicable
    scholar_confirmation: pending
```

Cobertura del inventario: 17 obras de investigación distintas (incluido el preprint UMA), cuatro pares candidatos de temas compartidos. Es un escaneo limitado, no una comparación exhaustiva de todos los pares ni una confirmación humana. Las interpretaciones propuestas quedan disponibles para revisión del usuario; no crean una nueva aprobación obligatoria.

## Prioridades y sensibilidad

| Nivel propuesto | Entrega | Dependencia y alternativa de menor intervención |
| --- | --- | --- |
| P0 | Extender evidencia emparejada y diagnósticos del ensemble/cero energético | Reutilizar harness y metadata actuales; analizar CSV con herramientas externas antes de añadir análisis al core |
| P1 compartida | Contrato de evaluación con identidad química, unidades, capacidades y resultados opcionales; adaptadores ASE | Usar calculators oficiales con proceso persistente antes de portar pesos/inferencia a Rust |
| P1 si agua/NPT | Presión restringida completa y/o MC molecular validado | Mantener NVT fijo o ejecutar NPT de referencia fuera de MODENA hasta pasar puertas |
| P2 por perfil | Scratch CPU, localidad/bloques/vecinos; control de error SPME | Retener camino serial/general como referencia; no reescribir sin cuello y error objetivo |
| P3 por caso y hardware | Paso GPU residente o MPI de producción | Elegir un backend según necesidad medida, conservar ejecución CPU o motor externo de referencia |
| Posterior al dominio | Física many-body/FF químico o adquisición ML | Elegir por material; integrar modelo validado puede costar menos que reimplementar familias completas |

La priorización cambia con el objetivo: agua eleva presión/MC; materiales ML eleva proveedor/ASE; LJ throughput eleva perfil CPU; producción multinodo eleva MPI tras cerrar trayectoria/identidad. No hay presupuesto/hardware/caso material del usuario que permita estimar plazos o un ranking de retorno universal.

## Vacíos y límites

Faltan snapshots comunes entre motores sobre unidades electrostáticas corregidas, propiedades de líquido con incertidumbre, benchmark de paso completo con recursos/errores iguales, evidencia de migración MPI real, modelo ML fijado y ejecutado, y campaña DFT convergida. No se mide retorno económico o adopción del producto. Las 284 pruebas registradas no cubren esos vacíos por sí mismas. Las conclusiones son recomendaciones de ingeniería respaldadas por mecanismos y código, sujetas a refutación por las puertas G0–G6 del método.

El [Checkpoint 2](CHECKPOINT-2.md) pasó sin Critical/Major. Se incorporan dos precisiones al reporte: iniciar el contrato con una sección vertical mínima de proveedor externo y paridad, generalizando solo capacidades demostradas; C4/C5 son familias de evaluación y requieren material, modelo y referencia fijados antes de afirmar eficacia química. No se solicita una nueva aprobación para este refinamiento metodológico.

Se usó IA en búsqueda, comprobación, síntesis y redacción. No se presupone revisión humana de fuentes ni réplica independiente de los papers.
