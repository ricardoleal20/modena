# Alcance y protocolo de auditoría

Fecha: 3 de octubre de 2026, America/Monterrey.
Snapshot inicial: `b68abd112339aa6970261f470a572da6633bc617`.
Solicitud: evaluar MODENA actual frente a LAMMPS en completud, velocidad y madurez para investigación.

## Pregunta y subpreguntas

¿Qué madurez científica y computacional demuestra MODENA actual frente a LAMMPS para investigación con dinámica molecular?

1. ¿Qué modelos físicos y flujos de trabajo están implementados, expuestos y validados?
2. ¿Qué evidencia permite confiar en los resultados y qué defectos limitan su uso científico?
3. ¿Qué rendimiento comparable está demostrado y qué cambios permitirían competir en cargas concretas?

Las tres subpreguntas heredan el mismo repositorio y snapshot. No se equipara dinámica molecular clásica con estructura electrónica.

FINER: factibilidad 4/5 (acceso local; sin cluster), interés 5/5, novedad 3/5 (auditoría aplicada), ética 5/5, relevancia 5/5; media 4.4/5.

## Método

Auditoría técnica comparativa: lectura de código, configuración, pruebas, CI y artefactos de validación; ejecución de verificaciones existentes; contraste con documentación oficial y publicaciones primarias de LAMMPS. No es una revisión sistemática PRISMA ni una medición exhaustiva de todas las combinaciones de hardware.

Jerarquía ajustada al dominio: código actual y ejecución reproducible para comportamiento de MODENA; documentación oficial para capacidades de LAMMPS; benchmarks conservados con metodología y condiciones para rendimiento histórico; README y roadmap como declaraciones que requieren corroboración.

Cada capacidad se clasifica por separado en implementada, accesible, validada y escalable. Ausencia de evidencia significa no verificado, salvo ausencia o defecto confirmado por código.

## Criterios fijados antes de la investigación

- Consistencia energía–fuerzas, continuidad pertinente y comprobaciones analíticas independientes.
- Convergencia al reducir paso temporal y aumentar precisión electrostática.
- Conservación NVE, momento y restricciones, con métricas físicamente interpretables.
- Distribuciones y fluctuaciones de los ensembles; presión/virial y grados de libertad correctos.
- Reinicio reproducible y estado completo del algoritmo cuando sea necesario.
- Interoperabilidad y trazabilidad del modelo, parámetros, unidades y resultados.
- Rendimiento a igual problema físico, precisión, hardware y recursos; distinguir wall time, CPU time y coste de inicialización/I/O.
- Robustez por clase de sistema, evitando extrapolar una validación de LJ a toda química molecular.

## Revisión crítica 1

Agente independiente `critical_review`: PASS, sin problemas críticos. Observación menor resuelta incorporando los criterios anteriores. Contrargumento central: un motor pequeño puede ser publicable y competitivo en un dominio validado sin paridad funcional general; una victoria local tampoco implica madurez general.

## Límites y transparencia

Se conservará la separación entre observado, medido, documentado e inferido. Los ensayos locales no sustituyen benchmarks HPC. El análisis utiliza herramientas de IA y revisión de agentes; sus conclusiones son juicios técnicos sustentados en evidencia, no una certificación.

El usuario confirmó el alcance en la conversación antes de iniciar la investigación detallada. La confirmación responde al requisito explícito de la skill `deep-research`: «User confirmation required after Phase 1 before proceeding».
