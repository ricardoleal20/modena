# Investigación del core de MODENA

Solicitud: investigar mejoras puntuales para acercar MODENA a sustituir flujos de simulación atomística actualmente cubiertos por LAMMPS, Quantum ESPRESSO y FAIR-Chem.

**Estado de la investigación: terminada, Fases 1–6.** Confirmación recibida: «Si, confirmo». Informe y backlog revisados. Los documentos de investigación conservan su base y sus conclusiones originales. Después, el usuario autorizó «Aplica el B06»; su implementación se documenta por separado en [B06](../../B06_ENERGY_DIAGNOSTICS.md) y su [verificación](../../B06_VERIFICATION.md).

Base local: `f166399070b2d6ebfd91f51e5934c132f04d46c0`. La fecha de la carpeta identifica la sesión UTC; no representa un corte de publicaciones ni una fecha de medición de resultados externos.

## Lectura recomendada

1. [Informe final](phase6_revision/REPORT.md): dictamen, ocho frentes de mejora, arquitectura propuesta y puertas de aceptación.
2. [Backlog por archivo](phase3_analysis/ACCEPTANCE-BACKLOG.md): 12 tareas B01–B12, riesgos, criterios refutables, alternativas y tamaños relativos.
3. [Fuentes verificadas](phase2_investigation/SOURCE-VERIFICATION.md): 38 registros, 37 unidades bibliográficas, condiciones por afirmación y acceso.
4. [Revisiones y cambios finales](phase6_revision/REVISION-LOG.md): tres checkpoints PASS, editorial Accept e integridad CLEARED; son dictámenes de agentes.

## Protocolo y evidencia

- [Pregunta y alcance](phase1_scope/RESEARCH-QUESTION.md): casos de uso verificables, tres subpreguntas y FINER 4.4/5.
- [Metodología](phase1_scope/METHODOLOGY.md): fuentes primarias, trazabilidad de evidencia, priorización y puertas G0–G6.
- [Inventario local](phase1_scope/CODE-INVENTORY.md): observaciones del código y preguntas a investigar; distingue correcciones ya implementadas de posibles siguientes mejoras.
- [Checkpoint crítico 1](phase1_scope/CHECKPOINT-1.md): PASS, sin Critical/Major. Las guardas sobre selección de casos y alternativa de integración se incorporaron al método.
- [Verificación local](phase1_scope/VERIFICATION.md): `cargo test --locked --offline --workspace` terminó con 284 pruebas aprobadas y una ignorada; no es evidencia de paridad externa.
- [Casos congelados](phase2_investigation/CASES.md): LJ, iones, agua, ML y DFT, con enmienda explícita de la alternativa MC molecular.
- [Bibliografía MD clásica](phase2_investigation/CLASSICAL-SOURCES.md), [ML](phase2_investigation/ML-SOURCES.md), [DFT/ASE](phase2_investigation/DFT-SOURCES.md) y [validación](phase2_investigation/VALIDATION-SOURCES.md), con sus registros de búsqueda.
- [Observación algebraica del diagnóstico](phase2_investigation/DIAGNOSTIC-OBSERVATION.md): dependencia del cero energético; no es una trayectoria ejecutada.
- [Síntesis y tensiones](phase3_analysis/SYNTHESIS.md): evidencia integrada, alternativas y límites de transferencia.

## Alcance del resultado

La confirmación exigida por la habilidad [deep-research](/Users/ricardoleal20/.agents/skills/deep-research/SKILL.md) se recibió expresamente. Esa confirmación abarcó la investigación; la implementación posterior de B06 tiene autorización independiente del usuario.

Las recomendaciones son propuestas de ingeniería. Durante la investigación no se ejecutaron nuevas campañas LAMMPS, ML, DFT, GPU o MPI; no se prometen aceleraciones ni paridad global. La elección del material, propiedades y recursos precede a una campaña de aceptación material. B06 añade diagnóstico y análisis estadístico; las otras recomendaciones siguen como propuestas.

Material elaborado con asistencia de IA. Los revisores automáticos no equivalen a revisión humana independiente. Las hipótesis de mejora no son mediciones de aceleración ni validaciones materiales.
