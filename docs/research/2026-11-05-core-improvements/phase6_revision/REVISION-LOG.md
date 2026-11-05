# Revisión final y respuesta a observaciones

Fase 6, iteración 1. Fuente preservada: [borrador de Fase 4](../phase4_composition/DRAFT.md). Resultado: [informe final](REPORT.md). No se modificó código de producción ni se implementaron recomendaciones.

## Dictámenes recibidos

- [Checkpoint 1](../phase1_scope/CHECKPOINT-1.md): PASS, sin Critical/Major.
- [Checkpoint 2](../phase3_analysis/CHECKPOINT-2.md): PASS, sin Critical/Major.
- [Revisión editorial](../phase5_review/EDITORIAL.md): Accept, puntuación de agentes 4.45/5, dos Minor.
- [Revisión de integridad](../phase5_review/INTEGRITY.md): CLEARED, sin correcciones obligatorias. Contraste documentado de 21/37 unidades con modalidad declarada.
- [Checkpoint 3](../phase5_review/CHECKPOINT-3.md): PASS, sin Critical/Major, un Minor.

Son revisiones asistidas por agentes de la misma sesión, algunos participantes también en bibliografía o alcance. No equivalen a arbitraje humano independiente ni certifican el motor.

## Cambios aplicados

| ID | Procedencia y gravedad | Observación | Acción en la versión final | Estado |
| --- | --- | --- | --- | --- |
| DA2-1 | Checkpoint 2, Minor | Evitar una generalización extensa como paso inevitable antes de probar un proveedor | El reporte empieza con snapshot, proveedor persistente, paridad E/F y ausencias explícitas; permite que el adaptador sea solución suficiente | Resuelto ya en el borrador y conservado |
| DA2-2 | Checkpoint 2, Minor | C4/C5 no son un benchmark químico ya fijado | Separa paridad del contrato de eficacia química; exige material/modelo/referencia y evaluación fijados antes de conclusiones del dominio | Resuelto ya en el borrador y conservado |
| ED-01 | Editorial, Minor | Antecedente ambiguo de «segunda obra» al hablar de error P3M | Nombra expresamente Deserno y Holm (1998b; CL10) | Resuelto |
| ED-02 | Editorial, Minor | Reproducibilidad de búsquedas podía sugerir idénticos hits | Precisa repetición de estrategia/auditoría de decisiones sin garantizar conjunto u orden de resultados | Resuelto |
| DA3-1 | Checkpoint 3, Minor | Fila de nueva química podía posponer selección del dominio | Antes de la tabla exige decidir flujo/material, propiedades y recursos; distingue esa elección de implementar física nueva | Resuelto |
| NAV-01 | Editorial, sugerencia | Navegación entre ocho áreas y doce tareas | Añade IDs B01–B12 a los títulos de las áreas correspondientes | Incorporado |
| TECH-01 | Revisión del orquestador, precisión | Declarar PBC no implementa automáticamente fronteras parciales | El vertical mínimo puede limitarse a PBC tridimensional y rechazar otros casos; también se precisa el alcance inicial del backlog | Incorporado |
| TECH-02 | Revisión del orquestador, precisión | Reutilizar scratch no elimina su memoria proporcional a N por partición | Se aclara en el cuerpo, consistente con B02 | Incorporado |

La corrección del localizador eSEN A.1 se aplicó en Fase 2 antes de la composición; queda documentada en [la verificación de fuentes](../phase2_investigation/SOURCE-VERIFICATION.md). No hubo una segunda iteración de revisión porque no quedaron observaciones Critical/Major ni condiciones obligatorias de integridad y las modificaciones finales fueron precisiones acotadas.

## Límites que permanecen declarados

No hay perfiles nuevos, speedups, validación de líquidos, migración MPI real, inferencia de un modelo ni referencia DFT convergida ejecutados en este estudio. Las 284 pruebas aprobadas y una ignorada pertenecen al comando de la [verificación local registrada](../phase1_scope/VERIFICATION.md). No se repitió esa suite después de cambiar solo documentación. Las fuentes mutables, acceso incompleto y ausencia de búsqueda sistemática de retractaciones permanecen en el informe.

La necesidad de fijar material/hardware para campañas futuras es una limitación de las propuestas, no una revisión pendiente de este informe. No se requiere una nueva confirmación del usuario para cerrar la investigación autorizada.

Elaborado con asistencia de IA.
