# Devil’s Advocate — Checkpoint 2

Fecha: 3 de octubre de 2026. Material revisado: `SYNTHESIS.md`, `PERFORMANCE-NOTES.md`, `probes.rs`, `PROBE-RESULTS.txt`, el corpus primario previamente verificado y los puntos pertinentes de `runner.rs`, `integrator.rs` y `thermostat.rs`. Revisión de solo lectura; no se modificó código de producción.

## Veredicto: PASS

No se identifican problemas críticos que invaliden la síntesis. La evidencia sustenta un motor existente y capaz, con defectos físicos concretos que impiden atribuir confiabilidad general a todas sus rutas disponibles. La síntesis también limita correctamente los resultados: no invalida NVE de pares por errores de otros algoritmos y no afirma un ranking universal de velocidad.

## Problemas menores que debe resolver el informe

1. **Muestras no equivalen a muestras independientes.** Los probes de Bussi y BAOAB usan una semilla y series temporales correlacionadas. El número de pasos refuerza la reproducibilidad, pero no puede presentarse como tamaño estadístico efectivo ni como estimación de incertidumbre. Recomendación: llamar a las cifras medias de series con la longitud indicada, no certificaciones estadísticas exhaustivas. La desviación grande, el control BAOAB y la coincidencia con predicciones analíticas hacen sólida la existencia de los defectos aun sin estimar ESS o intervalos de confianza. [Código del probe](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/probes.rs), [resultados](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/PROBE-RESULTS.txt), [metodología de validación física](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202764).

2. **La razón histórica de Ewald no identifica una brecha de rendimiento bajo el mismo Hamiltoniano.** La nota de rendimiento explicita diferencias en fuerzas Buckingham y condiciones de ejecución. Por ello, una razón aproximada de 8.9 es descripción del artefacto histórico, no una penalización aislada del solver ni un factor MODENA/LAMMPS físicamente equiparable. Recomendación: mantener la advertencia junto a esa cifra o reservarla al apéndice. La conclusión general defendible sigue siendo que no existe una ventaja demostrada bajo condiciones equivalentes. [Notas de rendimiento](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/PERFORMANCE-NOTES.md), [criterios oficiales de medición](https://docs.lammps.org/Speed_measure.html).

3. **“Robustez general refutada” necesita una definición explícita y localización por ruta.** Los contraejemplos rechazan la hipótesis de que todas las opciones actualmente disponibles satisfacen sus propiedades físicas nominales. No demuestran que toda investigación con MODENA sea inválida. Recomendación: formular el dictamen como falta de confiabilidad demostrada para uso general y precisar las rutas: CSVR/Bussi, BAOAB con constraints y estado de fuerzas/energía después de aplicar el barostato del runner. La síntesis ya hace buena parte de esta delimitación; conservarla en el resumen y la respuesta final. [Código del probe](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/probes.rs), [resultados](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/PROBE-RESULTS.txt).

## Revisión de contraevidencia

Las pruebas que pasan, el código de interacciones e integradores y las mejoras internas de rendimiento son contraevidencia relevante frente a la tesis exagerada de que MODENA es sólo una demo o es inutilizable. La síntesis las conserva. Por otra parte, pasar una suite de software no contradice un contraejemplo físico cuando la suite no comprueba esa propiedad. La existencia de A7 de Bussi, su corrección del autor y una implementación independiente que usa el factor 2 fortalece una triangulación que no depende de asumir infalibilidad de LAMMPS. [Fuentes CSVR](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/SOURCES.md), [corrección primaria del autor](https://sites.google.com/site/giovannibussi/Research/algorithms), [código LAMMPS](https://github.com/lammps/lammps/blob/stable/src/EXTRA-FIX/fix_temp_csvr.cpp).

La prueba de barostato usa acoplamiento deliberadamente fuerte: demuestra discrepancia del estado almacenado frente al estado físico presente, pero su error relativo no puede extrapolarse a parámetros típicos ni interpretarse como una desviación experimental universal. La síntesis declara esta limitación. [Configuración del probe](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/probes.rs), [salida registrada](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/PROBE-RESULTS.txt).

## Contrargumento más fuerte

MODENA podría ser útil y competitiva para investigación delimitada aunque dos métodos de control térmico y una ruta de barostato tengan errores. LAMMPS también tiene paquetes opcionales, restricciones y combinaciones que necesitan validación; exigir paridad de catálogo como condición para publicar investigación sería una vara equivocada. Un dominio pequeño, explícito y validado puede constituir una plataforma de investigación válida.

Este contrargumento limita el alcance del diagnóstico, pero no rehabilita las rutas defectuosas ni demuestra rendimiento equivalente. La consecuencia práctica es definir y validar el dominio que se quiere garantizar, no declarar inútil la plataforma ni proclamar madurez general prematuramente. [Capacidades y restricciones oficiales](https://docs.lammps.org/Intro_features.html), [validación física](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202764).

## Qué evidencia falta

- Réplicas independientes, incertidumbre y observables científicos de sistemas de producción; estos probes son contraejemplos controlados, no certificación de líquidos o proteínas reales.
- Comparación actual equivalente de paso MD completo, con Hamiltoniano, precisión, recursos, I/O y longitud explícitos; la nueva comparación exploratoria todavía no forma parte de esta síntesis.
- Pruebas HPC/GPU y migración distribuida de producción; no pueden deducirse de ejemplos estáticos ni de revisión de arquitectura.
- Trazabilidad completa de los comandos, herramientas y SHA del probe en el informe final, para que el lector reproduzca los artefactos.

## Pruebas de estrés

| Prueba | Resultado |
|---|---|
| Retirar el benchmark más favorable o desfavorable | La conclusión de falta de ventaja general demostrada se mantiene. |
| Retirar LAMMPS como autoridad de CSVR | La fórmula de normales cuadradas, la errata del autor y el contraejemplo de MODENA mantienen el diagnóstico. |
| Invertir la hipótesis: MODENA sí sirve para investigación delimitada | Es plausible y compatible con los datos; requiere validar el caso de uso concreto. |
| Generalizar a otros Hamiltonianos, condiciones o hardware | No está justificado y la síntesis lo restringe. |
| Utilidad práctica | Las evidencias permiten priorizar corrección física, validación independiente, benchmarks comparables e integración de backends. |

Análisis asistido por IA. No se identificó un conflicto de fuentes que obligue a bloquear esta fase; las limitaciones anteriores deben quedar visibles en la redacción final.
