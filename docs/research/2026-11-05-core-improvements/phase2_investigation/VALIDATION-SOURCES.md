# Fuentes de validación física y presión

Consulta: 2026-11-05 UTC. Fuentes primarias abiertas mediante búsqueda y lectura de contenido. La calificación es relativa al propósito: una documentación oficial demuestra una interfaz o algoritmo declarado; no demuestra aceleración de MODENA. No se aplican jerarquías clínicas a experimentos de software.

## V1 — Validez física

Merz, P. T., & Shirts, M. R. (2018). Testing for physical validity in molecular simulations. *PLOS ONE, 13*(9), e0202764. https://doi.org/10.1371/journal.pone.0202764

- Texto completo primario: [PLOS](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202764). Publicación revisada por pares; identidad bibliográfica y métodos visibles. Confianza alta para el método.
- Aporta pruebas de convergencia del integrador, equipartición y distribuciones termodinámicas, además de integración en validación de releases. Una trayectoria con energía aparentemente estable no prueba todas estas propiedades.
- Límite: casos concretos y métodos de prueba; no valida el checkout de MODENA ni certifica todos los sistemas que pasen controles simples. Autores desarrolladores de la biblioteca de validación; independencia algorítmica debe evaluarse al reutilizarla.

## V2 — Contraste de ensembles

Shirts, M. R. (2013). Simple quantitative tests to validate sampling from thermodynamic ensembles. *Journal of Chemical Theory and Computation, 9*(2), 909–926. https://doi.org/10.1021/ct300688p

- Identidad y referencia editorial corroboradas en [el registro del autor en arXiv](https://arxiv.org/abs/1208.0910); contenido de la editorial recuperado por búsqueda, apertura DOI con HTTP 403. Se usa el manuscrito primario [PDF](https://arxiv.org/pdf/1208.0910) para métodos; no se presume acceso íntegro al PDF editorial.
- El contraste de simulaciones emparejadas permite probar la distribución termodinámica sin conocer la densidad de estados. Relevante para aceptar NVT/NPT.
- Límite: exige datos estadísticos suficientes y solapamiento adecuado. No diagnostica por sí solo la causa de un fallo ni valida la superficie de energía.

## V3 — Equilibración y muestras efectivas

Chodera, J. D. (2016). A simple method for automated equilibration detection in molecular simulations. *Journal of Chemical Theory and Computation, 12*(4), 1799–1805. https://doi.org/10.1021/acs.jctc.5b00784

- Manuscrito completo y metadatos comprobados en [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4945107/), además del [registro del laboratorio autor](https://www.choderalab.org/publications/2015/6/30/a-simple-method-for-automated-equilibration-detection-in-molecular-simulations). La apertura del DOI editorial devolvió 403.
- Propone seleccionar descarte inicial maximizando muestras efectivamente independientes. Da una base para analizar incertidumbre considerando autocorrelación.
- Límite: diagnóstico de una serie observable; no garantiza haber cruzado barreras lentas no observadas ni asegura convergencia termodinámica global.

## V4 — Herramientas de series temporales

PyMBAR developers. (s. f.). *The timeseries module pymbar.timeseries*. [Documentación oficial](https://pymbar.readthedocs.io/en/stable/timeseries.html).

- Página consultada rotulada `4.2.1.dev0+gcf12100f8.d20251016`; no se infiere versión instalada. Documenta `detect_equilibration`, ineficiencia estadística y subsampling.
- Evidencia alta de la API documentada, no evidencia de que el core ya la utilice. El método se declara heurístico; las opciones aproximadas y búsqueda binaria tienen límites explícitos.

## V5 — Restricciones y virial

LAMMPS developers. (s. f.). *fix shake command; fix rattle command*. [Manual oficial](https://docs.lammps.org/fix_shake.html).

- Versión visible: release, git info `30Sep2026`. Confirma contribuciones de restricciones a presión y esfuerzo y sus condiciones de uso.
- Sirve como referencia operativa del término ausente en la presión restringida actual de MODENA. No se copia sin derivar la discretización compatible con VV/BAOAB y cuerpos rígidos del motor.

## V6 — Tensor de presión

LAMMPS developers. (s. f.). *compute pressure command*. [Manual oficial](https://docs.lammps.org/stable/compute_pressure.html).

- Versión visible: stable, git info `30Sep2026 update1`. Documenta presión escalar, tensor simétrico y contribuciones de interacciones/restricciones.
- Fundamenta separar virial configuracional, término cinético y convención de estrés. El resultado no puede reconstruirse siempre desde fuerza total y posiciones envueltas: hay dependencia explícita de celda y periodicidad.
- Límite: capacidad/documentación, no validación de un nuevo tensor en MODENA.

## V7 — NPT Monte Carlo molecular

OpenMM developers. (s. f.). *Standard forces: MonteCarloBarostat*. En *OpenMM User Guide 8.6*. [Algoritmo oficial](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#montecarlobarostat).

- Escala centros moleculares y acepta usando energía, presión externa, cambio de volumen y término entrópico dependiente del número de moléculas. Requiere control térmico compatible; no necesita virial instantáneo para la aceptación.
- Esto respalda una alternativa a implementar primero un barostato dinámico tensorial para agua rígida. No autoriza reusar el MC atómico de MODENA sin cambiar medida, propuesta y validación.
- Límite: algoritmo descrito y estimador molecular particular; no prueba equivalencia de todos los estimadores de presión o NPT restringido general.

## Registro de búsqueda y verificación

Cadenas ejecutadas:

1. `Merz Shirts 2018 Testing for physical validity in molecular simulations PLOS 0202764`
2. `Shirts 2013 Simple quantitative tests to validate sampling from thermodynamic ensembles ct300688p`
3. `Chodera 2016 simple method automated detection equilibration times molecular simulations 10.1021`
4. `site.docs.openmm.org Monte Carlo barostat molecule centers coordinates pressure energy volume acceptance`
5. `site.docs.lammps.org fix shake virial constraint pressure stress tensor`

Se abrieron páginas primarias de PLOS, arXiv, PMC, PyMBAR, LAMMPS y OpenMM; se comprobaron además título, autor y referencia editorial de V1–V3. Los resultados secundarios —Wikipedia, ResearchGate, Reddit y espejos de terceros— se excluyeron como apoyo técnico. No se informa un total de hits del buscador ni un flujo PRISMA: la interfaz devuelve resultados limitados. Dos URLs editoriales DOI devolvieron 403; se documenta el acceso alternativo primario, sin inventar lectura del texto editorial.

Los enlaces V1–V7 son siete entradas, de las cuales tres son artículos y cuatro documentación. Los registros arXiv/editorial de V2 y laboratorio/PMC de V3 representan la misma obra, no estudios adicionales. Elaborado con asistencia de IA.
