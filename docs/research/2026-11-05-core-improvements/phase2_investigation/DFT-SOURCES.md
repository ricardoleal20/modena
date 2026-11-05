# Bibliografía anotada: DFT, ASE e interoperabilidad

Fase 2 — investigación bibliográfica para C5 y la frontera de C4. Consulta: **2026-11-05 UTC**. Ocho registros bibliográficos de documentación y código primarios abiertos; el registro DFT-05 agrupa tres capítulos del mismo manual de ASE. No son ocho estudios independientes ni evidencia experimental de MODENA. Las anotaciones tienen 60–100 palabras; los localizadores posteriores separan el contenido comprobado de sus límites. Búsquedas y decisiones: [DFT-SEARCH-LOG.md](DFT-SEARCH-LOG.md).

La lectura se realizó sobre las páginas oficiales y el PDF original enlazados. No se descargó un archivo de cada fuente ni se declara lectura humana. No se generaron campos de un Material Passport ni identificadores bibliográficos no comprobados. Las páginas de ASE consultadas son documentación mutable, con pie de actualización **2026-11-02 16:21:16**, sin release o commit fijados por estas URL. Los documentos de QE identifican explícitamente **7.5**. Se usa «s. f.» cuando la fuente no ofrece fecha de publicación; copyright, fecha de rastreo y actualización del sitio no se convierten en año de publicación.

## DFT-01 — Alcance electrónico de PWscf y BOMD

**Referencia APA 7:** *PWscf User’s Guide (v.7.5).* (s. f.). Quantum ESPRESSO. Recuperado el 5 de noviembre de 2026, de https://www.quantum-espresso.org/Doc/user_guide_PDF/pw_user_guide.pdf

**Anotación:** El manual oficial delimita PWscf como un solver DFT con ondas planas y pseudopotenciales que obtiene energía, orbitales Kohn–Sham, fuerzas y estrés. Describe optimización, MD de Born–Oppenheimer y cálculos electrónicos posteriores, incluidos bandas y DOS. Su pertinencia para C5 es distinguir etiquetas de energía/fuerza de las necesidades de estructura electrónica. La guía también advierte sobre convergencia de base y consistencia energía–fuerza. Es una fuente primaria de capacidades, sin demostrar precisión material ni rendimiento para una configuración concreta de MODENA.

**Evidencia y localizadores:** §§1.1, 3.3, 3.4, 3.6 y 5; páginas impresas 2, 7–10, 12 y 26–27. §3.4 identifica la MD de `pw.x` como estado fundamental electrónico, distinta de Car–Parrinello. §5 vincula errores SCF con inconsistencia energía–fuerza y pérdida de invariancia traslacional por discretización. §3.4 advierte del cambio de base en la celda final de optimización variable. [Manual original](https://www.quantum-espresso.org/Doc/user_guide_PDF/pw_user_guide.pdf).

**Límite de transferencia:** inferencia acotada para C5: reproducir E/F/estrés mediante un proveedor no demuestra reproducir orbitales, densidad, bandas o DOS. Este registro no evalúa un MLIP ni prueba paridad electrónica.

## DFT-02 — Contrato de cálculo y convergencia de pw.x

**Referencia APA 7:** *pw.x: Input description* (Quantum ESPRESSO, versión 7.5). (s. f.). Quantum ESPRESSO. Recuperado el 5 de noviembre de 2026, de https://www.quantum-espresso.org/Doc/INPUT_PW.html

**Anotación:** La especificación oficial define entradas de estructura, especies, pseudopotenciales, funcional, estado electrónico, malla k y umbrales SCF. Explica que `conv_thr` limita un error energético estimado extensivo, mientras los cutoffs controlan aproximaciones diferentes y pueden introducir ruido en fuerzas y estrés. Es pertinente para registrar la identidad de una referencia DFT en C5 y distinguir convergencia electrónica de convergencia de discretización. Su limitación es deliberada: documenta parámetros y advertencias, pero no prescribe una tolerancia universal ni una malla suficiente para todo sistema.

**Evidencia y localizadores:** entradas `tstress`, `tprnfor`, `ecutwfc`, `ecutrho`, `conv_thr`, `scf_must_converge`, `input_dft`; tarjetas `ATOMIC_SPECIES`, `ATOMIC_POSITIONS`, `CELL_PARAMETERS`, `K_POINTS`. Los valores de entrada usan unidades QE, generalmente Ry/atómicas; densificar malla k y cutoffs es distinto de cerrar SCF. El funcional predeterminado procede de los pseudopotenciales; el archivo UPF se asigna por especie. La necesidad de barrer esos parámetros hasta estabilizar el observable es una inferencia de validación para C5, no un umbral proporcionado por el manual. [Especificación original](https://www.quantum-espresso.org/Doc/INPUT_PW.html).

**Límite:** `assume_isolated` describe correcciones para una supercelda tridimensional; no equivale automáticamente a los flags de periodicidad de un motor clásico. Que una ejecución termine no certifica convergencia del observable.

## DFT-03 — Dinámica electrónica de CP y unidades temporales

**Referencia APA 7:** *cp.x: Input description* (Quantum ESPRESSO, versión 7.5). (s. f.). Quantum ESPRESSO. Recuperado el 5 de noviembre de 2026, de https://www.quantum-espresso.org/Doc/INPUT_CP.html

**Anotación:** La entrada oficial de CP explicita la propagación de grados electrónicos y una masa electrónica efectiva en el Lagrangiano CP. Documenta además minimización electrónica mediante gradiente conjugado y la recuperación de la superficie Born–Oppenheimer cuando se pierde adiabaticidad. Esta fuente permite evitar equiparar cualquier MD ab initio con CPMD o asumir que el nombre del ejecutable fija todo el algoritmo. También evidencia una diferencia de unidades temporales frente a PW. No establece que una configuración CP cualquiera produzca las mismas etiquetas que un SCF convergido.

**Evidencia y localizadores:** `electron_dynamics`, `emass`, `dt`. CP permite propagar o minimizar grados electrónicos; su timestep usa unidades atómicas Hartree, mientras PW usa unidades atómicas Rydberg, con un factor dos entre las unidades temporales indicadas. [Entrada oficial CP](https://www.quantum-espresso.org/Doc/INPUT_CP.html).

**Límite:** esta especificación no sustituye una prueba de adiabaticidad, conservación o convergencia para una trayectoria; no se ejecutó CP en esta investigación.

## DFT-04 — Propiedades disponibles y estado del calculator ASE

**Referencia APA 7:** ASE developers. (s. f.). *ase.calculators.calculator* [Código fuente documentado]. ASE documentation. Recuperado el 5 de noviembre de 2026, de https://docs.ase-lib.org/_modules/ase/calculators/calculator.html

**Anotación:** El código primario define el contrato del calculator y comprueba la propiedad solicitada contra `implemented_properties`. Si la propiedad no está soportada o falta tras el cálculo, produce un error explícito. La comparación del estado incluye posiciones, números atómicos, celda, periodicidad y estados iniciales de carga y magnetización. Es evidencia directa para que una integración distinga capacidad anunciada y resultado realmente disponible. La base no valida por sí sola el dominio químico, la convergencia o la consistencia física de cada calculator externo.

**Evidencia y localizadores:** `BaseCalculator.get_property`, `Calculator`, `compare_atoms`, `all_changes`. Una propiedad de estrés ausente debe permanecer ausente o fallar; la API no promete que todos los calculators entreguen E/F/estrés. Un cambio físico registrado puede invalidar resultados almacenados. [Código original publicado por ASE](https://docs.ase-lib.org/_modules/ase/calculators/calculator.html).

**Límite:** estas son garantías de la clase base documentada; una implementación externa o envoltorio puede añadir requisitos y comportamientos. No se ha probado un adaptador de MODENA.

## DFT-05 — Manual ASE: identidad, unidades y proveedores ML

**Referencia APA 7:** ASE developers. (s. f.). *Atomic Simulation Environment* [Documentación; capítulos «The Atoms object», «Units» y «Calculators»]. Recuperado el 5 de noviembre de 2026, de https://docs.ase-lib.org/

**Anotación:** Estos capítulos del manual describen la identidad atómica mediante números o símbolos, una celda y periodicidad por dirección, además de las unidades internas de ASE. Presentan una interfaz compartida para códigos electrónicos y potenciales, con interfaces externas para MACE, CHGNet y DeePMD-kit, entre otros. El objeto Atoms advierte que la energía DFT puede depender de la convención térmica y que no todo calculator soporta todas las propiedades. La interoperabilidad documentada facilita comparar proveedores; no certifica equivalencia de Hamiltoniano ni cobertura química de modelos ML.

**Evidencia y localizadores:**

- [The Atoms object](https://docs.ase-lib.org/ase/atoms.html): «Unit cell and boundary conditions», «Adding a calculator», API de `get_stress`. Orden Voigt: xx, yy, zz, yz, xz, xy; la contribución cinética se añade solo bajo la opción correspondiente. La elección entre energía extrapolada y energía libre consistente con fuerzas depende del calculator.
- [Units](https://docs.ase-lib.org/ase/units.html): eV, Å, masa atómica y kelvin son las unidades base especificadas; el tiempo interno no es directamente fs. Fuerza y estrés derivados corresponden a eV/Å y eV/Å³.
- [Calculators](https://docs.ase-lib.org/ase/calculators/calculators.html): distingue envoltorios de software externo, implementaciones internas e interfaces de modelos. La disponibilidad de una interfaz no incorpora el motor externo al paquete ASE.

**Límite:** los números atómicos identifican elementos; no sustituyen tipos de un force field ni expresan por sí solos carga/spin o el dominio admitido de un modelo.

## DFT-06 — Calculator Espresso y restricciones de entrada

**Referencia APA 7:** ASE developers. (s. f.). *Espresso*. ASE documentation. Recuperado el 5 de noviembre de 2026, de https://docs.ase-lib.org/ase/calculators/espresso.html

**Anotación:** La página oficial muestra la preparación de `Espresso` con perfil de ejecución, diccionario de pseudopotenciales e `input_data`. Es pertinente para separar la conversión estructural de los parámetros físicos: estos últimos deben llegar en unidades QE, sin defaults adicionales de ASE. Explica muestreo k y restricciones compatibles, y reconoce límites en validación de entradas. Existe una discrepancia interna sobre listas de k-points, registrada abajo. Esta fuente describe una integración disponible; no acredita que sus valores de ejemplo sean referencias convergidas para C5.

**Evidencia y localizadores:** secciones «Parameters», «Parallelism», «Input/Output → pw.x». El mapeo de pseudopotenciales usa símbolos de especies; se convierten `FixAtoms` y `FixCartesian`. No se soporta `ibrav` distinto de cero en el escritor documentado. La API acepta matrices k de cuatro columnas. [Documentación oficial](https://docs.ase-lib.org/ase/calculators/espresso.html).

**Discrepancia comprobada:** la lista «Not implemented» afirma que listas de k-points no están implementadas, pero el mismo texto admite `np.ndarray` N×4; DFT-07 confirma la rama implementada. No se generaliza aquella lista a «k-points arbitrarios imposibles».

**Criterio de pin para reproducirla:** registrar la release instalada de ASE y, si procede de Git, su commit; conservar el código o hash de `write_espresso_in` y la entrada k utilizada. La capacidad se atribuirá a ese pin solo al comprobar que genera la tarjeta esperada con los pesos conservados. Aquí se verificó la documentación live, sin ejecución ni pin instalado: esa comprobación permanece pendiente.

## DFT-07 — Conversión real de energía, fuerza, estrés y periodicidad

**Referencia APA 7:** ASE developers. (s. f.). *ase.io.espresso* [Código fuente documentado]. ASE documentation. Recuperado el 5 de noviembre de 2026, de https://docs.ase-lib.org/_modules/ase/io/espresso.html

**Anotación:** El lector oficial permite verificar las conversiones efectivas y no solo su descripción. Usa constantes CODATA 2006, convierte energía y fuerzas y cambia el signo del tensor de estrés respecto al output de QE. Reconstruye configuraciones con periodicidad verdadera y conserva resultados mediante un calculator de punto único. Esto expone información que un transporte de etiquetas debe preservar fuera del parseo, incluida la periodicidad original. La inspección está limitada al código publicado en la documentación mutable; no es una prueba ejecutada del parser ni de todas sus variantes de archivo.

**Evidencia y localizadores:** `_get_energy`, `_get_forces`, `_get_stress`, `read_espresso_out`, `write_espresso_in`. Conversión de fuerzas por Ry/Bohr; estrés por **−Ry/Bohr³**, reordenado a Voigt ASE. El lector asigna `pbc=True`. El escritor contiene la rama de matrices k N×4. Los comentarios explican que con smearing la energía variacional puede ser energía libre; el lector asigna el mismo valor leído a energía y energía libre. [Código primario](https://docs.ase-lib.org/_modules/ase/io/espresso.html).

**Límite:** no asumir que el parseo conserva toda la intención de frontera, pseudopotencial o energía térmica del cálculo original; DFT-02 y DFT-05 delimitan esos metadatos.

## DFT-08 — Sockets i-PI como interfaz de cálculo

**Referencia APA 7:** ASE developers. (s. f.). *Communication with calculators over sockets*. ASE documentation. Recuperado el 5 de noviembre de 2026, de https://docs.ase-lib.org/ase/calculators/socketio/socketio.html

**Anotación:** La documentación explica un servidor ASE que mueve átomos y un cliente externo que devuelve energía, fuerzas y estrés mediante el protocolo i-PI. Incluye soporte para Quantum ESPRESSO y un cliente ASE, sin exigir instalar i-PI para usar el protocolo. La conexión persistente puede evitar arranques e I/O repetidos y permitir reutilizar estados electrónicos. Es pertinente para la frontera entre integración y solver. La fuente también revela un límite esencial: los números atómicos no se transmiten y deben inicializarse coherentemente en ambos extremos.

**Evidencia y localizadores:** «Introduction», tabla de compatibilidad, ejemplo GPAW y documentación de `SocketIOCalculator`. Se transmiten coordenadas y celda; el ejemplo debe leer la estructura inicial porque no recibe los elementos por socket. Los programas pueden requerir opciones o compilación específicas para modo driver. [Fuente primaria ASE](https://docs.ase-lib.org/ase/calculators/socketio/socketio.html).

**Límite:** eliminar sobrecoste de arranque no elimina el cálculo electrónico ni prueba aceleración para MODENA. El protocolo compartido no valida unidades, identidad física, convergencia o disponibilidad de estrés de un proveedor concreto.

## Límites de esta bibliografía

Selección dirigida de documentación primaria para contratos de software; no revisión sistemática ni benchmark. No se descargaron pseudopotenciales, ejecutaron cálculos QE/CP o midieron la precisión de modelos ML. La concentración en manuales y código es intencional para estas preguntas de capacidad e interfaz; no se presentan como publicaciones revisadas por pares. Las inferencias de validación señaladas son preguntas pendientes para C5, no resultados materiales. Investigación y redacción realizadas con asistencia de IA.
