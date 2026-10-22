# Evidencia técnica: física, ensembles y uso científico

Fecha: 2026-10-03. Auditoría de lectura del código actual. Investigación asistida por IA. Estas notas documentan evidencia local y deducciones verificables; no sustituyen una campaña de validación independiente. El agente coordinador ejecuta las pruebas: este agente no ejecutó la suite ni modificó código de producción.

## Método y niveles de evidencia

Se inspeccionaron physics/config/runner/restart/tests de `modena-core` y `python/modena`. Se contrastó cada capacidad con su ruta pública y pruebas. Los documentos son claims, no demostración.

- **Código confirmado**, **deducción matemática**, **prueba existente** y **madurez no demostrada** se distinguen de mediciones nuevas.

El coordinador ejecutó microprobes mediante API público y la suite default: 152 passed, 1 ignored, 19 suites, 176.83 s. [Resultados reproducibles](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/PROBE-RESULTS.txt), [código de probes](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/audits/2026-10-03-modena-vs-lammps/probes.rs).

## Capacidades existentes y accesibles

| Capacidad | Implementación/configuración | Evidencia de validación y límite |
|---|---|---|
| MD clásica con potenciales de pares | LJ, Buckingham, Morse, soft-sphere y tablas; selección `ModelSpec`, runner | Para los cuatro analíticos existen pruebas de fuerza=-gradiente, finitud y cutoff continuo. Tabulated queda excluido expresamente de la prueba de gradiente. |
| Integración NVE | VelocityVerlet; fuerzas iniciales cacheadas por runner | Pruebas de energía/momento y propiedades; es el camino con mejor evidencia local. |
| NVT Langevin sin restricciones | BAOAB con ruido gaussiano por masa y seed | Prueba de thermalization es smoke: una temperatura final dentro de ±0.5 del target, no validación de distribución. |
| NVT Nosé–Hoover chain | NHC con masas, cadena y propagación Suzuki–Yoshida | Conservación H* relativa <1%, media T dentro de30%; no prueba independiente de fluctuaciones o respuesta. DOF se fija incorrectamente para COM/remoción y restricciones. |
| Bonded | Bonds/angles armónicos y dihedrals Fourier; topology TOML explícita | FD, suma de fuerzas y conservación; no equivalencia de force field completo de proteína. |
| Distancias rígidas / agua | SHAKE/RATTLE y proyección Newton 3×3 por agua; agua genera restricciones/exclusiones | Pruebas geométricas y conservación de una molécula; fallos de convergencia no se devuelven, presión no incluye virial de constraint. |
| Electrostatics | Ewald directo y SPME con corrección erf para exclusiones; ambos accesibles en TOML | Madelung, FD y NVE iónico; SPME incorpora virial por strain pero directo no. SPME solo ortorrómbico. |
| NPT | Berendsen y supuesto Parrinello–Rahman accesibles | Los tests prueban que el volumen cambia/permanece positivo; no distribución NPT. Hay problemas de acoplamiento y coherencia de fuerzas. |
| Restart | Caja, x/v, contador y f64 exactos | Buen seam NVE comprobado por tests; no state de RNG/NHC/pistón ni identidad del modelo. |
| Python | TOML+CLI vía subprocess, API0.1 | Es capa de configuración/lectura de log, no bindings al motor. Tipado cubre solo subset del TOML. |

Referencias: [ModelSpec](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/io/input.rs:327), [conservation gates analíticos](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/property_consistency.rs:303), [VelocityVerlet](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/integrator.rs:71), [Langevin smoke](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/conservation.rs:125), [NHC test](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/ensembles.rs:61), [bonded](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/bonded.rs:22), [electrostatics tests](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/electrostatics.rs:56), [NVE restart seam](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/restart.rs:132), [Python API](/Users/ricardoleal20/Documents/Code/Misc/modena/python/modena/config.py:141).

## Hallazgos que bloquean una declaración amplia de robustez científica

### P1. Stochastic velocity rescaling usa una distribución Gamma de escala incorrecta

**Confirmado en código y por deducción matemática.** En [thermostat.rs:125](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/thermostat.rs:125), `r1=sample_gamma((Nf−1)/2)`; el sampler Marsaglia–Tsang devuelve Gamma(shape,scale=1), visible en [149](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/thermostat.rs:149). Se introduce `r1+r2²` en [128](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/thermostat.rs:128). La variable requerida para el método Bussi es chi-cuadrado con Nf−1 grados: `2*Gamma((Nf−1)/2,1)`.

La ecuación implementada tiene media `E[K'|K]=aK+(1−a)Ktarget*(Nf+1)/(2Nf)`. Para un gas libre en régimen estacionario (ignorando el clamp, irrelevante en el límite de decorrelación) obtiene `E[K]=Ktarget*(Nf+1)/(2Nf)`, casi la mitad del target. La prueba [159–180](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/thermostat.rs:159) se deriva de esta misma distribución incorrecta y espera explícitamente `(Nf+1)/2`, por lo que puede pasar y validar la ecuación equivocada.

Para el deck de64 aguas, Nf=576 y target298.15K, la predicción de temperatura **reportada** es149.3338107638889K. [WATER_MODEL.md, Measured](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/WATER_MODEL.md:48) declara130–160K al final del NVT: es consistente con esta falla. Es una hipótesis causal fundada para esa ejecución histórica, no una medición nueva ni prueba de que sea la única causa.

**Microprobe confirmado:** Nf=12, target T=1, 200k muestras: media T=0.54118396 frente a 13/24=0.54166667 predicho. **Cierre:** distribución canónica independiente, varios Nf/tau/dt.

### P1. La ruta Langevin con restricciones cambia también los modos libres

**Confirmado en código y deducción.** BAOAB actualiza `vO=c*v+s*G` en [integrator.rs:190–206](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/integrator.rs:190), pero al existir cualquier `Constraints` reemplaza **cada** velocidad por `(xnew−xold)/dt` en [211–219](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/integrator.rs:211). Como el drift es A½–O–A½, un átomo libre que no aparece en ningún constraint queda con `v'=(v+vO)/2`. RATTLE no lo vuelve a tocar.

Para fuerza cero, `v'=(1+c)v/2+sG/2` tiene varianza estacionaria `kT/m*(1+c)/(3+c)`, convergiendo a 0.5*kT/m cuando gamma*dt es pequeño. También afecta COM rígido. **Microprobe confirmado:** gamma=100, dt=0.01, T=1, 100k muestras: varianza 1.00073504 sin constraints y 0.40729533 con dimer constraint; predicción 0.40615452.

**Cierre:** Maxwell–Boltzmann en sistemas mixtos y traslación/rotación de agua; constraints consistentes con splits BAOAB.

### P1. DOF de temperatura y termostatos siempre es3N

**Confirmado.** [system.rs:405–418](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/system.rs:405) devuelve `3*n` y calcula `T=2K/(3N)`. No contiene contador de constraints ni tratamiento COM. Inicialización [runner.rs:680–694](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/runner.rs:680) remueve COM y rescalea usando ese3N. NHC [nose_hoover.rs:107](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/nose_hoover.rs:107) y Bussi [thermostat.rs:79](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/thermostat.rs:79) reciben3N sin resta adicional pese al comentario que afirma lo contrario.

Para Nw aguas rígidas no lineales hay6Nw grados cinéticos y9Nw coordenadas; si además se impone COM=0 hay6Nw−3. Usar9Nw subreporta la temperatura física por factor2/3 aproximadamente. Un thermostat que busca Ktarget=(9Nw/2)kT puede sobrecalentar los grados reales; las otras fallas de thermostat no deben confundirse con esta ni suponerse que la compensan.

### P1. Fuerzas, PE y virial quedan obsoletos después del barostat

**Confirmado.** Runner obtiene `e.potential`,virial y fuerza a la posición postintegrador [runner.rs:373–379](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/runner.rs:373), y después aplica barostat que cambia x y caja [394–400](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/runner.rs:394). Solo vuelve a construir `Energies` usando **el mismo** `current_pe/current_virial`; no llama a model.compute. Siguiente paso usa fuerzas cacheadas para primer half-kick [integrator.rs:82](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/integrator.rs:82).

Así el log de PE/virial no corresponde al frame/volumen posterior al barostat y el primer kick siguiente usa fuerza evaluada en otra geometría. La reconstrucción de neighbors en un compute posterior no corrige ese kick previo.

**Microprobe confirmado:** una etapa, 32 átomos, dt=0.005, tau=0.1, beta=1: PE reportada −172.821524 frente a PE recomputada −189.144900; error L2 relativo de fuerzas 0.51591483. Sin barostat ambos errores 0. Acoplamiento deliberadamente fuerte para visibilidad; no generalizar esa magnitud. **Cierre:** cache coherente y NPT consistente.

### P1. Presión electrostática: Ewald directo todavía usa el virial incorrecto

**Confirmado y reconocido en pruebas actuales.** [ewald.rs:381–390](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/ewald.rs:381) acumula `Σ ri·Fi` para reciprocal y suma self/background solo a energía. En periodic reciprocal-space, el virial correcto corresponde a derivative de strain de la celda y coordenadas: contiene dependencia explícita de V y k. [spme.rs:937–940](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/spme.rs:937) llama esto el defecto de presión documentado y prueba que viriales directo/SPME **no** deben coincidir.

SPME tiene el strain-form reciprocal y testFD [spme.rs:855](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/spme.rs:855), mejora real. Sin embargo el background `Ubg=−const/V` también necesita `Wbg=3Ubg`; [spme.rs:388–391](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/spme.rs:388) devuelve Ubg pero omite ese virial. El background está presente para net charge distinto de0; en sistema neutral ese término es0.

**Criterio de cierre:** FD de energía frente a strain total en cajas diferentes, sistemas neutral/no neutral y con exclusiones; equivalencia de presión con referencia externa, no solo energía/fuerzas.

### P1. Constraint virial no existe; no hay tensor de stress real

**Confirmado.** [Constraints::shake/rattle](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/constraints.rs:136) retornan `()`; no guardan impulso ni virial. Integrador llama a model.compute y luego RATTLE antes de `Energies::from`, sin añadir contribution de constraints [integrator.rs:110–124](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/integrator.rs:110). [thermo.rs:23](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/thermo.rs:23) usa2K+modelW. La presión de sistemas rígidos, incluido agua, no está completa.

[Virial](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/force_model.rs:40) solo tiene `w:f64`. Anisotropic Berendsen y PR usan W/3 en cada eje [barostat.rs:91–102](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:91), [164–175](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:164). No pueden representar anisotropía de stress configuracional o shear. Una fuerza de par solo en x exige Wxx distinto de Wyy/Wzz; el tensor actual no puede expresarlo.

### P1/P2. El barostat llamado Parrinello–Rahman no demuestra muestreo NPT correcto

**Faltantes confirmados; ensemble distribution no medida por este agente.** [barostat.rs:155–238](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:155) hace kicks a pistón y scalea caja/posiciones. No acopla velocidades/momentos de partículas al pistón, no thermostat de pistón ni correcciónMTK, y usa el mismo W antes/después de escala. No ofrece evolución full-cell. El comentario [215–217](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:215) dice que uniform scaling preserva pair distances: preserva coordenadas fraccionales, pero multiplica distancias físicas por mu.

El test [355–368](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:355) solo exige V>0; el test de integración [ensembles.rs:188–196](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/ensembles.rs:188) cubre Berendsen y volumen diferente. No respaldan el claim «samples NPT».

Rigid water + cualquier barostat scalea también distancias internas después de SHAKE/RATTLE [barostat.rs:113](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:113); no las reimpone en final state. Hace falta una formulación compatible con molecules/constraints y testensemble, no simplemente cambiar nombre del algoritmo.

## Faltantes y defectos de amplitud funcional

### Triclinic: objeto existe, pero no es soporte validado general

TOML [StructureSpec](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/io/input.rs:98) no expone h; [SPME](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/spme.rs:333) rechaza triclinic. Minimum-image Rust redondea fraccionales [system.rs:151](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/system.rs:151), sin garantizar imagen Euclídea más cercana.

**Contraejemplo matemático exacto:** a=(10,0,0),b=(4.9,10,0),c=(0,0,10); d=(6.37,3,0)=.49a+.30b. Redondeo deja d, norma7.041086563876346. Imagen d−a=(-3.63,3,0) tiene norma4.709235606762524. Con cutoff4.8, el modelo omite una interacción cercana aun usando tilt<mitad del edge. Tests [electrostatics.rs:127–209](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/electrostatics.rs:127) cubren equivalencia diagonal/wrapping/lattice translations, no este caso.

Cell-list usa bounding-box cartesiano/wrap [neighbor.rs:332](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/neighbor.rs:332): requiere validación de caras oblicuas. `scale_uniform` [system.rs:190](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/system.rs:190) no cambia vectors/inverse; no se valida triclinic+barostat restaurado por checkpoint.

### Imported chemistry y force fields

[XYZ reader](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/io/xyz.rs:63) carga posiciones y etiquetas pero species/masses/charges mantienen defaults de [System::resize](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/system.rs:264): especie0,masa1,carga0. RunnerXYZ solo admite **una** masa/carga global para todos los átomos [runner.rs:558](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/runner.rs:558). Las etiquetas de elemento son cosméticas [system.rs:240](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/system.rs:240). Una molécula heterogénea/proteína importada no recibe atom typing, cargas ni masas realistas vía este camino.

Topology manual [exclude_bonded](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/runner.rs:174) borra 1–2, 1–3 y 1–4; sin pesos diferenciados 1–4 LJ/Coulomb. No equivale a force fields AMBER/CHARMM/OPLS: faltan importación/tipado/excepciones en las rutas core.

### Tabulated: potencial/fuerza independientes y no conservativos

[tabulated.rs:131–133](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/potential/tabulated.rs:131) interpola linealmente U y F independientes. Incluso si tabla incluye U y −U' exactos de una función suave en nodos, la derivada de Ulineal es constante por intervalo mientras Flineal generalmente no lo es. Shifted-force no elimina ese desacuerdo. Los comentarios [tests/property_consistency.rs:246](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/property_consistency.rs:246) reconocen que se excluye la tabla del testgradiente [304–309](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/property_consistency.rs:304). Es limitación conocida, no cobertura científica completa por estar «registrado».

**Cierre:** interpolación de un único funcional (p.ej. Hermite con derivadas), fuerza derivada de él; out-of-range claro y FD/NVE representativos.

### Restart: buen NVE, estado incompleto para otros ensembles

[Checkpoint](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/io/restart.rs:43) carece de RNG/NHC/pistón/model hash; no reproduce estado extendido/estocástico continuo. [runner.rs:136](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/runner.rs:136) comprueba solo N, admite masas/cargas/modelos diferentes. [restart.rs:99](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/io/restart.rs:99) trunca/escribe sin tempfile+rename atómico.

### API adicionales que existen pero no equivalen a plataforma operativa

- **MonteCarloBarostat** existe, pero [build_barostat](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/runner.rs:1027) no lo construye y `apply` es no-op [barostat.rs:321](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:321). Su aceptación tiene entropy sign incorrecto [297–298](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:297):+N ln(Vnew/Vold) dentro de dh rechazado por exp(−dh). Para variableV convencional corresponde−Nln ratio; propuesta uniforme en logV [278–280](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/barostat.rs:278) necesita además Jacobian,Hastings resultando−(N+1)ln ratio para partículas libres. No interpretar esto como defecto de todos los runs NPT: es un API no conectado.
- **ReplicaExchange** [constructor](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/replica_exchange.rs:71) fija neighborscutoff2.5 aunque acepta ForceModel arbitrario. Itera replicas secuencialmente [99–104](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/replica_exchange.rs:99) pese al comentario de parallel. Swap [160–173](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/replica_exchange.rs:160) solo cambia posiciones/PE, mantiene fuerza de configuración vieja y no swapea velocities antes de rescalearlas. Si sloti permaneceTi y conserva sus velocidades deTi, multiplicar por sqrt(Tj/Ti) les daTj, contrario al objetivo. Debe swapear configuración mecánica compatible/recalcular fuerza y mantener o transformar correctamente velocities. Testactual solo comprueba aceptación>0/energías finitas [ensembles.rs:104](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/tests/ensembles.rs:104).
- **Python** [Config](/Users/ricardoleal20/Documents/Code/Misc/modena/python/modena/config.py:144): subset reduced sin electrostatics/constraints/NHC/barostats/agua/restart. `[[types]]` sin `model` [config.py:42](/Users/ricardoleal20/Documents/Code/Misc/modena/python/modena/config.py:42), strings sin escape y paths relativas contra TOML temporal necesitan corrección.

### Configuración y fail-fast

[validate](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/io/input.rs:670) comprueba cadencias nozero, initialT>=0 y algunos campos de agua; no revisa dt>0/finitud, masa>0, tau>0,gamma>=0, índices bonded/constraints, geometría válida de la caja o cutoff vs imagen mínima. Algunos constructores sí validan potencial/cutoff y runner speciesrange, pero no es validación física integral. Constraint algorithms pueden agotar iterations sin error [constraints.rs:144–180](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/constraints.rs:144), y aguaNewton abandona si solve es singular [281–288](/Users/ricardoleal20/Documents/Code/Misc/modena/crates/modena-core/src/constraints.rs:281). Reportar convergencia/failfast es necesario para no continuar silenciosamente con geometría inválida.

## Contraevidencia a interpretaciones optimistas de los documentos

[WATER_MODEL.md](/Users/ricardoleal20/Documents/Code/Misc/modena/docs/WATER_MODEL.md:48) informa driftNVE1.5% en500steps y4.9% en2000steps; atribuirlo a «structural relaxation, not integrator defect» no demuestra la causa. En NVE exacto una relajación estructural intercambia U y K preservando E; un estado inicial violento puede amplificar errores de timestep/constraints, pero no justifica por sí solo drift de total. El single-water test con1e−4 y seam exacto prueba aspectos locales diferentes del comportamiento del líquido.

El mismo documento reconoce la comparación estática de64 aguas como **partial**: referenciaLAMMPS corregida analíticamente≈−204kcal frente aMODENA≈−301kcal, residuo≈1.5kcal/mol·water y un bonded+SHAKE deck aún pendiente. Ese contraste mezcla convenciones y métodos y no es prueba concluyente de bug; tampoco es evidencia de equivalencia completa. No afirmar producción de líquidoSPC/E validado antes de repetir con mismo Hamiltonian, box y cutoff/constraints/electrostatics y observables convergidos.

## Prioridades sugeridas al coordinador

1. Cerrar Bussi, BAOAB restringido, DOF, cache NPT y stress de constraints/Ewald.
2. Validar LJ, iones SPME y agua con referencias independientes/distribuciones T/V/P.
3. Fail-fast, restarts completos, atom typing, 1–4 e importación de force fields.
4. Triclinic, amplitud química, enhanced sampling y Python nativo son desarrollo adicional.
