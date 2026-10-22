# Reproducciones históricas de la revisión final

Los doce archivos de `manifest.json` se preservan byte por byte. Reprodujeron
los defectos antes del fix, durante la revisión del HEAD `10bd9ee`, usando el
core release registrado en `library_provenance` (hash `794cca85…`). No son una
campaña de rendimiento ni salidas del código corregido.

`spme_com.rs` conserva su base relativa original:
`.superpowers/sdd/2026-10-03-scientific-hardening/final-probe`. Para ejecutar la
fuente sin editarla, desde la raíz de un checkout de esa versión histórica:

```sh
rtk proxy mkdir -p .superpowers/sdd/2026-10-03-scientific-hardening/final-probe
rtk proxy cp docs/scientific-hardening/evidence/final-review-probe/charges.xyz .superpowers/sdd/2026-10-03-scientific-hardening/final-probe/charges.xyz
rtk proxy rustc --edition=2021 -O -C lto=thin -C codegen-units=1 docs/scientific-hardening/evidence/final-review-probe/spme_com.rs --extern modena_core=target/release/libmodena_core.rlib -L dependency=target/release/deps -o target/spme-review-probe
rtk proxy env RAYON_NUM_THREADS=1 target/spme-review-probe
```

El archivo de evidencia se incorporó después del HEAD histórico: hay que copiar
este directorio desde la entrega final al checkout histórico antes de esos pasos.
`alloy.rs` se compila con la misma plantilla del manifiesto, cambiando fuente y
binario; conserva su input como texto embebido. Recompilar el core histórico
requiere `cargo build --release -p modena-core` y el toolchain apropiado. El hash
de una nueva compilación puede variar según ruta/toolchain; no se promete
reproducir los bytes del binario ni del log (por ejemplo, IDs de proceso).

El fixture `charges.xyz`, inicialmente omitido por una regla de Git, se añadió
con `git add --force` sin modificarlo. Su SHA256 es
`f69364f16fc3a62970cc8782efe69ce9fa1eaa3a5c4eb2351d7e4df5351b14a9`.
El manifiesto original no se reescribe; las doce comprobaciones de hash pasan.
