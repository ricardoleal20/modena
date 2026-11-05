# Observación algebraica del indicador de conservación

Tipo: inspección de fórmula y ejemplo aritmético, **no ejecución de trayectoria ni medición de rendimiento**. Base: `f166399070b2d6ebfd91f51e5934c132f04d46c0`, [conservation.rs](../../../../crates/modena-core/src/conservation.rs), métodos `ConservationReport::new`, `observe` y `ConservationVerdict::from_drift`.

La fórmula implementada es `(E_final-E0)/max(abs(E0),1e-12)`, con umbrales absolutos del indicador de 1e-4 y 1e-3. Si se suma una constante a la energía potencial, las fuerzas no cambian pero sí el denominador. La conclusión es algebraica y acotada al indicador; no demuestra errores en fuerzas, integradores o ensembles.

| Cero energético | E0 | E_final | ΔE aproximadamente | Indicador aproximadamente | Etiqueta de la fórmula |
| --- | --- | --- | --- | --- | --- |
| Original | 1 | 1.0002 | 0.0002 | 2e-4 | Warning |
| Mismos valores +999 | 1000 | 1000.0002 | 0.0002 | 2e-7 | Healthy |

Se evaluó la misma expresión y clasificación en Python estándar, sin invocar el motor. Los últimos bits de ΔE varían por representación flotante, sin alterar el ejemplo. Reproducción:

```sh
rtk proxy python3 -c 'for c in (0.0,999.0):
 e0=1.0+c; ef=1.0002+c; d=(ef-e0)/max(abs(e0),1e-12)
 print(c,d,"Healthy" if abs(d)<1e-4 else "Warning" if abs(d)<1e-3 else "Drifting")'
```

La documentación del módulo declara calibración LJ y uso consultivo. El [runner](../../../../crates/modena-core/src/runner.rs), `run_streaming`, observa energía total también al aplicar controladores térmicos y de volumen. El alcance de cada etiqueta requiere considerar el ensemble y su metadato. El hallazgo no invalida los tests científicos existentes ni demuestra que un usuario haya interpretado mal la salida.

Elaborado con asistencia de IA. Fuente matemática local y comando reproducible; no hay fuente bibliográfica externa para el ejemplo numérico.
