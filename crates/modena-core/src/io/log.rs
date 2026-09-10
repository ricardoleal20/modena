//! CSV energy/thermostat log writer.

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::Result;

pub struct LogWriter {
    writer: BufWriter<File>,
}

impl LogWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "step,time,kinetic,potential,total,temperature,pressure,virial"
        )?;
        writer.flush()?;
        Ok(Self { writer })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn log(
        &mut self,
        step: usize,
        time: f64,
        kinetic: f64,
        potential: f64,
        total: f64,
        temperature: f64,
        pressure: f64,
        virial: f64,
    ) -> Result<()> {
        writeln!(
            self.writer,
            "{step},{time:.6},{kinetic:.10},{potential:.10},{total:.10},{temperature:.10},{pressure:.10},{virial:.10}"
        )?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}
