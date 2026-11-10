//! Strict six-mode argument contract. Paths and authority are explicit.
use super::*;

pub struct Arguments {
    pub mode: String,
    pub values: BTreeMap<String, String>,
    pub extended: bool,
}
impl Arguments {
    pub fn parse(args: Vec<String>) -> Result<Self, CalibrationFailure> {
        let mode = args
            .first()
            .ok_or_else(|| invalid("expected prepare/verify/sweep/confirm/profile/report"))?
            .clone();
        let (required, optional): (&[&str], &[&str]) = match mode.as_str() {
            "prepare" => (&["--out", "--family", "--replicas", "--threads"], &[]),
            "verify" => (
                &["--protocol", "--snapshots"],
                &[
                    "--special",
                    "--special-sha256",
                    "--madelung",
                    "--madelung-sha256",
                    "--elapsed-ns",
                ],
            ),
            "sweep" => (
                &[
                    "--protocol",
                    "--snapshots",
                    "--verification-sha256",
                    "--out",
                ],
                &["--elapsed-ns"],
            ),
            "confirm" | "profile" => (
                &["--campaign", "--expected-sha256", "--out"],
                &["--elapsed-ns"],
            ),
            "report" => (
                &[
                    "--campaign",
                    "--expected-sha256",
                    "--memory",
                    "--memory-sha256",
                    "--out",
                ],
                &[
                    "--profiles",
                    "--profiles-sha256",
                    "--profile-receipt-sha256",
                ],
            ),
            _ => {
                return Err(invalid(
                    "expected prepare/verify/sweep/confirm/profile/report",
                ))
            }
        };
        let mut values = BTreeMap::new();
        let mut extended = false;
        let mut i = 1;
        while i < args.len() {
            let key = &args[i];
            if mode == "prepare" && key == "--reference-extended" && !extended {
                extended = true;
                i += 1;
                continue;
            }
            if !required.contains(&key.as_str()) && !optional.contains(&key.as_str()) {
                return Err(invalid(format!("unknown or duplicate option: {key}")));
            }
            let value = args
                .get(i + 1)
                .filter(|v| !v.starts_with("--"))
                .ok_or_else(|| invalid(format!("missing value for {key}")))?;
            if values.insert(key.clone(), value.clone()).is_some() {
                return Err(invalid(format!("duplicate option: {key}")));
            }
            i += 2;
        }
        for key in required {
            if !values.contains_key(*key) {
                return Err(invalid(format!("missing {key}")));
            }
        }
        for (path, hash) in [
            ("--special", "--special-sha256"),
            ("--madelung", "--madelung-sha256"),
            ("--profiles", "--profiles-sha256"),
        ] {
            if values.contains_key(path) != values.contains_key(hash) {
                return Err(invalid(format!("{path} and {hash} must appear together")));
            }
        }
        if values.contains_key("--profiles") != values.contains_key("--profile-receipt-sha256") {
            return Err(invalid(
                "profile import requires authoritative profile receipt hash",
            ));
        }
        Ok(Self {
            mode,
            values,
            extended,
        })
    }
    pub fn value(&self, key: &str) -> Result<&str, CalibrationFailure> {
        self.values
            .get(key)
            .map(String::as_str)
            .ok_or_else(|| invalid(format!("missing {key}")))
    }
    pub fn path(&self, key: &str) -> Result<&Path, CalibrationFailure> {
        Ok(Path::new(self.value(key)?))
    }
}
