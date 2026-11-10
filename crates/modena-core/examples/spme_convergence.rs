//! Immutable six-phase empirical SPME calibration helper.
mod spme_calibration;
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args == ["--help"] || args == ["-h"] {
        println!("SPME calibration: prepare | verify | sweep | confirm | profile | report\nSee docs/scientific-hardening/B05-SPME-PROTOCOL.md for exact flags, authority hashes and immutable artifact contract.");
        return;
    }
    if let Err(error) =
        spme_calibration::cli::Arguments::parse(args).and_then(spme_calibration::campaign::run)
    {
        eprintln!("{error}");
        std::process::exit(error.exit_code());
    }
}
