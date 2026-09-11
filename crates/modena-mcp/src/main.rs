//! `modena-mcp` — the MCP stdio server binary. Speaks newline-delimited
//! JSON-RPC 2.0 over stdin/stdout so any MCP-capable agent (Claude, Cherry
//! Studio, etc.) can drive MODENA simulations through typed tools.
//!
//! Run with `modena-mcp`. Logging goes to stderr (stdout is reserved for the
//! protocol). The working directory is the process cwd; configure it by
//! launching from the desired directory.

use std::io::{self, BufReader, BufWriter};

use modena_mcp::{server, RunManager, DEFAULT_MAX_CONCURRENT};

fn main() -> io::Result<()> {
    let cwd = std::env::current_dir()?;
    let mgr = RunManager::new(DEFAULT_MAX_CONCURRENT, cwd);

    let stdin = io::stdin();
    let stdout = io::stdout();

    // Lock once for the whole session. stdin is line-buffered by `BufReader`;
    // stdout is wrapped in a `BufWriter` and flushed after every response inside
    // `server::serve`.
    let input = BufReader::new(stdin.lock());
    let output = BufWriter::new(stdout.lock());

    if let Err(e) = server::serve(&mgr, input, output) {
        eprintln!("modena-mcp: server error: {e}");
        return Err(e);
    }
    Ok(())
}
