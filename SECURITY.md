# Security policy

## Scope and supported code

MODENA is in early development. Security fixes are prioritized for the current
code on `main` and the latest release; older releases have no separate
maintenance commitment. An issue report does not imply a guaranteed response
time or a support agreement.

Examples of security concerns include unintended file access, command
execution, exposure of private data, or unsafe processing of untrusted inputs.
Numerical disagreements without a security impact belong in the
[numerical discrepancy issue form](https://github.com/ricardoleal20/modena/issues/new?template=numerical_discrepancy.yml).

## Reporting a vulnerability

Do not disclose exploit details, credentials, or private data in a public issue
or pull request. Use GitHub's [private vulnerability report form](https://github.com/ricardoleal20/modena/security/advisories/new).
Private vulnerability reporting is enabled for this repository. If the form is
unavailable, use a private contact method listed on the maintainer's
[GitHub profile](https://github.com/ricardoleal20). If no private route is
available, open an issue asking for one without including vulnerability details.

After a private route is established, include:

- The affected version or commit, component, platform, and relevant configuration.
- A minimal reproducer or proof of concept and the steps needed to trigger it.
- The expected impact, prerequisites, and any suggested mitigation.
- Whether the issue has already been disclosed elsewhere.

Use only data you are authorized to share and remove secrets from evidence.
The maintainer will assess the report and coordinate a fix and disclosure
when possible. Agree on public disclosure timing through the private channel;
there is no fixed response or remediation deadline.

## Running MODENA

Input files can select output paths, and simulations can consume substantial
CPU time, memory, and disk space. Review configurations and run untrusted jobs
with appropriate operating-system isolation and resource limits. The MCP
interface does not make MODENA a sandbox for untrusted workloads.
