# Opt-in update diagnostics

`diagnostics` and its two copied kernel variants are compiled only for tests. Production
factorization and the Lanczos correctness fix are unchanged. The baseline is the existing
`products` + `subtract`; the candidate reproduces the reverted two-dimensional update.
Both originate in this repository, based on the parent module's Ng–Peyton/Rothberg–Gupta
references. Keep MR=8, NR=4 and ascending-k accumulation unchanged.

Run the ignored `cholesky::diagnostics::four_phase_update_experiment` test explicitly in
release. `FC_DIAG_SIZE=tiny` selects the capture smoke test; otherwise it runs four independent
2048-row, 1024-column, 512-term updates. `FC_DIAG_THREADS` selects workers;
`FC_DIAG_VARIANT=baseline|candidate|single|kernel` selects the path.
`single` serializes inner work, retaining outer parallelism; `kernel` changes only the
microkernel's compilation boundary. `FC_DIAG_PROFILE=0` disables fine-grained clocks/counters.
No environment variable affects production solving or target CPU selection.

`DIAG_SUMMARY.command_seconds` is the complete update experiment, not `solve.static` or
`solve.frequency`. Input creation and bitwise comparison are outside it; a process timer also
includes those. `DIAG_PHASE` records packing, dense arithmetic, scheduling preparation and
writeback, with call counts and modelled logical bytes. Bytes describe source-level operand
streams (including padded upper bounds), **not measured DRAM traffic**. Thread-span sums
include preemption and overlap and are **not CPU time or additive wall time**.
`DIAG_DISPATCH` reports inclusive region spans, task spans and ready-to-start spans;
ready-to-start includes normal finite-worker queueing, not just scheduler overhead.
The packing allocation/dispatch envelope also contributes to unassigned whole-command time.
Counters are thread-local and drained at task boundaries, outside inner dot products.

Always compare fine-timing on/off, and compare every output bit with production
`products` + `subtract`. Code-generation changes can affect even an unselected branch's
performance: do not attribute across-build gains solely to the isolated-kernel path.
The experiment does not adopt the reverted candidate, establish FE-model performance or
replace frozen results. No full solver/oracle rerun is needed for these test-only copies.
