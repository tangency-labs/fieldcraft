# Crash reports, the log file and Report a Bug

## Privacy

No crash report, log line or model is uploaded automatically. Crash reports and the log stay on
this computer (anonymous usage, if you agreed to share it, never carries them: see
[privacy.md](privacy.md)). The only way they leave it is **Help › Report a Bug** (or **Report It…** after a crash), which opens a GitHub issue
form in your browser; you read it there and submit it, or not.

**Help › Crash Reports…** lists every saved report (when, what, how many times), each with
**Report…**, **Copy Path** and **Delete**, and **Delete All** (`diagnostics.crash_reports`,
`diagnostics.crash_reports.delete`).

What a bug report carries:

| Always (you can edit it) | Only with **Share my session** | Only with **Include the model file** as well |
|---|---|---|
| the version and OS, what happened (the last error by default), for a crash its message, place and backtrace top; your home folder is shown as `~` | the commands that built the model (`fieldcraft-session.json`, replayable with `fieldcraft-cli run`); they can name your files | the model (`model.fieldcraft`), which you attach yourself |

The form shows exactly what will be filled in before you click **Open GitHub Issue**. GitHub
refuses long new-issue links, so a session that does not fit goes out as a file instead: it is
saved to `bug-reports/<time>/` in the data folder below and copied to the clipboard, and the issue
text says to attach or paste it. The model is always a file you attach.

## Where things are

The app's data folder holds `recovery/`, `crashes/`, `logs/` and `bug-reports/`:

| System | Folder |
|---|---|
| macOS | `~/Library/Application Support/FieldCraft` |
| Linux | `~/.local/share/fieldcraft` (or `$XDG_DATA_HOME/fieldcraft`) |
| Windows | `%LOCALAPPDATA%\FieldCraft` |

Runs driven by automation (`--control`, or `CI` set) use `fieldcraft-automation/` in the
temporary folder instead, so test runs never reach your folders. `FIELDCRAFT_CRASH_DIR` moves the
crash folder, `FIELDCRAFT_RECOVERY_DIR` the recovery folder. The browser build keeps no crash
folder and logs to the browser console.

## Solver diagnostics

A solve that fails, or stops before the end of its step, leaves a small record in the crash
folder (`kind: "solver_error"`, `format: "fieldcraft-diagnostic/1"`): when it does not converge,
its system is singular, it goes unstable or diverges, an element fails, or its device (GPU)
fails; and when a steady flow reaches its time cap without settling or an explicit step stops at
its increment cap. A refusal of the request leaves none: a work cap such as `max_dofs` or
`max_cells`, a parameter that is wrong, a feature not supported yet. Neither does a cancelled
solve, nor a transient that reaches its end time.

A record holds the command and its parameters (up to 16 KiB), the model's size, the step and the
progress reached, the backend asked for and the one that ran, and, for a GPU solve, what its
setup saw: the adapter, backend, driver, the limits that matter (storage buffers per stage,
buffer sizes) and the devices and instances the process has opened (the fluid solver's own;
the viewport's and the driver's are not counted). What was not observed stays `null`. Home
folders are shown as `~`. No model, replay, log or backtrace is copied, and a record never
brings up the "FieldCraft closed unexpectedly" question.

The same failure again counts in its record; at most five new records are written per process.
Records share the crash folder's limits (20 reports, 100 MB, 90 days) and are the first to go
when it is over them, before any crash report. Writing a record never changes the solve's
outcome: its error stays the same error, and a stopped result gains the record's path as
`diagnostic_report`.

Help › Crash Reports, `diagnostics.crash_reports` and `fieldcraft-cli crashes` list them;
`diagnostics.crash_report` reads one. The app, `fieldcraft-cli run` and `fieldcraft-cli mcp`
record into the crash folder; the CLI prints the path of the record on stderr and the MCP reply
of a failed command ends with it. A `Session` embedded in another program records nothing
unless its host names a folder (`Session::solver_diagnostics`).

### Headless workers

In a container or on a cloud instance, the crash folder goes with the machine. Set
`FIELDCRAFT_CRASH_DIR` to a persistent folder before starting the CLI or the MCP server, and keep
that folder as a job artifact:

```sh
FIELDCRAFT_CRASH_DIR=/work/fieldcraft-reports fieldcraft-cli run model.json
```

With confined file access (`--automation-write-root`), the folder must be inside the write root,
or nothing is recorded.

## The log file

`logs/fieldcraft.log`: one line per record, `time LEVEL [pid] target: message`. Past 2 MB it
becomes `fieldcraft.log.1` and a new one starts. Warnings and errors also go to the standard
error stream.

## Crash reports

`crashes/<YYYYMMDD-HHMMSS>-<signature>/`:

- `report.json` (`"format": "fieldcraft-crash/1"`): `kind` (`command_panic`: a command hit an
  internal error, caught, the model left unchanged; `recovered_panic`: one caught inside a command;
  `panic`: an internal error outside commands; `unclean_exit`: the app ended without a clean exit
  — an abort, a stack overflow, a crash in a graphics driver, the process being killed outright
  (SIGKILL) — found at the next start), `message`, `location`, `thread`, `backtrace`, `app`, `system`, `command`,
  `recent_commands` (up to 32), `last_activity` (the last control request, such as
  `control ui.capture`), `startup_phase` (`{step, ms}`: the step a run that never drew its first
  frame died in), `native` (Windows: what the exception filter caught), `base`, `count`,
  `first_seen`, `last_seen`, `notes`; `system.environment` names the renderer, GPU and driver
- `model.fieldcraft`: the model the replay starts from
- `replay.json`: `{"commands": [{command, params}, …]}`
- `log.txt`: the crashed run's lines of the log
- `minidump.dmp` (Windows, a crash in native code): open it with WinDbg or Visual Studio

### Startup steps and native crashes

While the app gets up, the journal keeps the step it is in: the app's own (`window`, `app setup`,
`viewport pipelines (shader compile)`, `first frame`, and `running` once a frame was presented)
and egui-wgpu's in between (`graphics instance`, `window surface`, `adapter request`, `adapter
chosen`, when the renderer is recorded, `device creation`, `egui renderer (shader compile)`,
`surface configure`, `first frame (present)`). A run that dies before `running` is reported "during startup (device creation, 1.2 s
after launch)", so a graphics driver crash names its step.

On Windows an unhandled-exception filter (`apps/fieldcraft/src/native_crash.rs`) catches a crash
in native code (an access violation in a graphics driver, a stack overflow): before the process
goes it writes the exception code, the faulting module and the offset into it, and a minidump,
next to the journal. The next start puts them into the report: "crashed with
EXCEPTION_ACCESS_VIOLATION (0xC0000005) in igvk64.dll+0x64EA72 during startup (device
creation…)", `backtrace` the code, the module and offset, and the thread, `native` the raw fields,
and `minidump.dmp`. After a stack overflow there is no stack left for a minidump; the rest is
still written.

The next start after a crash asks "FieldCraft closed unexpectedly": **Report It…** opens Report a
Bug for that report.

### For agents

| From | List | One report | Replay |
|---|---|---|---|
| Command (MCP `execute`, control channel, scripts) | `diagnostics.crash_reports` | `diagnostics.crash_report {id}` | |
| CLI | `fieldcraft-cli crashes [--json]` | `fieldcraft-cli crashes show <id>` | `fieldcraft-cli crashes replay <id>` |

`crashes replay` opens the report's model, runs its commands (files they write can only go to a
scratch folder) and says whether the same internal error happens again. `crashes dir` prints the
folder, `crashes clear [<id>]` deletes reports. `diagnostics.bug_report` returns the draft (link,
fields, attachments) without opening anything; over the control channel `ui.start
diagnostics.bug_report` opens the form, `ui.prompt {"answer": "set", …}` sets its fields, and
`ui.inspect` shows the draft. Only a click opens the browser.

**Check that reporting works:** with `FIELDCRAFT_CRASH_TEST=1`, `diagnostics.crash_test {mode}`
crashes on purpose: `panic`, `recovered` or `abort` (reported at the next start).

A request to terminate is no crash: on Linux and macOS, SIGTERM (`kill PID`) and SIGINT (Ctrl-C)
close the window as **Close** does (unsaved work is autosaved for recovery, the
journal ends, no report). A second signal ends the app at once.

### How it stays cheap and small

The app keeps a journal of the window's top-level commands (`crashes/sessions/*.jsonl`, at most
512 KB, two small writes per command that changes the model, none for commands that only read),
restarted from each autosave and each model opened. The journal also records external result imports
and display changes. An autosave embeds the imported data. Before that snapshot, replay reads the
original result files. A clean exit deletes the journal. Its lock is held while the app runs, so an
unlocked journal belongs to an app that died. A crash is identified by its place and message,
numbers ignored: repeats only raise `count`. At most 5 new reports per run, 20 reports, 100 MB,
90 days. A diagnostic model copy is limited to 20 MiB. FieldCraft checks the size of an in-memory
result series before encoding it for a report. The separate user recovery autosave can hold an
imported series up to the usual 512 MiB limit. Tests, the CLI and the browser build never install
the hook.

Adapted from SolidCraft's crash reports (MIT OR Apache-2.0, see ATTRIBUTION.md).
