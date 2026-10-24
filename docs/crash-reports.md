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

## Ordinary solver diagnostics

A returned solver error for a work/memory cap, failed convergence, instability or backend
failure writes a small local record with `kind: "solver_error"` and
`format: "fieldcraft-diagnostic/1"`. A configured steady solve that reaches its cap without
settling is recorded too; reaching the requested end of a physical transient is success.
Malformed command parameters and deliberate cancellation do not create solver records.
Warning-only truncated arc-length paths currently return a successful summary and do not
create these records; they remain a separate diagnostic coverage gap.

These records use the crash folder and its retention limits, but are not crash reports: no
model, replay, log or backtrace is copied, and they never produce a startup "app crashed"
question. Help › Crash Reports and `diagnostics.crash_reports` list them; a single record can
be read with `diagnostics.crash_report`. Nothing is uploaded or attached publicly.

A record names the command and parameters, model size, and the progress actually observed.
Backend setup can carry the selected adapter, driver, relevant adapter/device limits and
scoped process counters. Selection, a device request, device creation and cached-device reuse
are distinct phases. Refusal before allocation does not invent a device, time or increment;
unobserved data and unavailable checkpoint paths remain null. Counter scopes identify the
managed contexts counted and external/driver-internal contexts that are not counted.
Successful CPU fluid summaries identify their backend. A CPU fluid failure has no backend
identity observation today, so its actual backend remains null even if progress was reported;
the requested choice alone is not evidence that allocation or execution started.

Parameters are limited to 16 KiB, setup context to 64 KiB, with an explanation when omitted.
Home paths are redacted as values before JSON serialization. Repeats update one record;
at most five new solver records are written per process, independently of the crash-report
budget. The shared folder keeps at most 20 reports, 100 MB and 90 days. Failure to persist a
diagnostic never replaces the original solve result or error.

### Persistent reports on headless workers

CLI scripts and headless MCP sessions use the temporary automation folder. In containers or
cloud workers, set `FIELDCRAFT_CRASH_DIR` to a persistent mounted folder **before starting the
process**, then retain that folder as a job artifact. For example, in PowerShell:

```powershell
$env:FIELDCRAFT_CRASH_DIR = 'C:/work/fieldcraft-reports'
fieldcraft-cli run model.json
```

Every written solver diagnostic prints its actual path on stderr. A command failure carries
the same path in its error text, including the MCP reply; a returned stopped result adds
`diagnostic_report`. Workers can collect the record without searching a throwaway filesystem.
When MCP write roots are confined, the report directory must be inside the granted write
root (and the read root to list/read it). An inaccessible or unwritable directory leaves the
original solve outcome unchanged. Keep `FIELDCRAFT_CRASH_DIR` set when subsequently using
`fieldcraft-cli crashes` to inspect those records.

Embedded `Session`s default to no automatic persistence; a host may set
`Session::solver_diagnostics` or supply an explicit `FIELDCRAFT_CRASH_DIR`. The desktop's
installed reporter also enables it. Native compute reporting supports CPU and Vulkan; historical
GL, Metal and DX12 failures do not advertise those native backends as supported. The browser
still has no native report folder: shared progress/context types compile there, but filesystem
persistence and native adapter counters are unavailable. Existing browser renderer/bug-report
details remain separate from the solver's observed backend context.

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
restarted from each autosave and each model opened. A clean exit deletes it; it is locked while its
app runs, so a journal nobody holds belongs to an app that died. A crash is identified by its
place and message, numbers ignored: repeats only raise `count`. At most 5 new reports per run,
20 reports, 100 MB, 90 days; a model over 20 MB is not copied. Tests, the CLI and the browser build
never install the hook.

Adapted from SolidCraft's crash reports (MIT OR Apache-2.0, see ATTRIBUTION.md).
