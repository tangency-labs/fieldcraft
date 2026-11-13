# Privacy and anonymous usage

FieldCraft can share **anonymous usage** with its developers: which tools are used, where people
get stuck and how far they get from a design to results. It is **off until you agree**. The first
time the app starts it asks once:

> Help improve FieldCraft by sharing anonymous usage: which tools you use, nothing from your models.

**Share Anonymous Usage** turns it on, **Don't Share** (or closing the question) leaves it off.
Change it any time in **Help › Settings… › Privacy**, which also has **Include scripts and AI
agents (MCP)**: commands run by scripts, the control channel and AI agents are counted apart from
yours (`from`), and this leaves them out.

Checking for updates is separate from usage sharing: at start-up, at most once a day (and from
Help › Check for Updates), the desktop app asks GitHub for the latest release, which sends the
request itself (FieldCraft's version in the `User-Agent`; GitHub sees the IP address) and nothing
else. Turn it off in **Help › Settings… › Updates** (or with the first run's question): then no
request is made ([updates.md](updates.md)).

Crash reports and the log never leave your computer on their own; Help › Report a Bug sends only
what you submit yourself ([crash-reports.md](crash-reports.md)).

## Never sent

File names or paths, geometry, meshes, model contents, parameter values (sizes, materials,
loads…), names you give things, error messages, or anything else you typed. The events below are
built from FieldCraft's own lists (its command ids, sample names and requirement names) and
durations, so they cannot carry any of it; a test runs a scripted session and fails if a payload
holds anything else.

## When nothing is sent

- until you agree, and after you turn it off (what was waiting is dropped);
- in builds without an analytics key: everything built from source, development and fork builds.
  Only the release builds have one;
- with `FIELDCRAFT_TELEMETRY=0` (or `false`, `off`, `no`) in the environment;
- in runs driven by automation (`--control`) and in CI (`CI` set), unless
  `FIELDCRAFT_TELEMETRY=1` asks for it;
- `fieldcraft-cli` and the MCP server never send anything themselves; the app counts the
  commands they send it.

## Where it goes

[PostHog](https://posthog.com) Cloud, EU region (`https://eu.i.posthog.com`), a batch of events
every 30 seconds and when the window closes, on a background thread (the web build uses `fetch`).
A batch that cannot be sent (offline, a slow network) is dropped; at most 1000 events wait. Each
event goes with `$process_person_profile: false` (no person profile is built) and
`$geoip_disable: true` (no location is derived from the connection's IP address).
`FIELDCRAFT_TELEMETRY_HOST` sends to another host instead (to check what is sent against a local
listener).

## What is sent

Every event has these properties:

| Property | Value |
|---|---|
| `distinct_id` | a random id made for this install when you agree (a UUID), kept in the preferences; nothing about the computer or you |
| `app_version` | FieldCraft's version, `0.5.8` |
| `platform` | `desktop` or `web` |
| `os` | `macos`, `linux`, `windows`, or `web` |
| `arch` | `aarch64`, `x86_64`, or `wasm32` |
| `$process_person_profile` | `false` |
| `$geoip_disable` | `true` |

The events:

| Event | Properties | When |
|---|---|---|
| `app_started` | — | the app starts with sharing on, or you turn it on |
| `app_closed` | `duration_s`: seconds since `app_started` | the window closes (the page is closed or left) |
| `command_run` | `command`: the command's id (`mesh.tet.generate`); `kind`: `model` or `session`; `duration_ms`; `ok`: whether it succeeded; `error` (when it failed): the class of the error, one of `unknown_command`, `params`, `model`, `geometry`, `mesh`, `file`, `deck`, `solver`, `explicit`, `fluid`, `unavailable`, `internal`; `from`: `ui` (a click, a form, a menu), `palette` (the command search), `script` (a script or the control channel) or `mcp` (an AI agent) | a command ends; a command run by another (a script's steps inside it) is not counted on its own |
| `milestone` | `name`: `design_opened` (a design, mesh or model opened), `meshed` (a mesher made elements), `solved` (a solve finished), `results_viewed` (results shown) | the first time in a run |
| `sample_opened` | `name`: the sample's name (`cantilever`) | a sample opens |
| `form_cancelled` | `command`: the form's command id | a command's form closes without running it |
| `unavailable_clicked` | `command`: the command id; `missing`: what the model lacks, by requirement name (`material`, `step`, …) | a command is clicked before the model has what it needs |

A `command_run` as sent:

```json
{"event": "command_run", "distinct_id": "5f0c…-…", "properties": {"command": "solve.run", "kind": "model", "duration_ms": 412, "ok": true, "from": "ui", "app_version": "0.5.8", "platform": "desktop", "os": "macos", "arch": "aarch64", "$process_person_profile": false, "$geoip_disable": true}}
```

## For developers

The events and the policy are `fieldcraft_engine::usage`; the consent question and the settings
are `crates/ui-egui/src/privacy.rs`; the desktop sender is `apps/fieldcraft/src/usage.rs`, the
web one `apps/fieldcraft-web/src/usage.rs`. The PostHog project key is the build-time variable
`FIELDCRAFT_POSTHOG_KEY`, which the release workflow and the Pages build of `main` set from the `release` environment secret
`POSTHOG_KEY`; it is an ingest-only key. A new event or property is added here first.
