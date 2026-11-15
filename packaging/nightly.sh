#!/usr/bin/env bash
# Nightly builds (release.yml's `schedule`; docs/releasing.md, "Nightly builds"). Adapted from
# SolidCraft's packaging/nightly.sh (MIT OR Apache-2.0; see ATTRIBUTION.md).
#
# Usage:
#   packaging/nightly.sh decide          which commit tonight's nightly builds, if any
#   packaging/nightly.sh version BASE    the nightly version of BASE (MAJOR.MINOR.PATCH[-…])
#   packaging/nightly.sh prune KEEP      delete every published nightly but the newest KEEP
#
# Environment:
#   GITHUB_REPOSITORY, GH_TOKEN          the repository (gh)
#   GITHUB_EVENT_NAME, GITHUB_SHA        the run's event and commit (decide)
#   GATE_RUN_ID, GATE_SHA, GATE_URL, GATE_CONCLUSION
#                                        the CI run whose completion started a nightly (decide)
#   GITHUB_RUN_NUMBER, GITHUB_RUN_ATTEMPT  the run (version)
#   GITHUB_OUTPUT, GITHUB_STEP_SUMMARY   where decide writes its outputs and says why it skips
#
# `decide` writes build=true|false, sha=<the commit to build> and, for a nightly, gate=<the CI
# run that verified it>. Any other event builds the run's own commit. A nightly starts when CI's
# nightly run on develop completes (ci.yml, 03:23 UTC: the full gate, the heavy tests and the
# oracle) and builds that run's commit, so it reuses that run instead of verifying again; it skips
# the night when that run failed, or when that commit is the newest published nightly's.
set -euo pipefail

# The jobs of ci.yml's nightly run this reuses (their `name:`s), and its plan job, which only
# the weekly macOS-oracle run skips.
GATE_JOBS=("cargo xtask ci" "cargo xtask test --heavy" "cargo xtask oracle")
PLAN_JOB="Plan (affected crates)"

# A nightly's tag: v<MAJOR.MINOR.PATCH>-nightly.<YYYYMMDD>.<run>.<attempt> (`version`).
NIGHTLY_TAG='^v(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)-nightly\.[0-9]{8}\.[0-9]+\.[0-9]+$'

summary() {
  echo "$*"
  if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then echo "$*" >>"$GITHUB_STEP_SUMMARY"; fi
}

# The published nightly releases' tags, newest first. A tag counts only when it has the nightly
# form and its release is a published pre-release, so no other release is ever listed.
nightlies() {
  gh api --paginate --slurp "repos/$GITHUB_REPOSITORY/releases?per_page=100" | jq -r --arg re "$NIGHTLY_TAG" '
    [.[][] | select(.draft == false and .prerelease == true and .published_at != null and (.tag_name | test($re)))]
    | sort_by(.published_at) | reverse | .[].tag_name'
}

decide() {
  if [ "$GITHUB_EVENT_NAME" != workflow_run ]; then
    { echo build=true; echo "sha=$GITHUB_SHA"; } >>"$GITHUB_OUTPUT"
    return
  fi
  # The CI run that started this one: a scheduled run on develop (release.yml's version job runs
  # for no other), the daily nightly or the weekly macOS oracle.
  local sha="$GATE_SHA" url="$GATE_URL" jobs
  # Its jobs' conclusions, `<conclusion><TAB><name>` per line.
  jobs="$(gh api --paginate "repos/$GITHUB_REPOSITORY/actions/runs/$GATE_RUN_ID/jobs?per_page=100" \
    --jq '.jobs[] | "\(.conclusion)\t\(.name)"')"
  conclusion() { awk -F '\t' -v name="$1" '$2 == name { print $1 }' <<<"$jobs"; }
  if [ "$(conclusion "$PLAN_JOB")" = skipped ]; then
    echo build=false >>"$GITHUB_OUTPUT"
    summary "### No nightly: CI's run $url is the weekly macOS oracle, not the nightly"
    return
  fi
  local job result green="$GATE_CONCLUSION"
  for job in "${GATE_JOBS[@]}"; do
    result="$(conclusion "$job")"
    if [ "$result" != success ]; then green="$job: ${result:-missing}"; fi
  done
  if [ "$green" != success ]; then
    echo build=false >>"$GITHUB_OUTPUT"
    summary "### Nightly skipped: develop is red (CI's nightly run on \`${sha::9}\` did not pass: $green; $url)"
    return
  fi
  local last last_sha
  last="$(nightlies | head -n 1)"
  if [ -n "$last" ]; then
    # The tag's own commit (a release's target_commitish can name a moving branch).
    last_sha="$(gh api "repos/$GITHUB_REPOSITORY/commits/$last" --jq .sha)"
    if [ "$last_sha" = "$sha" ]; then
      echo build=false >>"$GITHUB_OUTPUT"
      summary "### Nightly skipped: develop has not moved since $last (\`${sha::9}\`)"
      return
    fi
  fi
  { echo build=true; echo "sha=$sha"; echo "gate=$url"; } >>"$GITHUB_OUTPUT"
  summary "Nightly of develop at \`${sha::9}\`, verified by CI's nightly run $url"
}

version() {
  local base="${1%%-*}"
  printf '%s-nightly.%s.%s.%s\n' "$base" "$(date -u +%Y%m%d)" "$GITHUB_RUN_NUMBER" "$GITHUB_RUN_ATTEMPT"
}

prune() {
  local keep="$1" tag
  case "$keep" in '' | *[!0-9]*) echo "error: KEEP must be a number, not '$keep'" >&2; exit 2 ;; esac
  nightlies | tail -n +"$((keep + 1))" | while read -r tag; do
    [ -n "$tag" ] || continue
    gh release delete "$tag" --repo "$GITHUB_REPOSITORY" --cleanup-tag --yes
    summary "Deleted the older nightly $tag"
  done
}

case "${1:-}" in
  decide) decide ;;
  version) version "${2:?usage: nightly.sh version BASE}" ;;
  prune) prune "${2:?usage: nightly.sh prune KEEP}" ;;
  *) echo "usage: packaging/nightly.sh decide | version BASE | prune KEEP" >&2; exit 2 ;;
esac
