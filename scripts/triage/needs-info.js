// `needs-info`: a short hand-written bug report without reproduction steps is asked for them once.
// The author's reply (or an edit that adds them) clears the label; the daily sweep closes what
// stays unanswered for INFO_DAYS. Reports filed from the app (Help › Report a Bug,
// crates/engine/src/bug_report.rs) carry their session, crash report and version, so they are
// never asked. Used by .github/workflows/triage.yml.
'use strict';

const { ensureLabel, isMaintainer, labelNames, section } = require('./triage.js');

const LABEL = 'needs-info';
const LABEL_COLOR = 'd4c5f9';
const LABEL_DESCRIPTION = 'waiting on the reporter for reproduction steps';
const INFO_DAYS = 14;
const DAY = 24 * 3600 * 1000;
const MAX_TEXT = 20000;

function repoArgs(context) {
  return { owner: context.repo.owner, repo: context.repo.repo };
}

// What only the app's own reports carry (crates/engine/src/bug_report.rs).
const FROM_APP = /Replay the attached session|^Crash report `|^Backtrace \(top\):|^Session \(session\.json\):/m;

// Requests, questions and plans, which aren't reproduced.
const NOT_A_BUG = /^\s*\[(feat|feature|change|enhancement)|feature request|suggestion|\?\s*$|^\s*(tracker|refactor|support|add|allow|how|why|request|missing)\b/i;
// A report this long describes what it saw; only short ones are asked. Counted without the
// form's headings and its `_No response_` placeholders.
const SHORT = 300;
const ATTACHMENT = /user-attachments\/|\.(fieldcraft|solidcraft|inp|bdf|k|vtu|step|stp|3mf|stl|json)\b/i;
const COMMAND = /["`][a-z_]+\.[a-z_][a-z_.]*["`]|\bfieldcraft(-cli)?\s+(--sample|run|snapshot)\b/;
const VERSION = /\bv?\d+\.\d+\.\d+\b|\b[0-9a-f]{7,40}\b/i;

function prose(body) {
  return body
    .split('\n')
    .filter((l) => !/^#{2,3}\s/.test(l) && l.trim() !== '_No response_')
    .join(' ')
    .replace(/\s+/g, ' ')
    .trim();
}

// What a bug report lacks: [] when it can be reproduced.
function missing(title, body) {
  const t = String(title || '');
  const b = String(body || '').slice(0, MAX_TEXT).replace(/<!--[\s\S]*?-->/g, '');
  if (NOT_A_BUG.test(t) || FROM_APP.test(b) || prose(b).length >= SHORT) return [];
  const steps = `${section(b, 'Reproduction')}\n${section(b, 'Steps to reproduce')}`.replace(/_No response_/g, '');
  const reproducible =
    steps.replace(/\s+/g, '').length >= 30 || ATTACHMENT.test(b) || COMMAND.test(b) || /^\s*1[.)]\s/m.test(b);
  if (reproducible) return [];
  const asks = ['steps'];
  if (!VERSION.test(`${t}\n${b}`)) asks.push('version');
  return asks;
}

function askBody(asks) {
  const lines = ['Thanks for the report! To reproduce it we need a little more:'];
  lines.push(
    '- **Steps to reproduce**: numbered steps from a sample (`fieldcraft --sample cantilever`) or a script ' +
      '(`fieldcraft-cli run script.json`), and the `.fieldcraft` or `.solidcraft` file that shows it if it is your own model. ' +
      '**Help › Report a Bug** in the app fills this in for you.',
  );
  if (asks.includes('version')) lines.push('- **FieldCraft version** (or commit) and your OS.');
  lines.push(
    '',
    `Reply here or edit the issue and it goes back to the triage queue. If we don't hear back within ${INFO_DAYS} days it is closed; reopen it any time.`,
  );
  return lines.join('\n');
}

// What to ask the reporter of `issue`, or [] (maintainers' and bots' issues are never asked).
async function asksFor(github, context, issue, cache) {
  if (issue.state !== 'open' || issue.pull_request || issue.locked) return [];
  const names = labelNames(issue);
  if (names.includes(LABEL) || names.includes('enhancement')) return []; // the feature form labels its issues
  const asks = missing(issue.title, issue.body);
  if (!asks.length || !issue.user || issue.user.type === 'Bot') return [];
  const login = issue.user.login;
  if (!cache.has(login)) cache.set(login, await isMaintainer(github, context, login));
  return cache.get(login) ? [] : asks;
}

async function ask(github, context, issue, asks) {
  await ensureLabel(github, context, LABEL, LABEL_COLOR, LABEL_DESCRIPTION);
  await github.rest.issues.addLabels({ ...repoArgs(context), issue_number: issue.number, labels: [LABEL] });
  await github.rest.issues.createComment({ ...repoArgs(context), issue_number: issue.number, body: askBody(asks) });
}

// issues.opened.
async function handleOpened({ github, context }) {
  const issue = context.payload.issue;
  if (!issue) return false;
  const asks = await asksFor(github, context, issue, new Map());
  if (!asks.length) return false;
  await ask(github, context, issue, asks);
  return true;
}

// Workflow dispatch: the same ask over the open issues filed before the bot. With `dryRun` it
// changes nothing and returns what it would ask.
async function backfill({ github, context, dryRun = true }) {
  const issues = await github.paginate(github.rest.issues.listForRepo, {
    ...repoArgs(context),
    state: 'open',
    per_page: 100,
  });
  const cache = new Map();
  const plan = [];
  for (const issue of issues) {
    const asks = await asksFor(github, context, issue, cache);
    if (!asks.length) continue;
    plan.push({ number: issue.number, title: issue.title, asks });
    if (!dryRun) await ask(github, context, issue, asks);
  }
  return plan;
}

async function clear(github, context, number) {
  try {
    await github.rest.issues.removeLabel({ ...repoArgs(context), issue_number: number, name: LABEL });
  } catch (e) {
    if (e.status !== 404) throw e;
  }
}

// issue_comment.created / issues.edited: the reporter answered.
async function handleReply({ github, context }) {
  const { issue, comment } = context.payload;
  if (!issue || !labelNames(issue).includes(LABEL)) return false;
  const author = issue.user && issue.user.login;
  const actor = comment ? comment.user.login : context.payload.sender && context.payload.sender.login;
  if (!author || actor !== author) return false;
  if (!comment && missing(issue.title, issue.body).length) return false; // an edit that still lacks steps
  await clear(github, context, issue.number);
  return true;
}

// When the label was last added, from the timeline; the issue's update time if it can't be found.
async function labeledAt(github, context, issue) {
  const events = await github.paginate(github.rest.issues.listEventsForTimeline, {
    ...repoArgs(context),
    issue_number: issue.number,
    per_page: 100,
  });
  let at = null;
  for (const ev of events) if (ev.event === 'labeled' && ev.label && ev.label.name === LABEL) at = ev.created_at;
  return new Date(at || issue.updated_at).getTime();
}

// The daily sweep: close what has waited INFO_DAYS for its reporter.
async function sweep({ github, context, now = Date.now() }) {
  const issues = await github.paginate(github.rest.issues.listForRepo, {
    ...repoArgs(context),
    state: 'open',
    labels: LABEL,
    per_page: 100,
  });
  const closed = [];
  for (const issue of issues) {
    if (issue.pull_request || issue.locked) continue;
    if (now - (await labeledAt(github, context, issue)) < INFO_DAYS * DAY) continue;
    await github.rest.issues.createComment({
      ...repoArgs(context),
      issue_number: issue.number,
      body: `Closing for now: we couldn't reproduce this without more detail. Reopen it with steps to reproduce and we'll pick it up.`,
    });
    await github.rest.issues.update({
      ...repoArgs(context),
      issue_number: issue.number,
      state: 'closed',
      state_reason: 'not_planned',
    });
    closed.push(issue.number);
  }
  return closed;
}

module.exports = { LABEL, missing, handleOpened, handleReply, sweep, backfill };
