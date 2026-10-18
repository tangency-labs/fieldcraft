// The triage bot's logic, used by .github/workflows/triage.yml through actions/github-script.
// New issues and PRs get `area:*` and type labels from fixed rules; an issue nobody has given a
// priority yet carries `needs-triage`. Maintainers set labels with `/triage`. The rules only ever
// add labels, so a person's correction is never undone.
// Everything takes the github-script `github` client and `context`, so the tests can pass fakes.
'use strict';

const NEEDS_TRIAGE = 'needs-triage';
const NEEDS_TRIAGE_COLOR = 'ededed';
const NEEDS_TRIAGE_DESCRIPTION = 'no priority yet; a maintainer sorts it with /triage';
const MAX_TOKENS = 10;
const MAX_TEXT = 20000;

// Changed paths of a PR → area. First match wins per file.
const PATH_AREAS = [
  [/^crates\/mesh\//, 'area:mesh'],
  [/^crates\/lbm[^/]*\//, 'area:cfd'],
  [/^crates\/(solver|fem|explicit)\//, 'area:solver'],
  [/^(crates\/ui-egui|apps\/[^/]+)\//, 'area:ui'],
  [/^crates\/(io|decks)\//, 'area:io'],
  [/^(docs\/|[^/]+\.md$)/, 'area:docs'],
  [/^(\.github\/|xtask\/|packaging\/|scripts\/|tools\/)/, 'area:infra'],
];

// Command id prefix (the in-app reports quote it: "Internal error in mesh.tet.generate") → area.
const COMMAND_AREAS = {
  mesh: 'area:mesh',
  geometry: 'area:mesh',
  model: 'area:solver',
  interaction: 'area:solver',
  solve: 'area:solver',
  file: 'area:io',
  deck: 'area:io',
  results: 'area:ui',
  view: 'area:ui',
};

// Words in an issue's title → area. Conservative: a miss leaves the issue for a person.
const WORD_AREAS = [
  [/\b(mesh(es|ing|er)?|tet(rahedral|s)?|hex(ahedral|es)?|octree|jacobian|sliver|element quality|remesh)\b/i, 'area:mesh'],
  [/\b(lbm|lattice boltzmann|cfd|fluid|flow|turbulen(ce|t)|velocity field)\b/i, 'area:cfd'],
  [/\b(solver?|solves?|solving|explicit|implicit|static|frequency|modal|buckling|contact|plasticity|hyperelastic|convergence|newton|material|thermal)\b/i, 'area:solver'],
  [/\b(abaqus|nastran|ls-?dyna|deck|vtk|vtu|paraview|gmsh|step|stl|3mf|csv|import|export)\b/i, 'area:io'],
  [/\b(toolbar|dialog|viewport|panel|button|shortcut|menu|window|legend|plot|theme)\b/i, 'area:ui'],
  [/\b(docs?|documentation|readme|tutorial|typo)\b/i, 'area:docs'],
  [/\b(ci|release|packaging|flatpak|installer|homebrew|winget)\b/i, 'area:infra'],
];

const CRASH = /^\s*crash\b|\bpanic(ked)?\b|^\s*internal error\b|\bclosed unexpectedly\b/i;
const FEATURE = /^\s*\[?feature( request)?\]?\b/i;
// The issue forms' own headings (.github/ISSUE_TEMPLATE): their labels are set by GitHub when the
// form is used, so these only matter for the backfill and for issues edited into shape.
const BUG_FORM = /^###\s*observed failure\s*$/im;
const FEATURE_FORM = /^###\s*what you need and why\s*$/im;

function repoArgs(context) {
  return { owner: context.repo.owner, repo: context.repo.repo };
}

// The text under a `## <name>` heading (`###` in issue forms), up to the next heading.
function section(body, name) {
  const lines = body.split('\n');
  const start = lines.findIndex((l) => new RegExp(`^#{2,3}\\s*${name}\\s*$`, 'i').test(l.trim()));
  if (start < 0) return '';
  const rest = lines.slice(start + 1);
  const end = rest.findIndex((l) => /^#{2,3}\s/.test(l));
  return (end < 0 ? rest : rest.slice(0, end)).join('\n');
}

function labelNames(item) {
  return (item.labels || []).map((l) => (typeof l === 'string' ? l : l.name));
}

function isPriority(name) {
  return /^priority\s*:/.test(name);
}

function hasPriority(names) {
  return names.some(isPriority);
}

// Labels for an issue from its title and body.
function classifyIssue(title, body) {
  const t = String(title || '').slice(0, 500);
  const b = String(body || '').slice(0, MAX_TEXT);
  const labels = new Set();
  if (CRASH.test(t) || BUG_FORM.test(b)) {
    labels.add('bug');
  } else if (FEATURE.test(t) || FEATURE_FORM.test(b)) {
    labels.add('enhancement');
  }
  // The title names the feature; a command id in the title (the in-app report's) is the fallback.
  // The body is not read: a report's session quotes every earlier command too.
  let area = null;
  for (const [re, a] of WORD_AREAS) {
    if (re.test(t)) {
      area = a;
      break;
    }
  }
  if (!area) {
    const m = /\b([a-z_]+)\.[a-z_]+(\.[a-z_]+)*\b/.exec(t);
    area = (m && COMMAND_AREAS[m[1]]) || null;
  }
  if (area) labels.add(area);
  return [...labels];
}

// Labels for a PR from its changed paths.
function classifyPaths(paths) {
  const labels = new Set();
  for (const p of paths) {
    for (const [re, area] of PATH_AREAS) {
      if (re.test(p)) {
        labels.add(area);
        break;
      }
    }
  }
  return [...labels].sort();
}

async function ensureLabel(github, context, name, color, description) {
  try {
    await github.rest.issues.getLabel({ ...repoArgs(context), name });
  } catch (e) {
    if (e.status !== 404) throw e;
    try {
      await github.rest.issues.createLabel({ ...repoArgs(context), name, color, description });
    } catch (e2) {
      if (e2.status !== 422) throw e2; // created meanwhile
    }
  }
}

async function pullPaths(github, context, number) {
  const files = await github.paginate(github.rest.pulls.listFiles, {
    ...repoArgs(context),
    pull_number: number,
    per_page: 100,
  });
  return files.map((f) => f.filename);
}

// The labels the rules would add to an issue or PR that it doesn't already carry.
async function proposal(github, context, item) {
  const current = labelNames(item);
  let wanted;
  if (item.pull_request) {
    wanted = classifyPaths(await pullPaths(github, context, item.number));
  } else {
    wanted = classifyIssue(item.title, item.body);
    if (!hasPriority(current) && !hasPriority(wanted)) wanted.push(NEEDS_TRIAGE);
  }
  return wanted.filter((l) => !current.includes(l));
}

async function apply(github, context, number, labels) {
  if (!labels.length) return;
  if (labels.includes(NEEDS_TRIAGE)) {
    await ensureLabel(github, context, NEEDS_TRIAGE, NEEDS_TRIAGE_COLOR, NEEDS_TRIAGE_DESCRIPTION);
  }
  await github.rest.issues.addLabels({ ...repoArgs(context), issue_number: number, labels });
}

// issues.opened / pull_request_target.opened.
async function handleOpened({ github, context }) {
  const issue = context.payload.issue;
  const pull = context.payload.pull_request;
  const item = issue || (pull && { ...pull, pull_request: {} });
  if (!item || item.state !== 'open') return [];
  const labels = await proposal(github, context, item);
  await apply(github, context, item.number, labels);
  return labels;
}

// issues.labeled: a priority settles triage.
async function handleLabeled({ github, context }) {
  const issue = context.payload.issue;
  const label = context.payload.label;
  if (!issue || !label || !isPriority(label.name)) return false;
  if (!labelNames(issue).includes(NEEDS_TRIAGE)) return false;
  try {
    await github.rest.issues.removeLabel({ ...repoArgs(context), issue_number: issue.number, name: NEEDS_TRIAGE });
  } catch (e) {
    if (e.status !== 404) throw e;
  }
  return true;
}

// `/triage area:mesh high S bug -enhancement`: tokens name existing labels, matched without
// regard to case, `-`/`_` for spaces or the space after a colon (`in-progress`, `size:S`); `high`
// is short for `priority: high` and `S`/`M`/`L` for `size: S`…, each replacing any other label of
// its kind; a leading `-` removes.
const ALIASES = { high: 'priority: high', s: 'size: S', m: 'size: M', l: 'size: L' };

function parseTriage(body) {
  const line = String(body || '').trimStart().split('\n')[0];
  const m = /^\/triage(\s.*)?$/.exec(line.trimEnd());
  if (!m) return null;
  const add = [];
  const remove = [];
  for (const raw of (m[1] || '').trim().split(/\s+/).filter(Boolean).slice(0, MAX_TOKENS)) {
    const minus = raw.startsWith('-');
    let name = minus ? raw.slice(1) : raw;
    name = ALIASES[name.toLowerCase()] || name;
    if (name) (minus ? remove : add).push(name);
  }
  return { add, remove };
}

function normalize(name) {
  return String(name).toLowerCase().replace(/[-_\s]+/g, ' ').replace(/\s*:\s*/g, ':').trim();
}

// The label kind that replaces its siblings: `priority` or `size`.
function kind(name) {
  const m = /^(priority|size)\s*:/.exec(name);
  return m ? m[1] : null;
}

async function isMaintainer(github, context, login) {
  if (login === context.repo.owner) return true;
  try {
    const { data } = await github.rest.repos.getCollaboratorPermissionLevel({ ...repoArgs(context), username: login });
    return ['admin', 'write', 'maintain'].includes(data.permission);
  } catch (e) {
    return false;
  }
}

async function react(github, context, content) {
  await github.rest.reactions.createForIssueComment({
    ...repoArgs(context),
    comment_id: context.payload.comment.id,
    content,
  });
}

// issue_comment.created.
async function handleComment({ github, context }) {
  const { issue, comment } = context.payload;
  if (!issue || !comment || comment.user.type === 'Bot' || comment.user.login.endsWith('[bot]')) return 'skip';
  const cmd = parseTriage(comment.body);
  if (!cmd) return 'skip';
  if (!(await isMaintainer(github, context, comment.user.login))) {
    await react(github, context, 'eyes');
    return 'denied';
  }
  const known = new Map(
    (await github.paginate(github.rest.issues.listLabelsForRepo, { ...repoArgs(context), per_page: 100 })).map((l) => [
      normalize(l.name),
      l.name,
    ]),
  );
  const resolve = (n) => known.get(normalize(n));
  const unknown = [...cmd.add, ...cmd.remove].filter((n) => !resolve(n));
  if (unknown.length) {
    await react(github, context, 'confused');
    await github.rest.issues.createComment({
      ...repoArgs(context),
      issue_number: issue.number,
      body: `No such label: ${unknown.map((n) => '`' + n + '`').join(', ')}. Nothing was changed.`,
    });
    return 'unknown';
  }
  const addNames = cmd.add.map(resolve);
  const current = labelNames(issue);
  const remove = new Set(cmd.remove.map(resolve));
  for (const k of new Set(addNames.map(kind).filter(Boolean))) {
    for (const n of current) if (kind(n) === k && !addNames.includes(n)) remove.add(n);
  }
  if (addNames.length || remove.size) remove.add(NEEDS_TRIAGE);
  for (const name of remove) {
    if (!current.includes(name)) continue;
    await github.rest.issues.removeLabel({ ...repoArgs(context), issue_number: issue.number, name });
  }
  const add = addNames.filter((n) => !current.includes(n));
  if (add.length) await github.rest.issues.addLabels({ ...repoArgs(context), issue_number: issue.number, labels: add });
  await react(github, context, '+1');
  return 'triage';
}

// Workflow dispatch: the rules over every open issue and PR with no labels at all. With `dryRun`
// it changes nothing and returns what it would do.
async function backfill({ github, context, dryRun = true }) {
  const items = await github.paginate(github.rest.issues.listForRepo, {
    ...repoArgs(context),
    state: 'open',
    per_page: 100,
  });
  const plan = [];
  for (const item of items) {
    if (labelNames(item).length || item.locked) continue;
    const labels = await proposal(github, context, item);
    plan.push({ number: item.number, pull: Boolean(item.pull_request), title: item.title, labels });
    if (!dryRun) await apply(github, context, item.number, labels);
  }
  return plan;
}

module.exports = {
  NEEDS_TRIAGE,
  classifyIssue,
  classifyPaths,
  ensureLabel,
  isMaintainer,
  labelNames,
  section,
  parseTriage,
  handleOpened,
  handleLabeled,
  handleComment,
  backfill,
};
