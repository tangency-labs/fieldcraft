// The /claim bot's logic, used by .github/workflows/claim.yml through actions/github-script.
// Everything takes the github-script `github` client and `context`, so the tests can pass fakes.
'use strict';

// CONTRIBUTING.md, "Before you start: claim an issue": a claimed issue is assigned and labelled
// `claimed`; `in progress` marks the core team's own work and can't be claimed; contributors pick
// issues labelled `help wanted` or `good first issue`; no visible progress for 14 days frees it.
const LABEL = 'claimed';
const LABEL_COLOR = '0052cc';
const LABEL_DESCRIPTION = 'Assigned to a contributor; see CONTRIBUTING.md';
const CORE = 'in progress';
const OPEN = ['help wanted', 'good first issue'];
const CLAIM_DAYS = 14;
const DAY = 24 * 3600 * 1000;

function parseCommand(body) {
  const text = String(body || '').trimStart();
  if (/^\/claim(\s|$)/.test(text)) return 'claim';
  if (/^\/unclaim(\s|$)/.test(text)) return 'unclaim';
  return null;
}

function fmt(date) {
  return new Date(date).toISOString().slice(0, 16).replace('T', ' ') + ' UTC';
}

function repoArgs(context) {
  return { owner: context.repo.owner, repo: context.repo.repo };
}

async function isMaintainer(github, context, login, cache) {
  if (login === context.repo.owner) return true;
  if (cache && cache.has(login)) return cache.get(login);
  let result = false;
  try {
    const { data } = await github.rest.repos.getCollaboratorPermissionLevel({
      ...repoArgs(context),
      username: login,
    });
    result =
      ['admin', 'write', 'maintain'].includes(data.permission) ||
      ['admin', 'maintain', 'write'].includes(data.role_name);
  } catch (e) {
    result = false;
  }
  if (cache) cache.set(login, result);
  return result;
}

async function ensureLabel(github, context) {
  try {
    await github.rest.issues.getLabel({ ...repoArgs(context), name: LABEL });
  } catch (e) {
    if (e.status !== 404) throw e;
    try {
      await github.rest.issues.createLabel({
        ...repoArgs(context),
        name: LABEL,
        color: LABEL_COLOR,
        description: LABEL_DESCRIPTION,
      });
    } catch (e2) {
      if (e2.status !== 422) throw e2; // created meanwhile
    }
  }
}

// When the label was last added, from the timeline; the issue's update time if it can't be found.
async function labeledAt(github, context, issue, timeline) {
  const events = timeline || (await getTimeline(github, context, issue.number));
  let at = null;
  for (const ev of events) {
    if (ev.event === 'labeled' && ev.label && ev.label.name === LABEL) at = ev.created_at;
  }
  return new Date(at || issue.updated_at).getTime();
}

function getTimeline(github, context, number) {
  return github.paginate(github.rest.issues.listEventsForTimeline, {
    ...repoArgs(context),
    issue_number: number,
    per_page: 100,
  });
}

// Open pull requests that mention the issue, from the timeline's cross-references.
function linkedPulls(timeline) {
  const seen = new Map();
  for (const ev of timeline) {
    const src = ev.event === 'cross-referenced' && ev.source && ev.source.issue;
    if (src && src.pull_request && src.state === 'open') seen.set(src.number, src);
  }
  return [...seen.values()];
}

async function lastPush(github, context, pr) {
  const { data } = await github.rest.pulls.get({ ...repoArgs(context), pull_number: pr.number });
  const { data: commit } = await github.rest.repos.getCommit({
    ...repoArgs(context),
    ref: data.head.sha,
  });
  return new Date(commit.commit.committer.date).getTime();
}

// Where a claimed issue stands: when it lapses, and whether it is free to release now.
async function claimStatus(github, context, issue, now) {
  const timeline = await getTimeline(github, context, issue.number);
  const since = await labeledAt(github, context, issue, timeline);
  const pulls = linkedPulls(timeline);
  if (pulls.length === 0) {
    const expires = since + CLAIM_DAYS * DAY;
    return { pulls, since, expires, release: now >= expires, reason: 'no PR' };
  }
  let newest = since;
  for (const pr of pulls) newest = Math.max(newest, await lastPush(github, context, pr));
  const expires = newest + CLAIM_DAYS * DAY;
  return { pulls, since, expires, release: now >= expires, reason: 'idle PR' };
}

async function react(github, context, content) {
  await github.rest.reactions.createForIssueComment({
    ...repoArgs(context),
    comment_id: context.payload.comment.id,
    content,
  });
}

async function reply(github, context, body) {
  await github.rest.issues.createComment({
    ...repoArgs(context),
    issue_number: context.payload.issue.number,
    body,
  });
}

function labelNames(issue) {
  return (issue.labels || []).map((l) => (typeof l === 'string' ? l : l.name));
}

function hasLabel(issue, name = LABEL) {
  return labelNames(issue).includes(name);
}

// Open to contributors: labelled `help wanted` or `good first issue`, and not the core team's.
function claimable(issue) {
  const names = labelNames(issue);
  return !names.includes(CORE) && OPEN.some((n) => names.includes(n));
}

async function release(github, context, issue, assignees) {
  if (assignees.length) {
    await github.rest.issues.removeAssignees({
      ...repoArgs(context),
      issue_number: issue.number,
      assignees,
    });
  }
  try {
    await github.rest.issues.removeLabel({
      ...repoArgs(context),
      issue_number: issue.number,
      name: LABEL,
    });
  } catch (e) {
    if (e.status !== 404) throw e;
  }
}

async function claim(github, context, now) {
  const issue = context.payload.issue;
  const user = context.payload.comment.user.login;
  const n = issue.number;
  const cache = new Map();
  const assignees = (issue.assignees || []).map((a) => a.login);

  if (assignees.includes(user)) {
    await react(github, context, '+1'); // already theirs
    return;
  }

  if (hasLabel(issue, CORE)) {
    await react(github, context, 'eyes');
    await reply(github, context, `The core team is working on this one (\`${CORE}\`), so it can't be claimed.`);
    return;
  }
  if (!claimable(issue) && !(await isMaintainer(github, context, user, cache))) {
    await react(github, context, 'eyes');
    await reply(
      github,
      context,
      `Only issues labelled \`help wanted\` or \`good first issue\` can be claimed, @${user}. ` +
        'Say here what you would like to do and a maintainer will open it up or assign it to you.',
    );
    return;
  }

  // Someone else is assigned: theirs unless the claim has lapsed (never for a maintainer).
  let held = assignees.length > 0 || hasLabel(issue);
  let status = null;
  if (held) {
    let maintained = false;
    for (const a of assignees) maintained = maintained || (await isMaintainer(github, context, a, cache));
    if (!maintained) {
      status = await claimStatus(github, context, issue, now);
      if (status.release) held = false;
    }
    if (held) {
      await react(github, context, 'eyes');
      const who = assignees.length ? assignees.map((a) => '@' + a).join(', ') + ' has' : 'Someone has';
      const when = maintained || !status ? '' : ` The claim lapses on ${fmt(status.expires)} unless there's progress.`;
      await reply(github, context, `${who} claimed this already.${when}`);
      return;
    }
  } else {
    // Not claimed, but someone's open PR may already cover it.
    const timeline = await getTimeline(github, context, n);
    const other = linkedPulls(timeline).find((p) => p.user && p.user.login !== user);
    if (other) {
      await react(github, context, 'eyes');
      await reply(github, context, `#${other.number} already covers this one.`);
      return;
    }
  }

  // Free (or lapsed): hand it to the commenter. Clear a lapsed claimant first.
  if (assignees.length) {
    await github.rest.issues.removeAssignees({ ...repoArgs(context), issue_number: n, assignees });
  }
  await ensureLabel(github, context);
  await github.rest.issues.addLabels({ ...repoArgs(context), issue_number: n, labels: [LABEL] });
  try {
    await github.rest.issues.addAssignees({
      ...repoArgs(context),
      issue_number: n,
      assignees: [user],
    });
  } catch (e) {
    // Assignment can be refused; the label and the reply naming them still stand.
  }
  await react(github, context, '+1');
  await reply(
    github,
    context,
    `Claimed by @${user}. Please open a draft PR early with \`Closes #${n}\`; the claim lapses after ${CLAIM_DAYS} days without visible progress. \`/unclaim\` releases it.`,
  );
}

async function unclaim(github, context) {
  const issue = context.payload.issue;
  const user = context.payload.comment.user.login;
  const assignees = (issue.assignees || []).map((a) => a.login);
  if (!hasLabel(issue) && assignees.length === 0) {
    await react(github, context, '+1'); // nothing to release
    return;
  }
  const allowed = assignees.includes(user) || (await isMaintainer(github, context, user));
  if (!allowed) {
    await react(github, context, 'eyes');
    await reply(github, context, `Only the person who claimed this, or a maintainer, can release it, @${user}.`);
    return;
  }
  await release(github, context, issue, assignees);
  await react(github, context, '+1');
  await reply(github, context, `Released. This one is free again, so anyone can \`/claim\` it.`);
}

// The comment event.
async function handleComment({ github, context, now = Date.now() }) {
  const { issue, comment } = context.payload;
  if (!issue || issue.pull_request || !comment) return 'skip';
  if (issue.state !== 'open' || issue.locked) return 'skip';
  if (comment.user.type === 'Bot' || comment.user.login.endsWith('[bot]')) return 'skip';
  const cmd = parseCommand(comment.body);
  if (!cmd) return 'skip';
  if (cmd === 'claim') await claim(github, context, now);
  else await unclaim(github, context);
  return cmd;
}

// Issue numbers a PR claims: `#N` in its title (`<summary> (#N)`) and closing keywords in its body.
const MAX_REFS = 5;
function referencedIssues(title, body) {
  const found = new Set();
  for (const m of String(title || '').matchAll(/(?:^|[\s(])#(\d{1,7})\b/g)) found.add(Number(m[1]));
  const closing = /\b(?:close[sd]?|fix(?:e[sd])?|resolve[sd]?)\s*:?\s+#(\d{1,7})\b/gi;
  for (const m of String(body || '').matchAll(closing)) found.add(Number(m[1]));
  return [...found].filter((n) => n > 0).slice(0, MAX_REFS);
}

// A pull request was opened: it claims each open-for-contributors issue it references that nobody
// holds yet. The core team's PRs claim nothing (their issues carry `in progress` instead).
async function handlePull({ github, context }) {
  const pull = context.payload.pull_request;
  if (!pull || pull.state !== 'open') return [];
  const author = pull.user && pull.user.login;
  if (!author || pull.user.type === 'Bot' || author.endsWith('[bot]')) return [];
  if (await isMaintainer(github, context, author)) return [];
  const claimed = [];
  for (const n of referencedIssues(pull.title, pull.body)) {
    if (n === pull.number) continue;
    let issue;
    try {
      ({ data: issue } = await github.rest.issues.get({ ...repoArgs(context), issue_number: n }));
    } catch (e) {
      if (e.status === 404 || e.status === 410) continue;
      throw e;
    }
    if (issue.pull_request || issue.state !== 'open' || issue.locked) continue;
    const assignees = (issue.assignees || []).map((a) => a.login);
    if (assignees.includes(author)) {
      if (!hasLabel(issue)) {
        await ensureLabel(github, context);
        await github.rest.issues.addLabels({ ...repoArgs(context), issue_number: n, labels: [LABEL] });
      }
      continue;
    }
    if (hasLabel(issue, CORE)) {
      await github.rest.issues.createComment({
        ...repoArgs(context),
        issue_number: pull.number,
        body: `#${n} is being worked on by the core team (\`${CORE}\`). Please ask on the issue before going further.`,
      });
      continue;
    }
    if (assignees.length || hasLabel(issue)) {
      const who = assignees.length ? assignees.map((a) => '@' + a).join(', ') : 'someone else';
      await github.rest.issues.createComment({
        ...repoArgs(context),
        issue_number: pull.number,
        body: `#${n} is already claimed by ${who}. If this PR duplicates their work, please close it or offer to help there.`,
      });
      continue;
    }
    // An older open PR for the issue wins (lowest number), even if it never recorded its claim.
    const timeline = await getTimeline(github, context, n);
    const older = linkedPulls(timeline).find((p) => p.number < pull.number && p.user && p.user.login !== author);
    if (older) {
      await github.rest.issues.createComment({
        ...repoArgs(context),
        issue_number: pull.number,
        body: `#${older.number} claimed #${n} first (the lowest-numbered open PR wins). Please close this one or offer to help there.`,
      });
      continue;
    }
    if (!claimable(issue)) {
      await github.rest.issues.createComment({
        ...repoArgs(context),
        issue_number: pull.number,
        body: `#${n} isn't open for contributors yet (no \`help wanted\` or \`good first issue\` label), so this PR doesn't claim it. A maintainer will reply on the issue.`,
      });
      continue;
    }
    await ensureLabel(github, context);
    await github.rest.issues.addLabels({ ...repoArgs(context), issue_number: n, labels: [LABEL] });
    try {
      await github.rest.issues.addAssignees({ ...repoArgs(context), issue_number: n, assignees: [author] });
    } catch (e) {
      // Assignment can be refused; the label and the comment naming them still stand.
    }
    await github.rest.issues.createComment({
      ...repoArgs(context),
      issue_number: n,
      body: `Claimed by @${author} through #${pull.number}. The claim lapses after ${CLAIM_DAYS} days without pushes.`,
    });
    claimed.push(n);
  }
  return claimed;
}

// The daily sweep over open issues labelled `claimed`.
async function sweep({ github, context, now = Date.now() }) {
  const issues = await github.paginate(github.rest.issues.listForRepo, {
    ...repoArgs(context),
    state: 'open',
    labels: LABEL,
    per_page: 100,
  });
  const cache = new Map();
  const released = [];
  for (const issue of issues) {
    if (issue.pull_request || issue.locked || hasLabel(issue, CORE)) continue;
    const assignees = (issue.assignees || []).map((a) => a.login);
    let skip = false;
    for (const a of assignees) skip = skip || (await isMaintainer(github, context, a, cache));
    if (skip) continue;
    const status = await claimStatus(github, context, issue, now);
    if (!status.release) continue;
    await release(github, context, issue, assignees);
    await github.rest.issues.createComment({
      ...repoArgs(context),
      issue_number: issue.number,
      body:
        status.reason === 'no PR'
          ? `No pull request showed up within ${CLAIM_DAYS} days, so this one is free again. Anyone can \`/claim\` it.`
          : `Nothing has been pushed for ${CLAIM_DAYS} days, so this one is free again. Anyone can \`/claim\` it.`,
    });
    released.push(issue.number);
  }
  return released;
}

module.exports = {
  LABEL,
  CORE,
  claimable,
  parseCommand,
  referencedIssues,
  isMaintainer,
  linkedPulls,
  claimStatus,
  handleComment,
  handlePull,
  sweep,
};
