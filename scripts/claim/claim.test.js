'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const claim = require('./claim.js');

const NOW = Date.parse('2026-10-16T12:00:00Z');
const H = 3600 * 1000;
const D = 24 * H;
const ago = (ms) => new Date(NOW - ms).toISOString();

// A fake github-script client that records writes. `world` describes the repo.
function fake(world = {}) {
  const calls = [];
  const rec = (name) => async (args) => {
    calls.push([name, args]);
    if (world.fail && world.fail[name]) throw Object.assign(new Error(name), { status: world.fail[name] });
    return { data: {} };
  };
  const github = {
    paginate: async (fn, args) => (await fn(args)).data,
    rest: {
      repos: {
        getCollaboratorPermissionLevel: async ({ username }) => {
          const permission = (world.perms || {})[username] || 'none';
          return { data: { permission } };
        },
        getCommit: async () => ({ data: { commit: { committer: { date: world.pushed || ago(H) } } } }),
      },
      pulls: { get: async () => ({ data: { head: { sha: 'abc' } } }) },
      reactions: { createForIssueComment: rec('react') },
      issues: {
        getLabel: async () => {
          if (world.labelMissing) throw Object.assign(new Error('nf'), { status: 404 });
          return { data: {} };
        },
        createLabel: rec('createLabel'),
        addLabels: rec('addLabels'),
        removeLabel: rec('removeLabel'),
        addAssignees: rec('addAssignees'),
        removeAssignees: rec('removeAssignees'),
        createComment: rec('comment'),
        listForRepo: async () => ({ data: world.issues || [] }),
        get: async ({ issue_number }) => {
          const found = (world.byNumber || {})[issue_number];
          if (!found) throw Object.assign(new Error('nf'), { status: 404 });
          return { data: found };
        },
        listEventsForTimeline: async () => ({ data: world.timeline || [] }),
      },
    },
  };
  return { github, calls };
}

const names = (calls) => calls.map((c) => c[0]);
const find = (calls, n) => calls.find((c) => c[0] === n);
const labeled = (at) => ({ event: 'labeled', label: { name: 'claimed' }, created_at: at });
const pr = (number, login = 'other', state = 'open') => ({
  event: 'cross-referenced',
  source: { issue: { number, state, pull_request: {}, user: { login } } },
});

function ctx(body, { assignees = [], labels = ['help wanted'], user = 'alice', extra = {} } = {}) {
  return {
    repo: { owner: 'tangency-labs', repo: 'fieldcraft' },
    payload: {
      issue: {
        number: 7,
        state: 'open',
        locked: false,
        updated_at: ago(1000 * H),
        assignees: assignees.map((login) => ({ login })),
        labels: labels.map((name) => ({ name })),
        ...extra,
      },
      comment: { id: 99, body, user: { login: user, type: 'User' } },
    },
  };
}

test('parseCommand', () => {
  assert.equal(claim.parseCommand('/claim'), 'claim');
  assert.equal(claim.parseCommand('  /claim please'), 'claim');
  assert.equal(claim.parseCommand('/unclaim'), 'unclaim');
  assert.equal(claim.parseCommand('/claimed'), null);
  assert.equal(claim.parseCommand('I would /claim this'), null);
  assert.equal(claim.parseCommand(undefined), null);
});

test('claim on a free issue assigns, labels, reacts and replies', async () => {
  const { github, calls } = fake({ labelMissing: true });
  assert.equal(await claim.handleComment({ github, context: ctx('/claim'), now: NOW }), 'claim');
  assert.deepEqual(find(calls, 'addAssignees')[1].assignees, ['alice']);
  assert.deepEqual(find(calls, 'addLabels')[1].labels, ['claimed']);
  assert.equal(find(calls, 'createLabel')[1].color, '0052cc');
  assert.equal(find(calls, 'react')[1].content, '+1');
  assert.match(find(calls, 'comment')[1].body, /^Claimed by @alice\. .*`Closes #7`.*14 days/);
});

test('claim still labels and replies when assignment is refused', async () => {
  const { github, calls } = fake({ fail: { addAssignees: 422 } });
  await claim.handleComment({ github, context: ctx('/claim'), now: NOW });
  assert.ok(find(calls, 'addLabels'));
  assert.match(find(calls, 'comment')[1].body, /@alice/);
});

test('claim by the current claimant is a quiet thumbs up', async () => {
  const { github, calls } = fake();
  await claim.handleComment({ github, context: ctx('/claim', { assignees: ['alice'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.deepEqual(names(calls), ['react']);
});

test('claim on an issue claimed by someone else names them and when it lapses', async () => {
  const { github, calls } = fake({ timeline: [labeled(ago(10 * H))] });
  await claim.handleComment({ github, context: ctx('/claim', { assignees: ['bob'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.equal(find(calls, 'react')[1].content, 'eyes');
  const body = find(calls, 'comment')[1].body;
  assert.match(body, /@bob/);
  assert.match(body, /2026-10-30 02:00 UTC/); // 10 h ago + 14 days
  assert.equal(find(calls, 'addAssignees'), undefined);
});

test('a lapsed claim is taken over', async () => {
  const { github, calls } = fake({ timeline: [labeled(ago(15 * D))] });
  await claim.handleComment({ github, context: ctx('/claim', { assignees: ['bob'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.deepEqual(find(calls, 'removeAssignees')[1].assignees, ['bob']);
  assert.deepEqual(find(calls, 'addAssignees')[1].assignees, ['alice']);
});

test('a claim with a recent PR is not taken over', async () => {
  const { github, calls } = fake({ timeline: [labeled(ago(30 * D)), pr(50, 'bob')], pushed: ago(2 * H) });
  await claim.handleComment({ github, context: ctx('/claim', { assignees: ['bob'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.equal(find(calls, 'addAssignees'), undefined);
});

test('a maintainer\'s claim is never taken over', async () => {
  const { github, calls } = fake({ timeline: [labeled(ago(500 * H))], perms: { carol: 'write' } });
  await claim.handleComment({ github, context: ctx('/claim', { assignees: ['carol'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.equal(find(calls, 'addAssignees'), undefined);
  assert.equal(find(calls, 'react')[1].content, 'eyes');
});

test('an unclaimed issue with someone else\'s open PR points at the PR', async () => {
  const { github, calls } = fake({ timeline: [pr(42, 'bob')] });
  await claim.handleComment({ github, context: ctx('/claim'), now: NOW });
  assert.equal(find(calls, 'addAssignees'), undefined);
  assert.match(find(calls, 'comment')[1].body, /#42/);
});

test('bots, locked issues, closed issues, PRs and other comments are ignored', async () => {
  const { github, calls } = fake();
  const bot = ctx('/claim');
  bot.payload.comment.user = { login: 'dependabot[bot]', type: 'Bot' };
  assert.equal(await claim.handleComment({ github, context: bot, now: NOW }), 'skip');
  assert.equal(await claim.handleComment({ github, context: ctx('/claim', { extra: { locked: true } }), now: NOW }), 'skip');
  assert.equal(await claim.handleComment({ github, context: ctx('/claim', { extra: { state: 'closed' } }), now: NOW }), 'skip');
  assert.equal(await claim.handleComment({ github, context: ctx('/claim', { extra: { pull_request: {} } }), now: NOW }), 'skip');
  assert.equal(await claim.handleComment({ github, context: ctx('hello'), now: NOW }), 'skip');
  assert.equal(calls.length, 0);
});

test('unclaim by the claimant releases', async () => {
  const { github, calls } = fake();
  await claim.handleComment({ github, context: ctx('/unclaim', { assignees: ['alice'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.deepEqual(find(calls, 'removeAssignees')[1].assignees, ['alice']);
  assert.equal(find(calls, 'removeLabel')[1].name, 'claimed');
});

test('unclaim by a maintainer releases, by a stranger does not', async () => {
  let w = fake({ perms: { carol: 'admin' } });
  await claim.handleComment({ github: w.github, context: ctx('/unclaim', { user: 'carol', assignees: ['bob'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.ok(find(w.calls, 'removeLabel'));
  w = fake();
  await claim.handleComment({ github: w.github, context: ctx('/unclaim', { user: 'mallory', assignees: ['bob'], labels: ['help wanted', 'claimed'] }), now: NOW });
  assert.equal(find(w.calls, 'removeLabel'), undefined);
  assert.equal(find(w.calls, 'react')[1].content, 'eyes');
});

const issue = (number, assignees = ['bob'], extra = []) => ({
  number,
  locked: false,
  updated_at: ago(1000 * H),
  assignees: assignees.map((login) => ({ login })),
  labels: [{ name: 'claimed' }, ...extra.map((name) => ({ name }))],
});
const repo = { owner: 'tangency-labs', repo: 'fieldcraft' };

test('sweep releases a claim older than 14 days with no PR', async () => {
  const { github, calls } = fake({ issues: [issue(7)], timeline: [labeled(ago(15 * D))] });
  assert.deepEqual(await claim.sweep({ github, context: { repo }, now: NOW }), [7]);
  assert.ok(find(calls, 'removeLabel'));
  assert.deepEqual(find(calls, 'removeAssignees')[1].assignees, ['bob']);
  assert.match(find(calls, 'comment')[1].body, /free again/);
});

test('sweep keeps a young claim and a claim with a fresh PR, using the latest label time', async () => {
  let w = fake({ issues: [issue(7)], timeline: [labeled(ago(30 * D)), labeled(ago(5 * D))] });
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), []);
  w = fake({ issues: [issue(7)], timeline: [labeled(ago(30 * D)), pr(9)], pushed: ago(3 * D) });
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), []);
});

test('sweep releases a claim whose PR has been idle for 14 days', async () => {
  const w = fake({ issues: [issue(7)], timeline: [labeled(ago(30 * D)), pr(9)], pushed: ago(15 * D) });
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), [7]);
  assert.match(find(w.calls, 'comment')[1].body, /14 days/);
});

test('sweep ignores a closed PR', async () => {
  const w = fake({ issues: [issue(7)], timeline: [labeled(ago(15 * D)), pr(9, 'other', 'closed')] });
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), [7]);
});

test('sweep never touches core-team or maintainer-claimed issues', async () => {
  const w = fake({
    issues: [issue(1, ['bob'], ['in progress']), issue(2, ['carol']), issue(3, ['bob', 'carol'])],
    timeline: [labeled(ago(900 * H))],
    perms: { carol: 'write' },
  });
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), []);
  assert.equal(w.calls.length, 0);
});

test('sweep is idempotent: a released issue is no longer listed', async () => {
  const world = { issues: [issue(7)], timeline: [labeled(ago(15 * D))] };
  const w = fake(world);
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), [7]);
  world.issues = [];
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), []);
});

test('referencedIssues reads the title and closing keywords, capped', () => {
  assert.deepEqual(claim.referencedIssues('Fillet: chain (#123)', 'Fixes #9, closes: #10, see #11'), [123, 9, 10]);
  assert.deepEqual(claim.referencedIssues('no refs', 'mentions #4 only'), []);
  assert.deepEqual(claim.referencedIssues('#1 #2 #3 #4 #5 #6 #7', ''), [1, 2, 3, 4, 5]);
  assert.deepEqual(claim.referencedIssues(undefined, null), []);
});

const openIssue = (number, assignees = [], labels = ['help wanted']) => ({
  number,
  state: 'open',
  locked: false,
  assignees: assignees.map((login) => ({ login })),
  labels: labels.map((name) => ({ name })),
});
const pullCtx = (number, title, body = '', login = 'alice') => ({
  repo,
  payload: { pull_request: { number, state: 'open', title, body, user: { login, type: 'User' } } },
});

test('a PR referencing a free issue claims it for the author', async () => {
  const { github, calls } = fake({ byNumber: { 7: openIssue(7) } });
  assert.deepEqual(await claim.handlePull({ github, context: pullCtx(50, 'Thing (#7)') }), [7]);
  assert.deepEqual(find(calls, 'addAssignees')[1], { ...repo, issue_number: 7, assignees: ['alice'] });
  assert.deepEqual(find(calls, 'addLabels')[1].labels, ['claimed']);
  assert.match(find(calls, 'comment')[1].body, /Claimed by @alice through #50/);
});

test('a PR on an issue someone else holds comments on the PR and changes nothing', async () => {
  const { github, calls } = fake({ byNumber: { 7: openIssue(7, ['bob'], ['help wanted', 'claimed']) } });
  assert.deepEqual(await claim.handlePull({ github, context: pullCtx(50, 'Thing (#7)') }), []);
  assert.equal(find(calls, 'addAssignees'), undefined);
  const c = find(calls, 'comment')[1];
  assert.equal(c.issue_number, 50);
  assert.match(c.body, /@bob/);
});

test('the lowest-numbered open PR wins an unrecorded claim', async () => {
  const { github, calls } = fake({ byNumber: { 7: openIssue(7) }, timeline: [pr(42, 'bob')] });
  assert.deepEqual(await claim.handlePull({ github, context: pullCtx(50, 'Thing (#7)') }), []);
  assert.match(find(calls, 'comment')[1].body, /#42 claimed #7 first/);
});

test('the claimant\'s own PR just makes sure the label is there', async () => {
  const { github, calls } = fake({ byNumber: { 7: openIssue(7, ['alice']) } });
  assert.deepEqual(await claim.handlePull({ github, context: pullCtx(50, 'Thing (#7)') }), []);
  assert.deepEqual(names(calls), ['addLabels']);
});

test('PRs from bots, references to PRs, closed and missing issues are ignored', async () => {
  const w = fake({
    byNumber: { 7: { ...openIssue(7), pull_request: {} }, 8: { ...openIssue(8), state: 'closed' } },
  });
  assert.deepEqual(await claim.handlePull({ github: w.github, context: pullCtx(50, 'x (#7) (#8) (#9)') }), []);
  const bot = pullCtx(51, 'x (#1)');
  bot.payload.pull_request.user = { login: 'dependabot[bot]', type: 'Bot' };
  assert.deepEqual(await claim.handlePull({ github: w.github, context: bot }), []);
  assert.equal(w.calls.length, 0);
});

test('an issue the core team holds (`in progress`) cannot be claimed', async () => {
  const { github, calls } = fake();
  await claim.handleComment({ github, context: ctx('/claim', { labels: ['help wanted', 'in progress'] }), now: NOW });
  assert.equal(find(calls, 'addAssignees'), undefined);
  assert.match(find(calls, 'comment')[1].body, /core team/);
});

test('only help wanted / good first issue issues can be claimed, except by a maintainer', async () => {
  let w = fake();
  await claim.handleComment({ github: w.github, context: ctx('/claim', { labels: ['bug'] }), now: NOW });
  assert.equal(find(w.calls, 'addAssignees'), undefined);
  assert.match(find(w.calls, 'comment')[1].body, /help wanted/);
  w = fake();
  await claim.handleComment({ github: w.github, context: ctx('/claim', { labels: ['good first issue'] }), now: NOW });
  assert.ok(find(w.calls, 'addAssignees'));
  w = fake({ perms: { carol: 'write' } });
  await claim.handleComment({ github: w.github, context: ctx('/claim', { user: 'carol', labels: [] }), now: NOW });
  assert.deepEqual(find(w.calls, 'addAssignees')[1].assignees, ['carol']);
});

test('a PR on a core-team or not-yet-open issue comments on the PR and claims nothing', async () => {
  const { github, calls } = fake({
    byNumber: { 7: openIssue(7, [], ['help wanted', 'in progress']), 8: openIssue(8, [], ['bug']) },
  });
  assert.deepEqual(await claim.handlePull({ github, context: pullCtx(50, 'Thing (#7) (#8)') }), []);
  assert.equal(find(calls, 'addLabels'), undefined);
  const bodies = calls.filter((c) => c[0] === 'comment').map((c) => c[1].body);
  assert.match(bodies[0], /core team/);
  assert.match(bodies[1], /isn't open for contributors/);
});

test('a maintainer\'s PR claims nothing', async () => {
  const { github, calls } = fake({ byNumber: { 7: openIssue(7) }, perms: { alice: 'admin' } });
  assert.deepEqual(await claim.handlePull({ github, context: pullCtx(50, 'Thing (#7)') }), []);
  assert.equal(calls.length, 0);
});

test('a claim whose PR is stale but which was claimed recently stands', async () => {
  const w = fake({ issues: [issue(7)], timeline: [labeled(ago(2 * D)), pr(9)], pushed: ago(30 * D) });
  assert.deepEqual(await claim.sweep({ github: w.github, context: { repo }, now: NOW }), []);
});
