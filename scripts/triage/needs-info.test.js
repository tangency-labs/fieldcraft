'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const info = require('./needs-info.js');

const repo = { owner: 'tangency-labs', repo: 'fieldcraft' };
const NOW = Date.parse('2026-10-18T12:00:00Z');
const DAY = 24 * 3600 * 1000;
const ago = (ms) => new Date(NOW - ms).toISOString();

function fake(world = {}) {
  const calls = [];
  const rec = (name) => async (args) => {
    calls.push([name, args]);
    return { data: {} };
  };
  const github = {
    paginate: async (fn, args) => (await fn(args)).data,
    rest: {
      repos: {
        getCollaboratorPermissionLevel: async ({ username }) => ({
          data: { permission: (world.perms || {})[username] || 'none' },
        }),
      },
      issues: {
        getLabel: async () => ({ data: {} }),
        createLabel: rec('createLabel'),
        addLabels: rec('addLabels'),
        removeLabel: rec('removeLabel'),
        createComment: rec('comment'),
        update: rec('update'),
        listForRepo: async () => ({ data: world.issues || [] }),
        listEventsForTimeline: async () => ({ data: world.timeline || [] }),
      },
    },
  };
  return { github, calls };
}
const find = (calls, n) => calls.find((c) => c[0] === n);
const labels = (names) => names.map((name) => ({ name }));
const long = 'x'.repeat(400);

const form = (observed, reproduction, version = '_No response_') =>
  `### Observed failure\n\n${observed}\n\n### Reproduction\n\n${reproduction}\n\n### Evidence\n\n_No response_\n\n### Expected behaviour\n\nworks\n\n### Version and platform\n\n${version}`;

test('missing: short reports without steps are asked, with the version when it is absent', () => {
  assert.deepEqual(info.missing('Solve does not work', 'The solver is not working.'), ['steps', 'version']);
  assert.deepEqual(info.missing('Solve does not work', 'Broken on v0.6.0.'), ['steps']);
  assert.deepEqual(info.missing('Hex mesh broken', form('it fails', 'mesh it')), ['steps', 'version']);
  assert.deepEqual(info.missing('Hex mesh broken', form('it fails', 'mesh it', 'v0.6.0, macOS')), ['steps']);
});

test('missing: reproducible, long, app-filed and non-bug issues are left alone', () => {
  assert.deepEqual(info.missing('Mesh fails', form('x', '1. fieldcraft --sample cantilever\n2. Mesh ▸ Tet')), []);
  assert.deepEqual(info.missing('Mesh fails', 'see bracket.fieldcraft attached'), []);
  assert.deepEqual(info.missing('Mesh fails', 'runs `mesh.tet.generate` and fails'), []);
  assert.deepEqual(info.missing('Mesh fails', long), []);
  assert.deepEqual(info.missing('Internal error in solve.static', form('x', 'Replay the attached session (3 command(s)): `fieldcraft-cli run session.json`.')), []);
  assert.deepEqual(info.missing('FieldCraft closed unexpectedly', form('x', '', 'v0.6.0') + '\nCrash report `abc` (panic), seen 1 time(s)'), []);
  assert.deepEqual(info.missing('Feature request: dark mode', 'please'), []);
  assert.deepEqual(info.missing('[FEAT] composites', 'please'), []);
  assert.deepEqual(info.missing('How do I apply a pressure?', ''), []);
});

const opened = (title, body, login = 'reporter') => ({
  repo,
  payload: { issue: { number: 4, state: 'open', title, body, labels: [], user: { login, type: 'User' } } },
});

test('an opened thin report is labelled and asked once', async () => {
  const { github, calls } = fake();
  assert.equal(await info.handleOpened({ github, context: opened('Solve broken', 'does not work') }), true);
  assert.deepEqual(find(calls, 'addLabels')[1].labels, ['needs-info']);
  const body = find(calls, 'comment')[1].body;
  assert.match(body, /Steps to reproduce/);
  assert.match(body, /FieldCraft version/);
  assert.match(body, /14 days/);
});

test('a maintainer\'s note is never asked', async () => {
  const w = fake({ perms: { jyviko: 'admin' } });
  assert.equal(await info.handleOpened({ github: w.github, context: opened('Solve broken', 'x', 'jyviko') }), false);
  assert.equal(w.calls.length, 0);
});

const reply = (login, { comment = true, body = 'x', issueLabels = ['needs-info'] } = {}) => ({
  repo,
  payload: {
    issue: { number: 4, title: 'Solve broken', body, labels: labels(issueLabels), user: { login: 'reporter' } },
    ...(comment ? { comment: { user: { login } } } : { sender: { login } }),
  },
});

test('the reporter\'s comment clears needs-info; someone else\'s does not', async () => {
  let w = fake();
  assert.equal(await info.handleReply({ github: w.github, context: reply('reporter') }), true);
  assert.equal(find(w.calls, 'removeLabel')[1].name, 'needs-info');
  w = fake();
  assert.equal(await info.handleReply({ github: w.github, context: reply('bystander') }), false);
  w = fake();
  assert.equal(await info.handleReply({ github: w.github, context: reply('reporter', { issueLabels: [] }) }), false);
  assert.equal(w.calls.length, 0);
});

test('an edit clears needs-info only once the steps are there', async () => {
  let w = fake();
  assert.equal(await info.handleReply({ github: w.github, context: reply('reporter', { comment: false, body: 'still broken' }) }), false);
  w = fake();
  const fixed = reply('reporter', { comment: false, body: '## Reproduction\n1. open the cantilever sample\n2. solve it' });
  assert.equal(await info.handleReply({ github: w.github, context: fixed }), true);
});

test('the sweep closes what waited 14 days, and keeps younger ones', async () => {
  const issue = { number: 4, locked: false, updated_at: ago(30 * DAY), labels: labels(['needs-info']) };
  const labeled = (at) => ({ event: 'labeled', label: { name: 'needs-info' }, created_at: at });
  let w = fake({ issues: [issue], timeline: [labeled(ago(15 * DAY))] });
  assert.deepEqual(await info.sweep({ github: w.github, context: { repo }, now: NOW }), [4]);
  assert.deepEqual(find(w.calls, 'update')[1], { ...repo, issue_number: 4, state: 'closed', state_reason: 'not_planned' });
  w = fake({ issues: [issue], timeline: [labeled(ago(20 * DAY)), labeled(ago(3 * DAY))] });
  assert.deepEqual(await info.sweep({ github: w.github, context: { repo }, now: NOW }), []);
  assert.equal(w.calls.length, 0);
});

test('backfill asks open thin reports once, and the dry run writes nothing', async () => {
  const issues = [
    { number: 1, state: 'open', title: 'Solve broken', body: 'nope', labels: [], user: { login: 'a', type: 'User' } },
    { number: 2, state: 'open', title: 'Solve broken', body: 'nope', labels: labels(['needs-info']), user: { login: 'a' } },
    { number: 3, state: 'open', title: 'x', body: 'nope', labels: [], pull_request: {}, user: { login: 'a' } },
  ];
  let w = fake({ issues });
  const plan = await info.backfill({ github: w.github, context: { repo }, dryRun: true });
  assert.deepEqual(plan.map((p) => p.number), [1]);
  assert.equal(w.calls.length, 0);
  w = fake({ issues });
  await info.backfill({ github: w.github, context: { repo }, dryRun: false });
  assert.equal(find(w.calls, 'addLabels')[1].issue_number, 1);
});

test('an issue from the feature form is never asked for steps', async () => {
  const { github, calls } = fake();
  const ctx = opened('Dark mode', '### What you need and why\nwork at night');
  ctx.payload.issue.labels = labels(['enhancement']);
  assert.equal(await info.handleOpened({ github, context: ctx }), false);
  assert.equal(calls.length, 0);
});
