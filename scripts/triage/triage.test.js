'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const triage = require('./triage.js');

const repo = { owner: 'tangency-labs', repo: 'fieldcraft' };

// A fake github-script client that records writes. `world` describes the repo.
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
      pulls: { listFiles: async () => ({ data: (world.files || []).map((filename) => ({ filename })) }) },
      reactions: { createForIssueComment: rec('react') },
      issues: {
        getLabel: async () => {
          if (world.labelMissing) throw Object.assign(new Error('nf'), { status: 404 });
          return { data: {} };
        },
        createLabel: rec('createLabel'),
        addLabels: rec('addLabels'),
        removeLabel: rec('removeLabel'),
        createComment: rec('comment'),
        listForRepo: async () => ({ data: world.items || [] }),
        listLabelsForRepo: async () => ({
          data: (world.labels || ['bug', 'enhancement', 'area:mesh', 'priority: high', 'size: S', 'size: M', 'in progress', 'needs-triage']).map(
            (name) => ({ name }),
          ),
        }),
      },
    },
  };
  return { github, calls };
}

const find = (calls, n) => calls.find((c) => c[0] === n);
const all = (calls, n) => calls.filter((c) => c[0] === n);
const labels = (names) => names.map((name) => ({ name }));

test('classifyIssue: in-app reports get bug and the area from the title or its command', () => {
  assert.deepEqual(triage.classifyIssue('Internal error in mesh.tet.generate', ''), ['bug', 'area:mesh']);
  assert.deepEqual(triage.classifyIssue('FieldCraft closed unexpectedly in solve.static', ''), ['bug', 'area:solver']);
  assert.deepEqual(triage.classifyIssue('Internal error in results.show', ''), ['bug', 'area:ui']);
  assert.deepEqual(triage.classifyIssue('Crash on startup', ''), ['bug']);
});

test('classifyIssue: the forms give the type, the title the area', () => {
  const bug = '### Observed failure\nwrong\n\n### Reproduction\n1. x\n';
  assert.deepEqual(triage.classifyIssue('Lattice Boltzmann diverges at Re 200', bug), ['bug', 'area:cfd']);
  const feat = '### Current behaviour\nnone\n\n### What you need and why\nVTU\n';
  assert.deepEqual(triage.classifyIssue('ParaView export', feat), ['enhancement', 'area:io']);
  assert.deepEqual(triage.classifyIssue('Feature request: dark theme', ''), ['enhancement', 'area:ui']);
  assert.deepEqual(triage.classifyIssue('Hex meshing of thin walls', ''), ['area:mesh']);
  assert.deepEqual(triage.classifyIssue('Typo in the docs', ''), ['area:docs']);
});

test('classifyIssue: the body\'s session never gives the area, unmatched text gets nothing', () => {
  assert.deepEqual(triage.classifyIssue('Wrong numbers', '"mesh.tet.generate"\n"solve.static"'), []);
  assert.deepEqual(triage.classifyIssue('Refactor with types', 'words'), []);
  assert.deepEqual(triage.classifyIssue(undefined, undefined), []);
});

test('classifyPaths maps changed files to areas', () => {
  assert.deepEqual(
    triage.classifyPaths([
      'crates/mesh/src/a.rs',
      'crates/lbm-gpu/b.rs',
      'crates/fem/c.rs',
      'apps/fieldcraft-cli/d.rs',
      'crates/decks/e.rs',
      'docs/oracle.md',
      'README.md',
      '.github/workflows/ci.yml',
      'Cargo.lock',
    ]),
    ['area:cfd', 'area:docs', 'area:infra', 'area:io', 'area:mesh', 'area:solver', 'area:ui'],
  );
  assert.deepEqual(triage.classifyPaths(['crates/engine/src/cmd/x.rs']), []);
});

test('parseTriage', () => {
  assert.deepEqual(triage.parseTriage('/triage area:mesh high S -enhancement'), {
    add: ['area:mesh', 'priority: high', 'size: S'],
    remove: ['enhancement'],
  });
  assert.deepEqual(triage.parseTriage('/triage'), { add: [], remove: [] });
  assert.equal(triage.parseTriage('/triaged'), null);
  assert.equal(triage.parseTriage('please /triage'), null);
  assert.equal(triage.parseTriage('/triage ' + 'bug '.repeat(50)).add.length, 10);
});

test('an opened issue gets its labels plus needs-triage, which is created if missing', async () => {
  const { github, calls } = fake({ labelMissing: true });
  const issue = { number: 5, state: 'open', title: 'Internal error in mesh.hex.octree', body: '', labels: [] };
  assert.deepEqual(await triage.handleOpened({ github, context: { repo, payload: { issue } } }), [
    'bug',
    'area:mesh',
    'needs-triage',
  ]);
  assert.ok(find(calls, 'createLabel'));
  assert.deepEqual(find(calls, 'addLabels')[1].labels, ['bug', 'area:mesh', 'needs-triage']);
});

test('an opened issue keeps the labels it was filed with and skips needs-triage with a priority', async () => {
  const { github, calls } = fake();
  const issue = { number: 5, state: 'open', title: 'Solver diverges', body: '', labels: labels(['bug', 'priority: high']) };
  assert.deepEqual(await triage.handleOpened({ github, context: { repo, payload: { issue } } }), ['area:solver']);
  assert.deepEqual(find(calls, 'addLabels')[1].labels, ['area:solver']);
});

test('an opened PR is labelled from its paths and never gets needs-triage', async () => {
  const { github, calls } = fake({ files: ['crates/mesh/src/tet.rs'] });
  const pull_request = { number: 9, state: 'open', title: 'Mesh: thing', body: '', labels: [] };
  assert.deepEqual(await triage.handleOpened({ github, context: { repo, payload: { pull_request } } }), ['area:mesh']);
  assert.deepEqual(find(calls, 'addLabels')[1], { ...repo, issue_number: 9, labels: ['area:mesh'] });
});

test('a priority label clears needs-triage', async () => {
  const { github, calls } = fake();
  const issue = { number: 5, labels: labels(['needs-triage', 'priority: high']) };
  assert.equal(await triage.handleLabeled({ github, context: { repo, payload: { issue, label: { name: 'priority: high' } } } }), true);
  assert.equal(find(calls, 'removeLabel')[1].name, 'needs-triage');
  const other = fake();
  assert.equal(
    await triage.handleLabeled({ github: other.github, context: { repo, payload: { issue, label: { name: 'bug' } } } }),
    false,
  );
  assert.equal(other.calls.length, 0);
});

const comment = (body, login = 'carol', issueLabels = []) => ({
  repo,
  payload: {
    issue: { number: 5, labels: labels(issueLabels) },
    comment: { id: 1, body, user: { login, type: 'User' } },
  },
});

test('/triage by a maintainer sets labels, replaces the size and clears needs-triage', async () => {
  const { github, calls } = fake({ perms: { jyviko: 'admin' } });
  const ctx = comment('/triage area:mesh HIGH S', 'jyviko', ['size: M', 'needs-triage', 'bug']);
  assert.equal(await triage.handleComment({ github, context: ctx }), 'triage');
  assert.deepEqual(all(calls, 'removeLabel').map((c) => c[1].name).sort(), ['needs-triage', 'size: M']);
  assert.deepEqual(find(calls, 'addLabels')[1].labels, ['area:mesh', 'priority: high', 'size: S']);
  assert.equal(find(calls, 'react')[1].content, '+1');
});

test('/triage matches labels with spaces written with dashes', async () => {
  const { github, calls } = fake({ perms: { jyviko: 'admin' } });
  assert.equal(await triage.handleComment({ github, context: comment('/triage -in-progress', 'jyviko', ['in progress']) }), 'triage');
  assert.equal(find(calls, 'removeLabel')[1].name, 'in progress');
});

test('/triage from a non-maintainer or with an unknown label changes nothing', async () => {
  let w = fake();
  assert.equal(await triage.handleComment({ github: w.github, context: comment('/triage bug', 'mallory') }), 'denied');
  assert.equal(find(w.calls, 'addLabels'), undefined);
  w = fake({ perms: { carol: 'admin' } });
  assert.equal(await triage.handleComment({ github: w.github, context: comment('/triage area:nope') }), 'unknown');
  assert.equal(find(w.calls, 'addLabels'), undefined);
  assert.match(find(w.calls, 'comment')[1].body, /area:nope/);
});

test('bots and other comments are ignored', async () => {
  const { github, calls } = fake();
  const bot = comment('/triage bug');
  bot.payload.comment.user = { login: 'x[bot]', type: 'Bot' };
  assert.equal(await triage.handleComment({ github, context: bot }), 'skip');
  assert.equal(await triage.handleComment({ github, context: comment('hello') }), 'skip');
  assert.equal(calls.length, 0);
});

test('backfill only touches unlabelled items, and the dry run writes nothing', async () => {
  const items = [
    { number: 1, title: 'Internal error in mesh.tet.generate', body: '', labels: [] },
    { number: 2, title: 'Already sorted', body: '', labels: labels(['bug']) },
    { number: 3, title: 'x', body: '', labels: [], pull_request: {} },
  ];
  let w = fake({ items, files: ['crates/io/src/step.rs'] });
  const plan = await triage.backfill({ github: w.github, context: { repo }, dryRun: true });
  assert.deepEqual(
    plan.map((p) => [p.number, p.labels]),
    [
      [1, ['bug', 'area:mesh', 'needs-triage']],
      [3, ['area:io']],
    ],
  );
  assert.equal(w.calls.length, 0);
  w = fake({ items, files: ['crates/io/src/step.rs'] });
  await triage.backfill({ github: w.github, context: { repo }, dryRun: false });
  assert.deepEqual(all(w.calls, 'addLabels').map((c) => c[1].issue_number), [1, 3]);
});
