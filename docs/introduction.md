# FieldCraft

FieldCraft is an open-source CAE application in pure Rust. It takes a parametric SolidCraft part
to a verified analysis: tetrahedral and hexahedral meshing, its own implicit, explicit, thermal
and lattice Boltzmann solvers, and decks for Abaqus, Nastran and LS-DYNA.

![The FieldCraft window: a clamped cantilever with its supports, pressure, gravity and a coupled reference point drawn in the viewport](media/setup-symbols.png)

It runs as:

- a **desktop app** for macOS, Windows and Linux, from the
  [releases](https://github.com/tangency-labs/fieldcraft/releases/latest);
- a **web app** in the browser, [try it here](app/);
- a **command line**, `fieldcraft-cli`, that runs, solves and exports models;
- an **MCP server** that AI agents analyse with ([AI agents](mcp.md)).

Everything the user can do is a command with JSON parameters. The toolbar, the dialogs, the
command line, the control channel and MCP all run the same commands ([features](features.md)).

Every capability is checked against closed forms and published benchmarks, with the results
frozen so they cannot drift ([verification and validation](verification.md)).

This book has five parts:
- the **user guide**, for people running analyses;
- the **solver manual**, the methods and their limits;
- **verification**, the evidence for each solver;
- the **developer guide**, for people working on FieldCraft;
- the **literature** the solvers and meshers are built from.

The source is on [GitHub](https://github.com/tangency-labs/fieldcraft).
