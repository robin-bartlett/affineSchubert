# Repository Guidelines

## Project Structure & Module Organization
This workspace is currently small and centered on Macaulay2 code.

- [`CLAUDE.md`](/Users/qp251849/Documents/grassmannian%20stability%20locus/CLAUDE.md) records tool-specific notes and the expected Macaulay2 workflow.


## Build, Test, and Development Commands
Use the Macaulay2 interpreter directly. For one of the X.m2 files in the folder:

- `M2` starts an interactive session.
- `M2 --script X.m2` runs the main file non-interactively.
- `M2 --script X.m2` exits on the first error, which is useful for quick verification.
- Inside `M2`, run `load "X.m2"` to reload the file while iterating.

There is no separate build step at the moment.


## Mistakes to avoid

Do not add simplifications of objects in the code unless explicitly asked. For example, in a previous version the command minimalPresentation was used without being asked. This caused certain computations to break.


## Testing Guidelines
There is no dedicated test suite yet. Until one exists, treat reproducible script execution as the baseline check.

- Run `M2 --stop --script grass.m2` before opening a pull request.
- For new helper files, add a small executable example at the bottom or in a separate `.m2` script.
- When behavior depends on a specific Grassmannian or endomorphism, document the input choice in comments.
