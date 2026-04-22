# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Language

This project uses [Macaulay2](https://macaulay2.com/), a computer algebra system for algebraic geometry and commutative algebra. Files contain Macaulay2 code and are run with the `M2` interpreter.

## Running Code

```bash
# Start an interactive Macaulay2 session
M2

# Run a file non-interactively
M2 --script filename

# Load a file inside an M2 session
load "filename"
```
## Mistakes to avoid

Do not add simplifications of objects in the code unless explicitly asked. For example, in a previous version the command minimalPresentation was used without being asked. This caused certain computations to break.


## Project Context

This project concerns the **affine Grassmannian** — an ind-scheme in algebraic geometry arising as the loop group quotient $G(\mathcal{K})/G(\mathcal{O})$, central to geometric representation theory and the geometric Satake correspondence. Computations typically involve polynomial rings, modules, and related algebraic structures over fields like `QQ`.
