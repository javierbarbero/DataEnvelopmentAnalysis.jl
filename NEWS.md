DataEnvelopmentAnalysis.jl release notes
========================================

Version 0.5.0 (May 5, 2021)
----------------------------
- Change default optimizer of Modified Directional Distance Function model to GLPK.

Version 0.4.0 (April 10, 2021)
----------------------------

- Add Hölder distance function model.
- Add Reverse directional distance function model.

Version 0.3.0 (March 18, 2021)
----------------------------

- Add Russell input, output, and graph models.
- Add Enhanced Russell Graph Slack Based Measure model.
- Add Modified Directional Distance Function model.
- Add an option to solve models using an optimizer specified by the user.
- Add an option to report economic efficiency in monetary terms.

Version 0.2.0 (January 3, 2021)
----------------------------

- Add input and output orientation in additive model.
- Add weak disposability in additive model.
- Cusotm weights in additive model are now specified with `rhoX` and `rhoY` instead of `wX` and `wY`.
- Add parameter to specify DMU names in all model functions.
- Add `peers` and `peersmatrix` functions to retur the peers.
- Add `targets` function to return optimal input and output targets.
- Add `normfactor` function to return normalization factor in economic models.
- Compatibility with GLPK 0.14.

Version 0.1.2 (Feb 25, 2020)
----------------------------

- Compatibility with JuMP 0.21.
- Weak disposability in radial, cost, and revenue models.
- Minor bugs fixed.

Version 0.1.1 (Oct 29, 2019)
----------------------------

- Compatibility with GLPK 0.12.

Version 0.1.0 (Oct 3, 2019)
---------------------------

- Initial release.
