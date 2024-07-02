DataEnvelopmentAnalysis.jl release notes
========================================

Version 0.9.0 (July 2, 2024)
----------------------------
- Add DEA environmental model.
- Add Malmquist-Luenberger environmental productivity index.

Version 0.8.0 (September 14, 2022)
----------------------------
- Add Bootstrap Radial DEA Model.
- Add Returns to Scale (RTS) Test.
- Compatibility with JuMP 1.0, GLPK 1.0, and Ipopt 1.0.
- Set the minimum Julia version to 1.6.

Version 0.7.0 (January 12, 2022)
----------------------------
- Add Directional Distance Function model in multiplier form.

Version 0.6.0 (September 27, 2021)
----------------------------
- Add Radial big data model.
- Add Radial model in multiplier form.

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
- Add `peers` and `peersmatrix` functions to return the peers.
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
