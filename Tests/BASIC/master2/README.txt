This is after program is set up to automatically add new constraints, on every repeat of master.
It is to test the addition of constraints and the integration of "main" and "master".
It still only tests one iteration since the subproblem is not included except at the start.

This agrees with results found with GAMS as in the previous test.

It is still signifigantly different than the final problem, for example by
not including a DC model of flow, and the subproblem still being an early example.

RuizC has several starting networks to test different starting conditions.