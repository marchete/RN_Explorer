# RN_Explorer
Efficient Solver for Number Shifiting puzzle (https://www.codingame.com/ide/puzzle/number-shifting)

About 300x-1000x performance improvement from previous Solvers.

**Solving time in milliseconds**
Level | Dancing Links | LAHC | LAHC + RN_Explorer
------------ | ------------ | ------------- | -------------
 225 | 6318000 | 490000 | 435
 336 | Unsolved | 14400000 | 43000
 
 Level 1000 can be solved in under 60 minutes. Due to the random nature of the Solver the solving time isn't fixed. To improve times it can be parallelized, and with more cores the solver tends to solve all levels in under 20 minutes. I had 1 Full RN_Explorer node + N Half RN_Explorers, sending good aproximations to that central RN_Explorer.
 
 
See https://github.com/marchete/RN_Explorer/blob/main/Efficient%20Solver%20for%20Number%20Shifiting%20puzzle.md for more info.

It was an evolution of an older code: https://github.com/marchete/Codingame/blob/master/Optimization/Number%20Shift/NumberShift_LAHC.cpp
This old code doesn't perform that well.
