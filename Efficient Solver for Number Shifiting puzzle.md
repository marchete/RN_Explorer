## Efficient Solver for Number Shifiting puzzle

#### Abstract

A new search endgame algorithm is presented for Number Shifting puzzle, based on an exhaustive search of all permutations of remaining numbers of the puzzle. Tests on levels 225 and 336 show that this new approach improves search time of other algorithms (Dancing Links, LAHC) by up to 30.000%.

#### Game Description and Rules

From the game IDE:
>"You are given a grid, filled with numbers. You can move a number horizontally or vertically by exactly as many cells as the value of the number. The number has to be pushed on another non-zero number. The moved number will then be added to the other number or subtracted. The absolute value will be taken on subtraction. The goal is to clear the board and not have any numbers remaining.
The top left corner has the coordinate (0,0). X increases to the right, y to the bottom."

![Example Level](https://github.com/marchete/RN_Explorer/raw/main/img/Levelb.jpg)


Highest levels have up to 1050 numbers on a 56x32 grid. On average each number have 50 neighbour numbers (they share X or Y coordinate). Any number can potentially be sent to 50 numbers, and each number can receive values from its neighbours. In addition to this, order of operations matter to the final value of a cell. The search space is huge, and this makes Number Shifting a challenging puzzle to solve.

#### Game properties and characteristics

##### Total Numbers are defined at start

Each move always removes a number. If the result of the operation is valued zero, then two numbers are removed. 

##### Operations are noncommutative

This is a key property of the puzzle. Order of the moves are important in the game. 
Number shifting operation is defined as ```destValue = abs(destValue <sign> sourceValue)```
Splitted in addition and subtraction:
Addition: ```destValue = destValue+srcValue```
Subtraction: ```destValue = abs(destValue-srcValue)```
Example of noncommutative operations:
```3+1-4 = 0```
```3-4+1 = 2```
This property leads to more complexity on the search. The solver not only needs to make a certain list of moves, but it needs to be executed in a correct order. As number of moves increase finding the correct order of movements become one of the main challenges of the puzzle.

##### Move definition

A move is the action to use a number against other. Each move will always affect 2 numbers.
A move can be defined as:
```
class Move{
  int sourceID;  
  int destinationID;
  int sign;
  int sourceValue;
  int destinationValue
};
```
Coordinates X,Y of each number are static, they can retrieved from lookup tables by ID. Direction of the move can be inferred or saved on a variable. Source and destination values can be retrieved from the current Game State, but they are important for coherence checks.

##### Subtraction is more frequent than addition

In general all levels have around 10% additions and 90% subtractions. This can be a collateral result of the level generation.

##### Movements can be grouped in lists of moves (Strategies)

As each Move have a source and a destination ID, it creates a directed graph. Replaying a list of moves must be done in the same order as they are stored.
 In the Solver a list of moves is defined as a class with a ```vector<Move>```, and it's called ```Strategy```

##### Movements inside the same Strategy are highly correlated. Strategies are independent.
There aren't any interaction between lists of moves. Any change on one Strategy only affects itself.
Unfortunately the graph is very sensitive to any change on the Move list. In general any change on a Move will invalidate all succesive moves, as it is a directed graph. Values are stored on the Move class to ensure that a change won't destroy more moves than necessary. If a Move is going to be applied the Solver first checks if values are coherent. If they don't the move isn't applied.
##### Complete/Incomplete Strategies.

A complete Strategy ends with a subtraction Move where srcValue == destValue, so it removes both numbers, and all numbers used in the Strategy are zero. 
An incomplete Strategy always have one single number different than zero. In the Solver that number is marked as a Remaining Number (RN). Unused numbers of the level can also be considered as incomplete Strategies with zero moves.

### RN_Explore Solver

RN\_Explore solver is based on a heavily modified LAHC Solver (Late Acceptance Hill Climbing). On late game, when the solver reaches a low amount of Remaining Numbers (RN) it changes the behaviour to RN Explore Algorithm. This Algorithm keeps all new candidates as long as  ```RN <= MAX_INSERT_NUMBERS```, and apply an Exhaustive Search on each one. Bitwise operations are performed to avoid supersets or copies of candidates with the same RN values. I.e. if candidate **RN={3}** exists, remove **RN={3,4}**, **RN={1,3,7}**

##### Basic steps of the Solver:

0. Read initial data, create N thread Workers for the LAHC algorithm
1. Create a new candidate solution, based on a previous accepted APX, changing little things (removing a number, truncating list of numbers, etc)
2. Score it, based on some objective. I picked reducing both numbers in the grid and total sum of these numbers
3. Based on the LAHC algorithm, accept it as a new accepted APX or discard it.
4. Save candidate on GLOBAL_RN.APX pool if it was a good APX.  This is unrelated to point 3.
5. Once you have enough good APX ( with RN <= MAX_INSERT_NUMBERS), change LAHC for RN Explore Search Algorithm
6. Goto 1, repeat until you have a real solution with 0 points and no numbers in the grid.


##### Merging Strategies.

Incomplete Strategies have only one RN, so finding a way to complete that number without causing coherence issues is critical.

One number can be joined to an existing Strategy if ```2*destValue = srcValue```

 ```
destValue = abs(destValue-srcValue)
destValue = abs(destValue-2*destValue)
destValue = abs(-destValue)
destValue = destValue
 ```
These special moves are Merges, they don't affect the destination value, so they haven't coherence issues. This is a key property on the Solver, and they are widely used to merge incomplete Strategies on another Strategies (completes or not).

This has one indesired side effect, and it's the tendence of the solver to create long Strategies with a lot of Merged numbers. Solutions found usually have 1 or 2  long Strategies with 90% of the numbers, and smaller ones with a couple of moves each.

##### Permute Incomplete Strategies

Any incomplete Strategy have a list of {1..N} Moves where destinationID is equal to the remaining number. By doing permute operations in it can achieve new candidates with lower RN count. Permutation is complex, and expensive. Move list is not only permuted on position, but also in operation type (subtraction and addition).

**_Self solving an Incomplete Strategy_**
 Doing a permutation on all Moves affecting the ending ID of the Strategy generates different ending Values. 
 If ending value == 0 is found, the Strategy is completed. We create a new Strategy and try to add it to the pool of good approximations. It's done inside the Exhaustive Search.

**_Exhaustive Search for merging Incomplete Strategies_**
At endgame the list of moves are tighly integrated, and any mutation on Moves tend to fail. At this point with each approximation the RN Explore Algorithm do an Exhaustive Search to find hidden merges.

*Steps of Exhaustive Search:*
 ```
  1-Get list of all Remaining Numbers (RN).
  2-For each RN, get all possible permutations for the ending value, create a hashmap<RN,pair<endValue,PermuteInfo>>, only if endValue== 0 or at endValues that have valid neighbours at distance=endValue. This limitation is to have a valid points to merge.
  3-For reach RN, iterate over possible pair<endValue,PermuteInfo>
     a) if any endValue == 0 then that Incomplete Strategy is Self-Solved. Apply and continue to next RN.
     b) if endValue!=0, permute all valid neighbours. On these neighbours we need to search for the mergeValue ( 2*destValue = endValue). Note that the search is for all intermediate values, it's irrelevant when the mergeValue happens, the critical point is to achieve that value at any time. 
     c) if endValue != 0 and exists a target neighbour with value=2*endValue at any point on time, recreate moves to achieve the merge.
```

This search is expensive, in some candidates it can take seconds to complete. For the sake of speed the operation permutations are limited, and the whole Exhaustive search has timer limits for each candidate.

##### Modifications to LAHC Algortihm

Previous versions of the Solver were based entirely on LAHC. The algorithm struggled to pass levels 330+, and according to metadata it missed a LOT of solutions that I was able to recover with these tweaks.

*Changes to LAHC algorithm:*
1- Keep a history queue of the last N best APX. With a timer I keep track of the last improvement time and I reset the worker if too much time passed without a new best APX.
2- LFA size is changed based on domain knowledge of the problem. The search space is so huge on early and mid game that the algorithm was wasting too much time on initial steps of the search. Early game in Number Shifting there is little interest on accepting a candidate with a RN much higher than last Accepted candidate. My approach was taking into account RN to calculate a coverage percentage. According to the coverage:
- Coverage < 90%: LFA= 150
- 90%..96%: LFA: linearly increase from 150 to 4600
- \>96%: LFA:4600

3- High LFA values doesn't explore enough the best candidates. With a lower LFA it does it much more often, but it's more prone to fall in local maximum.
4- Due to 3, An Exhaustive Search is done on LFA timeouts if RN is low. This is an expensive step.
5- Due to 3, LFA timeouts have a Backtrack/Flashback feature. It recovers a previous best APX from the history queue.
According to metadata many solutions needed up to 3 flashbacks to get the solution. This means LFA falls on a local optimum, it may be avoided with a bigger LFA size or mutators with more points, but convergence times went too high.
6- With enough good APX the LAHC algorithm is switched to RN Explorer. RN Explorer have better better chances to find a valid solution.


##### RN Explore Search Algorithm
All threads use and feed the GLOBAL RN.APX pool of APX. It's a vector<vector<LAHC_Node>>, splitted by RN count.

+Steps to accept a new Best approximation:
1. Only save LAHC\_Nodes with RN <= MAX\_INSERT\_NUMBERS. Score is irrelevant.
2. If there are a many APX with RN <= N, then remove all APX with RN > N and don't accept them.
This is to avoid overflowing the system with so many APX.
3. With remainingNumbers do bitset validations to remove worse superset candidates. I.E. if I have a RN={3}, remove RN={3,4}, RN={1,3,7}, etc.
Pseudocode: ```if (newBest.RN.isFullyContainedIn(GLOBAL_RN.APX[d][b].RN)) then remove(GLOBAL_RN.APX[d][b]);```
4. If ```newBest.Equals(GLOBAL_RN.APX[d][b])```, then keep the one with the Best Score.


+Steps to search new candidates:
1. -Change behaviour of LAHC Worker once RN_Explore have some approximations. With bad quality approximations (totalNumbers>=5) only change a couple of threads.
2. Once the RN Explore have high quality approximations (totalnumbers < 4) change more workers to the new search mode.
3. Clone RN Explore on the local Worker, to avoid collisions and race conditions.
4. Worker with ID=0 will run the function Mutate_Exhaustive_6() on those approximations that haven't done it, then it will mark this approximation as exhausted.
5. Use a weighted random selection to pick an approximation, those with lower remaining numbers should have higher chance.
Use this as lastAccepted on LAHC worker.
6. Mutate to generate a new candidate, if it's a good approximation add it on GLOBAL_RN.APX .
7. Goto 1

##### Experimental Results

RN Explore Solver were used on levels 200 to 600. Tests were performed on a single Corei7-8700K, using 10 Threads. Solving time per level is stored at https://github.com/marchete/RN_Explorer/blob/main/Solver_Performance_Data.txt . These times are C++ execution time, not including Python script time, connections to Codingame or other external processes. All 400 levels were solved in 10104 seconds (less than 3 hours of real execution time). Solving time for levels <300 are negligible, below 500ms. The final version of the solver needs several hundreds of milliseconds to just prepare the cache and the worker threads.
Certain levels take much more time to solve, this is a graph with solving times per level.
![Performance Times](https://github.com/marchete/RN_Explorer/raw/main/img/PerformanceTime_200to600.PNG)

Levels 515, 547 and 570 were orders of magnitude harder than others.

Levels above 800 were solved with RN Explore Solver too, but it was in a distributed search with 5 nodes (10 threads + 4x VM with 8vCPU each). Average solving time was below 20minutes per level.

##### Conclusions

In this document a new approach for Number Shifting puzzle is presented. The proposed RN Explore Algorithm solved all 1000 levels of the puzzle in reasonable time, with a single computer. Solving time is 300 times better than previous approaches.

##### References
E.K.Burke and Y.Bykov, . "The Late Acceptance Hill-Climbing Heuristic".European Journal of Operational Research.
https://pdfs.semanticscholar.org/28cf/19a305a3b02ba2f34497bebace9f74340ade.pdf?_ga=2.123170187.1062017964.1581936670-1790746444.1580199348
