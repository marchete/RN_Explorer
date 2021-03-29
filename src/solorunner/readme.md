# Single CPU Solver for Number Shifting.

***This code is not for running on Codingame directly, but for using it locally***

Solutions are saved on files, that can be submitted to Codingame with PHP code. Submittable solutions are at the end of SOLUTION_* files, it starts with a password_level + list of moves. Tested on Ubuntu 18.04 LTS. 

## Prerrequisites:

- Python 3 (Tested on Python 3.6.9)
- Clang++-9
- Codingame account (user + password, not a github/google linked account. If you have one you can force a password change to have a password)
- cg_email.txt file with your user credentials. **NEVER SHARE THIS FILE WITH ANYONE**
- cg_password.txt file with your Codingames' password. **NEVER SHARE THIS FILE WITH ANYONE**

## Running the solver

python3 submit.py RN_ExploDiv_7 <THREADS> <LAHC_TIME_LIMIT> <LFA_SIZE> <K_A> <K_B> <K_C> <K_D> <INC_TIME> <INC_LFA> <RN_COUNT>

- \<THREADS>: On a physical machine up to 2x CPU core count. On Intel® Core™ i7-8700K I used 10.
- <LAHC_TIME_LIMIT>: Max time in seconds for LAHC Search. When it's timeout the thread resets and restart.
- <LFA_SIZE>:
- <K_A>: Remaining Numbers score
- <K_B>: Number of X Rows and Y Rows with numbers. I didn't use it on my final run.
- <K_C>: Remaining Points score
- <K_D>: Remaining Squared Points score.
- <INC_TIME>: Increase Search time on milliseconds, after restarts. Like 50ms to allow more time if the level is hard. The code limit the max search time to 70sec.
- <INC_LFA>: Increase LFA array on restarts. The code limit the Max LFA size to LFA_SIZE + 40*INC_LFA.
- <RN_COUNT>: Defines how many workers will change from LAHC mode to RN_Explore mode. A value <= <THREADS>. I always use THREADS-1 or THREADS-2, leaving a worker always in LAHC mode.

**Note:** There are a lot more of options and parameters inside the .cpp code. You have more parameters to tweak (probably more important than K_A..K_D) on lines 1491-1498 and 138-169.
Constant value *MIN_DEGREE_TOINSERT* defines at what amount of Remaining Numbers a Worker changes from LAHC to RN_Explorer. It can be anything from 3 to MAX_NUMBERS. I think a value 
between 6 and 9 is good for level 1000. A high value will remove LAHC and always use RN_Explorer. It will eventually solve it, but RN_Explorer has no reset timers, so you can end with
a local minimum that will take some time to go out.
 

## Features of submit.py

It's an evolved version of the recommended submit.py: https://github.com/eulerscheZahl/NumberShifting/tree/master/solver

- It's resilient to common errors: Solver crashes, wrong solutions, wrong replays, login errors
- It recompiles the code to the correct W,H and MAX_NUMBERS. This ensures maximum performance.
- On solver crashes (due to bugs on the code) it tries to recover the solution from SOLUTION\_*.txt file.
- It allows a _parameters.txt_ file. Useful when deploying multiples nodes, for changing parameters on the fly.
- It allows a _runningprocess.txt_ file.  Useful when deploying multiples nodes, for changing the running code on the fly.

## Performance

With an Intel® Core™ i7-8700K it solved levels 200 to 500 in 2hrs 41min (using 8/10 THREADS depending on the time), while I was doing other things in that PC.

## Compiling
CPP code was tested on CLANG and Visual Studio. It needs these compiler options ```/GS /GL /W3 /Gy- /Zc:wchar_t /Zi /Gm- /Ox /Ob2 /sdl- /Zc:inline /fp:precise /D "NDEBUG" /D "_CONSOLE" /D "_CRT_SECURE_NO_WARNINGS" /D "_UNICODE" /D "UNICODE" /errorReport:prompt /WX- /Zc:forScope /arch:AVX2 /Gd /Oy /Oi /MD /std:c++17 /FC /Fa"x64\Release\" /EHsc /nologo /Fo"x64\Release\" /Ot /diagnostics:classic ```
