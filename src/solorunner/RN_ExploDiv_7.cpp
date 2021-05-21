/**
Solver for Number Shifting https://www.codingame.com/ide/puzzle/number-shifting

Abbreviations:
RN : Remaining Numbers on a candidate solution. These are numbers on the grid that its value is different from zero. Once the solver reaches RN==0 the level is solved.
APX: Approximation solution, a good candidate with low RN. I consider these with potential to be near the real solution.

***This code is not for running on Codingame directly, but for using it locally, it works better with the submit.py script***
There are small parts omitted on the code, but simple enough to recreate it in 10 minutes.They are marked with a TODO:
The idea is not only Copy&Paste it, but learn from it. As I reached the max level I no longer care about improvements.

A Corei7-8700K CPU is able to solve levels 900+ in about 15-30 minutes on average, using 10 threads.
This code is highly parallelizable and restartable. The more threads, the better. But it has diminished returns.
At highest levels it's hard to go below 10 minutes because LAHC algorithm takes some time to reach a good APX.
Nevertheless having more CPU's reduces the worst time on all levels.
Also it seems that raw power per thread is very important too, on last levels my Corei7-8700K solved 90% of the levels, vs 4x8vCPUs on cloud computing.
My local Corei7 achieved around 100k candidates per thread per second vs 50k candidates on cloud vCPU.


Loosely based on Late Acceptance Hill-Climbing:
//E.K.Burke and Y.Bykov, . "The Late Acceptance Hill-Climbing Heuristic".European Journal of Operational Research.
// https://pdfs.semanticscholar.org/28cf/19a305a3b02ba2f34497bebace9f74340ade.pdf?_ga=2.123170187.1062017964.1581936670-1790746444.1580199348
With a lot of changes. LAHC by itself struggled to pass levels 330+, and according to metadata it missed a LOT of solutions that I was able to recover
with these tweaks.
Tweaks to LAHC:
1-I keep an history queue of the last N best APX. With a timer I keep track of the last improvement time and I reset the Worker if too much time passed without a new best APX.
1-LFA size is changed based on domain knowledge of the problem. The search space is so huge on early and mid game that the algorithm was wasting
too much time on initial steps of the search. In this game there is little interest on accepting a candidate with a RN much higher than last Accepted candidate.
My approach was taking into account RN to calculate a coverage percentage. According to the coverage:
Coverage < 90%: LFA= 150
90%..96%: LFA: linearly increase from 150 to 4600
>96%: LFA:4600
2-I think high LFA values doesn't explore enough the best candidates. With a lower LFA it does it much more often, but it's more prone to fall in local maximum.
3-Due to 2, I try to do an Exhaustive (==expensive) on LFA timeouts if RN is low.
4-Due to 2, I changed LFA timeouts to do a Backtrack/Flashback. I recover a previous best APX from the history queue.
According to metadata many solutions needed up to 3 flashbacks to get the solution. This means LFA falls on a local optimum, it may be avoided with a bigger LFA
size, but convergence times went too high.
5-With enough good APX I just disable LAHC and no longer use it. In highest levels at that point I usually get a solution in less than 10 minutes.


Basic steps of the Solver:
0- Read initial data, create N thread Workers for the LAHC algorithm
1- Create a new candidate solution, based on a previous accepted APX, changing little things (removing a number, truncating list of numbers, etc)
2- Score it, based on some objective. I picked reducing both numbers in the grid and total sum of these numbers
3- Based on the LAHC algorithm, accept it as a new accepted APX or discard it.
4- Save candidate on GLOBAL_RN.APX pool if it was a good APX.  This is unrelated to point 3.
5- Once you have enough good APX ( with RN <= MAX_INSERT_NUMBERS), change LAHC for RN Explore Search Algorithm
6- Goto 1, repeat until you have a real solution with 0 points and no numbers in the grid.

*/

/********************  RN Explore Search Algorithm ************************************
All threads use and feed the GLOBAL_RN.APX pool of APX. It's a vector<vector<LAHC_Node>>, splitted by RN count.
I.e. all APX are stored in GLOBAL_RN.APX[APX.totalNumbers.count()-1] vector

+Steps to accept a new Best approximation:
1-Only save LAHC_Nodes with RN <= MAX_INSERT_NUMBERS. Score is irrelevant.
2-If there are a lot APX with RN <= N, then remove all APX with RN > N and don't accept them.
This is to avoid overflowing the system with so many APX.
3-With remainingNumbers do bitset validations to remove worse superset candidates. I.E. if I have a RN={3}, remove RN={3,4}, RN={1,3,7}, etc.
Pseudocode: if (newBest.RN.isFullyContainedIn(GLOBAL_RN.APX[d][b].RN)) then remove(GLOBAL_RN.APX[d][b]);
4-If newBest.Equals(GLOBAL_RN.APX[d][b]), then keep the one with the Best Score.


+Steps to search new candidates:
1-Change behaviour of LAHC Worker once RN_Explore have some approximations. With bad quality approximations (totalNumbers>=5) only change a couple of threads.
2-Once the RN Explore have high quality approximations (totalnumbers < 4) change all workers to the new search mode.
3-Clone RN Explore on the local Worker, to avoid collisions and race conditions.
4-Worker with ID=0 will run the function Mutate_Exhaustive_6() on those approximations that haven't done it, then it will mark this approximation as exhausted.
5-Use a weighted random selection to pick an approximation, those with lower remaining numbers should have higher chance.
Use this as lastAccepted on LAHC worker.
6-Mutate to generate a new candidate, if it's a good approximation add it on GLOBAL_RN.APX .
7-Rinse and repeat
***************************************************************************************/

//No pragmas, this is offline
#include <immintrin.h> //AVX SSE Extensions
#include <algorithm>
#include <atomic>
#include <bitset>
#include <cassert>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <string>
#include <thread> 
#include <unordered_set>
#include <unordered_map>
#ifndef _MSC_VER
#include <bits/stdc++.h> //All main STD libraries
#endif

using namespace std;

#ifdef _MSC_VER
//Windows stack increase, on Linux it's done with stacksize()
#pragma comment(linker, "/STACK:33554432")
#endif
#ifdef _MSC_VER
#define ALIGN
//#define USE_AVX
//#define ALIGN __declspec(align(32))
#else 
//#define USE_AVX
#define ALIGN
//#define ALIGN __attribute__((aligned(32)))
#endif

class Random;

struct Grid; //Current GameState
struct Move; //A single Move
struct ALIGN Strategy; //A list of related Moves
struct LAHC_Node; //Lists of Strategies + current Grid+remainingNumbers bitset+some metadata

				  //I adjust these constants on each level to  W,H and COUNT_NUMBERS+1, and I recompile aftewards. This improves speed on low levels
const int MAX_W = 46;
const int MAX_H = 26;
const int MAX_NUMBERS = 666;


//I had two pseudo-random number generators, that switch from one to another
#define ALTERNATE_RANDOM
#define SHUFFLE_CODE //Shuffling moves on Strategies to try to get different intermediate values on numbers.
//#define DEBG_MODE //This enables some asserts


const uint16_t LIM_P = 8;
const uint16_t LIM_SQ = 8;

int K_B_NUM = 5;
int CNT_RECOMB = 0;

const int SHUFFLE_TRIES = 2;
const int K_ALLOWSUM = 10;
const int K_MAX_POINTS_REMOVE = 4;
const int K_FLASHBACK_NUMBERS = 4;
const int K_FLASHBACK_POINTS = 17;
const int FLASHBACK_COUNT = 2; //ERA 3

							   //I use an small LFA (300) until I get a coverage equals to K_LFA_PERC_MIN percent.
							   //Then I linearly increase it until coverage is K_LFA_PERC_MAX
							   //This is domain knowledge, early game there is too much room for move generation to allow much worse movements.
							   //That tweak improves convergence time without a big impact on endgame diversity.
const double K_LFA_PERC_MIN = 90.0;//85.0;
const double K_LFA_PERC_MAX = 96.0;

/***** STRINGS_4 PARAMETERS *****/
const int K_LOW_POINTS_REMOVE = 2;//2;
const int K_CHANCE_POINTS_REMOVE = 4;
const int K5_MERGE0 = 750;
const int K5_MERGE1 = 920;

const int K_CHANCE_SIGN = 5;
const int K_CHANCE_IGNORE = 10;
const int K_CHANCE_MERGE = 920;
const int K_CHANCE_EXPLODE = 20;
const int K_PREFER_ENDINGS = 600;
const int K_TAIL_PERCENT = 75;
const int K_SEL_VAL_RANDOM = 70;
const int K_SEL_VAL_TAIL = 300;
const int K_SEL_VAL_OPT = 300;
const int K_SEL_VAL_ORPHAN = 300;
const int K_SEL_VAL_CROSS = 100;
const int K_SEL_VAL_TOTAL = K_SEL_VAL_RANDOM + K_SEL_VAL_TAIL + K_SEL_VAL_OPT + K_SEL_VAL_ORPHAN + K_SEL_VAL_CROSS;
/***** STRINGS_4 PARAMETERS *****/

//Real values of the level
int W, H;
int COUNT_NUMBERS = MAX_NUMBERS;
int level = -1;
string passwordLevel = "first_level";

//Solutions vars
atomic<bool> solved;

//Node vars, to keep track what computer solved the level
string PROGRAM_NAME = "NONAMED";
string COMPUTER_NAME = "NONAMED";
int THREADS = 4;

//LAHC Global variables.
int SIZE_LFA = 3600;
int CAN_INCREASE_LFA = 0;
int CAN_INCREASE_TIME = 0;
long long LIMIT_TIME_IMPROVEMENT = 59 * 1000; //In ms

											  //For scoring purposes
double K_A = 60000.0;
double K_B = 0.0;
double K_C = 20000.0;
double K_D = 40000.0;
// to avoid divisions
double INV_COUNT_NUMBERS = 1.0;
double INV_GROUPS = 1.0;
double INV_POINTS = 1.0;
double INV_TOTALNUMBERS = 1.0;
double INV_SQUAREPOINTS = 1.0;
double INV_RESX = 1.0;
double INV_RESY = 1.0;


int ALLOW_GLOBAL_BEST = 999; //Numbers of threads that can go from LAHC to Recombinate

							 //Thread mutexes, to avoid problems between threads working on the same data
std::mutex mutexRN_Explorer;
std::mutex mutexSOL;
std::mutex mutexSAVE;
std::mutex mutexAPROX;

//For Aproximations save, to avoid rewriting when not needed
long filesize = 0;

//Performance counters
atomic<uint64_t> newFound;
atomic<uint64_t> JACKPOT_OK;
atomic<uint64_t> SUPERCOMBINATOR;
atomic<uint64_t> MAX_PERMUTATIONS;
atomic<uint64_t> JACKPOT_TOTAL;
int64_t performance;
atomic<uint64_t> SimCount;
#ifdef SHUFFLE_CODE
atomic<uint64_t> correctShuffles;
atomic<uint64_t> countShuffles;
#endif

//Not very useful
int spawns = -1;
int GROUPS = -1;

//***************** Number Shift constants **************/
const int U = 0;
const int D = 1;
const int L = 2;
const int R = 3;
const int DIR_X[] = { 0,0,-1,1 };
const int DIR_Y[] = { -1,1,0,0 };
const string strDIRS[] = { "U","D","L","R" };

const uint16_t INVALID_ID = (uint16_t)(MAX_NUMBERS + 2);

#ifdef DEBG_MODE
#define ASSERT(x) assert(x)
#else 
#define ASSERT(x) 
#endif

typedef int32_t I;

#ifndef _MSC_VER
#include <sys/resource.h>
#include <unistd.h>
#endif


//***************** Timing Control - precise timing is not critical in this case **************/
#include <chrono>
#define Now() chrono::high_resolution_clock::now()
struct Stopwatch {
	chrono::high_resolution_clock::time_point c_time, c_timeout;
	void Start(uint64_t us) { c_time = Now(); c_timeout = c_time + chrono::microseconds(us); }
	void setTimeout(uint64_t us) { c_timeout = c_time + chrono::microseconds(us); }
	inline bool Timeout() {
		return Now() > c_timeout;
	}
	long long EllapsedMicroseconds() { return chrono::duration_cast<chrono::microseconds>(Now() - c_time).count(); }
	long long EllapsedMilliseconds() { return chrono::duration_cast<chrono::milliseconds>(Now() - c_time).count(); }
} stopwatch;


//***************** Random Numbers - There are two alternatives, switched with the #define ALTERNATE_RANDOM **************/
#include <random>
uint64_t randomKey = 0;
#ifdef ALTERNATE_RANDOM
class Random {
public:
	uint64_t _seed1;
	uint64_t _seed2;

	Random(uint64_t SEED) {
		std::mt19937_64 gen(SEED);
		std::uniform_int_distribution<uint64_t> dis;
		_seed1 = dis(gen);
		_seed2 = dis(gen);
		if (SEED == randomKey)
		{
			cerr << "**************************************************************" << endl;
			cerr << "**************************************************************" << endl;
			cerr << "*******WARNIIIIIIINGGG USING PRESEED NUMBERS!!!!!!!!!!!*******" << endl;
			cerr << "**************************************************************" << endl;
			cerr << "**************************************************************" << endl;
			cerr << "Seed1:" << bitset<64>(_seed1) << endl;
			cerr << "Seed2:" << bitset<64>(_seed2) << endl;
		}
	}
	Random() {
		std::random_device rd;
		std::seed_seq seedseq1{ rd(), rd(), rd() , rd() };
		std::mt19937_64 gen(seedseq1);
		std::uniform_int_distribution<uint64_t> dis;
		_seed1 = dis(gen);
		_seed1 ^= Now().time_since_epoch().count();
		_seed2 = dis(gen);
		_seed2 ^= (Now().time_since_epoch().count() * 7);
	}
	inline uint64_t xrandom() {
		auto s0 = _seed1;
		auto s1 = _seed2;
		auto result = s0 + s1;
		s1 ^= s0;
		_seed1 = ((s0 << 55) | (s0 >> (64 - 55))) ^ s1 ^ (s1 << 14); // a, b
		_seed2 = ((s1 << 36) | (s1 >> (64 - 36))); // c
		return result;
	}
	inline bool NextBool() {
		return (xrandom() & 4) == 4;
	}

	//	inline uint32_t Next1024() {return (uint32_t)xrandom() & 1023;}

	inline uint32_t NextInt(const uint32_t& range) {
		if (range == 0)
			return 0;
		return (uint32_t)xrandom() % range;
	}
	inline int32_t NextInt(const int32_t& a, const int32_t&  b) {
		return  (int32_t)NextInt((uint32_t)(b - a + 1)) + a;
	}
	inline float NextFloat() {
		uint32_t xr = (uint32_t)xrandom();
		if (xr == 0U) return 0.0f;
		union
		{
			float f;
			uint32_t i;
		} pun = { (float)xr };
		pun.i -= 0x10000000U;
		return  pun.f;
	}
	inline float NextFloat(const float& a, const float& b) {
		return NextFloat()*(b - a) + a;
	}

};
#else

class Random {
public:
	const uint64_t K_m = 0x9b60933458e17d7d;
	const uint64_t K_a = 0xd737232eeccdf7ed;
	mt19937 E4;
	uint64_t seed;

	Random(int SEED) {
		std::mt19937_64 gen(SEED);
		std::uniform_int_distribution<uint64_t> dis;
		seed = dis(gen);
		if (SEED == randomKey)
		{
			cerr << "**************************************************************" << endl;
			cerr << "**************************************************************" << endl;
			cerr << "*******WARNIIIIIIINGGG USING PRESEED NUMBERS!!!!!!!!!!!*******" << endl;
			cerr << "**************************************************************" << endl;
			cerr << "**************************************************************" << endl;

			cerr << "Seed:" << bitset<64>(seed) << endl;
		}

	}
	Random() {
		std::random_device rd;
		std::seed_seq seedseq1{ rd(), rd(), rd() , rd() }; // is there an optimal number of rd() to use?
		std::mt19937_64 gen(seedseq1);
		std::uniform_int_distribution<uint64_t> dis;
		seed = dis(gen);
		seed ^= Now().time_since_epoch().count();
		cerr << "Seed:" << bitset<64>(seed) << endl;
	}
	inline uint32_t xrandom() {
		//PCG 
		seed = seed * K_m + K_a;
		return (uint32_t)(seed >> (29 - (seed >> 61)));
	}
	inline bool NextBool() {
		return (xrandom() & 4) == 4;
	}
	inline uint32_t NextInt(const uint32_t& range) {
		return xrandom() % range;
		/*uniform_int_distribution<uint32_t> dist(0, range-1);
		return dist(E4);*/
	}
	inline int32_t NextInt(const int32_t& a, const int32_t&  b) {
		return  (int32_t)NextInt((uint32_t)(b - a + 1)) + a;
	}
	inline float NextFloat() {
		uint32_t xr = xrandom();
		if (xr == 0U) return 0.0f;
		union
		{
			float f;
			uint32_t i;
		} pun = { (float)xr };
		pun.i -= 0x10000000U;
		return  pun.f;
	}
	inline float NextFloat(const float& a, const float& b) {
		return NextFloat()*(b - a) + a;
	}

};
#endif
//Random rnd; //Main random number generator


//***************** Helper Functions  **************/

//Fisher–Yates shuffle. I use this function a lot. I shuffle lists of moves to give diversity
template <class T>
inline void do_shuffle(T& smallDLX, Random& rnd)
{
	for (int i = (int)smallDLX.size() - 1; i > 0; --i) {
		int r = rnd.NextInt(i + 1);
		if (i != r) {
			swap(smallDLX[i], smallDLX[r]);
		}
	}
}
inline bool file_exists(const std::string& name) {
	ifstream f(name.c_str());
	return f.good();
}

void Add(vector<size_t>&A, const vector<size_t>& B)
{
	for (auto&b : B)
		A.push_back(b);
}

long GetFileSize(std::string filename)
{
	ifstream in_file(filename, ios::binary);
	in_file.seekg(0, ios::end);
	return (long)in_file.tellg();
}
//Weighted random search
int rndLowBound(const vector<size_t>& ListSizes, const int& valPoint)
{
	int selList = (int)(lower_bound(ListSizes.begin(), ListSizes.end(), valPoint, [](auto &a, auto &b) { return a <= b; }) - ListSizes.begin());
	return min(max(0, selList), (int)ListSizes.size() - 1);
}

vector<string> split(const string& str, const string& delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == string::npos) pos = str.length();
		string token = str.substr(prev, pos - prev);
		if (!token.empty()) tokens.push_back(token);
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}
//This comes from the referee, not very useful
void loadValues(int _level) {
	level = _level;
	spawns = 3 + level / 2;
	if (level > 150) spawns = 3 + level - 75;
	H = 5;
	W = H * 16 / 9;
	while (W * H < spawns * 2) {
		spawns -= 2;
		H++;
		W = H * 16 / 9;
	}
};

struct Coord { int X; int Y; uint16_t ID; };
vector<Coord> NMB;
//Lookup table to convert number ID to coords.
uint16_t IDX[MAX_W][MAX_H];


#ifdef USE_AVX
inline bool all_zeros(__m256i x)
{
	return _mm256_testz_si256(x, x) == 1;
}
#endif

#ifdef _MSC_VER
#define __builtin_popcountll _mm_popcnt_u64
#endif

#ifdef USE_AVX
const int AVX_W = ((MAX_NUMBERS / 256) + (MAX_NUMBERS % 256 ? 1 : 0));
#endif


//***************** Bitset-like class **************/

const int INT64_W = ((MAX_NUMBERS / 64) + (MAX_NUMBERS % 64 ? 1 : 0));
struct ALIGN NumberSet {
	static ALIGN NumberSet Mask;
	union ALIGN {
#ifdef USE_AVX
		__m256i  X[AVX_W];
		uint64_t W[4 * AVX_W];
#else
		uint64_t W[INT64_W];
#endif	
		//		uint8_t  B[32 * AVX_W];
	};

	NumberSet() { clear(); }
	NumberSet(vector<int> V) {
		clear();
		for (auto&v : V)
			set((size_t)v);
	}
	NumberSet(vector<size_t> V) {
		clear();
		for (auto&v : V)
			set(v);
	}
	inline void set(const size_t& N) {
		uint64_t element = N / 64;
		uint64_t pos = 63 - (N % 64);
		W[element] |= (1ULL << pos);
	}
	inline bool get(const int& N)const {
		uint64_t element = N / 64;
		uint64_t pos = 63 - (N % 64);
		return (W[element] & (1ULL << pos)) != 0;
	}
	void unset(const int& N) {
		uint64_t element = N / 64;
		uint64_t pos = 63 - (N % 64);
		W[element] &= ~(1ULL << pos);
	}

	void negate() {
		for (int i = 0; i < INT64_W; ++i)
		{
			W[i] = ~W[i];
		}
		this->_and(Mask);
	}
	bool Equals(const NumberSet& DST) const
	{
		for (int i = 0; i < INT64_W; ++i)
		{
			if (W[i] != DST.W[i])
				return false;
		}
		return true;
	}
	void clear() {
#ifdef USE_AVX
		for (int i = 0; i < AVX_W; ++i)
		{
			X[i] = _mm256_setzero_si256();
		}
#else 
		for (int i = 0; i < INT64_W; ++i)
		{
			W[i] = 0;
		}
#endif

	}
	void setones() {
		*this = Mask;
	}

	inline int count() const {
		int result = 0;

		for (int i = 0; i < INT64_W; ++i)
		{
			result += (int)__builtin_popcountll(W[i]);
		}
		return result;
	}

	bool Intersects(const NumberSet& B)
	{
#ifdef USE_AVX

		for (int i = 0; i < AVX_W; ++i)
		{
			if (!all_zeros(_mm256_and_si256(X[i], B.X[i])))
				return true;
		}
#else

		for (int i = 0; i < INT64_W; ++i)
		{
			if ((W[i] & B.W[i]) != 0)
				return true;
		}

#endif
		return false;
	}
	bool Contains(int N) {
		uint64_t element = N / 64;
		uint64_t pos = 63 - (N % 64);
		return (((W[element] >> pos) & 1) == 1);
	}


	void _and(const NumberSet& b) {

#ifdef USE_AVX
		for (int IIN = 0; IIN < AVX_W; ++IIN)
		{
			X[IIN] = _mm256_and_si256(X[IIN], b.X[IIN]);
		}
#else
		for (int IIN = 0; IIN < INT64_W; ++IIN)
		{
			W[IIN] &= b.W[IIN];
		}
#endif
	}

	void _or(const NumberSet& b) {
#ifdef USE_AVX

		for (int IIN = 0; IIN < AVX_W; ++IIN)
		{
			X[IIN] = _mm256_or_si256(X[IIN], b.X[IIN]);
		}
#else

		for (int IIN = 0; IIN < INT64_W; ++IIN)
		{
			W[IIN] |= b.W[IIN];
		}
#endif

	}

	int CountIntersects(const NumberSet& B)
	{
		int count = 0;
		union {
			__m256i cross;
			uint64_t cW[4];
		};

#ifdef USE_AVX
		for (int i = 0; i < AVX_W; ++i)
		{
			cross = _mm256_and_si256(X[i], B.X[i]);
			count += (int)(__builtin_popcountll(cW[0]) + __builtin_popcountll(cW[1]) + __builtin_popcountll(cW[2]) + __builtin_popcountll(cW[3]));
		}
#else
		for (int i = 0; i < INT64_W; ++i)
		{
			count += (int)__builtin_popcountll(W[i] & B.W[i]);
		}
#endif
		return count;
	}
	inline bool Disjoint(const NumberSet& B) { return CountIntersects(B) == 0; }
	bool isFullyContainedIn(const NumberSet& B)const {
		for (int u = 0; u < INT64_W; ++u)
		{
			if ((W[u] & B.W[u]) != W[u])
				return false;
		}
		return true;
	}
	static NumberSet CreateMask(int CountNumbers) {
		NumberSet n;
		for (size_t i = 0; i < CountNumbers; ++i)
			n.set(i); //Ugly but it works...
		return n;
	}
}
;
ALIGN NumberSet NumberSet::Mask;
ostream& operator<<(ostream& os, const NumberSet& m)
{
	for (int i = 0; i < INT64_W; ++i)
	{
		os << bitset<64>(m.W[i]);
	}
	return os;
};

//Lookup of valid targets of each move. Used on Orphans numbers mutator
vector<NumberSet> NEIGHBOURS;

//On Move class I store the srcVal and destVal, I use it to keep coherence on moves.
//If I try to make a move and srcVal or destVal is not the same on the current Grid I don't do the move.
//This ensures that a move Mutation doesn't alter further moves
struct Move {
	int16_t srcIDX;
	int16_t destIDX;
	int16_t srcVal;
	int16_t destVal;
	uint8_t dir;
	uint8_t sign;
	Move() {}

	Move(int _srcID, int _destID, int _dir, int _sign, int _srcVal, int _destVal) {
		srcIDX = _srcID;
		destIDX = _destID;
		dir = _dir;
		sign = _sign;
		srcVal = _srcVal;
		destVal = _destVal;
	}

	//Hash isn't stored on the class
	size_t CreateHash() const {
		size_t hash = (size_t)srcIDX; hash <<= 14;
		hash += (size_t)destIDX; hash <<= 14;
		hash += (size_t)dir; hash <<= 1;
		hash += (size_t)sign; hash <<= 3;
		hash += (size_t)srcVal; hash <<= 14;
		hash += (size_t)destVal;
		return hash;
	}
	/*inline bool operator==(const Move& rhs) {
	return (srcIDX == rhs.srcIDX) &&
	(destIDX == rhs.destIDX) &&
	(srcVal == rhs.srcVal) &&
	(destVal == rhs.destVal) &&
	(dir == rhs.dir) &&
	(sign == rhs.sign);
	}
	inline bool operator!=(const Move& rhs) {
	return (srcIDX != rhs.srcIDX) ||
	(destIDX != rhs.destIDX) ||
	(srcVal != rhs.srcVal) ||
	(destVal != rhs.destVal) ||
	(dir != rhs.dir) ||
	(sign != rhs.sign);
	}*/

};
ostream& operator<<(ostream& os, const Move& m)
{
	os << (int)NMB[m.srcIDX].X << ' ' << (int)NMB[m.srcIDX].Y << ' ' << strDIRS[m.dir] << " " << (m.sign <= 0 ? "-" : "+");
	return os;
};

//Lookup table for doMerge function. if srcIDX and destIDX can't combine, you get a nullptr.
Move* CACHE_MERGE[INVALID_ID + 1][INVALID_ID + 1];
inline void Add(vector<Move>&A, const vector<Move>& B)
{
	for (auto&b : B)
		A.push_back(b);
	//A.insert(A.end(), B.begin(),B.end());
}

#ifdef SHUFFLE_CODE
struct ALIGN ShuffleStrat {
	list<Move> linear;
	inline void clear() {
		linear.resize(0);
	}
};
#endif

//***************** An Strategy it's just a list of related Moves, a set of Moves independent to another Strategies. **************/
struct ALIGN Strategy {
	vector<Move> linear;
	inline void addMove(const Move& m) {
		linear.push_back(m);
	}
	inline void resize(const int& N) {
		linear.resize(N);
	}
	inline void addPlan(const Strategy& S) {
		Add(linear, S.linear);
	}
	//The list of moves ends with 0. Then it's not an Incomplete list of moves.
	inline bool isZeroSumSet() {
		if (linear.size() == 0)
			return false;
		auto& m = linear.back();
		return (m.srcVal == m.destVal && m.sign == 0);
	}
	inline void clear() {
		linear.resize(0);
	}
	//Do all moves, it also tries to Merge unused Numbers with some probability
	void ApplyMoves(Grid& g, const uint32_t& prob1000, const uint32_t& mark, bool dontPlayEndpoint, vector<int>& UnusedNumbers, Random& rnd);
};
Strategy SolvedArray;//Keeps solution

					 //***************** Current Game State. I don't store coordinates but just values **************/
struct Grid {
	uint16_t Val[MAX_NUMBERS];
	Strategy Plan_NMB[MAX_NUMBERS]; //vector<> was slower than that. Maybe that wasn't true on higher levels, but it was on 300-400 levels
	int totalPoints = 0; //val
	int squaredPoints = 0; //val*val
	int totalNumbers = 0; //val==0?1:0
	void clear() {
		for (int i = 0; i < COUNT_NUMBERS; ++i)
		{
			Plan_NMB[i].clear();
		}
		totalPoints = MAX_NUMBERS - 1;
		squaredPoints = 0;
		totalNumbers = MAX_NUMBERS - 1;
	}

	//I've had many problems with wrong solutions that was trying to do moves with zero srcVal or destVal, I double check that now...
	bool ValidateMove(const Move& m) {
		if (m.srcIDX == m.destIDX || m.srcIDX < 0 || m.destIDX < 0 || m.srcIDX >= COUNT_NUMBERS || m.destIDX >= COUNT_NUMBERS)
			return false;
		if (NMB[m.srcIDX].X != NMB[m.destIDX].X && NMB[m.srcIDX].Y != NMB[m.destIDX].Y)
			return false;

		/*		if (pMOVE == nullptr || pMOVE->srcVal != m.srcVal)
		return false;
		*/
		int dist = max(abs((int)NMB[m.srcIDX].X - (int)NMB[m.destIDX].X), abs((int)NMB[m.srcIDX].Y - (int)NMB[m.destIDX].Y));
		if (dist != m.srcVal || dist == 0 || dist > W)
			return false;
		if (Val[m.srcIDX] == 0 || Val[m.destIDX] == 0 || m.srcVal == 0 || m.srcVal >= W || m.destVal == 0)
		{
			return false;
		}
		return true;
	}

	//Apply Move keeps coherence. If the Grid gamestate doesn't match with both srcVal and dstVal it's ignored.
	//My solver is very performance dependent, the faster the move generation, the better.
	inline bool ApplyMove(const Move& m)
	{
		//Coherence checks.This check shouldn't be needed but I might have some bug on move Mutations
		if (m.srcIDX == m.destIDX || m.srcIDX < 0 || m.destIDX < 0 || m.srcIDX >= COUNT_NUMBERS || m.destIDX >= COUNT_NUMBERS)
		{
			return false;
		}
		if (Val[m.srcIDX] == 0 || Val[m.destIDX] == 0)
			return false;
		auto srcScore = Val[m.srcIDX];
		if (m.srcVal != srcScore)
			return false;
		auto destScore = Val[m.destIDX];
		if (destScore != m.destVal)
			return false;

		//If it's OK, I apply the move
		--totalNumbers;
		Val[m.srcIDX] = 0;
		auto Dest = (m.sign == 0 ? abs(destScore - srcScore) : destScore + srcScore);
		Val[m.destIDX] = Dest;
		if (Dest == 0)
		{
			--totalNumbers;
		}
		//Strategy update. I transfer all moves from src Strategy to dest Strategy.
		//This is the fastest way I get to keep track movelists.
		auto& SRC = Plan_NMB[m.srcIDX];
		auto& DST = Plan_NMB[m.destIDX];
		if (SRC.linear.size() != 0) {
			DST.addPlan(SRC);
			SRC.clear();
		}
		DST.addMove(m);
		return true;
	}

#ifdef SHUFFLE_CODE
	//This is an alternative ApplyMove, it does the move but it tries to alter the order of the moves without breaking the final result.
	//The goal of that is to have different intermediate values on each node.
	inline void ShuffleMove(const Move& m, ShuffleStrat Shuffle_NMB[MAX_NUMBERS], Random& rnd)
	{
		if (m.srcIDX == m.destIDX)
			return;
		if (m.srcIDX == INVALID_ID || m.destIDX == INVALID_ID)
			return;
		auto srcScore = Val[m.srcIDX];
		if (m.srcVal != srcScore) return; //incoherent
		auto destScore = Val[m.destIDX];
		if (destScore == 0 || destScore != m.destVal) return; //Incoherent

		--totalNumbers;
		Val[m.srcIDX] = 0;
		auto Dest = (m.sign == 0 ? abs(destScore - srcScore) : destScore + srcScore);
		Val[m.destIDX] = Dest;
		if (Dest == 0)
			--totalNumbers;
		{
			auto& SRC = Shuffle_NMB[m.srcIDX];
			auto& DST = Shuffle_NMB[m.destIDX];
			if (SRC.linear.size() == 0)
			{
			}
			else if (DST.linear.size() == 0)
			{
				swap(Shuffle_NMB[m.srcIDX], Shuffle_NMB[m.destIDX]);
			}
			else {
				int inject = rnd.NextInt(max(1, (int)DST.linear.size() / 2));
				auto pos = DST.linear.begin();
				for (auto&l : SRC.linear)
				{
					for (int j = 0; j < inject; ++j)
						if (pos != DST.linear.end())
							pos++;
					DST.linear.insert(pos, l);
					inject = rnd.NextInt(1, 3);
				}
				SRC.linear.clear();
			}
			//To do: go back in time to find a proper value.
			if (Dest == 0)
			{
				//End, it must be the last
				DST.linear.push_back(m);
			}
			else {
				int countDown = 0;
				int changeValue = -1;
				for (auto it = DST.linear.crbegin(); it != DST.linear.crend(); ++it)
				{ //Backtrack moves to see how far I can send that move
				  //Independent moves, so I can pass it.
					if (it->srcIDX != m.srcIDX && it->destIDX != m.srcIDX && it->srcIDX != m.destIDX && it->destIDX != m.destIDX)
					{
						++countDown;
						continue;
					}

					//That will affect result
					if (it->srcIDX == m.srcIDX || it->destIDX == m.srcIDX || it->srcIDX == m.destIDX)
					{
						break;
					}
					if (it->destIDX == m.destIDX && (m.srcVal >= m.destVal || it->srcVal >= it->destVal))
					{
						break;
					}
					//Dependent moves, I may pass it or not.
					if (m.sign == 1)
					{

						++countDown;
						changeValue = countDown;
						continue;
					}
					if (it->sign == 1)
					{
						if (it->destVal <= m.srcVal)
							break;
						else {
							++countDown;
							changeValue = countDown;
							continue;
						}
					}
					//both signs 0
					if (it->destVal <= m.srcVal || it->srcVal <= m.destVal)
						break;

					++countDown;
				}
				if (countDown == 0)
				{
					DST.linear.push_back(m);
				}
				else {
					//Shuffling a move
					int randomStop = rnd.NextInt(countDown);
					int nStart = changeValue + 1;
					int nEnd = countDown - 1;
					if (changeValue >= 0 && nStart <= nEnd && rnd.NextInt(1000) < 600) {
						randomStop = rnd.NextInt(nStart, nEnd);
					}
					auto pos = DST.linear.end();
					for (int i = 0; i < randomStop; ++i)
					{
						pos--;
						if (pos == DST.linear.begin())
							break;
					}
					DST.linear.insert(pos, m);
				}

			}
		}
	}
#endif	

	/*Domain specific property of Numbershift's subtract
	If you have srcVal = 2*destVal, the subtract is
	abs(destVal-srcVal)=abs(destVal-2*destVal)=abs(-destVal)=destVal
	, so destVal is unaffected and this won't invalidate any move due to coherence
	I abuse it to inject an Incomplete Strategy on another one, and it doesn't affect anything
	*/
	inline bool doMerge(const int& srcID, const int& destID, Random& rnd) {
		if (srcID == destID || srcID >= COUNT_NUMBERS || destID >= COUNT_NUMBERS || srcID < 0 || destID < 0)
			return false;
		//Lookup to check if a merge can be done or not
		Move* lm = CACHE_MERGE[srcID][destID];
		if (lm == nullptr)return false;
		if (Val[srcID] != lm->srcVal || Val[srcID] == 0) return false;
		auto targetVAL = Val[destID];
		if ((targetVAL == 0) || (2 * targetVAL != Val[srcID])) return false;
		Move m = *lm; //Move(srcID, destID, d, 0, srcVal, targetVAL);
		m.destVal = targetVAL;//To pass coherence.
		ApplyMove(m);
		return true;
	}
	//This is the same function but with a hash check to avoid some valid merges.
	inline bool doMerge(const int& srcID, const int& destID, Random& rnd, unordered_set<size_t>& Forbidden) {
		if (srcID == destID || srcID >= COUNT_NUMBERS || destID >= COUNT_NUMBERS || srcID < 0 || destID < 0)
			return false;
		Move* lm = CACHE_MERGE[srcID][destID];
		if (lm == nullptr)return false;
		auto srcVal = Val[srcID];
		if (srcVal != lm->srcVal || srcVal == 0) return false;
		auto targetVAL = Val[destID];
		if ((targetVAL == 0) || (2 * targetVAL != srcVal)) return false;

		Move m = *lm; //Move(srcID, destID, d, 0, srcVal, targetVAL);
		m.destVal = targetVAL;
		if (Forbidden.find(m.CreateHash()) == Forbidden.end())
		{
			ApplyMove(m);
			return true;
		}
		return false;
	}
	//Destructive join, this will change destValue, so coherence will break further moves of that Strategy
	inline bool doForceJoin(const int& srcID, const int& destID, uint32_t randomChance, Random& rnd) {
		if (srcID == destID || srcID >= COUNT_NUMBERS || destID >= COUNT_NUMBERS || srcID < 0 || destID < 0)
			return false;
		Move* lm = CACHE_MERGE[srcID][destID];
		if (lm == nullptr)return false;
		if (Val[srcID] != lm->srcVal || Val[srcID] == 0) return false;
		auto targetVAL = Val[destID];
		if (targetVAL == 0) return false;
		if (rnd.NextInt(1000) >= randomChance)
			return false;
		Move m = *lm; //Move(srcID, destID, d, 0, srcVal, targetVAL);
		m.destVal = targetVAL;
		ApplyMove(m);
		return true;
	}

	//Create a list of Moves for a number. All subtract moves and some random Add moves
	inline void ExplodeMoves(const int& srcID, vector<Move>& moves, Random* prnd)
	{
		if (srcID < 0 || srcID >= COUNT_NUMBERS)
			return;
		auto srcVal = Val[srcID];
		if (srcVal > W || srcVal == 0) return;

		auto srcX = NMB[srcID].X;
		auto srcY = NMB[srcID].Y;

		for (int d = 0; d < 4; ++d)
		{
			auto destX = srcX + DIR_X[d] * srcVal;
			if (destX < 0 || destX >= W) continue;
			auto destY = srcY + DIR_Y[d] * srcVal;
			if (destY < 0 || destY >= H) continue;
			auto destID = IDX[destX][destY];
			if (destID >= COUNT_NUMBERS || destID < 0)
				continue;
			auto destVal = Val[destID];
			if (destVal == 0) continue;

			Move m = Move(srcID, destID, d, 0, srcVal, destVal);
			moves.push_back(m);
			if (prnd != nullptr && prnd->NextInt(100) < K_ALLOWSUM)
			{
				m = Move(srcID, destID, d, 1, srcVal, destVal);
				moves.push_back(m);
			}
		}
		return;
	}

	//Generate a random valid move. This function should be fast too.
	//It didn't have any "clever" search of moves. On my tests simcount is much better than smart moves.
	//At levels 400-600 I tried to do lookups, I didn't see any real improvement on performance.
	inline bool doRandomMove(Random& rnd) {
		int srcID = rnd.NextInt(COUNT_NUMBERS);
		uint32_t ren = (uint32_t)rnd.xrandom();
		for (int NN = 0; NN < COUNT_NUMBERS; ++NN)
		{
			++srcID;
			if (srcID >= COUNT_NUMBERS)
				srcID = 0;
			auto& srcVal = Val[srcID];
			if (srcVal >= W || srcVal == 0)
				continue;

			auto srcX = NMB[srcID].X;
			auto srcY = NMB[srcID].Y;

			int d = ren & 3;
			ren >>= 2;
			if (ren == 0) //Yeah ren check to make another random isn't a good practice, but it works enough.
				ren = (uint32_t)rnd.xrandom();
			for (int stard = 0; stard < 4; ++stard)
			{
				++d;
				if (d >= 4)
					d = 0;
				auto destX = srcX + DIR_X[d] * srcVal;
				if (destX < 0 || destX >= W) continue;
				auto destY = srcY + DIR_Y[d] * srcVal;
				if (destY < 0 || destY >= H) continue;
				auto destID = IDX[destX][destY];
				if (destID == INVALID_ID)
					continue;
				auto& destVal = Val[destID];
				if (destVal == 0) continue;
				int signo = rnd.NextInt(100) < K_ALLOWSUM ? 1 : 0;
				Move m = Move(srcID, destID, d, signo, srcVal, destVal);
				ApplyMove(m);
				return true;
			}
		}
		return false;
	}
} ALIGN;
ostream& operator<<(ostream& os, const Grid& g)
{
	int width;
	if (g.totalPoints >= 1000) width = 4;
	else if (g.totalPoints >= 100) width = 3;
	else if (g.totalPoints >= 10) width = 2;
	else width = 1;


	os << W << " " << H << endl;
	for (int y = 0; y < H; y++) {
		for (int x = 0; x < W; x++) {
			int idx = IDX[x][y];
			os << std::setw(width) << (idx == INVALID_ID || g.Val[idx] == 0 ? "." : to_string(g.Val[idx])) << " ";
		}
		os << endl;
	}
	return os;
};
ALIGN Grid MAP;//Here I keep the original Grid. It's used on clear() to restore the initial value.

void Strategy::ApplyMoves(Grid& g, const uint32_t& prob1000, const uint32_t& mark, bool dontPlayEndpoint, vector<int>& UnusedNumbers, Random& rnd) {
	if (linear.size() == 0)
		return;
	int endIDX = linear.back().destIDX;
	for (int iii = 0; iii < linear.size(); ++iii)
	{
		auto&m = linear[iii];
		//Try to add optionals as srcID
		if (UnusedNumbers.size() > 0 && rnd.NextInt(1000) < prob1000)
			for (auto& unused : UnusedNumbers)
			{
				g.doMerge(unused, m.srcIDX, rnd);
			}
		if ((m.destIDX == endIDX) && dontPlayEndpoint)
			continue;
		g.ApplyMove(m);
		if (UnusedNumbers.size() > 0 && rnd.NextInt(1000) < prob1000)
			for (auto& unused : UnusedNumbers)
			{
				g.doMerge(unused, m.destIDX, rnd);
			}
	}
}

/************** LAHC_Node keeps a Grid, a list of Strategies, RN NumberSet and a Score *****************************/
struct ALIGN LAHC_Node {
	NumberSet remainingNumbers; //Numbers != 0. Critical for GLOBAL_RN.APX algorithm
	NumberSet usedInMovesNumbers; //All srcIDX and destIDX, may be 0 or not. Not that important
	NumberSet untouchedNumbers; //Numbers not used on any move. Not that important
	NumberSet completedNumbers; //Numbers from complete Lists. Not that important

	double BestScore = 0.0;
	double A, B, C, D, ExtraC, ExtraD;
	Grid grid;
	//Indexes related to grid.Plan_NMB[]
	vector<int> Movelists; //All moves
	vector<int> ZerosumSets; //I end on val=X-X=0, it's a complete set of moves
	vector<int> IncompleteSets; //Ending on val!=0

	int MutateType = 9; //Helps calculating probabilities of what Mutator is working better. 9 means pure random
	bool exhausted = false; //For Mutator_Exhaustive_6. It's not random and it's costly, so I did only once.

	void CalcStats() {
		//Recalculate lists, totalPoints, used rows/columns, etc.
		exhausted = false;
		Movelists.resize(0);
		ZerosumSets.resize(0);
		IncompleteSets.resize(0);
		BestScore = 0.0;
		A = 0.0; B = 0.0; C = 0.0; D = 0.0; ExtraC = 0.0; ExtraD = 0.0;

		remainingNumbers.clear();
		usedInMovesNumbers.clear();
		untouchedNumbers.clear();
		completedNumbers.clear();

		int extraPoints = 0;
		int extraSQ = 0;
		grid.totalNumbers = 0;
		grid.totalPoints = 0;
		grid.squaredPoints = 0;

		uint8_t X_CELLS[MAX_W] = { 0 };
		uint8_t Y_CELLS[MAX_H] = { 0 };

		for (auto& n : NMB)
		{
			auto idx = n.ID;
			auto& V = grid.Val[n.ID];
			if (V > 0)
			{
				remainingNumbers.set(n.ID);

				++X_CELLS[n.X];
				++Y_CELLS[n.Y];
				++grid.totalNumbers;
				grid.totalPoints += V;
				grid.squaredPoints += (V*V);
				if (V <= LIM_P)
				{
					extraPoints += 1;
				}
				else extraPoints += V;
				if (V <= LIM_SQ)
				{
					extraSQ += 1;
				}
				else extraSQ += (V*V);
			}

			if (grid.Plan_NMB[idx].linear.size() > 0)
			{

				Movelists.push_back(idx);
				//int Index = (int)Movelists.size() - 1;
				auto& plan = grid.Plan_NMB[idx];
				for (auto& P : plan.linear)
				{
					usedInMovesNumbers.set(P.srcIDX);
				}
				usedInMovesNumbers.set(plan.linear[plan.linear.size() - 1].destIDX);
				if (plan.isZeroSumSet())
				{
					ZerosumSets.push_back(idx);
					for (auto& P : plan.linear)
					{
						completedNumbers.set(P.srcIDX);
					}
					completedNumbers.set(plan.linear[plan.linear.size() - 1].destIDX);
				}
				else {
					IncompleteSets.push_back(idx);
				}
			}
		}

		if (grid.totalNumbers > K_B_NUM)
		{
			extraPoints = grid.totalPoints;
			extraSQ = grid.squaredPoints;
		}
		untouchedNumbers = usedInMovesNumbers;
		untouchedNumbers.negate();
		NumberSet check = untouchedNumbers;
		check._or(usedInMovesNumbers);

		//Simple simcount. It's a sum of all moves.
		uint64_t S = 0;
		for (auto& l : Movelists)
		{
			S += grid.Plan_NMB[l].linear.size();
		}
		SimCount += S;

		//Score Calculation
		//A is remaining number percentage
		A = K_A * (double)(COUNT_NUMBERS - grid.totalNumbers)*INV_COUNT_NUMBERS;
		//B counts 
		if (K_B != 0.0 && grid.totalNumbers <= K_B_NUM)
		{
			int resX = 0;
			int resY = 0;
			for (int x = 0; x < W; ++x)
			{
				resX += (X_CELLS[x] == 2 ? 4 : (X_CELLS[x] > 0 ? 10 : 0));
			}
			for (int y = 0; y < H; ++y)
			{
				resY += (Y_CELLS[y] == 2 ? 4 : (Y_CELLS[y] > 0 ? 10 : 0));
			}

			B = K_B * (1.0 - (((double)resX)* INV_RESX + ((double)resY)* INV_RESY));
		}
		else B = 0.0;
		//C is remaining totalPoints percentage
		C = K_C * 0.0; //TODO: (MAP.??? - grid.???) * normalization;
		//ExtraC is to allow that all numbers <= LIM_P scores the same in C. I think val=1 must score the same as val=3, because on low values it's not true the C scoring.
		ExtraC = K_C * (double)(MAP.totalPoints - extraPoints)* INV_POINTS - C;

		//D is remaining squaredPoints points percentage
		D = K_D *  0.0; //TODO: (MAP.??? - grid.???) * normalization;
		//Similar idea than ExtraC
		ExtraD = K_D * (double)(MAP.squaredPoints - extraSQ)*INV_SQUAREPOINTS - D;
		
		//TODO: This is incorrect
		BestScore = A*B*C*D + max(0.0, ExtraC) + max(0.0, ExtraD);
	}
	LAHC_Node() {}

	//Resets to starting position. 
	void clear() {
		exhausted = false;
		grid.clear();
		//Copy MAP's Val array, but not MAP.Plan_NMB[]
		copy(MAP.Val, MAP.Val + MAX_NUMBERS, grid.Val);

		//I avoid clear(), because it frees memory
		Movelists.resize(0);
		ZerosumSets.resize(0);
		IncompleteSets.resize(0);
		A = 0.0;
		B = 0.0;
		C = 0.0;
		D = 0.0;
		ExtraC = 0.0;
		ExtraD = 0.0;
		BestScore = 0.0;
		MutateType = 9;
	}

	//Saving a good approximation to a file. This is used to restart the work later, also to send it to the cloud for Recombination purposes.
	void SaveApprox(string s)const {
		//mutexAPROX.lock();
		std::ofstream outfile;
		outfile.open(s, std::ios_base::app); // append instead of overwrite
		for (auto& L : Movelists)
		{
			for (auto&m : grid.Plan_NMB[L].linear)
			{
				//srcID, destID, dir, sign, srcVal, destScore
				outfile << (int)m.srcIDX << " " << (int)m.destIDX << " " << (int)m.dir << " " << (int)m.sign << " " << (int)m.srcVal << " " << (int)m.destVal << "|";
			}
		}
		outfile << endl;
		outfile.close();
		filesize = GetFileSize(s);
		//mutexAPROX.unlock();
	}

	//Saving to a file with a more human friendly format. If it was solved it also writes a solution text at the end of the file to have a backup of solution.txt
	//It also writes some metadata, like the computer it solves it, some ExtraTest, Mutator's Probabilities Array, etc..
	void SaveToFile(string s, long long maxTimeImprovement = -1, int currLFA = -1, double PROB[] = nullptr, string ExtraText = "")
	{
		std::sort(std::begin(ZerosumSets), std::end(ZerosumSets),
			[=](int a, int b) {return grid.Plan_NMB[a].linear.size() > grid.Plan_NMB[b].linear.size(); });
		std::sort(std::begin(IncompleteSets), std::end(IncompleteSets),
			[=](int a, int b) {return grid.Plan_NMB[a].linear.size() > grid.Plan_NMB[b].linear.size(); });

		double points = (double)(MAP.totalPoints - grid.totalPoints)*100.0 / (double)MAP.totalPoints;
		double coverage = 100.0*(double)(COUNT_NUMBERS - grid.totalNumbers)*INV_COUNT_NUMBERS;
		//This mutex avoids two threads writing a file at the same time, thus corruption both solutions.
		mutexSAVE.lock();
		std::ofstream outfile;
		outfile.open(s, std::ios_base::app); // append instead of overwrite
											 //Print metadata - Useful for analysis.
		outfile << "********** Binary:" << PROGRAM_NAME << " HOST:" << COMPUTER_NAME << " Ellapsed Time:" << stopwatch.EllapsedMilliseconds() << "ms C/I" << ZerosumSets.size() << "/" << IncompleteSets.size() << endl;
		outfile << ExtraText << endl;
		outfile << " CurrLfa/SIZE_LFA:" << currLFA << "/" << SIZE_LFA << " K_ABCD:" << (int)K_A << "," << (int)K_B << "," << (int)K_C << "," << (int)K_D << " LIMIT_TIME_IMPROVEMENT:" << LIMIT_TIME_IMPROVEMENT;
		if (maxTimeImprovement > 0) outfile << " maxTimeImprovement:" << (maxTimeImprovement / 1000) << "s ";

		outfile << " Performance:" << performance << "s/found(" << newFound << ")";
		if (PROB != nullptr) outfile << "PROB:" << (int)PROB[0] << "," << (int)PROB[1] << "," << (int)PROB[2] << "," << (int)PROB[3] << "," << (int)PROB[4] << "," << (int)PROB[5];
		outfile << endl;
		outfile << " Best:" << (int)A << "," << (int)B << "," << (int)C << "," << (int)D << (int)ExtraC << "," << (int)ExtraD;
		outfile << " Points:" << grid.totalPoints << "/" << MAP.totalPoints << "=" << (double)(MAP.totalPoints - grid.totalPoints)*100.0 / (double)MAP.totalPoints;
		outfile << "% Lists C/U:" << ZerosumSets.size() << "/" << Movelists.size();
		outfile << " Numbers R / T:" << (grid.totalNumbers) << " / " << MAP.totalNumbers;
		outfile << " BC/U/T:" << completedNumbers.count() << "/" << usedInMovesNumbers.count() << "/" << COUNT_NUMBERS << " " << 100.0*(double)completedNumbers.count() / (double)COUNT_NUMBERS << "% ";

		//Grid state for incomplete APX
		if (grid.totalPoints != 0)
		{
			outfile << " GRID IS:" << endl;
			outfile << grid;
		}
		else outfile << endl;
		//Move Subsets, first completes then incompletes.
		for (auto& b : ZerosumSets)
		{
			outfile << "Complete:";
			outfile << grid.Plan_NMB[b].linear.size() << " |";
			for (auto&m : grid.Plan_NMB[b].linear)
			{
				outfile << m << "|";
			}
			outfile << endl;
		}
		for (auto& b : IncompleteSets)
		{
			outfile << "Incomplete:";
			outfile << grid.Plan_NMB[b].linear.size() << " |";
			for (auto&m : grid.Plan_NMB[b].linear)
			{
				outfile << m << "|";
			}
			outfile << endl;
		}
		if (grid.totalPoints == 0 && Movelists.size() > 0)
		{
			//Print the solution in Codingame's format.
			outfile << passwordLevel << endl;
			for (auto&b : Movelists)
				for (auto& m : grid.Plan_NMB[b].linear)
				{
					outfile << m << endl;
				}
		}
		outfile.close();
		mutexSAVE.unlock();

	}
};

//SOLVING checker. If totalNumbers==0, redo all moves to ensure they are all valid, then save the previous GameState + Solution to a file
// Previous GameState is just for analysis purposes, to know how it evolves on the last Mutation
bool IS_SOLVED(const LAHC_Node& previous, const LAHC_Node& candidate, string TEXT_)
{
	if (!solved && candidate.Movelists.size() > 0 && candidate.grid.totalPoints == 0 && candidate.grid.totalNumbers == 0)
	{
		cerr << TEXT_ << " SOLVED!" << endl;

		mutexSOL.lock();
		//Validate Solution
		LAHC_Node checkSol;
		checkSol.clear();
		SolvedArray.clear();
		for (auto&N : candidate.Movelists)
		{
			for (auto&m : candidate.grid.Plan_NMB[N].linear)
			{
				if (checkSol.grid.ValidateMove(m) //Double check, I've had wrong Moves at some point. ApplyMove ignored them but not CG referee
					&& checkSol.grid.ApplyMove(m))
				{
					SolvedArray.addMove(m); //Only AAA quality moves here
				}
			}
		}
		checkSol.MutateType = 9;
		checkSol.CalcStats(); //Needed to resfresh totalnumbers, totalPoints and BestScore
		if (checkSol.grid.totalNumbers == 0 && checkSol.grid.totalPoints == 0 && checkSol.BestScore > 0.0)
		{ //All looking good, save 
			string filename = "SOLUTION_" + to_string(level) + "_" + passwordLevel + ".txt";
			LAHC_Node p = previous;
			p.SaveToFile(filename, 0, 0, nullptr, TEXT_);
			checkSol.SaveToFile(filename, 0, 0, nullptr, TEXT_);
			solved = true;
			/*			cerr << "SOLUTION had " << candidate.Movelists.size() << " GROUPS" << endl;
			cerr << passwordLevel << endl;
			for (auto&m : SolvedArray.linear)
			{
			cerr << m << endl;
			}*/
			mutexSOL.unlock();
			return true;
		}
		else {
			cerr << "UPS! Wrong solution! Numbers:" << candidate.grid.totalNumbers << "<->" << checkSol.grid.totalNumbers << " Points" << candidate.grid.totalPoints << "<->" << checkSol.grid.totalPoints << endl;
			SolvedArray.clear();
		}
		mutexSOL.unlock();
	}
	return false;
}


//************** See GLOBAL_RN.APX Search Algorithm at the start of the code*************************************
const int LIMITE_GLOBAL = 1000;

//GLOBAL BEST vars
int swapDone = 0;
const int MAX_INSERT_NUMBERS = 7;
const int MAX_DEGREE_TOSAVE = 4; //To save in a file
const int MIN_DEGREE_TOINSERT = 7;//To start using RN_Explorer Algorithm
const int K_EXHAUSTIVE_MIN_NUMBERS = MIN_DEGREE_TOINSERT;

atomic<int> totalGLOBAL;
struct RN_Explorer {
	vector<vector<LAHC_Node>> APX;
	vector<int> Degree;
	vector<int> SizeAcum;
	int totalRN;
	bool isGLOBAL; //To use mutexes
	RN_Explorer(bool _isGLOBAL) {
		isGLOBAL = _isGLOBAL;
		clear();
	}
	void Calc_Sizes() {
		if (APX.size() == 0)
		{
			totalRN = 0;
			if (isGLOBAL)
				totalGLOBAL = 0;
			Degree.resize(1);
			Degree[0] = COUNT_NUMBERS + 1;
			SizeAcum.resize(1);
			SizeAcum[0] = 1;
			return;
		}

		SizeAcum.resize(APX.size());
		Degree.resize(APX.size());
		for (int i = 0; i< APX.size(); ++i)
		{
			int S = (int)APX[i].size();
			if (S > 0)
				Degree[i] = APX[i][0].grid.totalNumbers;
			else Degree[i] = COUNT_NUMBERS + 1;
			SizeAcum[i] = S;
			if (i > 0)
				SizeAcum[i] += SizeAcum[i - 1];
		}
		totalRN = SizeAcum.back();
		if (isGLOBAL)
			totalGLOBAL = totalRN;
	}

	void clear() {
		APX.resize(0);
		totalRN = 0;
		if (isGLOBAL)
			totalGLOBAL = 0;
		Degree.resize(1);
		Degree[0] = COUNT_NUMBERS + 1;
		SizeAcum.resize(1);
		SizeAcum[0] = 0;
	}

	void packDegrees() {
		if (APX.size()>0)
			for (int i = (int)APX.size() - 1; i >= 0; --i)
			{
				if (APX[i].empty())
				{
					APX.erase(APX.begin() + i);
				}
			}
	}

	bool RN_Insert(LAHC_Node& newBest, bool saveToFile)
	{
		if (newBest.grid.totalNumbers == 0)
		{ //Solved? Someone will took care of it
			return false;
		}
		if (isGLOBAL)
			mutexRN_Explorer.lock();
		{
			int mLimit = MIN_DEGREE_TOINSERT;
			if (APX.size()>0 && APX[0].size()>0)
			{
				mLimit = max(MIN_DEGREE_TOINSERT, APX[0][0].grid.totalNumbers + 2);
			}
			if (newBest.grid.totalNumbers > mLimit)
			{ //Not enough
				if (isGLOBAL)
					mutexRN_Explorer.unlock();
				return false;
			}
		}

		if (newBest.grid.totalNumbers > Degree.back() && APX.size() >= MAX_INSERT_NUMBERS)
		{
			if (isGLOBAL)
				mutexRN_Explorer.unlock();
			return false;
		}

		//Equals. Doesn't change anything else
		for (auto&A : APX)
			if (A.size() > 0 && A[0].grid.totalNumbers == newBest.grid.totalNumbers)
			{
				for (auto&b : A)
				{
					if (b.remainingNumbers.Equals(newBest.remainingNumbers))
					{
						if (newBest.BestScore > b.BestScore)
						{
							bool keepExhausted = b.exhausted;
							b = newBest;
							b.exhausted = keepExhausted;
							if (newBest.grid.totalNumbers <= MAX_DEGREE_TOSAVE)
								++swapDone;
						}
						if (isGLOBAL)
							mutexRN_Explorer.unlock();
						return false;
					}
				}
				break;
			}

		for (auto&A : APX)
		{
			if (A.size() > 0)
			{
				//Lower Quality - Has a subset.
				if (A[0].grid.totalNumbers < newBest.grid.totalNumbers)
				{
					for (auto&b : A)
					{
						if (b.remainingNumbers.isFullyContainedIn(newBest.remainingNumbers))
						{
							//Exists a better option
							if (isGLOBAL)
								mutexRN_Explorer.unlock();
							return false;
						}
					}
				}
				else //Superset
					if (A[0].grid.totalNumbers > newBest.grid.totalNumbers)
					{
						for (int b = (int)A.size() - 1; b >= 0; --b)
						{
							if (newBest.remainingNumbers.isFullyContainedIn(A[b].remainingNumbers))
							{
								if (b != A.size() - 1)
									A[b] = A.back();
								A.pop_back();
							}
						}
					}
			}
		}
		packDegrees();
		//Then we should add it, search 
		int indexDegree = -1; int insertAt = -1;
		if (APX.size()>0)
			for (int i = 0; i < APX.size(); ++i)
			{
				if (APX[i].size() > 0)
				{
					if (APX[i][0].grid.totalNumbers == newBest.grid.totalNumbers)
					{
						indexDegree = i;
						break;
					}
					else if (APX[i][0].grid.totalNumbers > newBest.grid.totalNumbers)
					{
						insertAt = i;
						break;
					}
				}
			}
		if (indexDegree == -1 && insertAt == -1)
		{
			if (APX.size() >= MAX_INSERT_NUMBERS)
			{
				if (isGLOBAL)
					mutexRN_Explorer.unlock();
				return false;
			}
			else {
				APX.resize(APX.size() + 1);
				indexDegree = APX.size() - 1;
			}
		}

		if (insertAt >= 0)
		{
			APX.insert(APX.begin() + insertAt, vector<LAHC_Node>());
			indexDegree = insertAt;
		}



		// new Approx
		++totalGLOBAL;
		APX[indexDegree].push_back(newBest);
		if (saveToFile && newBest.grid.totalNumbers <= MAX_DEGREE_TOSAVE)
		{
			if (PROGRAM_NAME != "./test")
			{
				if (swapDone == 0)
				{
					++swapDone;
				}
			}
		}
		packDegrees();
		Calc_Sizes();
		if (totalRN > LIMITE_GLOBAL && APX.size() > 5)
		{
			APX.pop_back();
			Calc_Sizes();
		}
		if (APX.size() > MAX_INSERT_NUMBERS)
		{
			APX.resize(MAX_INSERT_NUMBERS);
			Calc_Sizes();
		}
		if (isGLOBAL)
			mutexRN_Explorer.unlock();
		return true;
	}

	LAHC_Node* pickWeightedRandom(Random& rnd) {
		if (totalRN == 0 || APX.empty())
			return nullptr;
		//Weighted random select. Prefer lower RN.
		vector<size_t> rndWeight;
		rndWeight.resize(APX.size());
		size_t initialWeight = 350;
		for (int i = 0; i < APX.size(); ++i)
		{
			rndWeight[i] = 0;
			if (i > 0)
			{
				rndWeight[i] += rndWeight[i - 1];
			}
			if (APX[i].size() > 0)
			{
				rndWeight[i] += initialWeight;
				initialWeight -= 70;
				if (initialWeight < 50)
					initialWeight = 50;
			}
		}
		int rndList = -1;
		int rndIndex = -1;

		if (rndWeight.back() > 0)
		{
			int bigRecomb = rnd.NextInt((uint32_t)rndWeight.back());
			rndList = rndLowBound(rndWeight, bigRecomb);
			for (int i = rndList; i >= 0; --i)
				if (APX[rndList].size() >0)
				{
					rndIndex = rnd.NextInt((uint32_t)APX[rndList].size());
					break;
				}
		}
		if (rndIndex >= 0 && rndList >= 0)
			return &APX[rndList][rndIndex];
		return nullptr;
	}
};
RN_Explorer GLOBAL_RN(true);





#ifdef SHUFFLE_CODE
//Mutator that doesn't improve score, but adds diversity by creating different intermediate values.
//That's used to improve Merge chances.
bool Mutate_Shuffle(LAHC_Node& currentBest, LAHC_Node& nodeShuff, Random& rnd) {
	ShuffleStrat Shuffle_NMB[MAX_NUMBERS];
	++countShuffles;
	for (int TRIES = 0; TRIES < SHUFFLE_TRIES; ++TRIES)
	{
		nodeShuff.clear();
		nodeShuff.MutateType = currentBest.MutateType;
		for (int i = 0; i < COUNT_NUMBERS; ++i)
			Shuffle_NMB[i].clear();

		vector<int> ML = currentBest.Movelists;
		do_shuffle<vector<int>>(ML, rnd);
		for (auto& L : ML)
		{
			auto& ST = currentBest.grid.Plan_NMB[L];
			for (auto& m : ST.linear) {
				nodeShuff.grid.ShuffleMove(m, Shuffle_NMB, rnd);
			}
		}

		nodeShuff.clear();
		for (int i = 0; i < COUNT_NUMBERS; ++i)
			if (Shuffle_NMB[i].linear.size() > 0)
			{
				for (auto&m : Shuffle_NMB[i].linear) {
					m.destVal = nodeShuff.grid.Val[m.destIDX];
					nodeShuff.grid.ApplyMove(m);
				}
			}
		nodeShuff.MutateType = currentBest.MutateType;
		nodeShuff.CalcStats();
		if (nodeShuff.BestScore >= currentBest.BestScore
			&& nodeShuff.grid.totalNumbers <= currentBest.grid.totalNumbers
			&& nodeShuff.grid.totalPoints == currentBest.grid.totalPoints
			&& nodeShuff.grid.squaredPoints == currentBest.grid.squaredPoints
			)
		{
			//ACCEPTED! Now pass from list to vector
			for (int i = 0; i < COUNT_NUMBERS; ++i)
			{
				if (nodeShuff.grid.Plan_NMB[i].linear.size() > 0)
				{
					nodeShuff.grid.Plan_NMB[i].linear.resize(0);
				}
				if (Shuffle_NMB[i].linear.size() > 0)
				{
					for (auto&m : Shuffle_NMB[i].linear) {
						nodeShuff.grid.Plan_NMB[i].linear.push_back(m);
					}
				}
			}
			nodeShuff.MutateType = currentBest.MutateType;
			nodeShuff.CalcStats(); //Simcount are valid because Shuffle didn't count sim.
			if (nodeShuff.BestScore >= currentBest.BestScore && nodeShuff.grid.totalNumbers <= currentBest.grid.totalNumbers)
			{
				++correctShuffles;
				return true;
			}
		}
	}
	return false;
}
#endif
/************************************************************************************/
/************************     MUTATORS - CHANGE/REPLACES MOVES **********************/
/************************************************************************************/

//One of the first mutators I did, I think it's not really useful
inline void Mutate_Completes_0(LAHC_Node& currentBest, LAHC_Node& candidate, Random& rnd)
{
	if (currentBest.ZerosumSets.size() > 0) {
		candidate.MutateType = 0;

		bool tryMerge = currentBest.Movelists.size() > 8 && (currentBest.usedInMovesNumbers.count() > COUNT_NUMBERS * 70 / 100) && rnd.NextInt(100) < 60;
		vector<int> UnusedNumbers;
		uint64_t maskX = 0;
		uint64_t maskY = 0;
		auto Z = currentBest.ZerosumSets;

		NumberSet UsedNumbers;
		for (auto&z : Z)
		{
			if (currentBest.grid.Plan_NMB[z].linear.size() > 1 || (rnd.NextInt(1000) < 840))
				for (auto&l : currentBest.grid.Plan_NMB[z].linear)
				{
					UsedNumbers.set(l.srcIDX);
					UsedNumbers.set(l.destIDX);
				}
		}
		bool withoutEndpoint = rnd.NextBool();
		uint32_t chanceMerge = (uint32_t)rnd.NextInt(800, 950);
		uint32_t mark = (uint32_t)rnd.xrandom();

		//Some chance to use incompletes
		for (auto&i : currentBest.IncompleteSets)
			if (rnd.NextInt(1000) < 860)
			{
				Z.push_back(i);
			}
		do_shuffle<vector<int>>(Z, rnd);


		if (tryMerge)
		{
			for (int i = 0; i < COUNT_NUMBERS; ++i)
			{
				if (!UsedNumbers.get(i))
				{
					UnusedNumbers.push_back(i);
					maskX = 1ULL << (NMB[i].X);
					maskY = 1ULL << (NMB[i].Y);
				}
			}
			for (auto&i : currentBest.IncompleteSets)
				if (rnd.NextInt(1000) < 800)
				{
					auto& ST = currentBest.grid.Plan_NMB[i];
					auto& nn = ST.linear[ST.linear.size() - 1].destIDX;
					UnusedNumbers.push_back(nn);
					maskX = 1ULL << (NMB[nn].X);
					maskY = 1ULL << (NMB[nn].Y);
				}
		}


		for (auto&z : Z)
		{
			currentBest.grid.Plan_NMB[z].ApplyMoves(candidate.grid, chanceMerge, mark, withoutEndpoint, UnusedNumbers, rnd);
		}
	}
}

//Just pick random points and remove/ignore them
inline void Mutate_Points_1(LAHC_Node& currentBest, LAHC_Node& candidate, Random& rnd)
{
	candidate.MutateType = 1;
	vector<Move> linearMoves;
	vector<int> K;
	NumberSet UsedNumbers;
	for (auto&z : currentBest.ZerosumSets)
		if (rnd.NextInt(1000) < 970)
		{
			if (currentBest.grid.Plan_NMB[z].linear.size() > 1 || (rnd.NextInt(1000) < 700)) {
				K.push_back(z);
				for (auto&l : currentBest.grid.Plan_NMB[z].linear)
				{
					//linearMoves.push_back(l);
					UsedNumbers.set(l.srcIDX);
					UsedNumbers.set(l.destIDX);
				}
			}
		}
	for (auto&I : currentBest.IncompleteSets)
		if (rnd.NextInt(1000) < 970)
		{
			if (currentBest.grid.Plan_NMB[I].linear.size() > 1 || (rnd.NextInt(1000) < 700)) {
				K.push_back(I);
				for (auto&l : currentBest.grid.Plan_NMB[I].linear)
				{
					UsedNumbers.set(l.srcIDX);
				}
			}
		}
	do_shuffle<vector<int>>(K, rnd);
	for (auto&k : K)
		Add(linearMoves, currentBest.grid.Plan_NMB[k].linear);
	if (linearMoves.size() > 0)
	{
		uint64_t maskX = 0;
		uint64_t maskY = 0;
		bool tryMerge = (currentBest.usedInMovesNumbers.count() > COUNT_NUMBERS * 50 / 100) && rnd.NextInt(100) < 80;
		vector<int> UnusedNumbers;
		if (tryMerge)
		{
			for (int i = 0; i < COUNT_NUMBERS; ++i)
			{
				if (!UsedNumbers.get(i))
				{
					UnusedNumbers.push_back(i);
					maskX = 1ULL << (NMB[i].X);
					maskY = 1ULL << (NMB[i].Y);

				}
			}
		}
		int countRemove = rnd.NextInt(1, min((currentBest.grid.totalNumbers < 5 ? K_MAX_POINTS_REMOVE : 3), (int)linearMoves.size()));

		unordered_set<int> toRemove;
		for (int i = 0; i < countRemove; ++i)
		{
			int mutateIndex;
			do {
				mutateIndex = (int)rnd.NextInt((uint32_t)linearMoves.size());
			} while (toRemove.find(mutateIndex) != toRemove.end());
			toRemove.insert(mutateIndex);
		}
		for (int i = 0; i < linearMoves.size(); ++i)
			if (toRemove.find(i) == toRemove.end())
			{
				bool BU = (((1ULL << NMB[linearMoves[i].srcIDX].X) & maskX) != 0 || ((1ULL << NMB[linearMoves[i].srcIDX].Y) & maskY) != 0) && tryMerge && (rnd.NextInt(1000) < 900);
				if (BU)
				{
					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, linearMoves[i].srcIDX, rnd);
					}
				}
				candidate.grid.ApplyMove(linearMoves[i]);
				BU = (((1ULL << NMB[linearMoves[i].destIDX].X) & maskX) != 0 || ((1ULL << NMB[linearMoves[i].destIDX].Y) & maskY) != 0) && tryMerge && (rnd.NextInt(1000) < 900);
				if (BU)
				{
					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, linearMoves[i].destIDX, rnd);
					}
				}

			}
	}

}

//Remove tails on Lists
inline void Mutate_Lists_2(LAHC_Node& currentBest, LAHC_Node& candidate, Random& rnd)
{
	if (currentBest.Movelists.size() > 0)
	{
		uint64_t maskX = 0;
		uint64_t maskY = 0;
		candidate.MutateType = 2;
		unordered_set<int> usedPlan;
		bool tryMerge = /*currentBest.Movelists.size() > 5 && */(currentBest.usedInMovesNumbers.count() > COUNT_NUMBERS * 50 / 100) && rnd.NextInt(100) < 85;
		vector<int> UnusedNumbers;
		if (tryMerge)
		{
			for (int i = 0; i < COUNT_NUMBERS; ++i)
			{
				if (!currentBest.usedInMovesNumbers.get(i))
				{
					UnusedNumbers.push_back(i);
					maskX = 1ULL << (NMB[i].X);
					maskY = 1ULL << (NMB[i].Y);
				}
			}
		}

		if (rnd.NextInt(1000) < 400) //Do incomplete first, and try to add to the list of unused numbers
			for (auto&i : currentBest.IncompleteSets)
				if (rnd.NextInt(1000) < 800)
				{
					auto& LMOV = currentBest.grid.Plan_NMB[i].linear;
					//Do moves
					for (auto&m : LMOV)
					{
						candidate.grid.ApplyMove(m);
					}
					//Add the incomplete end to unused
					UnusedNumbers.push_back(LMOV[LMOV.size() - 1].destIDX);
					maskX = 1ULL << (NMB[LMOV[LMOV.size() - 1].destIDX].X);
					maskY = 1ULL << (NMB[LMOV[LMOV.size() - 1].destIDX].Y);
				}

		int ListsToMutate = rnd.NextInt(1, min(4, (int)currentBest.Movelists.size()));

		for (int i = 0; i < ListsToMutate; ++i)
		{
			int mutateIndex;
			do {
				mutateIndex = (int)rnd.NextInt((uint32_t)currentBest.Movelists.size());
			} while (usedPlan.find(mutateIndex) != usedPlan.end());
			usedPlan.insert(mutateIndex);
		}

		int i = rnd.NextInt((int)currentBest.Movelists.size());
		//for (int Ni = 0; Ni < (int)currentBest.Movelists.size(); ++Ni)
		for (auto& bbg : currentBest.Movelists)
		{
			++i;
			if (i >= currentBest.Movelists.size())
			{
				i = 0;
			}
			auto& UEL = currentBest.grid.Plan_NMB[currentBest.Movelists[i]];
			uint32_t limitCopy = (uint32_t)UEL.linear.size();
			if (usedPlan.find(i) != usedPlan.end())
			{ //To mutate, we'll cut these
				limitCopy = (int)rnd.NextInt(limitCopy);// max(1, (int)rnd.NextInt(limitCopy));
			}
			for (uint32_t j = 0; j < limitCopy; ++j)
			{
				//auto& n = NMB[UEL.linear[j].srcIDX];
				bool BU = (((1ULL << NMB[UEL.linear[j].srcIDX].X) & maskX) != 0 || ((1ULL << NMB[UEL.linear[j].srcIDX].Y) & maskY) != 0)
					&& tryMerge && (rnd.NextInt(1000) < 800);
				if (BU)
				{
					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, UEL.linear[j].srcIDX, rnd);
					}
				}
				candidate.grid.ApplyMove(UEL.linear[j]);
				BU = (((1ULL << NMB[UEL.linear[j].destIDX].X) & maskX) != 0 || ((1ULL << NMB[UEL.linear[j].destIDX].Y) & maskY) != 0)
					&& tryMerge && (rnd.NextInt(1000) < 800);
				if (BU)
				{
					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, UEL.linear[j].destIDX, rnd);
					}
				}

			}
		}
	}
}

//Similar to Points_1, but tries to get adjacent points
inline void Mutate_AdjPoints_3(LAHC_Node& currentBest, LAHC_Node& candidate, Random& rnd)
{
	candidate.MutateType = 3;
	vector<Move> linearMoves;
	vector<int> K;
	NumberSet UsedNumbers;
	for (auto&z : currentBest.ZerosumSets)
		if (rnd.NextInt(1000) < 970)
		{
			if (currentBest.grid.Plan_NMB[z].linear.size() > 1 || (rnd.NextInt(1000) < 900)) {
				K.push_back(z);
				for (auto&l : currentBest.grid.Plan_NMB[z].linear)
				{
					UsedNumbers.set(l.srcIDX);
					UsedNumbers.set(l.destIDX);
				}
			}
		}
	NumberSet CloseIDX;
	CloseIDX.clear();
	for (auto&I : currentBest.IncompleteSets)
		if (rnd.NextInt(1000) < 970)
		{
			CloseIDX.set(currentBest.grid.Plan_NMB[I].linear[currentBest.grid.Plan_NMB[I].linear.size() - 1].destIDX);
			if (currentBest.grid.Plan_NMB[I].linear.size() > 1 || (rnd.NextInt(1000) < 900)) {
				K.push_back(I);
				for (auto&l : currentBest.grid.Plan_NMB[I].linear)
				{
					UsedNumbers.set(l.srcIDX);
				}
			}
		}
	do_shuffle<vector<int>>(K, rnd);

	for (auto&k : K)
		Add(linearMoves, currentBest.grid.Plan_NMB[k].linear);
	if (linearMoves.size() > 0)
	{

		bool tryMerge = true;// rnd.NextInt(1000) < 960;
		vector<int> UnusedNumbers;
		uint64_t maskX = 0;
		uint64_t maskY = 0;
		if (tryMerge)
		{
			for (int i = 0; i < COUNT_NUMBERS; ++i)
			{
				if (!UsedNumbers.get(i))
				{
					UnusedNumbers.push_back(i);
					maskX = 1ULL << (NMB[i].X);
					maskY = 1ULL << (NMB[i].Y);
				}
			}
		}
		int countRemove = rnd.NextInt(1, min((currentBest.grid.totalNumbers < 5 ? 3 : 2), (int)linearMoves.size()));


		unordered_set<int> toRemove;
		for (int i = 0; i < countRemove; ++i)
		{
			int mutateIndex;
			mutateIndex = (int)rnd.NextInt((uint32_t)linearMoves.size());
			if (toRemove.find(mutateIndex) == toRemove.end())
				toRemove.insert(mutateIndex);

			//Extra remove
			int extraRem = rnd.NextInt(2);
			if (extraRem > 0)
				for (int ex = 1; ex <= extraRem; ++ex)
				{
					int ule = ex + mutateIndex;
					if (ule >= linearMoves.size())
						break;
					/*		if (linearMoves[ule - 1].destIDX != linearMoves[ule].srcIDX)
					break;*/
					if (toRemove.find(ule) == toRemove.end())
					{
						toRemove.insert(ule);
					}
				}
		}
		vector<pair<int, int>> explotar;
		for (auto&i : toRemove)
		{
			explotar.push_back(pair<int, int>{2 + i, linearMoves[i].destIDX});
			explotar.push_back(pair<int, int>{2 + i, linearMoves[i].srcIDX});
		}
		do_shuffle< vector<pair<int, int>>>(explotar, rnd);
		do_shuffle<vector<int>>(UnusedNumbers, rnd);
		vector<Move> tmpMoves;

		for (int i = 0; i < linearMoves.size(); ++i)
			if (toRemove.find(i) == toRemove.end())
			{
				bool BU = (((1ULL << NMB[linearMoves[i].srcIDX].X) & maskX) != 0 || ((1ULL << NMB[linearMoves[i].srcIDX].Y) & maskY) != 0)
					&& tryMerge && (rnd.NextInt(1000) < 900);
				if (BU)
				{

					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, linearMoves[i].srcIDX, rnd);
					}
				}
				if (tryMerge && (rnd.NextInt(1000) < 900))
				{
					for (auto& U : explotar)
						if (U.first > i + 1)
						{
							candidate.grid.doMerge(U.second, linearMoves[i].srcIDX, rnd);
						}
				}

				if (CloseIDX.get(linearMoves[i].destIDX)) {
					for (auto& U : explotar)
						if (U.first > i)
						{
							tmpMoves.resize(0);
							candidate.grid.ExplodeMoves(U.second, tmpMoves, &rnd);
							for (auto&mm : tmpMoves)
							{
								if (mm.destIDX == linearMoves[i].destIDX && (rnd.NextInt(1000) < 850))
								{
									candidate.grid.ApplyMove(mm);
								}
							}
						}
				}

				candidate.grid.ApplyMove(linearMoves[i]);
				BU = (((1ULL << NMB[linearMoves[i].destIDX].X) & maskX) != 0 || ((1ULL << NMB[linearMoves[i].destIDX].Y) & maskY) != 0)
					&& tryMerge && (rnd.NextInt(1000) < 900);
				if (BU)
				{
					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, linearMoves[i].destIDX, rnd);
					}
				}
				if (tryMerge && (rnd.NextInt(1000) < 900))
				{
					for (auto& U : explotar)
						if (U.first > i)
						{
							candidate.grid.doMerge(U.second, linearMoves[i].destIDX, rnd);
						}
				}
				/*if (CloseIDX.get(linearMoves[i].destIDX)) {
				for (auto& U : explotar)
				if (U.first > i)
				{
				candidate.grid.doForceJoin(U.second, linearMoves[i].destIDX, 850);
				}
				}*/
			}
	}

}

/*************** Mutate_Complex_4 *******************/
/*
It will remove Points, but based on different selectors:
-On a pure random, like Mutate_Points_1
-On tails of Strategies, like Mutate_Lists_2
-On orphan numbers, these are the numbers that doesn't have any valid neighbour so they can't be removed at all
-On Cross moves, these are the neighbours of incomplete numbers.

This is an evolution of previous Mutators, to be able to mix them so they are more dynamic.
*/
/* Random move, but from endings*/

//Pick a random move
void selectorRandom(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves, Random& rnd)
{
	int rndTotal = rnd.NextInt((int)totalMoves);
	int selList = rndLowBound(ListSizes, rndTotal);
	if (selList >= 0)
	{
		auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
		int selPoint = rnd.NextInt((int)ML.linear.size());
		int Code = selList * 10000 + selPoint;
		if (tmpPointsToRemove.find(Code) == tmpPointsToRemove.end())
		{
			tmpPointsToRemove.emplace(Code);
		}
	}
}
void selectorRandomTail(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves, Random& rnd)
{
	int rndTotal = rnd.NextInt((int)totalMoves);
	int selList = rndLowBound(ListSizes, rndTotal);
	if (selList >= 0)
	{
		auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
		int Code;
		if (ML.linear.size() >= 5)
		{
			int start = min((int)ML.linear.size() - 3, (int)ML.linear.size()*K_TAIL_PERCENT / 100);
			int startEND = (int)ML.linear.size() - 1;
			if (startEND == start && rnd.NextInt(1000) < 800)
				start = startEND - 1;
			start = min(max(0, start), startEND);
			int selPoint = rnd.NextInt(start, startEND);
			Code = selList * 10000 + selPoint;
		}
		else Code = selList * 10000 + rnd.NextInt((int)ML.linear.size());

		if (tmpPointsToRemove.find(Code) == tmpPointsToRemove.end())
		{
			tmpPointsToRemove.emplace(Code);
		}
	}
}
//Get a move that is 2*destVal = srcVal and subtract
void selectorOptional(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves, Random& rnd)
{
	int rndTotal = rnd.NextInt((int)totalMoves);
	int selList = rndLowBound(ListSizes, rndTotal);
	if (selList >= 0)
	{
		for (int uuv = 0; uuv < currentBest.Movelists.size(); ++uuv)
		{
			auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
			int searchOpt = rnd.NextInt((int)ML.linear.size());
			for (int uuz = 0; uuz < ML.linear.size(); ++uuz)
			{
				++searchOpt;
				if (searchOpt >= ML.linear.size())
					searchOpt = 0;
				auto& m = ML.linear[searchOpt];
				if (2 * m.destVal == m.srcVal && m.sign == 0 && m.sign == 0 && m.destVal != 0)
				{
					int Code = selList * 10000 + searchOpt;
					if (tmpPointsToRemove.find(Code) == tmpPointsToRemove.end()) {
						tmpPointsToRemove.emplace(Code);
						return;
					}
				}
			}

			++selList;
			if (selList >= currentBest.Movelists.size())
				selList = 0;
		}
	}
}
//Without neighbours
void selectorOrphanMoves(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves, const vector<int>& orphanNumbers, Random& rnd)
{

	if (orphanNumbers.size() == 1 && rnd.NextInt(1000) < 150)
	{ //Single number, try to change all points affecting it
		auto& orphan = NMB[orphanNumbers[0]];
		//for (auto& N : currentBest.Movelists)
		vector<int> affecting;
		for (int selList = 0; selList < ListSizes.size(); ++selList)
		{
			auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
			for (int searchOpt = 0; searchOpt < ML.linear.size(); ++searchOpt)
			{
				auto& m = ML.linear[searchOpt];
				if (m.srcIDX == orphan.ID || m.destIDX == orphan.ID)
				{
					int Code = selList * 10000 + searchOpt;
					affecting.push_back(Code);
				}
			}
		}


		if (affecting.size() > 0)
		{
			if (affecting.size() > 1 && rnd.NextInt(1000) < 400)
				affecting.resize(rnd.NextInt(1, (int)affecting.size()));
			for (auto& Code : affecting)
				if (tmpPointsToRemove.find(Code) == tmpPointsToRemove.end()) {
					tmpPointsToRemove.emplace(Code);
				}
			return;
		}
	}

	if (orphanNumbers.size() != 0)
	{
		auto& orphan = NMB[orphanNumbers[rnd.NextInt((int)orphanNumbers.size())]];
		auto neigh = NEIGHBOURS[orphan.ID];

		int rndTotal = rnd.NextInt((int)totalMoves);
		int selList = rndLowBound(ListSizes, rndTotal);
		if (selList >= 0)
		{
			for (int uuv = 0; uuv < currentBest.Movelists.size(); ++uuv)
			{
				auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
				int searchOpt = rnd.NextInt((int)ML.linear.size());
				for (int uuz = 0; uuz < ML.linear.size(); ++uuz)
				{
					++searchOpt;
					if (searchOpt >= ML.linear.size())
						searchOpt = 0;
					auto& m = ML.linear[searchOpt];
					if (neigh.get(m.srcIDX) || (rnd.NextInt(1000) < 400 && neigh.get(m.destIDX)))
					{
						int Code = selList * 10000 + searchOpt;
						if (tmpPointsToRemove.find(Code) == tmpPointsToRemove.end()) {
							tmpPointsToRemove.emplace(Code);
							return;
						}
					}
				}
				++selList;
				if (selList >= currentBest.Movelists.size())
					selList = 0;
			}
		}

	}
}
//Any move in the X,Y of a remaining number
void selectorCrossMoves(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves, const vector<int>& pendingNumbers, Random& rnd)
{
	if (pendingNumbers.size() > 0)
	{
		int rndID = pendingNumbers[rnd.NextInt((int)pendingNumbers.size())];
		auto& crossNumber = NMB[rndID];
		auto neigh = NEIGHBOURS[crossNumber.ID];

		int rndTotal = rnd.NextInt((int)totalMoves);
		int selList = rndLowBound(ListSizes, rndTotal);
		if (selList >= 0)
		{
			for (int uuv = 0; uuv < currentBest.Movelists.size(); ++uuv)
			{
				auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
				int searchOpt;
				//Prefer endings
				if (ML.linear.size() > 1)
				{
					if ((ML.linear.size() > 4) && rnd.NextInt(1000) < K_PREFER_ENDINGS)
					{
						int start = (int)ML.linear.size() * 80 / 100;
						int startEND = (int)ML.linear.size() - 1;
						start = min(max(0, start), startEND);
						searchOpt = rnd.NextInt(start, startEND);
					}
					else searchOpt = rnd.NextInt((int)ML.linear.size());
				}
				else searchOpt = 0;
				for (int uuz = 0; uuz < ML.linear.size(); ++uuz)
				{
					++searchOpt;
					if (searchOpt >= ML.linear.size())
						searchOpt = 0;
					auto& m = ML.linear[searchOpt];
					if (neigh.get(m.srcIDX) || neigh.get(m.destIDX))
					{
						int Code = selList * 10000 + searchOpt;
						if (tmpPointsToRemove.find(Code) == tmpPointsToRemove.end()) {
							tmpPointsToRemove.emplace(Code);
							return;
						}
					}
				}
				++selList;
				if (selList >= currentBest.Movelists.size())
					selList = 0;
			}
		}

	}
}
void calcOrphanMoves(const LAHC_Node& currentBest, vector<int>& orphanNumbers) {

	for (int i = 0; i < COUNT_NUMBERS; ++i)
	{
		if (currentBest.grid.Val[i] != 0)
		{
			//Check Neighbours
			if (NEIGHBOURS[i].Disjoint(currentBest.remainingNumbers))
			{
				orphanNumbers.push_back(i);
			}
		}
	}


}
inline void Mutate_Complex_4(const LAHC_Node& currentBest, LAHC_Node& candidate, Random& rnd)
{
	/*	NumberSet consumed;
	consumed.clear();*/
	vector<int> orphanmoves;
	calcOrphanMoves(currentBest, orphanmoves);

	unordered_set<int> tmpPointsToRemove;
	vector<int> pendingNumbers;
	uint64_t maskX = 0;
	uint64_t maskY = 0;

	for (int i = 0; i < COUNT_NUMBERS; ++i)
	{
		if (currentBest.grid.Val[i] != 0)
		{
			pendingNumbers.push_back(i);
		}
	}
	vector<pair<int, int>> cutPoint;
	{ //Preparing the cutpoints

		vector<size_t> ListSizes;
		ListSizes.resize(currentBest.Movelists.size());
		size_t totalMoves = 0;
		//for (auto& L : currentBest.Movelists)
		for (int n = 0; n < currentBest.Movelists.size(); ++n)
		{
			totalMoves += currentBest.grid.Plan_NMB[currentBest.Movelists[n]].linear.size();
			ListSizes[n] = totalMoves;
		}
		if (totalMoves < 2 * K_CHANCE_POINTS_REMOVE)
			return; //Don't waste my time
		int countRemove = rnd.NextInt(1, max(1, min((currentBest.grid.totalNumbers < 5 ? K_CHANCE_POINTS_REMOVE : K_LOW_POINTS_REMOVE), (int)totalMoves)));

		for (int T = 0; T < 3 * countRemove; ++T)
		{
			int chance = rnd.NextInt(K_SEL_VAL_RANDOM + K_SEL_VAL_TAIL + K_SEL_VAL_OPT + K_SEL_VAL_ORPHAN + K_SEL_VAL_CROSS);
			//Temporal
			//	selectorRandomTail(tmpPointsToRemove, currentBest, ListSizes, totalMoves);
			if (chance < K_SEL_VAL_ORPHAN && (orphanmoves.size() != 0))
				selectorOrphanMoves(tmpPointsToRemove, currentBest, ListSizes, totalMoves, orphanmoves, rnd);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL)
				selectorRandomTail(tmpPointsToRemove, currentBest, ListSizes, totalMoves, rnd);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL + K_SEL_VAL_OPT)
				selectorOptional(tmpPointsToRemove, currentBest, ListSizes, totalMoves, rnd);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL + K_SEL_VAL_OPT + K_SEL_VAL_RANDOM)
				selectorRandom(tmpPointsToRemove, currentBest, ListSizes, totalMoves, rnd);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL + K_SEL_VAL_OPT + K_SEL_VAL_RANDOM + K_SEL_VAL_CROSS)
				selectorCrossMoves(tmpPointsToRemove, currentBest, ListSizes, totalMoves, pendingNumbers, rnd);
			else selectorRandom(tmpPointsToRemove, currentBest, ListSizes, totalMoves, rnd);
			if (tmpPointsToRemove.size() >= countRemove)
				break;
		}
		if (tmpPointsToRemove.size() == 0)
			return; //Task failed succesfully
		NumberSet related;

		unordered_set<size_t> forbiddenMoves;

		vector<int> endings;
		for (auto&t : tmpPointsToRemove)
		{
			int selList = t / 10000;
			int selPoint = t % 10000;
			ASSERT(selList < currentBest.Movelists.size());
			auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
			ASSERT(ML.linear.size() > 0);
			ASSERT(selPoint < ML.linear.size());
			Move moveToCut = ML.linear[selPoint];

			forbiddenMoves.emplace(moveToCut.CreateHash());

			related.clear();
			related.set(moveToCut.srcIDX);
			related.set(moveToCut.destIDX);
			vector<int> repeatMoves;
			if (selPoint > 0)
				for (int i = selPoint - 1; i >= 0; --i) //Backtrack to Find all related moves
				{
					if ((related.get(ML.linear[i].srcIDX) || related.get(ML.linear[i].destIDX)))
					{
						repeatMoves.push_back(i);
						related.set(ML.linear[i].srcIDX);
						related.set(ML.linear[i].destIDX);
					}
				}
			//Apply those strings
			//for (auto& i : repeatMoves)
			if (repeatMoves.size() > 0)
				for (int uuj = (int)repeatMoves.size() - 1; uuj >= 0; --uuj)
				{
					int i = repeatMoves[uuj];
					if (ML.linear[i].destIDX == moveToCut.srcIDX || ML.linear[i].destIDX == moveToCut.destIDX)
					{

						Move m = ML.linear[i];
						int chances = rnd.NextInt(1000);
						if (chances < K_CHANCE_SIGN) //Change sign
						{
							m.sign = 1 - m.sign;
							pendingNumbers.push_back(m.destIDX);
							candidate.grid.ApplyMove(m);

						}
						else if (chances < K_CHANCE_SIGN + K_CHANCE_IGNORE) {
							//Just ignore the move
							pendingNumbers.push_back(m.srcIDX);
							if (rnd.NextInt(1000) < 300)
								pendingNumbers.push_back(m.destIDX);
						}
						else if (candidate.grid.ApplyMove(m)) {/*consumed.set(m.srcIDX);*/ };
					}
					else {
						if (candidate.grid.ApplyMove(ML.linear[i])) {/*consumed.set(ML.linear[i].destIDX);*/ }
					}
				}
			repeatMoves.resize(0);
			//Now repeat all, with forbidden moves.
		}



		do_shuffle<vector<int>>(pendingNumbers, rnd);
		vector<Move> tmpMoves;
		//Now repeat all, trying to Merge / explode but without forbidden moves.
		vector<int> Listas = currentBest.IncompleteSets;// currentBest.Movelists;
		if (rnd.NextInt(1000) < 500)
		{
			do_shuffle<vector<int>>(Listas, rnd);
			auto HJU = currentBest.ZerosumSets;
			do_shuffle<vector<int>>(HJU, rnd);
			for (auto&I : HJU)
				Listas.push_back(I);

		}
		else {
			for (auto&I : currentBest.ZerosumSets)
				Listas.push_back(I);
			do_shuffle<vector<int>>(Listas, rnd);
		}




		for (auto&i : pendingNumbers) {
			maskX = 1ULL << (NMB[i].X);
			maskY = 1ULL << (NMB[i].Y);
		}
		//speed up.
		//Limit pending number merges
		/*	if (pendingNumbers.size() > 5)
		pendingNumbers.resize(5);*/
		for (auto&L : Listas)
		{
			auto& ML = currentBest.grid.Plan_NMB[L];
			for (auto&uwm : ML.linear) {
				if (forbiddenMoves.find(uwm.CreateHash()) == forbiddenMoves.end())
				{
					bool BU = (rnd.NextInt(1000) < K_CHANCE_MERGE);
					if (BU)
					{
						if (((1ULL << NMB[uwm.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[uwm.srcIDX].Y) & maskY) != 0)
							for (auto& t : pendingNumbers)
								if (t != uwm.srcIDX && t != uwm.destIDX)
								{
									candidate.grid.doMerge(t, uwm.srcIDX, rnd);
								}
					}
					bool forced = false;
					if (rnd.NextInt(1000) < K_CHANCE_EXPLODE)
					{
						if (((1ULL << NMB[uwm.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[uwm.srcIDX].Y) & maskY) != 0)
							for (auto& unused : pendingNumbers)
							{
								if (candidate.grid.doForceJoin(unused, uwm.destIDX, 1000, rnd)) {
									forced = true;
								}

							}

					}
					if (!forced)
						candidate.grid.ApplyMove(uwm);
					else {
						Move m = uwm;
						m.destVal = candidate.grid.Val[m.destIDX];
						candidate.grid.ApplyMove(m);
					}
					BU = (rnd.NextInt(1000) < K_CHANCE_MERGE);
					if (BU)
					{
						if (((1ULL << NMB[uwm.destIDX].X) & maskX) != 0 || ((1ULL << NMB[uwm.destIDX].Y) & maskY) != 0)
							for (auto& t : pendingNumbers)
								if (t != uwm.srcIDX && t != uwm.destIDX)
								{
									candidate.grid.doMerge(t, uwm.destIDX, rnd);
								}

					}

				}
			}

			//If incomplete add ending number as pending
			if (candidate.grid.Val[ML.linear.back().destIDX] != 0)
			{
				int i = ML.linear.back().destIDX;
				pendingNumbers.push_back(i);
				maskX = 1ULL << (NMB[i].X);
				maskY = 1ULL << (NMB[i].Y);
			}
		}
	}
	candidate.MutateType = 4;
}

//Will remove a % of moves affecting a X,Y of remaining numbers.
inline bool Mutate_Cross_4(LAHC_Node& currentBest, LAHC_Node& candidate, Random& rnd) {
	vector<int> UnusedNumbers;
	uint64_t maskX = 0;
	uint64_t maskY = 0;

	int countRemove = rnd.NextInt(1, (currentBest.grid.totalNumbers < 5 ? K_MAX_POINTS_REMOVE : 3));
	vector<size_t> OPT;
	for (auto&M : currentBest.Movelists)
	{
		auto& ST = currentBest.grid.Plan_NMB[M];
		//for (auto&m : ST.linear)
		for (int i = 0; i < ST.linear.size(); ++i)
		{
			auto&m = ST.linear[i];
			if ((m.sign == 0) && (2 * m.destVal == m.srcVal)) //an optional
			{
				OPT.push_back(m.CreateHash());
				if (i > 1 && ST.linear[i].destIDX == m.srcVal && rnd.NextInt(1000) < 300)
				{
					OPT.push_back(ST.linear[i].CreateHash());
					if (rnd.NextInt(1000) < 800)
						UnusedNumbers.push_back(ST.linear[i].srcIDX);

				}
				if (rnd.NextInt(1000) < 800)
					UnusedNumbers.push_back(m.srcIDX);
			}
		}
	}

	if (OPT.size() == 0)
		return false;
	do_shuffle<vector<size_t>>(OPT, rnd);
	if (OPT.size() > countRemove)
		OPT.resize(countRemove);

	vector<Move> linearMoves;
	vector<int> K;
	/*	auto Z = currentBest.ZerosumSets;
	auto I = currentBest.IncompleteSets;*/
	/*	do_shuffle<vector<int>>(Z);
	do_shuffle<vector<int>>(I);*/
	NumberSet UsedNumbers;
	for (auto&z : currentBest.ZerosumSets)
		if (rnd.NextInt(1000) < 970)
		{
			if (currentBest.grid.Plan_NMB[z].linear.size() > 1 || (rnd.NextInt(1000) < 700)) {
				K.push_back(z);
				for (auto&l : currentBest.grid.Plan_NMB[z].linear)
				{
					//linearMoves.push_back(l);
					UsedNumbers.set(l.srcIDX);
					UsedNumbers.set(l.destIDX);
				}
			}
		}
	for (auto&I : currentBest.IncompleteSets)
		if (rnd.NextInt(1000) < 970)
		{
			if (currentBest.grid.Plan_NMB[I].linear.size() > 1 || (rnd.NextInt(1000) < 700)) {
				K.push_back(I);
				for (auto&l : currentBest.grid.Plan_NMB[I].linear)
				{
					//linearMoves.push_back(l);
					UsedNumbers.set(l.srcIDX);
					//UsedNumbers.set(l.destIDX);
				}
			}
		}
	do_shuffle<vector<int>>(K, rnd);
	for (auto&k : K)
		Add(linearMoves, currentBest.grid.Plan_NMB[k].linear);
	if (linearMoves.size() > 0)
	{

		bool tryMerge = /*currentBest.Movelists.size() > 8 &&*/ (currentBest.usedInMovesNumbers.count() > COUNT_NUMBERS * 50 / 100) && rnd.NextInt(100) < 80;

		if (tryMerge)
		{
			for (int i = 0; i < COUNT_NUMBERS; ++i)
			{
				if (!UsedNumbers.get(i))
				{
					UnusedNumbers.push_back(i);
				}
			}
		}
		for (auto&U : UnusedNumbers)
		{
			maskX = 1ULL << (NMB[U].X);
			maskY = 1ULL << (NMB[U].Y);
		}

		unordered_set<size_t> toRemove;
		for (auto& O : OPT)
			toRemove.insert(O);
		size_t mutateHash;
		countRemove = min(countRemove, (int)linearMoves.size());
		while (toRemove.size() < countRemove)
		{
			mutateHash = linearMoves[rnd.NextInt((uint32_t)linearMoves.size())].CreateHash();
			if (toRemove.find(mutateHash) == toRemove.end()) {
				toRemove.insert(mutateHash);
			}
		}

		for (int i = 0; i < linearMoves.size(); ++i)
		{
			auto& lm = linearMoves[i];
			size_t hashMove = lm.CreateHash();
			if (toRemove.find(hashMove) != toRemove.end())
				continue;
			bool BU = (((1ULL << NMB[lm.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[lm.srcIDX].Y) & maskY) != 0)
				&& tryMerge && (rnd.NextInt(1000) < 900);
			if (BU)
			{
				for (auto& unused : UnusedNumbers)
				{
					candidate.grid.doMerge(unused, lm.srcIDX, rnd);
				}
			}
			candidate.grid.ApplyMove(lm);
			BU = (((1ULL << NMB[lm.destIDX].X) & maskX) != 0 || ((1ULL << NMB[lm.destIDX].Y) & maskY) != 0)
				&& tryMerge && (rnd.NextInt(1000) < 900);
			if (BU)
			{
				for (auto& unused : UnusedNumbers)
				{
					candidate.grid.doMerge(unused, lm.destIDX, rnd);
				}
			}

		}
	}
	candidate.MutateType = 4;
	return true;
}

//Mutating ending values of Incomplete Strategies.
//According to stats, more than 40% of accepted candidates come from that Mutator, so it's high quality. I don't know why
bool Mutate_Shuffle_Incomplete_End_5(LAHC_Node& currentBest, LAHC_Node& candidate, Random& rnd)
{

	vector<int> INC;
	uint64_t maskX = 0;
	uint64_t maskY = 0;

	for (auto &I : currentBest.IncompleteSets)
	{
		auto& ST = currentBest.grid.Plan_NMB[I];
		if (ST.linear.size() > 1)
		{
			int lastIDX = ST.linear[ST.linear.size() - 1].destIDX;
			int cuenta = 0;
			for (auto&m : ST.linear)
			{
				if (m.destIDX == lastIDX)
					++cuenta;
			}
			if (cuenta > 0 && cuenta != ST.linear.size())
				INC.push_back(I);
		}
	}
	if (INC.size() == 0)
		return false;
	candidate.MutateType = 5;
	do_shuffle<vector<int>>(INC, rnd);
	int resizeI = min((int)INC.size(), (int)rnd.NextInt(1, 3));
	if (INC.size() > resizeI)
	{
		INC.resize(resizeI);
	}
	unordered_set<int> usedI;
	vector<int> TryMerge;
	for (int i = 0; i < COUNT_NUMBERS; ++i)
	{
		if (!currentBest.usedInMovesNumbers.get(i))
		{
			TryMerge.push_back(i);
			maskX = 1ULL << (NMB[i].X);
			maskY = 1ULL << (NMB[i].Y);
		}
	}
	do_shuffle<vector<int>>(TryMerge, rnd);
	for (auto&I : INC)
	{
		usedI.emplace(I);
		auto& ST = currentBest.grid.Plan_NMB[I];
		int lastIDX = ST.linear[ST.linear.size() - 1].destIDX;
		auto lastX = NMB[lastIDX].X;
		auto lastY = NMB[lastIDX].Y;
		vector<Move> lastMoves;
		for (auto& m : ST.linear)
		{
			if (m.destIDX == lastIDX)
			{
				lastMoves.push_back(m);
			}
			else
			{ //Do all moves except ending
				if (rnd.NextInt(1000) < 500)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.srcIDX, rnd);
				candidate.grid.ApplyMove(m);
				if (rnd.NextInt(1000) < 500)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.destIDX, rnd);
			}
		}
		do_shuffle<vector<Move>>(lastMoves, rnd);
		for (auto& m : lastMoves)
		{
			m.sign = rnd.NextInt(100) < K_ALLOWSUM ? 1 : 0;
			ASSERT(NMB[m.destIDX].X == lastX);
			ASSERT(NMB[m.destIDX].Y == lastY);
			m.destVal = candidate.grid.Val[lastIDX];
			if (((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
				if (rnd.NextInt(1000) < K5_MERGE0)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.srcIDX, rnd);
			candidate.grid.ApplyMove(m);
			if (((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
				if (rnd.NextInt(1000) < K5_MERGE0)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.destIDX, rnd);
		}
		//if (candidate.grid.Val[lastIDX] != 0) 
		{
			TryMerge.push_back(lastIDX);
			maskX = 1ULL << (NMB[lastIDX].X);
			maskY = 1ULL << (NMB[lastIDX].Y);
		}
	}
	do_shuffle<vector<int>>(TryMerge, rnd);
	vector<int> NMI = currentBest.IncompleteSets;
	if (rnd.NextInt(1000) < 500)
	{
		do_shuffle<vector<int>>(NMI, rnd);
		auto HJU = currentBest.ZerosumSets;
		do_shuffle<vector<int>>(HJU, rnd);
		for (auto&I : HJU)
			NMI.push_back(I);

	}
	else {
		for (auto&I : currentBest.ZerosumSets)
			NMI.push_back(I);
		do_shuffle<vector<int>>(NMI, rnd);
	}



	for (auto&I : NMI)
		if (usedI.find(I) == usedI.end())
		{
			auto& ST = currentBest.grid.Plan_NMB[I];
			int lastIDX = ST.linear[ST.linear.size() - 1].destIDX;
			for (auto& m : ST.linear)
			{
				if ((((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
					&& rnd.NextInt(1000) < K5_MERGE1)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.srcIDX, rnd);
				candidate.grid.ApplyMove(m);
				if ((((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
					&& rnd.NextInt(1000) < K5_MERGE1)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.destIDX, rnd);
			}
			if (candidate.grid.Val[lastIDX] != 0)
			{
				TryMerge.push_back(lastIDX);
				maskX = 1ULL << (NMB[lastIDX].X);
				maskY = 1ULL << (NMB[lastIDX].Y);
			}
		}
	return true;
}

//***************** Mutate_Exhaustive_6 - Exhaustive Mutator **************/
//Picks incomplete strategies. Permutates ending moves to get a list of all possible end values.
//Then it tries to merge these ending values with any neighbour.
//It's not completely exhaustive. It aborts if it takes too long (>N seconds), specially because ending moves also permutates by sign.
//If you have an ending number that comes from 8 moves, you have 8! permutations, and then inside each one (1<<8)-1 signs. That's a lot.
//I don't keep stats of this Mutator because it's a one time use.
struct EXHSender {
	int indexMoveList = -1;
	int senderIDX;
	int endValue;
	vector<int> movePos;
	vector<Move> moves;
	vector<Move> targetMoves;
	void create(const LAHC_Node& currentBest, int _destIDX, int _indexMoveList) {
		indexMoveList = _indexMoveList;
		auto& LM = currentBest.grid.Plan_NMB[indexMoveList];
		senderIDX = _destIDX;
		endValue = currentBest.grid.Val[senderIDX];
		//for (auto&m : LM.linear)
		movePos.resize(0);
		moves.resize(0);
		targetMoves.resize(0);
		//for (auto&m:LM.linear)
		for (int i = 0; i < LM.linear.size(); ++i)
		{
			auto&m = LM.linear[i];
			if (m.destIDX == senderIDX)
			{
				movePos.push_back(i);
			}
		}
	}
};
struct EXHMerger {
	int indexMoveList = -1;
	int mergerIDX;
	int endValue;
	vector<int> interValue;
	uint64_t maskValues;
	vector<Move> moves;
	vector<int> movePos;
	int lastPosMove;
	void create(const LAHC_Node& currentBest, int _destIDX) {
		mergerIDX = _destIDX;
		maskValues = 0;
		const Strategy* LM = nullptr;
		indexMoveList = -1;
		{
			//Search the list
			for (auto& M : currentBest.Movelists)
			{
				for (auto&m : currentBest.grid.Plan_NMB[M].linear)
				{
					if (m.srcIDX == mergerIDX || m.destIDX == mergerIDX) {
						indexMoveList = M;
						break;
					}
				}
				if (indexMoveList >= 0)
					break;
			}
		}
		if (indexMoveList >= 0)
			LM = &currentBest.grid.Plan_NMB[indexMoveList];
		moves.resize(0);
		movePos.resize(0);
		interValue.resize(0);
		lastPosMove = 0;
		if (LM != nullptr)
			for (int llm = 0; llm < LM->linear.size(); ++llm)
			{
				auto& m = LM->linear[llm];
				if (m.destIDX == mergerIDX)
				{
					movePos.push_back(llm);
					lastPosMove = max(lastPosMove, llm);
				}
			}
	}
};
void permuteEXHSender(const LAHC_Node& currentBest, vector<EXHSender>& src, unordered_map<int, vector<EXHSender>>& dest) {
	dest.clear();
	Grid G;
	vector<Move>tmpMoves;
	Stopwatch muchoTiempo;
	if (currentBest.grid.totalNumbers == 1)
		muchoTiempo.Start(10 * 1000 * 1000);
	else if (currentBest.grid.totalNumbers > 8)
		muchoTiempo.Start(2 * 100 * 1000);
	else muchoTiempo.Start(4 * 1000 * 1000);
	for (auto&E : src)
	{
		dest[E.senderIDX].resize(0);
		vector<EXHSender>& SS = dest[E.senderIDX];
		if (E.indexMoveList < 0)
		{ //Direct pass, no permutations
			EXHSender nuevaPerm = E;
			nuevaPerm.endValue = MAP.Val[E.senderIDX];
			SS.push_back(nuevaPerm);
			continue;
		}
		bool JackPot0 = false;
		do {
			uint64_t signos;
			if (currentBest.grid.totalNumbers == 1)
				signos = (1ULL << (int)E.movePos.size());
			else if (currentBest.grid.totalNumbers > 8)
				signos = (1ULL << min((E.movePos.size()>6 ? 3 : 4), (int)E.movePos.size()));
			else signos = (1ULL << min(6, (int)E.movePos.size()));
			if (muchoTiempo.Timeout())
				return;

			for (uint64_t s = 0; s < signos; ++s) {
				if (solved)
					return;
				EXHSender nuevaPerm;
				nuevaPerm.indexMoveList = E.indexMoveList;
				nuevaPerm.senderIDX = E.senderIDX;
				nuevaPerm.endValue = MAP.Val[E.senderIDX];
				for (int i = 0; i < E.movePos.size(); ++i) {
					Move m = currentBest.grid.Plan_NMB[E.indexMoveList].linear[E.movePos[i]];
					m.sign = (((1ULL << i) & s) == 0) ? 0 : 1;
					m.destVal = nuevaPerm.endValue;
					nuevaPerm.movePos.push_back(E.movePos[i]);
					nuevaPerm.moves.push_back(m);
					//Apply new value to get intermediate values and endValue
					nuevaPerm.endValue = (m.sign == 0 ? abs(m.destVal - m.srcVal) : m.destVal + m.srcVal);

					if (nuevaPerm.endValue == 0)
						break;
				}
				if (nuevaPerm.endValue == 0)
				{
					if (nuevaPerm.moves.size() != E.movePos.size())
						continue;
					//JACKPOT! Auto solved
					SS.resize(0);
					SS.push_back(nuevaPerm);
					JackPot0 = true;
					break;
				}
				if (nuevaPerm.endValue >= W) //Impossible to join
					continue;

				bool duplicado = false;
				for (auto&P : SS)
				{
					if (P.endValue == nuevaPerm.endValue)
					{
						duplicado = true;
						break;
					}
				}
				if (duplicado)
					continue;
				//With endValue check possible targets;
				nuevaPerm.targetMoves.resize(0);
				G = MAP;
				for (auto& N : currentBest.grid.Plan_NMB[nuevaPerm.indexMoveList].linear)
				{
					G.Val[N.srcIDX] = 0;
					G.Val[N.destIDX] = 0;
				}
				G.Val[E.senderIDX] = nuevaPerm.endValue;
				tmpMoves.resize(0);
				Random rr;
				G.ExplodeMoves(E.senderIDX, tmpMoves, nullptr);
				nuevaPerm.targetMoves.resize(0);
				for (auto&m : tmpMoves)
					if (m.sign == 0)
					{
						nuevaPerm.targetMoves.push_back(m);
					}
				if (nuevaPerm.targetMoves.size() > 0) //Possible connection because I can target a destination IDX
				{
					SS.push_back(nuevaPerm);
				}
			}

			if (JackPot0) break;
		} while (next_permutation(E.movePos.begin(), E.movePos.end()));

		//Validar no duplicados
		/*		if (SS.size() > 1)
		for (int i = 0; i < SS.size() - 1; ++i)
		{
		for (int j = i + 1; j < SS.size(); ++j) {
		if (SS[i].endValue == SS[j].endValue)
		{
		cerr << "error duplicado" << endl;
		}
		}
		}*/

	}
}
void permuteEXHMerger(const LAHC_Node& currentBest, vector<EXHMerger>& src, vector<EXHMerger>& dest, int searchMergeValue, int forceENDVALUE = -1)
{
	Stopwatch muchoTiempo;
	if (currentBest.grid.totalNumbers == 1)
		muchoTiempo.Start(10 * 1000 * 1000);
	else if (currentBest.grid.totalNumbers > 8)
		muchoTiempo.Start(2 * 100 * 1000);
	else muchoTiempo.Start(4 * 1000 * 1000);
	dest.resize(0);
	for (auto&E : src)
	{
		do {
			uint64_t signos;
			if (currentBest.grid.totalNumbers == 1)
				signos = (1ULL << (int)E.movePos.size());
			else if (currentBest.grid.totalNumbers > 8)
				signos = (1ULL << min((E.movePos.size()>6 ? 3 : 4), (int)E.movePos.size()));
			else signos = (1ULL << min(6, (int)E.movePos.size()));
			if (muchoTiempo.Timeout())
				return;
			for (uint64_t s = 0; s < signos; ++s)
			{
				if (solved)
					return;
				EXHMerger nuevaPerm;
				nuevaPerm.indexMoveList = E.indexMoveList;
				nuevaPerm.mergerIDX = E.mergerIDX;
				nuevaPerm.endValue = MAP.Val[E.mergerIDX];
				nuevaPerm.maskValues = 0;
				nuevaPerm.lastPosMove = E.lastPosMove;
				for (int i = 0; i < E.movePos.size(); ++i) {
					Move m = currentBest.grid.Plan_NMB[E.indexMoveList].linear[E.movePos[i]];
					m.sign = (((1ULL << i) & s) == 0) ? 0 : 1;
					m.destVal = nuevaPerm.endValue;
					nuevaPerm.movePos.push_back(E.movePos[i]);
					nuevaPerm.moves.push_back(m);
					//Apply new value to get intermediate values and endValue
					nuevaPerm.endValue = (m.sign == 0 ? abs(m.destVal - m.srcVal) : m.destVal + m.srcVal);
					nuevaPerm.interValue.push_back(nuevaPerm.endValue);
					if (nuevaPerm.endValue == 0)
						break;
				}

				if (nuevaPerm.endValue == 0 && nuevaPerm.moves.size() < E.movePos.size())
				{ //Invalid merge
					continue;
				}
				if (forceENDVALUE >= 0)
				{
					if (nuevaPerm.endValue == forceENDVALUE)
					{
						for (auto&I : nuevaPerm.interValue)
							if (searchMergeValue == 2 * I && searchMergeValue != 0)
							{
								dest.push_back(nuevaPerm); //Jackpot!
								return;
							}
					}
				}
				else
				{
					for (auto&I : nuevaPerm.interValue)
					{
						if (searchMergeValue == 2 * I && searchMergeValue != 0)
						{
							dest.push_back(nuevaPerm);
							return;
						}
					}

				}

			}
		} while (next_permutation(E.movePos.begin(), E.movePos.end()));

	}
}
bool Mutate_Exhaustive_6(const LAHC_Node& currentBest, LAHC_Node& candidate, LAHC_Node& tmpCandidate, Random& rnd)
{
	bool mejorado = false;
	vector<EXHSender> incompletes;
	//Untouched numbers
	for (auto&N : NMB)
	{
		if (currentBest.untouchedNumbers.get(N.ID) && currentBest.grid.Val[N.ID] > 0) {
			EXHSender ad;
			ad.senderIDX = N.ID;
			ad.endValue = currentBest.grid.Val[ad.senderIDX];
			incompletes.push_back(ad);
		}
	}
	//Incomplete Sets
	for (int i = 0; i < currentBest.IncompleteSets.size(); ++i)
	{
		const Move& lastMove = currentBest.grid.Plan_NMB[currentBest.IncompleteSets[i]].linear.back();
		EXHSender ad;
		ad.create(currentBest, lastMove.destIDX, currentBest.IncompleteSets[i]);
		incompletes.push_back(ad);
	}

	unordered_map<int, vector<EXHSender>> permIncompletes;
	permuteEXHSender(currentBest, incompletes, permIncompletes);
	int totalPermutations = 0;

	vector<int> numeros;
	for (auto&I : permIncompletes)
	{
		numeros.push_back(I.first);
		totalPermutations += (int)I.second.size();
		std::sort(I.second.begin(), I.second.end(), [=](auto& a, auto& b) {return b.endValue > a.endValue; });
	}

	if (numeros.size() > 1)
		for (int parteA = 0; parteA < numeros.size() - 1; ++parteA)
		{
			for (int parteB = parteA + 1; parteB < numeros.size(); ++parteB)
			{
				int IDX_A = numeros[parteA];
				int IDX_B = numeros[parteB];
				if (IDX_A == IDX_B)
					continue;
				if (NMB[IDX_A].X != NMB[IDX_B].X && NMB[IDX_A].Y != NMB[IDX_B].Y)
					continue;
				if (CACHE_MERGE[IDX_A][IDX_B] == nullptr || CACHE_MERGE[IDX_B][IDX_A] == nullptr)
					continue;
				if (permIncompletes[IDX_A].size() == 0 || permIncompletes[IDX_B].size() == 0)
					continue;
				//Neighbours
				int distance = CACHE_MERGE[IDX_A][IDX_B]->srcVal;
				int bestValue_A = permIncompletes[IDX_A][0].endValue;
				int bestIndex_A = 0;

				for (int i = 0; i < permIncompletes[IDX_A].size(); ++i)
				{
					if (permIncompletes[IDX_A][i].endValue == distance)
					{ //JACKPOT!
						bestValue_A = permIncompletes[IDX_A][i].endValue;
						bestIndex_A = i;
						break;
					}
					//nearest, to have the lowest end number
					if (abs(permIncompletes[IDX_A][i].endValue - distance) < abs(bestValue_A - distance))
					{
						bestValue_A = permIncompletes[IDX_A][i].endValue;
						bestIndex_A = i;
					}
				}

				int bestValue_B = permIncompletes[IDX_B][0].endValue;
				int bestIndex_B = 0;
				for (int i = 0; i < permIncompletes[IDX_B].size(); ++i)
				{
					if (permIncompletes[IDX_B][i].endValue == distance)
					{ //JACKPOT!
						bestValue_B = permIncompletes[IDX_B][i].endValue;
						bestIndex_B = i;
						break;
					}
					//nearest, to have the lowest end number
					if (abs(permIncompletes[IDX_B][i].endValue - distance) < abs(bestValue_B - distance))
					{
						bestValue_B = permIncompletes[IDX_B][i].endValue;
						bestIndex_B = i;
					}
				}
				if (bestValue_A == distance || bestValue_B == distance)
				{ //JACKPOT!
					EXHSender JJ[2] = { permIncompletes[IDX_A][bestIndex_A],permIncompletes[IDX_B][bestIndex_B] };
					tmpCandidate.clear();
					for (int nn = 0; nn < 2; ++nn)
					{
						if (JJ[nn].indexMoveList > 0)
						{
							auto& LM = currentBest.grid.Plan_NMB[JJ[nn].indexMoveList];
							//Redo all moves on the incomplete.
							for (auto&m : LM.linear)
							{
								if (m.srcIDX == JJ[nn].senderIDX || m.destIDX == JJ[nn].senderIDX)
									continue;
								tmpCandidate.grid.ApplyMove(m);
							}
						}
						//Change the ending to the permutations
						for (auto&m : JJ[nn].moves)
							if (m.destIDX < COUNT_NUMBERS)
							{
								m.destVal = tmpCandidate.grid.Val[m.destIDX];
								if (tmpCandidate.grid.ValidateMove(m))
									tmpCandidate.grid.ApplyMove(m);
							}
					}
					//Now mix them
					Move UUm = (bestValue_A == distance ? *CACHE_MERGE[IDX_A][IDX_B] : *CACHE_MERGE[IDX_B][IDX_A]);
					UUm.srcVal = tmpCandidate.grid.Val[UUm.srcIDX];
					UUm.destVal = tmpCandidate.grid.Val[UUm.destIDX];
					UUm.sign = 0;
					if (tmpCandidate.grid.ValidateMove(UUm))
						tmpCandidate.grid.ApplyMove(UUm);
					//Do all other moves
					for (auto&N : currentBest.Movelists)
						if (N != JJ[0].indexMoveList && N != JJ[1].indexMoveList)
						{
							auto& LM = currentBest.grid.Plan_NMB[N];
							for (auto&m : LM.linear)
							{
								tmpCandidate.grid.doMerge(UUm.destIDX, m.srcIDX, rnd);
								tmpCandidate.grid.doMerge(UUm.destIDX, m.destIDX, rnd);
								tmpCandidate.grid.ApplyMove(m);
								tmpCandidate.grid.doMerge(UUm.destIDX, m.srcIDX, rnd);
								tmpCandidate.grid.doMerge(UUm.destIDX, m.destIDX, rnd);
							}
						}
					while (tmpCandidate.grid.doRandomMove(rnd)) {}
					tmpCandidate.MutateType = 6;
					tmpCandidate.CalcStats();
					if (tmpCandidate.Movelists.size() == 0)
						continue;
					++SUPERCOMBINATOR;
					GLOBAL_RN.RN_Insert(tmpCandidate, false);
					if (tmpCandidate.BestScore > candidate.BestScore)
					{
						swap(candidate, tmpCandidate);
					}
					if (candidate.BestScore > currentBest.BestScore)
					{
						mejorado = true;
						++JACKPOT_OK;
					}
					if (IS_SOLVED(currentBest, candidate, "JACKPOT Inc-Inc"))
						return true;
				}

			}
		}
	/*if (currentBest.grid.totalNumbers > K_EXHAUSTIVE_MIN_NUMBERS + 4)
	return false;*/
	int LIMIT_PERM = 1500;

	if (MAX_PERMUTATIONS < totalPermutations)
		MAX_PERMUTATIONS = totalPermutations;
	if (totalPermutations > LIMIT_PERM)
	{  //what to do?
		double factor = (double)LIMIT_PERM / (double)totalPermutations;
		for (auto&I : permIncompletes)
		{
			int targetSize = max(1, (int)((double)I.second.size()*factor));
			if (targetSize < I.second.size())
				I.second.resize(targetSize);
		}
	}
	++JACKPOT_TOTAL;
	for (auto& AP : permIncompletes)
		//for (EXHSender& P : AP.second)
		for (int esI = 0; esI<AP.second.size(); ++esI)
		{
			EXHSender& P = AP.second[esI];
			if (solved)
				return false;
			//**********Self Jackpot************************
			if ((esI == 0 || P.endValue == 0) && P.indexMoveList >= 0)
			{

				tmpCandidate.clear();
				vector<int> copyLists = currentBest.Movelists;
				auto& LM = currentBest.grid.Plan_NMB[P.indexMoveList];
				//Redo all moves on the incomplete.
				for (auto&m : LM.linear)
				{
					if (m.srcIDX == P.senderIDX || m.destIDX == P.senderIDX)
						continue;
					tmpCandidate.grid.ApplyMove(m);
				}
				//Change the ending to the permutations
				for (auto&m : P.moves)
					if (m.destIDX < COUNT_NUMBERS)
					{
						m.destVal = tmpCandidate.grid.Val[m.destIDX];
						tmpCandidate.grid.ApplyMove(m);
					}
				for (auto&N : currentBest.Movelists)
					if (N != P.indexMoveList)
					{
						auto& LISTA = currentBest.grid.Plan_NMB[N];
						for (auto&m : LISTA.linear)
						{
							tmpCandidate.grid.ApplyMove(m);
						}
					}
				while (tmpCandidate.grid.doRandomMove(rnd)) {}
				tmpCandidate.MutateType = 6;
				tmpCandidate.CalcStats();
				if (tmpCandidate.Movelists.size() == 0)
					continue;
				GLOBAL_RN.RN_Insert(tmpCandidate, false);
				if (tmpCandidate.BestScore > candidate.BestScore)
				{
					swap(candidate, tmpCandidate);
				}
				if (candidate.BestScore > currentBest.BestScore)
				{
					mejorado = true;
					if (P.endValue == 0)
						++JACKPOT_OK;
				}
				if (IS_SOLVED(currentBest, candidate, "JACKPOT Inc-Self0"))
					return true;
			}
			for (auto&t : P.targetMoves)
			{
				//Potential point to merge, outside my IncompleteList and value is valid to merge.
				vector<EXHMerger> targetMerge;
				targetMerge.resize(1);
				targetMerge[0].create(currentBest, t.destIDX);
				vector<EXHMerger> explodeTargetMerge;
				if (targetMerge[0].indexMoveList > 0)
				{
					//Search end value of potential merge point
					int expectedEndValue = -1;
					for (auto&uum : currentBest.grid.Plan_NMB[targetMerge[0].indexMoveList].linear)
					{
						if (uum.srcIDX == t.destIDX)
						{
							expectedEndValue = uum.srcVal;
							break;
						}
					}
					permuteEXHMerger(currentBest, targetMerge, explodeTargetMerge, P.endValue, expectedEndValue);
					if (explodeTargetMerge.size() > 0)
					{
						++LIMIT_PERM;
						--LIMIT_PERM;
					}
				}
				if (explodeTargetMerge.size() > 0)
				{ //JACKPOT
					tmpCandidate.clear();
					vector<int> copyLists = currentBest.Movelists;
					if (P.indexMoveList > 0)
					{
						auto& incompleteLM = currentBest.grid.Plan_NMB[P.indexMoveList];
						//Redo all moves on the incomplete.
						for (auto&m : incompleteLM.linear)
						{
							if (m.srcIDX == P.senderIDX || m.destIDX == P.senderIDX)
								continue;
							tmpCandidate.grid.ApplyMove(m);
						}
					}
					//Change the ending to the permutations
					for (auto&m : P.moves)
						if (m.destIDX < COUNT_NUMBERS)
						{
							m.destVal = tmpCandidate.grid.Val[m.destIDX];
							tmpCandidate.grid.ApplyMove(m);
						}
					//mark incomplete list as done
					for (int i = 0; i < copyLists.size(); ++i)
					{
						if (copyLists[i] == P.indexMoveList)
						{
							copyLists[i] = -1;
							break;
						}
					}
					//Backtrack on the targetMerge.
					EXHMerger& dest = explodeTargetMerge[0];
					if (dest.indexMoveList >= 0)
					{
						//for (auto&m :currentBest.grid.Plan_NMB[dest.indexMoveList].linear)
						for (int i = 0; i < dest.lastPosMove; ++i)
						{
							auto&m = currentBest.grid.Plan_NMB[dest.indexMoveList].linear[i];
							for (auto&ignore : dest.moves)
								if (ignore.srcIDX == m.srcIDX || ignore.srcIDX == dest.mergerIDX || ignore.destIDX == dest.mergerIDX)
									continue;
							tmpCandidate.grid.ApplyMove(m);
						}
						for (auto&m : dest.moves)
							if (m.destIDX < COUNT_NUMBERS)
							{
								m.destVal = tmpCandidate.grid.Val[m.destIDX];
								tmpCandidate.grid.doMerge(P.senderIDX, m.srcIDX, rnd);
								tmpCandidate.grid.doMerge(P.senderIDX, m.destIDX, rnd);
								tmpCandidate.grid.ApplyMove(m);
								tmpCandidate.grid.doMerge(P.senderIDX, m.srcIDX, rnd);
								tmpCandidate.grid.doMerge(P.senderIDX, m.destIDX, rnd);
							}
						for (auto&m : currentBest.grid.Plan_NMB[dest.indexMoveList].linear) {
							tmpCandidate.grid.ApplyMove(m);
						}
						//mark incomplete list as done
						for (int i = 0; i < copyLists.size(); ++i)
						{
							if (copyLists[i] == dest.indexMoveList)
							{
								copyLists[i] = -1;
								break;
							}
						}
					}
					//Merge
					//Complete targetMerge
					//Do all other moveLists
					for (auto& N : copyLists)
						if (N >= 0)
						{
							auto& doL = currentBest.grid.Plan_NMB[N];
							for (auto&m : doL.linear)
							{
								tmpCandidate.grid.ApplyMove(m);
							}
						}
					while (tmpCandidate.grid.doRandomMove(rnd)) {
					}
					tmpCandidate.MutateType = 6;
					tmpCandidate.CalcStats();
					if (tmpCandidate.Movelists.size() == 0)
						continue;
					GLOBAL_RN.RN_Insert(tmpCandidate, false);
					if (tmpCandidate.BestScore > candidate.BestScore)
					{
						swap(candidate, tmpCandidate);
					}
					if (candidate.BestScore > currentBest.BestScore)
					{
						mejorado = true;
						++JACKPOT_OK;
					}
					if (IS_SOLVED(currentBest, candidate, "JACKPOT Inc-Merge"))
						return true;
				}
			}

		}
	return mejorado;
}


//***************** On Some Cases Mutators left a bad candidate, and in some cases I can save them by Applying Incompletes first and then completes.**************/
//Imo it's an ugly fix.
inline bool reArmCandidate(LAHC_Node& currentBest, LAHC_Node& tmpCandidate, Random& rnd) {
	vector<int> UnusedNumbers;
	//for (auto& n : NMB)
	for (int n = 0; n < NMB.size(); ++n)
	{
		if (currentBest.grid.Val[n] > 0)
			UnusedNumbers.push_back(n);
	}
	tmpCandidate.clear();
	vector<int> INCop = currentBest.IncompleteSets;
	do_shuffle<vector<int>>(INCop, rnd);
	for (auto& I : INCop)
	{
		auto& UEL = currentBest.grid.Plan_NMB[I];
		for (auto&m : UEL.linear)
		{
			bool intentaMerge = rnd.NextInt(1000) < 850;
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.srcIDX, rnd);
			tmpCandidate.grid.ApplyMove(m);
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.destIDX, rnd);
		}
	}
	INCop.resize(0);
	vector<int> COMP = currentBest.ZerosumSets;
	do_shuffle<vector<int>>(COMP, rnd);
	for (auto& I : COMP)
	{
		auto& UEL = currentBest.grid.Plan_NMB[I];
		for (auto&m : UEL.linear)
		{
			bool intentaMerge = rnd.NextInt(1000) < 850;
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.srcIDX, rnd);
			tmpCandidate.grid.ApplyMove(m);
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.destIDX, rnd);

		}
	}
	COMP.resize(0);
	while (tmpCandidate.grid.doRandomMove(rnd)) {}
	tmpCandidate.MutateType = currentBest.MutateType;
	tmpCandidate.CalcStats();
	return true;
}


/***********************************************************************************/
// Candidate generation. It applies some mutator to the currentBest Node to create a new candidate
// Each mutator n have a PROB[n] value of how good it was in the past, so Mutators that had more success
// will be used more often
/***********************************************************************************/
inline void Mutate(LAHC_Node& currentBest, LAHC_Node& candidate, double PROB[], Random& rnd)
{
	candidate.clear();
	if (currentBest.Movelists.size() > 0)
	{
		//Calculate sum of probabilites
		double P0 = (currentBest.ZerosumSets.size() > 0 ? PROB[0] : 0.0);
		double P2 = PROB[2] + P0;
		double P3 = PROB[3] + P2;
		double P4 = PROB[4] + P3;
		double P5 = (currentBest.IncompleteSets.size() > 0 ? PROB[5] : 0.0) + P4;
		double PTOTAL = PROB[1] + P5;

		double r = rnd.NextFloat()*PTOTAL;//  (double)rnd.NextInt((int)(PTOTAL*10000.0)) / 10000.0;
										  //Weighted Random
		if (r < P0) { Mutate_Completes_0(currentBest, candidate, rnd); }
		else if (r < P2) { Mutate_Lists_2(currentBest, candidate, rnd); }
		else if (r < P3) { Mutate_AdjPoints_3(currentBest, candidate, rnd); }
		else if (r < P4) { Mutate_Complex_4(currentBest, candidate, rnd); }
		else if (r < P5) { Mutate_Shuffle_Incomplete_End_5(currentBest, candidate, rnd); }
		//If nothing was done then defaults to simplest mutator.
		if (candidate.MutateType == 9)
			Mutate_Points_1(currentBest, candidate, rnd);
	}
	//Do any possible random Move it left.
	while (candidate.grid.doRandomMove(rnd)) {}
	//Calculate score, lists of moves, remaining values,etc.
	candidate.CalcStats();
}

//Old mode loader. It can load a SOLUTION_ file and add APX to Global_Best
void LoadRecombinators_LAHC(string s)
{
	vector<LAHC_Node> Recombinators;
	mutexSAVE.lock();
	std::ifstream infile(s);
	string line;
	if (infile.is_open())
	{
		Recombinators.resize(0);
		LAHC_Node newCandidate;
		newCandidate.clear();
		bool StartFound = false;
		while (getline(infile, line))
		{
			if (StartFound == false) //Seek
			{
				if (line.size() >= 12 &&
					(((line[0] == 'C') && (line[1] == 'o') && (line[2] == 'm'))
						|| ((line[0] == 'I') && (line[1] == 'n') && (line[2] == 'c'))
						)
					) {
					StartFound = true;
				}
			}
			if (StartFound)
			{
				if (line.size() >= 12 &&
					(((line[0] == 'C') && (line[1] == 'o') && (line[2] == 'm'))
						|| ((line[0] == 'I') && (line[1] == 'n') && (line[2] == 'c'))
						)
					) {
					//Add a new Movelist.
					vector<string> moveList = split(line, "|");
					if (moveList.size() > 1)
						for (int i = 1; i < moveList.size(); ++i)
						{
							vector<string> ACT = split(moveList[i], " "); //20 18 L -

							if (ACT.size() == 4 && ACT[0][0] != '*') {

								auto srcX = stoi(ACT[0]);
								auto srcY = stoi(ACT[1]);
								if (srcX >= W || srcY >= H || srcX < 0 || srcY < 0)
									continue;
								auto srcID = IDX[srcX][srcY];
								if (srcID >= COUNT_NUMBERS)
									continue;
								int dir = 0;
								for (int d = 0; d < 4; ++d)
									if (strDIRS[d][0] == ACT[2][0])
									{
										dir = d;
										break;
									}
								int sign = (ACT[3][0] == '-' ? 0 : 1);
								auto srcVal = newCandidate.grid.Val[srcID];
								int destX = srcX + DIR_X[dir] * srcVal;
								int destY = srcY + DIR_Y[dir] * srcVal;
								if (destX >= W || destY >= H || destX < 0 || destY < 0)
									continue;
								auto destID = IDX[destX][destY];
								auto destScore = newCandidate.grid.Val[destID];
								Move m = Move(srcID, destID, dir, sign, srcVal, destScore);
								if (newCandidate.grid.ValidateMove(m))
									newCandidate.grid.ApplyMove(m);
							}
						}
				}
				else {
					StartFound = false;
					//Close a current Candidate
					if (newCandidate.grid.totalPoints != MAP.totalPoints && newCandidate.grid.totalNumbers < 5)
					{
						newCandidate.CalcStats();
						Recombinators.push_back(newCandidate);
						newCandidate.clear();
					}
				}
			}
		}
		infile.close();

	}
	mutexSAVE.unlock();
	//Prefer less points
	std::sort(std::begin(Recombinators), std::end(Recombinators), [](auto& a, auto& b) {return a.grid.totalNumbers < b.grid.totalNumbers; });

	for (auto& R : Recombinators)
	{
		GLOBAL_RN.RN_Insert(R, true);
	}
	Recombinators.resize(0);
}

//APX Loader. It's used to restart a partial search or to inject APX values from other nodes.
void LoadRecombinators_APROX(string s, bool firstStart)
{
	mutexSAVE.lock();
	long tmpFilesize = GetFileSize(s);
	if (tmpFilesize <= 0 || filesize == tmpFilesize)
	{
		mutexSAVE.unlock();
		return;
	}
	filesize = tmpFilesize;
	std::ifstream infile(s);
	string line;
	int total = 0;
	int success = 0;
	if (infile.is_open())
	{
		LAHC_Node newCandidate;
		newCandidate.clear();
		bool StartFound = false;
		while (getline(infile, line))
		{
			vector<string> moveList = split(line, "|");
			if (moveList.size() > 1)
				for (int i = 0; i < moveList.size(); ++i)
				{
					vector<string> ACT = split(moveList[i], " ");
					if (ACT.size() >= 6)
					{
						Move m = Move(stoi("0" + ACT[0]), stoi("0" + ACT[1]), stoi("0" + ACT[2]), stoi("0" + ACT[3]), stoi("0" + ACT[4]), stoi("0" + ACT[5]));
						if (newCandidate.grid.ValidateMove(m))
							newCandidate.grid.ApplyMove(m);
					}
				}
			newCandidate.CalcStats();
			if (newCandidate.grid.totalPoints != MAP.totalPoints && newCandidate.grid.totalNumbers < 5)
			{
				++total;

				if (GLOBAL_RN.RN_Insert(newCandidate, false))
					++success;
				newCandidate.clear();
			}
		}
		infile.close();
	}
	mutexSAVE.unlock();
	if (success > 0)
		cerr << " Load Recombinators APROX ok/total:" << success << "/" << total << endl;
	if (firstStart && (success > 4) && (success < total - 10) && swapDone == 0)
	{
		swapDone = 1;
	}
}

//Similar to Mutate(), but with fixed probabilities: 25% Complex, 35% Incomplete, 40% other Mutators
void Recombinate(LAHC_Node currentBest, LAHC_Node& candidate, Random& rnd) {


	candidate.clear();
	int REN = rnd.NextInt(1000);
	if (REN < 250)
	{
		Mutate_Complex_4(currentBest, candidate, rnd);
	}
	else if (REN < 600)
	{
		Mutate_Shuffle_Incomplete_End_5(currentBest, candidate, rnd);
	}
	else
	{

		int Cycle = rnd.NextInt(5);
		switch (Cycle)
		{
		case 0:Mutate_Completes_0(currentBest, candidate, rnd);	break;
		case 1:Mutate_Points_1(currentBest, candidate, rnd);	break;
		case 2:Mutate_Lists_2(currentBest, candidate, rnd);	break;
		case 3:Mutate_AdjPoints_3(currentBest, candidate, rnd);	break;
		case 4:Mutate_Cross_4(currentBest, candidate, rnd);	break;
		}
	}
	while (candidate.grid.doRandomMove(rnd)) {}
	candidate.CalcStats();
}


//Heavily modified version of a LAHC Worker.
void Worker_LAHC(int ID)
{
	//I'll keep a queue of 25 APX. This is for Flashback purposes.
	const int TAM_HIS = 25;//15;
	ALIGN deque<LAHC_Node> historyBEST;
	int realHistoryLength = 0;
	int TOTAL_FLASH = 0;
	int SUCCESS_FLASH = 0;
	int Flashbacks = -1;
	Random rnd;
	if (randomKey != 0)
		rnd = Random(randomKey + 50 + ID);

	//Global Best Search Algorithm variables
	bool canBeGlobalBest = (ID < ALLOW_GLOBAL_BEST);
	bool isGlobalBestActive = false;
	RN_Explorer LOCAL_RN(false);

	//Variables for statitical success of each mutator Type.
	double MUTATE_STATS[10] = { 0 };
	double MUTATE_BESTSTATS[10] = { 0 };
	double PROB[10] = { 10.0,15.0,15.0,10.0,30.0,30.0,0.0,0.0,0.0,0.0 };

	//Dynamic LFA, variable with coverage percent.
	const int SMALL_LFA = 150;
	int FINAL_LFA = SIZE_LFA;
	int EXPAND_LFA = FINAL_LFA + 20 * CAN_INCREASE_LFA;
	int CurrLFA = SMALL_LFA; //Fast convergence, low diversity

							 //	cerr << (canBeGlobalBest ? "RECOMBINATOR " : "LAHC ") << ID << " Lfa size:" << FINAL_LFA << " EXP:" << EXPAND_LFA << endl;

							 //LAHC related candidates
	LAHC_Node lastAccepted;
	LAHC_Node candidate;
	uint64_t v = 0; //First iteration I=0;
	uint64_t I = 0; //First iteration I=0;
	vector<double> Fitness;

	//These aren't from standard LAHC
	LAHC_Node tmpCandidate;
	LAHC_Node bestStrategy;
	uint64_t previousI = 0;

	//Timers
	long long localTimeLimit = LIMIT_TIME_IMPROVEMENT;
	Stopwatch syncClock;
	syncClock.Start(0); //Force timeout at start
	Stopwatch notifyClock;
	notifyClock.Start(30 * 1000 * 1000);
	Stopwatch convergenceClock;
	convergenceClock.Start(10 * 1000 * 1000);
	Stopwatch lastImprovement;
	lastImprovement.Start(0);
	Stopwatch timeFromLastReset;
	timeFromLastReset.Start(0);
	long long maxTimeImprovement = 0;

	bool copiedFromGlobal = false;
	{ //Produce an initial solution s
		Flashbacks = -1;
		historyBEST.clear();
		lastAccepted.clear();
		Mutate(lastAccepted, candidate, PROB, rnd);
		lastAccepted = candidate;
		Fitness.resize(EXPAND_LFA);//We have an expanded size of the array from the start.
		for (int k = 0; k < FINAL_LFA; ++k)
		{
			Fitness[k] = lastAccepted.BestScore; //For all k € {0...Lfa-1} fk:=C(s)
		}
		CurrLFA = SMALL_LFA; //Fast convergence
		bestStrategy = lastAccepted;
		historyBEST.push_back(bestStrategy);
		++realHistoryLength;
	}


	while (!solved) //Do until a chosen stopping condition
	{
		//Global Best Search Algorithm - Sync Best APX each 10 secs
		if (
			(ID == 0 || canBeGlobalBest)
			&& (ID == 0
				|| bestStrategy.grid.totalNumbers <= MIN_DEGREE_TOINSERT
				|| GLOBAL_RN.Degree[0] < 3
				|| (ID == 1 && GLOBAL_RN.Degree[0] < 6)
				|| (ID == 2 && GLOBAL_RN.Degree[0] < 5)
				|| (ID == 3 && GLOBAL_RN.Degree[0] < 4)
				)
			&& (((I % 1000) == 0) && syncClock.Timeout()))
		{
			syncClock.Start(30 * 1000 * 1000);

			//Avoid race conditions
			mutexRN_Explorer.lock();

			if (LOCAL_RN.totalRN>0 && abs(LOCAL_RN.Degree[0] - GLOBAL_RN.Degree[0])<MAX_DEGREE_TOSAVE + 2)
			{
				for (auto& V : GLOBAL_RN.APX)
					for (auto& A : V)
					{
						LOCAL_RN.RN_Insert(A, false);
					}
			}
			else LOCAL_RN = GLOBAL_RN;
			LOCAL_RN.isGLOBAL = false;
			mutexRN_Explorer.unlock();
			//Change LAHC worker to Global Best:
			// ID=0. When there is any Global APX available
			// ID=1. Exists 1 APX of quality 1-5
			// ID=2. Exists 1 APX of quality 1-4
			// ID>2. Quality 1-4 count is greater than COUNT_RECOMB
			if (!isGlobalBestActive && canBeGlobalBest)
			{
				if (ID == 0) isGlobalBestActive = (LOCAL_RN.totalRN > 0);
				else if (ID == 1) {
					if (MAX_INSERT_NUMBERS >= 6)
						isGlobalBestActive = (LOCAL_RN.totalRN > 0 && LOCAL_RN.Degree[0] < 6);
					else isGlobalBestActive = (LOCAL_RN.SizeAcum.size() > 0 && LOCAL_RN.SizeAcum[min(4, (int)LOCAL_RN.SizeAcum.size() - 1)] > CNT_RECOMB);
				}
				else if (ID == 2) {
					if (MAX_INSERT_NUMBERS >= 5)
						isGlobalBestActive = (LOCAL_RN.totalRN > 0 && LOCAL_RN.Degree[0] < 5);
					else isGlobalBestActive = (LOCAL_RN.SizeAcum.size() > 0 && LOCAL_RN.SizeAcum[min(4, (int)LOCAL_RN.SizeAcum.size() - 1)]> CNT_RECOMB);
				}
				else isGlobalBestActive = (LOCAL_RN.SizeAcum.size() > 0 && LOCAL_RN.SizeAcum[min(4, (int)LOCAL_RN.SizeAcum.size() - 1)] > CNT_RECOMB);
			}


			//ID 0 does the Mutate_Exhaustive_6. Done locally, then translated to GLOBAL,
			if (isGlobalBestActive && ID == 0) {
				syncClock.Start(30 * 1000 * 1000); //Restart it just to have some timer for how long it takes the mutator.
				int runTotal = 0;
				int runSuccess = 0;
				candidate.clear();
				candidate.grid.doRandomMove(rnd);
				candidate.CalcStats();
				tmpCandidate.clear();
				for (int d = 0; d < LOCAL_RN.APX.size(); ++d)
					if (LOCAL_RN.APX[d].size() > 0)
					{
						if (LOCAL_RN.APX[d][0].grid.totalNumbers > 10 && syncClock.EllapsedMilliseconds() > 1000)
							break;
						for (int i = 0; i < LOCAL_RN.APX[d].size(); ++i)
						{
							if (LOCAL_RN.APX[d][i].grid.totalNumbers > 10 && syncClock.EllapsedMilliseconds() > 1000)
								break;
							if (!LOCAL_RN.APX[d][i].exhausted)
							{
								++runTotal;
								//Mutator done in LOCAL, to free GLOBAL in that time.
								if (Mutate_Exhaustive_6(LOCAL_RN.APX[d][i], candidate, tmpCandidate, rnd))
								{
									++runSuccess;
									if (solved)
										return; //We won. yay!
								}
								LOCAL_RN.APX[d][i].exhausted = true;
								//Search on GLOBAL to exhaust it. This way the mutex time is minimal. 
								//If we tried to do the Mutate_Exhaustive_6 on GLOBAL it will hang all other threads.
								mutexRN_Explorer.lock();
								for (int j = 0; j < GLOBAL_RN.APX[d].size(); ++j)
								{
									if (!GLOBAL_RN.APX[d][j].exhausted && GLOBAL_RN.APX[d][j].remainingNumbers.Equals(LOCAL_RN.APX[d][i].remainingNumbers))
									{
										GLOBAL_RN.APX[d][j].exhausted = true;
										break;
									}
								}
								mutexRN_Explorer.unlock();
							}
						}
					}
				auto elapTime = syncClock.EllapsedMilliseconds();
				//	if (runTotal > 0 && elapTime > 200) //Notify when Keeper took too much to end.
				//cerr << "      Keeper Exhaustive:" << runSuccess << "/" << runTotal << " Time:" << elapTime << "ms" << endl;
			}

			//On Global Best Search we fill the Fitness[k] with random scores of worst RN's
			if (isGlobalBestActive && LOCAL_RN.totalRN > 0)
			{
				bool completo = false;
				for (int d = (int)LOCAL_RN.APX.size() - 1; d >= 0; --d)
				{
					if (LOCAL_RN.APX[d].size() > 0)
					{
						for (int k = 0; k < FINAL_LFA; ++k)
						{
							Fitness[k] = LOCAL_RN.APX[d][rnd.NextInt((int)LOCAL_RN.APX[d].size())].BestScore;
						}
						break;
					}
				}
			}
		}

		//Global Best Algorithm
		if (isGlobalBestActive)
		{

			LAHC_Node* future = nullptr;
			if (rnd.NextInt(1000) < 200)
			{ //pick 
				int indexUse = rnd.NextInt(max(0, (int)historyBEST.size() - 8), max(1, (int)historyBEST.size() - 1));
				future = &historyBEST[indexUse];
			}
			else future = LOCAL_RN.pickWeightedRandom(rnd);


			if (future != nullptr)
			{
#ifdef SHUFFLE_CODE
				if (future->grid.totalNumbers < 10 && rnd.NextInt(1000) < 20 && Mutate_Shuffle(*future, tmpCandidate, rnd))
				{
					if (rnd.NextInt(1000) < 200)
						*future = tmpCandidate;
					lastAccepted = tmpCandidate;
				}
				else lastAccepted = *future;
#else
				lastAccepted = *future;
#endif
			}

			//Random pick Mutate or Recombinate mode, they are really similar.
			int whatToPick = rnd.NextInt(1000);
			if ((lastAccepted.grid.totalNumbers > 6 && whatToPick < 50)
				|| (whatToPick < 200))
				Recombinate(lastAccepted, candidate, rnd);
			else Mutate(lastAccepted, candidate, PROB, rnd);
			//Optional reArm. Some candidates can be fixed if Incomplete sets are played first.
			if (candidate.IncompleteSets.size() > 0 && (candidate.BestScore < lastAccepted.BestScore))
			{
				bool bb = reArmCandidate(candidate, tmpCandidate, rnd);
				if (bb && tmpCandidate.BestScore >= candidate.BestScore)
				{
					swap(tmpCandidate, candidate);
				}
			}
		}
		else
		{  //LAHC mutate from lastAccepted
			Mutate(lastAccepted, candidate, PROB, rnd); //Construct a candidate solution s*  and Calculate its cost function C(s*)
														//uint64_t candHash = candidate.usedNumbers.simpleHash();

														//Optional reArm. Some candidates can be fixed if Incomplete sets are played first.
			if (candidate.IncompleteSets.size() > 0 //We need incomplete sets to reArm
				&& candidate.MutateType != 5 //Not useful on 5
				&& candidate.MutateType != 0 //Not very useful on 0
				&& ((candidate.BestScore < lastAccepted.BestScore && candidate.BestScore < Fitness[v]) //Candidate will be rejected
					|| rnd.NextInt(1000) < 100) //On Accepted candidates just some random chance
				)
			{
				bool bb = reArmCandidate(candidate, tmpCandidate, rnd);
				if (bb && tmpCandidate.BestScore > candidate.BestScore)
				{
					swap(tmpCandidate, candidate);
				}
			}

		}
		MUTATE_STATS[candidate.MutateType] += 1.0;


		//Global Best Search Algorithm. Any new candidate will try to be inserted as a Global Best APX
		if (candidate.grid.totalPoints != 0
			&& (!(candidate.remainingNumbers.Equals(lastAccepted.remainingNumbers) && candidate.BestScore <= lastAccepted.BestScore))
			)
		{
			LOCAL_RN.RN_Insert(candidate, false);
			if (GLOBAL_RN.RN_Insert(candidate, true))
			{
				if (candidate.grid.totalNumbers <= MAX_DEGREE_TOSAVE)
					++newFound; //Perfomance counter
			}
		}

		//Flashback 1 - Rescue dying jobs. If the next LAHC will fail and it's a Timeout I'll do some exhaustive to save it.
		//I'll try to use Mutate_Exhaustive_6 on the last 5 best APX on historyBest. 
		if (((I + 1) % 1000 == 0) && convergenceClock.Timeout())
		{

			if (lastImprovement.EllapsedMilliseconds() > localTimeLimit - 10)
			{
				if ((candidate.BestScore <= bestStrategy.BestScore)
					&& bestStrategy.grid.totalNumbers <= K_EXHAUSTIVE_MIN_NUMBERS) //As it's exhaustive I can't do it with a big amount of RN
				{
					bool savingJob = false;
					//Try last 5 from history
					if (historyBEST.size() > 1)
						for (int i = (int)historyBEST.size() - 1; i >= 0 && i >= (int)historyBEST.size() - 6; --i) {
							Mutate_Exhaustive_6(historyBEST[i], candidate, tmpCandidate, rnd);
							if (candidate.BestScore > bestStrategy.BestScore)
							{
								savingJob = true;
								break;
							}
						}
					if (savingJob)
						cerr << "    Saved ID:" << ID << " from dying! BestScore:" << bestStrategy.BestScore << "<" << candidate.BestScore << " Numbers:" << bestStrategy.grid.totalNumbers << "<" << candidate.grid.totalNumbers << endl;
				}
			}

		}

		bool wasAccepted = false;
		if (
			(candidate.BestScore >= lastAccepted.BestScore || candidate.BestScore >= Fitness[v]) //If C(s*)<=fv or C(s*)<=C(s)
			)
		{ //Accept the candidate
			wasAccepted = true;
			if (candidate.grid.totalPoints == 0) //Stopping condition
			{ //Solve
				if (!solved)
				{
					mutexSOL.lock();
					//Validate Solution
					LAHC_Node checkSol;
					checkSol.clear();
					SolvedArray.clear();
					for (auto&N : candidate.Movelists)
					{
						for (auto&m : candidate.grid.Plan_NMB[N].linear)
						{
							if (checkSol.grid.ValidateMove(m) && checkSol.grid.ApplyMove(m))
							{
								SolvedArray.addMove(m);
							}
						}
					}
					checkSol.MutateType = 9;
					checkSol.CalcStats();
					if (checkSol.grid.totalNumbers == 0 && checkSol.grid.totalPoints == 0 && checkSol.BestScore > 0.0)
					{

						mutexSAVE.lock();
						string filename = "SOLUTION_" + to_string(level) + "_" + passwordLevel + ".txt";
						std::ofstream outfile;
						outfile.open(filename, std::ios_base::app); // append instead of overwrite
						outfile << (isGlobalBestActive ? "RECOMBINATOR " : "LAHC ") << "Stats: Lfa:" << FINAL_LFA << " Time From Reset:" << timeFromLastReset.EllapsedMilliseconds() << " maxTimeImprovement:" << maxTimeImprovement / 1000 << "s Last Improvement:" << lastImprovement.EllapsedMilliseconds() / 1000 << "s" << endl;
						outfile.close();
						mutexSAVE.unlock();
						lastAccepted.SaveToFile(filename, maxTimeImprovement, FINAL_LFA, PROB, " Last Accepted -  Flash:" + to_string(SUCCESS_FLASH) + "/" + to_string(TOTAL_FLASH));
						if (isGlobalBestActive)
							checkSol.SaveToFile(filename, maxTimeImprovement, FINAL_LFA, PROB, "RECOMBINATED!! Flash:" + to_string(SUCCESS_FLASH) + " / " + to_string(TOTAL_FLASH));
						else checkSol.SaveToFile(filename, maxTimeImprovement, FINAL_LFA, PROB, " Flash:" + to_string(SUCCESS_FLASH) + "/" + to_string(TOTAL_FLASH));
						solved = true;
						/*SolvedArray.clear();
						for (auto& c : candidate.Movelists)
						SolvedArray.addPlan(candidate.grid.Plan_NMB[c]);*/
						cerr << "SOLUTION had " << candidate.Movelists.size() << " GROUPS" << endl;
					}
					else {
						cerr << "UPS! Wrong solution! Numbers:" << candidate.grid.totalNumbers << "<->" << checkSol.grid.totalNumbers << " Points" << candidate.grid.totalPoints << "<->" << checkSol.grid.totalPoints << endl;
						SolvedArray.clear();
					}
					mutexSOL.unlock();
				}
				return;
			}

			if (candidate.BestScore > bestStrategy.BestScore) //Keep a best candidate too, not only the last Accepted
			{
				MUTATE_BESTSTATS[candidate.MutateType] += 1.0;
				if (Flashbacks >= 0)
				{
					++SUCCESS_FLASH;
					Flashbacks = -1; //Rearm flashbacks, it was good.
				}
				bestStrategy = candidate;
				++realHistoryLength;
				if (historyBEST.size() >= TAM_HIS) { //To reduce size of memory data...
					historyBEST.pop_front();
				}
				historyBEST.push_back(LAHC_Node(candidate));
				maxTimeImprovement = max(maxTimeImprovement, lastImprovement.EllapsedMilliseconds());
				lastImprovement.Start(0);
			}
			lastAccepted = candidate;
		}

		Fitness[v] = lastAccepted.BestScore; //Insert the current cost into the fitness array fv:=C(s)
		++I;//Increment the iteration number I : = I + 1
		++v;//v := I mod Lfa
		if (v >= CurrLFA)
		{
			v = 0;
		}

		//Some tweak to reuse historyBest[random] as lastAccepted.
		if (!wasAccepted && !isGlobalBestActive && historyBEST.size() > 1 && rnd.NextInt(10000) < 3)
		{
			int indexUse = rnd.NextInt(max(0, (int)historyBEST.size() - 5), max(1, (int)historyBEST.size() - 1));
			lastAccepted = historyBEST[indexUse];
		}

		if ((I % 1000 == 0) && convergenceClock.Timeout())
		{
			convergenceClock.Start(10 * 1000 * 1000);

#ifdef SHUFFLE_CODE			
			if (Mutate_Shuffle(lastAccepted, candidate, rnd)) //Shuffling idea
			{
				swap(lastAccepted, candidate);
			}
#endif		

			//Tweak. Dynamic LFA
			double coverage = 100.0*(double)(COUNT_NUMBERS - bestStrategy.grid.totalNumbers)*INV_COUNT_NUMBERS;
			if (coverage > K_LFA_PERC_MIN)
			{
				int TargetLFA = CurrLFA;
				double media = coverage;
				//High coverage, use Final LFA
				if (media > K_LFA_PERC_MAX)
				{
					TargetLFA = FINAL_LFA;
				}
				else {
					//Linear value.
					TargetLFA = SMALL_LFA + max(0, (int)((FINAL_LFA - SMALL_LFA)*(media - K_LFA_PERC_MIN) / max(K_LFA_PERC_MAX - K_LFA_PERC_MIN, 1.0)));
				}

				if (TargetLFA > CurrLFA) //Expand Fitness, use randomized current LFA values
				{
					for (int k = CurrLFA; k < TargetLFA; ++k)
					{
						Fitness[k] = Fitness[rnd.NextInt(CurrLFA)];
					}
					CurrLFA = TargetLFA;
				}
			}


			//Probabilities recalc
			double RAT[] = {
				(MUTATE_STATS[0] == 0 ? 0.0 : (double)MUTATE_BESTSTATS[0] * 100.0 / MUTATE_STATS[0]),
				(MUTATE_STATS[1] == 0 ? 0.0 : (double)MUTATE_BESTSTATS[1] * 100.0 / MUTATE_STATS[1]),
				(MUTATE_STATS[2] == 0 ? 0.0 : (double)MUTATE_BESTSTATS[2] * 100.0 / MUTATE_STATS[2]),
				(MUTATE_STATS[3] == 0 ? 0.0 : (double)MUTATE_BESTSTATS[3] * 100.0 / MUTATE_STATS[3]),
				(MUTATE_STATS[4] == 0 ? 0.0 : (double)MUTATE_BESTSTATS[4] * 100.0 / MUTATE_STATS[4]),
				(MUTATE_STATS[5] == 0 ? 0.0 : (double)MUTATE_BESTSTATS[5] * 100.0 / MUTATE_STATS[5]),
			};
			double SUMRAT = RAT[0] + RAT[1] + RAT[2] + RAT[3] + RAT[4] + RAT[5];
			if (SUMRAT == 0.0)
			{
				PROB[0] = 10.0;
				PROB[1] = 20.0;
				PROB[2] = 20.0;
				PROB[3] = 10.0;
				PROB[4] = 20.0;
				PROB[5] = 20.0;
			}
			else {
				PROB[0] = max(2.0, (100.0*RAT[0] / SUMRAT));
				PROB[1] = max(2.0, (100.0*RAT[1] / SUMRAT));
				PROB[2] = max(2.0, (100.0*RAT[2] / SUMRAT));
				PROB[3] = max(2.0, (100.0*RAT[3] / SUMRAT));
				PROB[4] = max(2.0, (100.0*RAT[4] / SUMRAT));
				PROB[5] = max(2.0, (100.0*RAT[5] / SUMRAT));
			}


			for (int k = 0; k < 10; ++k)//Decay, first hits are less relevant now
			{
				MUTATE_STATS[k] = MUTATE_STATS[k] * 0.90;
				MUTATE_BESTSTATS[k] = MUTATE_BESTSTATS[k] * 0.90;
			}
			if (!isGlobalBestActive)
				if (lastImprovement.EllapsedMilliseconds() > localTimeLimit)
				{
					lastImprovement.Start(0);

					//Increase limit time
					if (CAN_INCREASE_TIME > 0)
					{
						if (localTimeLimit < 80 * 1000)
						{
							localTimeLimit += CAN_INCREASE_TIME;
						}
					}
					//Increase LFA on resets
					if (CAN_INCREASE_LFA > 0)
					{
						FINAL_LFA += CAN_INCREASE_LFA;
						if (FINAL_LFA > EXPAND_LFA)
							FINAL_LFA = EXPAND_LFA;
					}

					//Flashbacks 2 - Recovering a previous Best State on timeout.
					//I also recreate Fitness[] vector with lower Score values from history.
					if ((historyBEST.size() > 8) && (Flashbacks != 0)
						&& (bestStrategy.grid.totalNumbers <= K_FLASHBACK_NUMBERS || bestStrategy.grid.totalPoints <= K_FLASHBACK_POINTS))
					{
						++TOTAL_FLASH;
						if (Flashbacks == -1)
						{
							Flashbacks = FLASHBACK_COUNT;
						}
						else
							--Flashbacks;
						//	cerr << "   ===>FLASHBACK LAHC ID:" << ID << " F:" << Flashbacks << ".Shrink:" << realHistoryLength << " to ";
						int newIndex;
						int lowLimit;
						if (bestStrategy.BestScore != historyBEST.back().BestScore)
						{
							newIndex = max(1, (int)historyBEST.size() - (int)rnd.NextInt(1, 2));
							lowLimit = max(0, newIndex - 3);
						}
						else {
							newIndex = max(1, (int)historyBEST.size() - (int)rnd.NextInt(1, 2));
							lowLimit = max(0, newIndex - (int)rnd.NextInt(3, 4));
						}
						for (int k = 0; k < FINAL_LFA; ++k)
						{
							int NNe = max(0, rnd.NextInt(lowLimit, newIndex));
							Fitness[k] = historyBEST[NNe].BestScore;
						}
						if (historyBEST.size() > newIndex + 1)
							historyBEST.resize(newIndex);

						if (historyBEST.size() > 0)
						{
							lastAccepted = historyBEST.back();
							//Try to Mutate Shuffle here
#ifdef SHUFFLE_CODE							
							{
								if (Mutate_Shuffle(lastAccepted, candidate, rnd)) //Shuffling idea
								{
									swap(lastAccepted, candidate);
								}
							}
#endif							
						}
					}
					else { //RESET WORKER - Dead End - we need to restart and produce an initial solution s
						cerr << "    **RESET LAHC ID:" << ID << " Lfa:" << FINAL_LFA;
						cerr << " Limit Time Improvement:" << localTimeLimit / 1000 << "s" << endl;
						CurrLFA = SMALL_LFA;
						SUCCESS_FLASH = 0;
						TOTAL_FLASH = 0;
						Flashbacks = -1;
						copiedFromGlobal = false;
						LOCAL_RN.clear();
						timeFromLastReset.Start(0);
						lastAccepted.clear();
						bestStrategy.clear();

						for (int k = 0; k < 10; ++k)
						{
							MUTATE_STATS[k] = 0.0;
							MUTATE_BESTSTATS[k] = 0.0;
							PROB[k] = 0.0;
						}
						maxTimeImprovement = 0;
						PROB[0] = 10.0;
						PROB[1] = 20.0;
						PROB[2] = 20.0;
						PROB[3] = 10.0;
						PROB[4] = 20.0;
						PROB[5] = 20.0;


						Mutate(lastAccepted, candidate, PROB, rnd);
						lastAccepted = candidate;
						bestStrategy = candidate;
						realHistoryLength = 1;
						//historyBEST.clear();
						for (int JG = 0; JG < historyBEST.size(); ++JG) {
							historyBEST[JG].clear();
						}
						historyBEST.resize(1);
						historyBEST[0] = candidate;
						for (int k = 0; k < FINAL_LFA; ++k)
						{
							Fitness[k] = lastAccepted.BestScore; //For all k € {0...Lfa-1} fk:=C(s)
						}
						I = 1;
						previousI = 0;
					}
				}

		}

		//Notification on screen
		if ((I % 1000 == 0) && notifyClock.Timeout())
		{
			notifyClock.Start(30 * 1000 * 1000);
			auto tiempo = lastImprovement.EllapsedMilliseconds();
			auto secs = tiempo / 1000;
			auto ms = tiempo % 1000;

			//Visual mark for high quality APROX, different by RN
			if (bestStrategy.grid.totalNumbers == 1) cerr << "++++";
			else if (bestStrategy.grid.totalNumbers == 2) cerr << "++| ";
			else if (bestStrategy.grid.totalNumbers == 3) cerr << "+|  ";
			else if (bestStrategy.grid.totalNumbers == 4) cerr << "|   ";
			else										  cerr << "    ";

			auto& show = !isGlobalBestActive ? bestStrategy : lastAccepted;
			if (isGlobalBestActive)
				cerr << "RC:" << ID;
			else cerr << "ID:" << ID;
			cerr << " Pt:" << setw(3) << show.grid.totalPoints;
			cerr << " No:" << setw(2) << (show.grid.totalNumbers);
			/*			cerr << " K:" << (int)show.A << "," << fixed << setprecision(2) << show.B << "," << (int)show.C << "," << (int)show.D;
			cerr << "," << setw(2) << (int)show.ExtraC << "," << setw(2) << (int)show.ExtraD;*/
			cerr << " max:" << setw(2) << (maxTimeImprovement / 1000) << "s";
			cerr << " P:[" << setw(2) << (int)PROB[0] << "," << setw(2) << (int)PROB[1] << "," << setw(2) << (int)PROB[2] << "," << setw(2) << (int)PROB[3] << "," << setw(2) << (int)PROB[4] << "," << setw(2) << (int)PROB[5] << "]";
			if (isGlobalBestActive)
			{
				cerr << " LCBEST:";
				if (LOCAL_RN.APX.size()>0)
					for (int lc = 0; lc < LOCAL_RN.APX.size(); ++lc)
						cerr << (LOCAL_RN.APX[lc].size()>0 ? LOCAL_RN.APX[lc][0].grid.totalNumbers : 0) << "=" << LOCAL_RN.APX[lc].size() << ",";
			}

			/*
			cerr << " B:[" << setw(2) << MUTATE_BESTSTATS[0] << "," << setw(2) << MUTATE_BESTSTATS[1] << "," << setw(2) << MUTATE_BESTSTATS[2] << "," << setw(2) << MUTATE_BESTSTATS[3] << "," << setw(2) << MUTATE_BESTSTATS[4] << "," << setw(2) << MUTATE_BESTSTATS[5] << "]";
			cerr << " S:[" << setw(2) << (int)MUTATE_STATS[0] << "," << setw(2) << (int)MUTATE_STATS[1] << "," << setw(2) << (int)MUTATE_STATS[2] << "," << setw(2) << (int)MUTATE_STATS[3] << "," << setw(2) << (int)MUTATE_STATS[4] << "," << setw(2) << (int)MUTATE_STATS[5] << "]";
			*/

			if (show.grid.totalNumbers < 4 || (isGlobalBestActive && show.grid.totalNumbers < 6))
			{ //Print them
				cerr << " Numbers:(";
				for (auto& n : NMB)
				{
					if (show.grid.Val[n.ID] > 0)
					{
						cerr << (int)n.X << "," << (int)n.Y << "=" << (int)show.grid.Val[n.ID] << "|";
					}
				}
				cerr << ")";
			}
			//cerr << " Lfa:" << setw(4) << CurrLFA << "/" << FINAL_LFA;
			if (isGlobalBestActive)
				cerr << " LOCAL:" << (LOCAL_RN.totalRN) << " GLOBAL:" << totalGLOBAL;
			cerr << " I:" << (I - previousI) << endl;
			previousI = I;
		}
	}
}

/*
void initialExhaustive() {
swapDone = 0;
//for (auto& G : GLOBAL_RN.APX)
for (int d = 0; d < 4; ++d)
for (int i = 0; i < GLOBAL_RN.APX[d].size(); ++i)
{
auto G = GLOBAL_RN.APX[d][i];
cerr << "N:" << (int)G.grid.totalNumbers << " P:" << (int)G.grid.totalPoints << " :[";
for (auto&N : NMB)
{
if (G.grid.Val[N.ID] > 0)
{
cerr << (int)N.X << "," << (int)N.Y << "=" << (int)G.grid.Val[N.ID] << " |";
}
}
cerr << "]";
cerr << endl;
LAHC_Node candidate, tmpCandidate;
candidate.clear();
tmpCandidate.clear();
Mutate_Exhaustive_6(G, candidate, tmpCandidate);
GLOBAL_RN.APX[d][i].exhausted = true;
if (IS_SOLVED(G, candidate, "EXHAUSTIVE"))
return;
else if (candidate.Movelists.size() > 0 && candidate.grid.totalNumbers < G.grid.totalNumbers) {
cerr << "IMPROVED!!!!!:" << G.grid.totalNumbers << " <->" << GLOBAL_RN.APX[d][i].grid.totalNumbers << " -> " << candidate.grid.totalNumbers << endl;
++swapDone;
}
}
if (swapDone > 0)
{
cerr << "Recreating APROX file due to " << swapDone << " swaps..." << endl;
string ss = "APROX_" + passwordLevel + ".txt";
remove(ss.c_str());
for (int d = 0; d < 4; ++d)
for (auto&B : GLOBAL_RN.APX[d])
{
B.SaveApprox(ss);
}
swapDone = 0;
}
}
*/


/*
Restore a saved state for Global Insert Search
Create worker Threads.
Prints some global performance stats on cerr
Once solved it waits all threads and print the solution to cout.
*/
void ParallelWork() {
	if (PROGRAM_NAME != "./test")
	{
		//Load APX stored in files. Restarts GLOBAL_Insert to a previous state.
		//It will need to recalculate the Mutate_Exhaustive_6, to be done on Thread ID=0
		LoadRecombinators_APROX("APROX_" + passwordLevel + ".txt", true);
		LoadRecombinators_LAHC("LAHC_" + passwordLevel + ".txt"); //Deprecated
	}
	//initialExhaustive();
	vector<thread> threads(max(1, THREADS));

	if (!solved)
	{
		cerr << "Creating " << THREADS << " Threads" << endl;
		//TODO: Create threads of type Worker_LAHC, with an ID starting at 0
	}

	int C = 0;
	while (!solved)
	{
		if (++C > 750) //Each 75 seconds
		{
			//Import Decentralized APX from other nodes. 
			//Move EXTERN_*.txt file to toImport.txt, process it and remove it afterwards
			if (ALLOW_GLOBAL_BEST > 0 && (PROGRAM_NAME != "./test"))
			{
				string externFile = "EXTERN_" + passwordLevel + ".txt";
				string importFile = "toImport.txt";
				remove(importFile.c_str());
				rename(externFile.c_str(), importFile.c_str());
				LoadRecombinators_APROX("toImport.txt", true); //Reload APROX
				remove(importFile.c_str());
			}

			//Recreate APROX file. Do with a temporary file first, then a system MOVE call to avoid risks of corruption on abnormal terminations.
			if (swapDone > 0 && (PROGRAM_NAME != "./test")
				&& (stopwatch.EllapsedMilliseconds() >= 130 * 1000) //Delay saving to file
				)
			{
				cerr << "Recreating APROX file due to " << swapDone << " swaps...";
				string safeFile = "SAFE_" + passwordLevel + ".txt";
				string finalFile = "APROX_" + passwordLevel + ".txt";
				mutexAPROX.lock();
				remove(safeFile.c_str());
				std::ofstream outfile;
				outfile.open(safeFile, std::ios_base::app); // append instead of overwrite
				for (int d = 0; d < GLOBAL_RN.APX.size(); ++d)
					if (GLOBAL_RN.APX[d].size()>0 && GLOBAL_RN.Degree[d] <= MAX_DEGREE_TOSAVE)
						for (auto&B : GLOBAL_RN.APX[d])
						{
							for (auto& L : B.Movelists)
								if (L >= 0 && L<COUNT_NUMBERS)
									for (auto&m : B.grid.Plan_NMB[L].linear)
									{
										//srcID, destID, dir, sign, srcVal, destScore
										outfile << (int)m.srcIDX << " " << (int)m.destIDX << " " << (int)m.dir << " " << (int)m.sign << " " << (int)m.srcVal << " " << (int)m.destVal << "|";
									}
							outfile << endl;
						}

				outfile.close();
				filesize = GetFileSize(safeFile);
				swapDone = 0;
				mutexAPROX.unlock();
				cerr << "done." << endl;
				//Swap files. This reduces chance of corruption
				rename(safeFile.c_str(), finalFile.c_str());
				cerr << "All OK" << endl;
			}


			//Just some performance 
			cerr << "LEVEL:" << level << " P:" << setw(4) << MAP.totalPoints << " Numbers:" << setw(4) << MAP.totalNumbers;
			cerr << " SIM/Sec:" << SimCount / (stopwatch.EllapsedMilliseconds() / 1000);


			if (newFound == 0)
				performance = stopwatch.EllapsedMilliseconds() / 1000 * THREADS;
			else performance = stopwatch.EllapsedMilliseconds() / 1000 * THREADS / newFound;
			cerr << " Performance:" << performance << "s/found(" << newFound << ")";
			cerr << " JACKPOTS:" << JACKPOT_OK << "/" << JACKPOT_TOTAL << " ";
			cerr << "SUPERCOMBINATOR:" << SUPERCOMBINATOR << " ";
			cerr << "MAX_PERMUTATIONS:" << MAX_PERMUTATIONS;

#ifdef SHUFFLE_CODE			
			//	cerr << " Valid Shuffles:" << correctShuffles << "/" << countShuffles << " : " << (countShuffles == 0 ? 0 : correctShuffles * 100 / countShuffles) << "%";
#endif			

			cerr << " GLOBAL_RN.APX:";
			mutexRN_Explorer.lock();
			for (int lc = 0; lc < GLOBAL_RN.APX.size(); ++lc)
				cerr << (GLOBAL_RN.APX[lc].size()>0 ? GLOBAL_RN.APX[lc][0].grid.totalNumbers : 0) << "=" << GLOBAL_RN.APX[lc].size() << ",";
			cerr << "=" << totalGLOBAL << endl;
			mutexRN_Explorer.unlock();
			C = 0;
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}

	cerr << "  -- SOLUTION -- SimCount:" << SimCount << " Took " << stopwatch.EllapsedMilliseconds() << "ms" << endl;
	if (SolvedArray.linear.size() == 0)
		cerr << "NO SOLUTION!!!" << endl;
	else {
		cerr << passwordLevel << endl;
		for (auto&m : SolvedArray.linear)
		{
			cerr << m << endl;
			cout << m << endl;
		}

	}

	for (auto&t : threads) //who cares? I win!
	{
		t.join();
	}

}

bool isPartOf(const char *a, const char *b) {
	if (std::strstr(b, a) != NULL) {    //Strstr says does b contain a
		return true;
	}
	return false;
}
void stacksize()
{
#ifndef _MSC_VER
	const rlim_t kStackSize = 256 * 1024 * 1024;   // min stack size = 16 MB
	struct rlimit rl;
	int result;

	result = getrlimit(RLIMIT_STACK, &rl);
	if (result == 0)
	{
		if (rl.rlim_cur < kStackSize)
		{
			rl.rlim_cur = kStackSize;
			result = setrlimit(RLIMIT_STACK, &rl);
			if (result != 0)
			{
				fprintf(stderr, "setrlimit returned result = %d\n", result);
			}
		}
	}
#endif
}


void PopulateCacheMoves()
{
	Grid ValidateCache;
	for (int srcID = 0; srcID < COUNT_NUMBERS; ++srcID)
	{
		for (int destID = 0; destID < COUNT_NUMBERS; ++destID)
		{
			if (srcID == destID) continue;
			if (srcID >= COUNT_NUMBERS || destID >= COUNT_NUMBERS) continue;
			int srcX = (int)NMB[srcID].X;
			int srcY = (int)NMB[srcID].Y;
			int targetX = (int)NMB[destID].X;
			int targetY = (int)NMB[destID].Y;

			if (srcX != targetX && srcY != targetY)
			{
				continue;
			}
			int dist = max(abs(srcX - targetX), abs(srcY - targetY));
			for (int d = 0; d < 4; ++d)
			{
				int destX = srcX + DIR_X[d] * dist;
				if (destX != targetX) continue;
				int destY = srcY + DIR_Y[d] * dist;
				if (destY != targetY) continue;
				CACHE_MERGE[srcID][destID] = new Move(srcID, destID, d, 0, dist, MAP.Val[destID]);
				ValidateCache = MAP;
				ValidateCache.Val[srcID] = dist;
				ValidateCache.ApplyMove(*CACHE_MERGE[srcID][destID]);
				if (ValidateCache.Val[srcID] != 0)
				{
					cerr << "ERROR DE CACHE src:" << (int)srcID << " dest:" << (int)destID << " dir;" << (int)d << " dist:" << (int)dist << " " << (int)srcX << "," << (int)srcY << " --> " << (int)destX << "," << (int)destY << endl;
					CACHE_MERGE[srcID][destID] = nullptr;
				}
			}
		}
	}

}

void CreateNeighbours() {
	NEIGHBOURS.resize(INVALID_ID + 2);
	for (int i = 0; i < INVALID_ID + 1; ++i)
		NEIGHBOURS[i].clear();

	//for (int i = 0; i < MAX_NUMBERS; ++i)
	for (auto& N : NMB)
	{
		for (int x = 0; x < MAX_W; ++x)
		{
			if (IDX[x][N.Y] < MAX_NUMBERS && (int)x != (int)N.X) {
				NEIGHBOURS[N.ID].set(IDX[x][N.Y]);
			}
		}
		for (int y = 0; y < MAX_H; ++y)
		{
			if (IDX[N.X][y] < MAX_NUMBERS && (int)y != (int)N.Y) {
				NEIGHBOURS[N.ID].set(IDX[N.X][y]);
			}
		}
	}


}

int main(int argc, char *argv[])
{
	Random rnd;
	{
		PROGRAM_NAME = argv[0];
		system("hostname > hostname.txt");
		ifstream f("hostname.txt");
		if (f.good()) { getline(f, COMPUTER_NAME); f.close(); }

		//if (PROGRAM_NAME == "./test") //For paramtuning
		{
			ifstream fn("randomKey.txt");

			if (fn.good()) {
				fn >> randomKey; fn.close();
				cerr << "Random KEY is:" << randomKey << endl;
				if (randomKey != 0)
					rnd = Random(randomKey);
			}
		}
	}
#ifdef SHUFFLE_CODE
	correctShuffles = 0;
	countShuffles = 0;
#endif	
	JACKPOT_TOTAL = 0; JACKPOT_OK = 0; SUPERCOMBINATOR = 0; MAX_PERMUTATIONS = 0;
	newFound = 0;

	stacksize();
	stopwatch.Start(0);

	if (argc > 1) {
		THREADS = max(1, atoi(argv[1]));
	}

	if (argc > 2) {
		LIMIT_TIME_IMPROVEMENT = atoi(argv[2]) * 1000;
	}
	if (argc > 3) {
		SIZE_LFA = atoi(argv[3]);
	}
	if (argc > 4) {
		K_A = (double)atof(argv[4]);
		if (K_A < 0)
			K_A = rnd.NextInt(100, 900);
	}
	if (argc > 5) {
		K_B = (double)atof(argv[5]);
		if (K_B < 0)
			K_B = rnd.NextInt(0, 100);

	}
	if (argc > 6) {
		K_C = (double)atof(argv[6]);
		if (K_C < 0)
			K_C = rnd.NextInt(100, 900);

	}
	if (argc > 7) {
		K_D = (double)atof(argv[7]);
		if (K_D < 0)
			K_D = rnd.NextInt(100, 900);

	}
	if (argc > 8) {
		CAN_INCREASE_TIME = atoi(argv[8]);
	}
	if (argc > 9) {
		CAN_INCREASE_LFA = atoi(argv[9]);
	}
	if (argc > 10) {
		ALLOW_GLOBAL_BEST = atoi(argv[10]);
	}
	if (argc > 11) {
		CNT_RECOMB = atoi(argv[11]);
		cerr << "CNT_RECOMB:" << CNT_RECOMB << endl;
	}

	cerr << "Program Name:" << argv[0] << endl;
	cerr << "LIMIT TIME IMPROVEMENT IS " << LIMIT_TIME_IMPROVEMENT << endl;
	cerr << "SIZE_LFA IS" << SIZE_LFA << " CAN_INCREASE_TIME:" << CAN_INCREASE_TIME << "ms CAN_INCREASE_LFA:" << CAN_INCREASE_LFA << endl;
	cerr << " K_A .. K_D = " << (int)K_A << "," << (int)K_B << "," << (int)K_C << "," << (int)K_D << endl;

	SimCount = 0;
	solved = false;
	srand(unsigned(time(0)));
	SolvedArray.linear.reserve(50);
	passwordLevel = "first_level";

	{
		ifstream f("level_password.txt");
		if (f.good()) { std::getline(f, passwordLevel); f.close(); }
		ifstream fn("number_level.txt");
		if (fn.good()) { fn >> level; loadValues(level); fn.close(); }
	}

	//#ifdef _MSC_VER
	if (file_exists("level.txt"))
	{
		auto r = freopen("level.txt", "r", stdin);
		r = freopen("solution.txt", "w", stdout);
	}
	//#endif
	//	cout << passwordLevel << endl;
	//	cout << passwordLevel << endl;
	cerr << "***********LEVEL:" << passwordLevel << endl;

	cin >> W >> H; cin.ignore();

	//cerr << W << " " << H << endl;

	if (W > MAX_W || H > MAX_H)
	{
		cerr << "Width/Height overflow: W>MAX_W" << W << " > " << MAX_W << endl;
		cerr << "Width/Height overflow: H>MAX_H" << H << " > " << MAX_H << endl;
		assert(W <= MAX_W);
		assert(H <= MAX_H);
	}
	//Load Initial State on MAP
	MAP.totalNumbers = 0;
	MAP.totalPoints = 0;
	MAP.squaredPoints = 0;
	for (int i = 0; i < MAX_NUMBERS; ++i)
		MAP.Val[i] = 0;
	for (int y = 0; y < MAX_H; y++) {
		for (int x = 0; x < MAX_W; x++) {
			IDX[x][y] = INVALID_ID;
		}
	}
	for (int y = 0; y < H; y++) {
		for (int x = 0; x < W; x++) {
			int val = 0;
			cin >> val; cin.ignore();
			if (val != 0)
			{
				IDX[x][y] = (uint16_t)MAP.totalNumbers;
				MAP.Val[IDX[x][y]] = val;
				NMB.push_back(Coord{ (uint8_t)x,(uint8_t)y, (uint16_t)IDX[x][y] });
				++MAP.totalNumbers;
				MAP.totalPoints += val;
				MAP.squaredPoints += val * val;
				ASSERT(NMB[IDX[x][y]].X == (uint8_t)x);
				ASSERT(NMB[IDX[x][y]].Y == (uint8_t)y);
			}
		}
	}
	COUNT_NUMBERS = MAP.totalNumbers;
	assert((int)NMB.size() == (int)COUNT_NUMBERS);

	NumberSet::Mask = NumberSet::CreateMask(COUNT_NUMBERS);
	if (COUNT_NUMBERS >= MAX_NUMBERS)
	{
		cerr << "INCREASE MAX_NUMBERS!!!!" << COUNT_NUMBERS << " > " << MAX_NUMBERS << endl;
		abort();
	}
	GROUPS = max(2, abs(COUNT_NUMBERS - spawns));
	GROUPS = (GROUPS % 4 != 0 ? 1 : 0) + GROUPS / 4; //It seems there are less

													 //Calculate divisions to use them as multiplications in score calculation.
	INV_COUNT_NUMBERS = 1.0 / (double)COUNT_NUMBERS;
	INV_GROUPS = 1.0 / (double)GROUPS;
	INV_POINTS = 1.0 / (double)MAP.totalPoints;
	INV_SQUAREPOINTS = 1.0 / (double)MAP.squaredPoints;
	INV_RESX = 0.5 / (10.0*(double)W);
	INV_RESY = 0.5 / (10.0*(double)H);

	cerr << "LEVEL:" << level << " NUMBERS:" << COUNT_NUMBERS << " SPAWNS:" << spawns << " GROUPS:" << GROUPS << " W,H:" << W << "," << H << endl;
	PopulateCacheMoves();
	CreateNeighbours();
	ParallelWork();

	return 0;
}
