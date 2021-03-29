// LEGACY CODE: //Code is messy, but I was trying to solve the puzzle, not doing production-level code. This isn't representative of production-level code.
//This was the real code I used for reaching level 1000. There are a lot of segfaults, mainly due to unfixed race conditions on multithread.
//As the whole solver+python script was resilient to errors, I don't care a lot as long as I don't waste a lot of time.
//There are small parts omitted on the code, but simple enough to recreate it in 10 minutes.They are marked with a TODO:
//The idea is not only Copy&Paste it, but learn from it. As I reached the max level I no longer care about improvements.


/*#pragma GCC optimize("O3","unroll-loops","omit-frame-pointer","inline")
#pragma GCC option("arch=native","tune=native","no-zeroupper")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2")*/
#include <immintrin.h> //SSE Extensions
#include <bits/stdc++.h> //All main STD libraries
#include <atomic>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <thread> 
#include <mutex>
#include <unordered_set>
#pragma comment(linker, "/STACK:33554432")

#define ALTERNATE_RANDOM
//#include <execution>

//Para usar caché de movimientos o no
//#define USE_CACHE_MOVES

//Menos atomics #define LOG_REARM
#define SHUFFLE_CODE

#define DISABLE_OPTIONALS
const uint16_t LIM_P = 8;
const uint16_t LIM_SQ = 8;
const int K_ENABLE_4_LIMIT = 2000; //EL _4 calcula muy poco //30;
int K_B_NUM = 5;
int CNT_RECOMB = 0;
const int K_EXHAUSTIVE_MIN_NUMBERS = 7;
const int SHUFFLE_TRIES = 2;
const int K_ALLOWSUM = 10;
const int K_MAX_POINTS_REMOVE = 4;
const int K_FLASHBACK_NUMBERS = 4;
const int K_FLASHBACK_POINTS = 17;
const int FLASHBACK_COUNT = 2; //ERA 3
const double K_LFA_PERC_MIN = 90.0;//85.0;
const double K_LFA_PERC_MAX = 96.0;

/***** STRINGS_4 PARAMETERS *****/
const int K_LOW_POINTS_REMOVE = 2;//2;
const int K_CHANCE_POINTS_REMOVE = 5;
const int K5_MERGE0 = 750;
const int K5_MERGE1 = 850;

const int K_CHANCE_SIGN = 5;
const int K_CHANCE_IGNORE = 15;
const int K_CHANCE_MERGE = 920;
const int K_CHANCE_EXPLODE = 50;
const int K_PREFER_ENDINGS = 600;
const int K_TAIL_PERCENT = 75;
const int K_SEL_VAL_RANDOM = 70;
const int K_SEL_VAL_TAIL = 300;
const int K_SEL_VAL_OPT = 300;
const int K_SEL_VAL_ORPHAN = 300;
const int K_SEL_VAL_CROSS = 100;
const int K_SEL_VAL_TOTAL = K_SEL_VAL_RANDOM + K_SEL_VAL_TAIL + K_SEL_VAL_OPT + K_SEL_VAL_ORPHAN + K_SEL_VAL_CROSS;
/***** STRINGS_4 PARAMETERS *****/

#ifndef DISABLE_OPTIONALS
#include <unordered_map>
#endif
/*************************************************************/
#define STACK_SIZE 64000
bool USE_SA = false;

double K_A = 700.0;
double K_B = 10.0;
double K_C = 100.0;
double K_D = 400.0;
int THREADS = 4;




//#define DEBG_MODE
#ifdef DEBG_MODE
#define ASSERT(x) assert(x)
#else 
#define ASSERT(x) 
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
typedef int32_t I;
using namespace std;
long long LIMIT_TIME_IMPROVEMENT = 110 * 1000;
#ifndef _MSC_VER
#include <sys/resource.h>
#include <unistd.h>
#endif

int SIZE_LFA = -1;
int CAN_INCREASE_LFA = 0;
int CAN_INCREASE_TIME = 0;


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


bool blockSearch = false;


string passwordLevel = "first_level";
string PROGRAM_NAME = "NONAMED";
string COMPUTER_NAME = "NONAMED";


#include <random>
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
		cerr << "**************************************************************" << endl;
		cerr << "**************************************************************" << endl;
		cerr << "*******WARNIIIIIIINGGG USING PRESEED NUMBERS!!!!!!!!!!!*******" << endl;
		cerr << "**************************************************************" << endl;
		cerr << "**************************************************************" << endl;
		cerr << "Seed1:" << bitset<64>(_seed1) << endl;
		cerr << "Seed2:" << bitset<64>(_seed2) << endl;
	}
	Random() {
		std::random_device rd;
		std::seed_seq seedseq1{ rd(), rd(), rd() , rd() }; // is there an optimal number of rd() to use?
		std::mt19937_64 gen(seedseq1);
		std::uniform_int_distribution<uint64_t> dis;
		_seed1 = dis(gen);
		_seed1 ^= Now().time_since_epoch().count();
		_seed2 = dis(gen);
		_seed2 ^= (Now().time_since_epoch().count() * 7);
		cerr << "Seed1:" << bitset<64>(_seed1) << endl;
		cerr << "Seed2:" << bitset<64>(_seed2) << endl;
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
		cerr << "**************************************************************" << endl;
		cerr << "**************************************************************" << endl;
		cerr << "*******WARNIIIIIIINGGG USING PRESEED NUMBERS!!!!!!!!!!!*******" << endl;
		cerr << "**************************************************************" << endl;
		cerr << "**************************************************************" << endl;

		cerr << "Seed:" << bitset<64>(seed) << endl;
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
Random rnd;

const int MAX_W = 56;
const int MAX_H = 32;
const int MAX_NUMBERS = 1055;

int COUNT_NUMBERS = MAX_NUMBERS;
struct Coord { int X; int Y; uint16_t ID; };
vector<Coord> NMB;
uint16_t IDX[MAX_W][MAX_H];
const uint16_t INVALID_ID = (uint16_t)(MAX_NUMBERS + 2);

int W, H;
int level = -1;
int spawns = -1;
int GROUPS = -1;
double INV_COUNT_NUMBERS = 1.0;
double INV_GROUPS = 1.0;
double INV_POINTS = 1.0;
double INV_TOTALNUMBERS = 1.0;
double INV_SQUAREPOINTS = 1.0;
double INV_MINIPUNTO = 1.0;
double INV_RESX = 1.0;
double INV_RESY = 1.0;

atomic<uint64_t> newFound;
atomic<uint64_t> JACKPOT_OK;
atomic<uint64_t> SUPERCOMBINATOR;
atomic<uint64_t> MAX_PERMUTATIONS;

atomic<uint64_t> JACKPOT_TOTAL;
int64_t performance;


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

const int U = 0;
const int D = 1;
const int L = 2;
const int R = 3;
const int DIR_X[] = { 0,0,-1,1 };
const int DIR_Y[] = { -1,1,0,0 };
const string strDIRS[] = { "U","D","L","R" };


inline bool file_exists(const std::string& name) {
	ifstream f(name.c_str());
	return f.good();
}


void Add(vector<size_t>&A, const vector<size_t>& B)
{
	for (auto&b : B)
		A.push_back(b);
}

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

	uint64_t simpleHash() {

		uint64_t H = 0;
		for (int i = 0; i < INT64_W; ++i)
		{
			H ^= W[i];
		}
		return H;
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
		/*for (int i = 0; i < AVX_W; ++i)
		X[i] = _mm256_set1_epi8(0xFF);*/
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
long filesize = 0;

ostream& operator<<(ostream& os, const NumberSet& m)
{
	for (int i = 0; i < INT64_W; ++i)
	{
		os << bitset<64>(m.W[i]);
	}
	return os;
};

template <class T>
void do_shuffle(T& smallDLX)
{
	for (int i = (int)smallDLX.size() - 1; i > 0; --i) {
		int r = rnd.NextInt(i + 1);
		if (i != r) {
			swap(smallDLX[i], smallDLX[r]);
		}
	}
}


//#define low_bits
struct Move {
	/*
	uint32_t srcIDX;
	uint32_t destIDX;
	uint32_t dir;
	uint32_t sign;
	uint32_t srcVal;
	uint32_t destVal;
	*/
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



	size_t CreateHash() const {
		size_t hash = (size_t)srcIDX; hash <<= 14;
		hash += (size_t)destIDX; hash <<= 14;
		hash += (size_t)dir; hash <<= 1;
		hash += (size_t)sign; hash <<= 3;
		hash += (size_t)srcVal; hash <<= 14;
		hash += (size_t)destVal;
		return hash;
	}
	inline bool operator==(const Move& rhs) {
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
	}

};
inline void Add(vector<Move>&A, const vector<Move>& B)
{
	for (auto&b : B)
		A.push_back(b);
	//A.insert(A.end(), B.begin(),B.end());
}

vector<NumberSet> NEIGHBOURS;

Move* CACHE_MERGE[INVALID_ID + 1][INVALID_ID + 1];
#ifdef USE_CACHE_MOVES
vector<Move>* CACHE_MOVES[INVALID_ID + 1][MAX_W + 1];
#endif

ostream& operator<<(ostream& os, const Move& m)
{
	os << (int)NMB[m.srcIDX].X << ' ' << (int)NMB[m.srcIDX].Y << ' ' << strDIRS[m.dir] << " " << (m.sign <= 0 ? "-" : "+");
	return os;
};
struct Grid;
struct ALIGN Strategy;
#ifndef DISABLE_OPTIONALS
inline void Add(unordered_map<int, vector<Strategy>>&A, const unordered_map<int, vector<Strategy>>&B);
#endif
#ifdef SHUFFLE_CODE
struct ALIGN ShuffleStrat {
	list<Move> linear;
	inline void clear() {
		linear.resize(0);
	}
};
#endif
struct ALIGN Strategy {
	vector<Move> linear;
	uint8_t Endpoint = 0;

#ifndef  DISABLE_OPTIONALS
	vector<Move> criticalPath;
	unordered_map<int, vector<Strategy>> optionals;
	int optionalCount = 0;
	int directOptCount = 0;
	uint16_t LinkValue = 40000;
	uint32_t markTrack = 0;
#endif // ! DISABLE_OPTIONALS
	inline void addMove(const Move& m) {
#ifndef DISABLE_OPTIONALS
		criticalPath.push_back(m);
#endif
		linear.push_back(m);
	}
	inline void resize(const int& N) {
#ifndef DISABLE_OPTIONALS
		criticalPath.resize(N);
#endif
		linear.resize(N);

	}
	inline void addPlan(const Strategy& S) {
		Add(linear, S.linear);
#ifndef  DISABLE_OPTIONALS
		Add(criticalPath, S.criticalPath);
		Add(optionals, S.optionals);
		optionalCount += S.optionalCount;
#endif
	}
#ifndef  DISABLE_OPTIONALS
	inline void addAsOptional(const Strategy& S, uint32_t linkidx, uint32_t linkvalue) {
		optionalCount = 1 + S.optionalCount;
		optionals[linkidx].push_back(S);
		Add(linear, S.linear);
		auto& o = optionals[linkidx][optionals[linkidx].size() - 1];
		o.LinkValue = linkvalue;
		++directOptCount;
	}
#endif
	inline void setEndPoint(int NumberIDX) { Endpoint = 1; }
	inline void clear() {
		linear.resize(0), Endpoint = 0;
#ifndef  DISABLE_OPTIONALS
		markTrack = 0;
		LinkValue = 40000;
		criticalPath.clear();
		optionals.clear(); directOptCount = 0; optionalCount = 0;
#endif
	}
	~Strategy() {
		/*optionals.clear();
		moves.clear();
		linear.clear();*/
	}


#ifndef  DISABLE_OPTIONALS
	string PrintOptionals(int linkID, int indent) {
		string s;
		if (linkID == -1)
			s = string(indent, ' ') + (Endpoint > 0 ? "C:" : "I:") + to_string(criticalPath.size()) + "/" + to_string(linear.size()) + "|";
		else s = string(indent, ' ') + "Opt:ID:" + to_string(linkID) + " V:" + to_string(LinkValue) + "->" + to_string(criticalPath.size()) + "/" + to_string(linear.size()) + "|";
		indent += 3;
		for (auto& m : criticalPath)
		{
			s += to_string(m.srcIDX) + "," + to_string(m.destIDX) + " " + to_string(m.srcVal) + (m.sign == 0 ? "-" : "+") + to_string(m.destVal) + "|";
		}
		s += "\r\n";
		for (auto& o : optionals)
		{
			for (auto& SEC : o.second)
			{
				s += SEC.PrintOptionals(o.first, indent);
			}
		}
		return s;
	}
#endif
	void enumerateOptionals(Grid&g, vector<pair<Strategy*, Move>>& _optionals);

	void ApplyMoves(Grid& g, const int& prob1000, const uint32_t& mark, bool dontPlayEndpoint, vector<int>& UnusedNumbers);
};

long GetFileSize(std::string filename)
{

	ifstream in_file(filename, ios::binary);
	in_file.seekg(0, ios::end);
	return (long)in_file.tellg();
}
static void print_Sol(vector<vector<size_t>>& P, vector<size_t>& idx) {
	cerr << "Solution:";
	for (auto& n : idx)
	{
		cerr << n << ",";
	}
	cerr << endl;

	vector<int> count;

	for (auto& n : idx) {
		for (auto&v : P[n])
		{
			if (count.size() < v + 1)
				count.resize(v + 1);
			++count[v];
		}
	}
	cerr << "  Item count:";
	for (auto& N : count) {
		cerr << N << ",";
	}
	cerr << endl;
}

#ifndef DISABLE_OPTIONALS
inline void Add(unordered_map<int, vector<Strategy>>&A, const unordered_map<int, vector<Strategy>>&B) {
	for (auto&b : B)
	{
		auto& Q = A[b.first];
		auto& G = b.second;
		Q.insert(
			Q.end(),
			std::make_move_iterator(G.begin()),
			std::make_move_iterator(G.end())
		);

		/*for (auto& s : b.second)
		{
		A[b.first].push_back(s);
		}*/
	}
}
#endif

atomic<uint64_t> SimCount;

#ifdef LOG_REARM
atomic<uint64_t> validRearm;
atomic<uint64_t> ATTEMPTRearm;
atomic<uint64_t> prevRearm;
#endif

#ifdef SHUFFLE_CODE
atomic<uint64_t> correctShuffles;
atomic<uint64_t> countShuffles;
#endif
struct Grid {
	uint16_t Val[MAX_NUMBERS];
	Strategy Plan_NMB[MAX_NUMBERS];
	int totalPoints = 0;
	int squaredPoints = 0;
	int totalNumbers = 0;
	void clear() {
		for (int i = 0; i < COUNT_NUMBERS; ++i)
		{
			Plan_NMB[i].clear();
		}
		totalPoints = MAX_NUMBERS - 1;
		squaredPoints = 0;
		totalNumbers = MAX_NUMBERS - 1;
	}
	bool ValidateMove(const Move& m) {
		if (m.srcIDX == m.destIDX || m.srcIDX <0 || m.destIDX <0 || m.srcIDX >= COUNT_NUMBERS || m.destIDX >= COUNT_NUMBERS)
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
	inline bool ApplyMove(const Move& m)
	{
		if (m.srcIDX == m.destIDX || m.srcIDX <0 || m.destIDX <0 || m.srcIDX >= COUNT_NUMBERS || m.destIDX >= COUNT_NUMBERS)
		{
			return false;
		}
#ifdef _MSC_VER
		if (NMB[m.srcIDX].X != NMB[m.destIDX].X && NMB[m.srcIDX].Y != NMB[m.destIDX].Y)
		{
			return false;
		}
		{
			int dist = max(abs((int)NMB[m.srcIDX].X - (int)NMB[m.destIDX].X), abs((int)NMB[m.srcIDX].Y - (int)NMB[m.destIDX].Y));
			if (dist != m.srcVal || dist == 0 || dist > W)
				return false;
		}
#endif
		if (Val[m.srcIDX] == 0 || Val[m.destIDX] == 0)
			return false;

		auto srcScore = Val[m.srcIDX];
		if (m.srcVal != srcScore)
			return false; //incoherent
		auto destScore = Val[m.destIDX];
		if (destScore != m.destVal)
			return false; //Incoherent

		--totalNumbers;
		Val[m.srcIDX] = 0;
		auto Dest = (m.sign == 0 ? abs(destScore - srcScore) : destScore + srcScore);
		Val[m.destIDX] = Dest;

		auto& SRC = Plan_NMB[m.srcIDX];
		auto& DST = Plan_NMB[m.destIDX];

		if (Dest == 0)
		{
			--totalNumbers;
			DST.setEndPoint(m.destIDX);
		}
#ifndef  DISABLE_OPTIONALS
		if (2 * destScore == srcScore && Dest != 0 && m.sign == 0 && DST.criticalPath.size() > 0)// (2 * destScore == srcScore && m.sign == 0)
		{
			SRC.addMove(m);
			DST.addAsOptional(SRC, m.destIDX, destScore);
			SRC.clear();
		}
		else
#endif
		{
			if (SRC.linear.size() != 0) {
				DST.addPlan(SRC);
				SRC.clear();
			}
			DST.addMove(m);
		}
		return true;
	}

#ifdef SHUFFLE_CODE
	inline void ShuffleMove(const Move& m, ShuffleStrat Shuffle_NMB[MAX_NUMBERS])
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
				int cuentaAtras = 0;
				int changeValue = -1;
				for (auto it = DST.linear.crbegin(); it != DST.linear.crend(); ++it)
				{
					if (it->srcIDX != m.srcIDX && it->destIDX != m.srcIDX && it->srcIDX != m.destIDX && it->destIDX != m.destIDX)
					{
						++cuentaAtras;
						continue;
					}
					if (it->srcIDX == m.srcIDX || it->destIDX == m.srcIDX || it->srcIDX == m.destIDX)
					{
						//That will affect result
						break;
					}
					if (it->destIDX == m.destIDX && (m.srcVal >= m.destVal || it->srcVal >= it->destVal))
					{
						break;
					}
					if (m.sign == 1)
					{

						++cuentaAtras;
						changeValue = cuentaAtras;
						continue;
					}
					if (it->sign == 1)
					{
						if (it->destVal <= m.srcVal)
							break;
						else {
							++cuentaAtras;
							changeValue = cuentaAtras;
							continue;
						}
					}
					//both signs 0
					if (it->destVal <= m.srcVal || it->srcVal <= m.destVal)
						break;

					++cuentaAtras;
				}
				if (cuentaAtras == 0)
				{
					DST.linear.push_back(m);
				}
				else {
					int randomStop = rnd.NextInt(cuentaAtras);
					int nIni = changeValue + 1;
					int nFin = cuentaAtras - 1;
					if (changeValue >= 0 && nIni <= nFin && rnd.NextInt(1000) < 600) {
						randomStop = rnd.NextInt(nIni, nFin);
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
	inline bool doMerge(const int& srcID, const int& destID, unordered_set<size_t>& Forbidden) {
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
	inline bool doMerge(const int& srcID, const int& destID) {
		if (srcID == destID || srcID >= COUNT_NUMBERS || destID >= COUNT_NUMBERS || srcID < 0 || destID < 0)
			return false;
		Move* lm = CACHE_MERGE[srcID][destID];
		if (lm == nullptr)return false;
		if (Val[srcID] != lm->srcVal || Val[srcID] == 0) return false;
		auto targetVAL = Val[destID];
		if ((targetVAL == 0) || (2 * targetVAL != Val[srcID])) return false;
		Move m = *lm; //Move(srcID, destID, d, 0, srcVal, targetVAL);
		m.destVal = targetVAL;
		ApplyMove(m);
		return true;
	}

	inline bool doForceJoin(const int& srcID, const int& destID, uint32_t randomChance) {
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

	inline void ExplodeMoves(const int& srcID, vector<Move>& moves)
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
			if (rnd.NextInt(100) < K_ALLOWSUM)
			{
				m = Move(srcID, destID, d, 1, srcVal, destVal);
				moves.push_back(m);
			}
		}
		return;
	}
#ifdef USE_CACHE_MOVES
	inline bool doRandomMove() {
		int srcID = rnd.NextInt(COUNT_NUMBERS);
		uint32_t ren = rnd.xrandom();
		for (int NN = 0; NN < COUNT_NUMBERS; ++NN)
		{
			++srcID;
			if (srcID >= COUNT_NUMBERS)
				srcID = 0;
			auto srcVal = Val[srcID];
			if (srcVal >= W || srcVal == 0)
				continue;
			auto mvl = CACHE_MOVES[srcID][srcVal];
			if (mvl == nullptr) continue;
			vector<Move>& moves = *mvl;
			int d = ren & 3;
			ren >>= 2;
			if (ren == 0)
				ren = rnd.xrandom();
			int L = (int)moves.size();
			for (int eje = 0; eje < L; ++eje)
			{
				++d;
				if (d >= L)
					d = 0;
				if (Val[moves[d].destIDX] != 0)
				{
					Move m = moves[d];
					m.destVal = Val[m.destIDX];
					m.sign = rnd.NextInt(100) < K_ALLOWSUM ? 1 : 0;
					ApplyMove(m);
					return true;
				}
			}
		}
		return false;
	}
#else 
	inline bool doRandomMove() {
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
			if (ren == 0)
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
#endif
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

ALIGN Grid MAP;
#ifndef DISABLE_OPTIONALS
void Strategy::enumerateOptionals(Grid&g, vector<pair<Strategy*, Move>>& _optionals) {
	vector<Move> MMM;
	for (auto&o : optionals)
	{
		for (auto&s : o.second)
			if (s.linear.size() > 1) //No single moves...
			{
				//Get the srcID from last move
				auto* m = &s.linear[s.linear.size() - 1];


				MMM.clear();
				auto BCKval = g.Val[m->srcIDX];
				g.Val[m->srcIDX] = m->srcVal;
				g.ExplodeMoves(m->srcIDX, MMM);
				g.Val[m->srcIDX] = BCKval;
				for (auto& posM : MMM)
					if (m->destIDX != posM.destIDX)
					{
						_optionals.push_back(pair<Strategy*, Move>{&s, posM});
					}
				//The LinkValue
				s.enumerateOptionals(g, _optionals);
			}

	}
}
#endif
#ifndef DISABLE_OPTIONALS
void Strategy::ApplyMoves(Grid& g, const int& prob1000, const uint32_t& mark, bool dontPlayEndpoint, vector<int>& UnusedNumbers) {
	if (markTrack == mark) //Already done
		return;
	markTrack = mark;
	for (int iii = 0; iii < criticalPath.size(); ++iii)
	{
		auto&m = criticalPath[iii];
		//Try to add optionals as srcID
#ifndef DISABLE_OPTIONALS
		int mergeIDX = m.srcIDX;
		if (optionals.count(mergeIDX))
		{
			for (auto&o : optionals[mergeIDX])
				if (o.markTrack != mark && o.LinkValue == g.Val[mergeIDX])
				{

					int chanceDo = rnd.NextInt(1000);
					if (chanceDo < prob1000)
					{
						o.ApplyMoves(g, prob1000, mark, false, UnusedNumbers);
					}
					else markTrack = mark; //marks it as used
				}
		}
#endif

		if (UnusedNumbers.size() > 0 && rnd.NextInt(1000) < 850)
			for (auto& unused : UnusedNumbers)
			{
				g.doMerge(unused, m.srcIDX);
			}
		if ((iii == criticalPath.size() - 1) && dontPlayEndpoint)
			continue;
		g.ApplyMove(m);
		if (UnusedNumbers.size() > 0 && rnd.NextInt(1000) < 850)
			for (auto& unused : UnusedNumbers)
			{
				g.doMerge(unused, m.destIDX);
			}
#ifndef DISABLE_OPTIONALS
		mergeIDX = m.destIDX;
		if (optionals.count(mergeIDX))
		{
			for (auto&o : optionals[mergeIDX])
				if (o.markTrack != mark && o.LinkValue == g.Val[mergeIDX])
				{

					int chanceDo = rnd.NextInt(1000);
					if (chanceDo < prob1000)
					{
						o.ApplyMoves(g, prob1000, mark, false, UnusedNumbers);
					}
					else markTrack = mark; //marks it as used
				}
		}
#endif
	}
}
#else

void Strategy::ApplyMoves(Grid& g, const int& prob1000, const uint32_t& mark, bool dontPlayEndpoint, vector<int>& UnusedNumbers) {
	for (int iii = 0; iii < linear.size(); ++iii)
	{
		auto&m = linear[iii];
		//Try to add optionals as srcID
		if (UnusedNumbers.size() > 0 && rnd.NextInt(1000) < 850)
			for (auto& unused : UnusedNumbers)
			{
				g.doMerge(unused, m.srcIDX);
			}
		if ((iii == linear.size() - 1) && dontPlayEndpoint)
			continue;
		g.ApplyMove(m);
		if (UnusedNumbers.size() > 0 && rnd.NextInt(1000) < 850)
			for (auto& unused : UnusedNumbers)
			{
				g.doMerge(unused, m.destIDX);
			}
	}
}
#endif
std::mutex mutexSOL;
atomic<bool> solved;
Strategy SolvedArray;
std::mutex mutexSAVE;
std::mutex mutexAPROX;

struct ALIGN LAHC_Node {
	NumberSet remainingNumbers; //Numbers != 0
	NumberSet usedInMovesNumbers; //All srcIDX and destIDX, may be 0 or not
	NumberSet untouchedNumbers; //Numbers not used on any move
	NumberSet completedNumbers; //Numbers from complete Lists

	double BestScore = 0.0;
	double A, B, C, D, ExtraC, ExtraD;
	//Strategy allMoves;
	Grid grid;
	vector<int> Movelists;
	vector<int> ZerosumSets; //I have a zero
	vector<int> IncompleteSets;



	//completeNumbers.clear();
	//int Version;
	int MutateType = 9;
	bool exhausted = false;
	void CalcStats() {
		//MutateType = _MutateType;
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
				if (plan.Endpoint > 0)
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

		untouchedNumbers = usedInMovesNumbers;
		untouchedNumbers.negate();
		NumberSet check = untouchedNumbers;
		check._or(usedInMovesNumbers);


		if (grid.totalNumbers > K_B_NUM)
		{
			extraPoints = grid.totalPoints;
			extraSQ = grid.squaredPoints;
		}


//>>>> TODO Score Calculation
		A = K_A * (double)(COUNT_NUMBERS - grid.TODOTODOTODO)*INV_COUNT_NUMBERS;
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
			/*	B = 1.0 - ((double)abs(GROUPS - (int)Movelists.size()))*INV_GROUPS;
			B -= (0.5*(double)IncompleteSets.size())*INV_GROUPS;
			B = K_B * B;*/
		}
		else B = 0.0;
		/*if (grid.totalNumbers <= K_B_NUM)
		B += K_B * (double)completeNumbers.count() *INV_COUNT_NUMBERS;*/
		C = K_C * (double)(MAP.TODOTODOTODO - grid.TODOTODOTODO)* INV_POINTS;
		ExtraC = K_C * (double)(MAP.totalPoints - extraPoints)* INV_POINTS - C;
		//		C += ExtraC;
		D = K_D * (double)(MAP.TODOTODOTODO - grid.TODOTODOTODO)*INV_SQUAREPOINTS;
		ExtraD = K_D * (double)(MAP.squaredPoints - extraSQ)*INV_SQUAREPOINTS - D;
//<<<<< TODO Score Calculation
		BestScore = A + B + C + D + max(0.0, ExtraC) + max(0.0, ExtraD);

		uint64_t S = 0;
		for (auto& l : Movelists)
		{
			S += grid.Plan_NMB[l].linear.size();
		}
		SimCount += S;
	}
	LAHC_Node() {}
	void clear() {
		exhausted = false;
		grid.clear();
		//copy(begin(MAP.Val), std::end(MAP.Val), begin(grid.Val));
		copy(MAP.Val, MAP.Val + MAX_NUMBERS, grid.Val);

		//allMoves.clear();
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

	void SaveToFile(string s, long long maxTimeImprovement = -1, int currLFA = -1, double PROB[] = nullptr, string ExtraText = "")
	{
		std::sort(std::begin(ZerosumSets), std::end(ZerosumSets),
			[=](int a, int b) {return grid.Plan_NMB[a].linear.size() > grid.Plan_NMB[b].linear.size(); });
		std::sort(std::begin(IncompleteSets), std::end(IncompleteSets),
			[=](int a, int b) {return grid.Plan_NMB[a].linear.size() > grid.Plan_NMB[b].linear.size(); });

		double points = (double)(MAP.totalPoints - grid.totalPoints)*100.0 / (double)MAP.totalPoints;
		double coverage = 100.0*(double)(COUNT_NUMBERS - grid.totalNumbers)*INV_COUNT_NUMBERS;
		mutexSAVE.lock();
		std::ofstream outfile;
		outfile.open(s, std::ios_base::app); // append instead of overwrite
		outfile << "********** Binary:" << PROGRAM_NAME << " HOST:" << COMPUTER_NAME << " C/I" << ZerosumSets.size() << "/" << IncompleteSets.size() << endl;
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
		outfile << ExtraText << endl << " GRID IS:";
		outfile << endl;
		if (grid.totalPoints != 0)
			outfile << grid;
#ifndef DISABLE_OPTIONALS		
		for (auto& b : ZerosumSets)
		{
			string ss = grid.Plan_NMB[b].PrintOptionals(-1, 0);
			outfile << ss;
		}
		for (auto& b : IncompleteSets)
		{
			string ss = grid.Plan_NMB[b].PrintOptionals(-1, 0);
			outfile << ss;
		}
#endif
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

		if (grid.totalPoints == 0)
		{
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

bool IS_SOLVED(const LAHC_Node& candidate, string TEXT_)
{
	if (!solved && candidate.Movelists.size()>0 && candidate.grid.totalPoints == 0 && candidate.grid.totalNumbers == 0)
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
			outfile << TEXT_ << endl;
			outfile.close();
			mutexSAVE.unlock();
			checkSol.SaveToFile(filename, 0, 0, nullptr, TEXT_);
			solved = true;
			cerr << "SOLUTION had " << candidate.Movelists.size() << " GROUPS" << endl;
			cerr << passwordLevel << endl;
			for (auto&m : SolvedArray.linear)
			{
				cerr << m << endl;
				cout << m << endl;
			}
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
int swapDone = 0;
const int MAX_INSERT_NUMBERS = 6;
std::mutex mutexGLOBAL_Best;
vector<LAHC_Node> GLOBAL_Best[MAX_INSERT_NUMBERS];
atomic<int> totalGLOBAL;
bool GLOBAL_Best_Insert(LAHC_Node& newBest, bool saveToFile)
{
	//if (IS_SOLVED(newBest, "INSERTING"))
	//return true;
	if (newBest.grid.totalNumbers > MAX_INSERT_NUMBERS)
		return false;
	mutexGLOBAL_Best.lock();
	const int LIMITE_GLOBAL = 250;
	if (MAX_INSERT_NUMBERS == 6 && GLOBAL_Best[0].size() + GLOBAL_Best[1].size() + GLOBAL_Best[2].size() + GLOBAL_Best[3].size() + GLOBAL_Best[4].size()> LIMITE_GLOBAL)
	{
		if (GLOBAL_Best[5].size() > 0)
			GLOBAL_Best[5].resize(0);
		if (newBest.grid.totalNumbers >= 6)
		{
			mutexGLOBAL_Best.unlock();
			return false;
		}
	}

	if (MAX_INSERT_NUMBERS >= 5 && GLOBAL_Best[0].size() + GLOBAL_Best[1].size() + GLOBAL_Best[2].size() + GLOBAL_Best[3].size()> LIMITE_GLOBAL)
	{
		if (GLOBAL_Best[4].size() > 0)
			GLOBAL_Best[4].resize(0);
		if (newBest.grid.totalNumbers >= 5)
		{
			mutexGLOBAL_Best.unlock();
			return false;
		}
	}
	if (MAX_INSERT_NUMBERS >= 4 && GLOBAL_Best[0].size() + GLOBAL_Best[1].size() + GLOBAL_Best[2].size() > LIMITE_GLOBAL)
	{
		if (GLOBAL_Best[3].size() > 0)
			GLOBAL_Best[3].resize(0);
		if (newBest.grid.totalNumbers >= 4)
		{
			mutexGLOBAL_Best.unlock();
			return false;
		}
	}
	if (GLOBAL_Best[0].size() + GLOBAL_Best[1].size() >= LIMITE_GLOBAL)
	{
		if (GLOBAL_Best[2].size() > 0)
			GLOBAL_Best[2].resize(0);
		if (newBest.grid.totalNumbers >= 3)
		{
			mutexGLOBAL_Best.unlock();
			return false;
		}
	}
	totalGLOBAL = (int)(GLOBAL_Best[0].size() + GLOBAL_Best[1].size() + GLOBAL_Best[2].size() + GLOBAL_Best[3].size() +
		(MAX_INSERT_NUMBERS >= 5 ? GLOBAL_Best[4].size() : 0) + (MAX_INSERT_NUMBERS >= 6 ? GLOBAL_Best[5].size() : 0)
		);


	if (newBest.grid.totalNumbers <= MAX_INSERT_NUMBERS)
	{
		int indexN = max(0, newBest.grid.totalNumbers - 1);
		if (newBest.grid.totalNumbers > 0)
		{
			if (newBest.grid.totalNumbers < MAX_INSERT_NUMBERS)
				for (int d = newBest.grid.totalNumbers; d < MAX_INSERT_NUMBERS; ++d)
					if (GLOBAL_Best[d].size() > 0)
					{ //Remove superset
						for (int b = (int)GLOBAL_Best[d].size() - 1; b >= 0; --b)
						{
							if (newBest.remainingNumbers.isFullyContainedIn(GLOBAL_Best[d][b].remainingNumbers))
							{
								if (b != GLOBAL_Best[d].size() - 1)
									GLOBAL_Best[d][b] = GLOBAL_Best[d][GLOBAL_Best[d].size() - 1];
								GLOBAL_Best[d].pop_back();
							}
						}
					}
			//Already exists

			for (auto&b : GLOBAL_Best[indexN])
			{
				if (b.remainingNumbers.Equals(newBest.remainingNumbers))
				{
					if (newBest.BestScore > b.BestScore)
					{
						bool keepExhausted = b.exhausted;
						b = newBest;
						b.exhausted = keepExhausted;
						if (newBest.grid.totalNumbers <= 4)
							++swapDone;
					}
					mutexGLOBAL_Best.unlock();
					return false;
				}
			}
			//Lower quality
			if (indexN > 0)
				for (int d = 0; d < indexN; ++d)
				{
					for (int b = 0; b < (int)GLOBAL_Best[d].size(); ++b)
					{
						if (GLOBAL_Best[d][b].remainingNumbers.isFullyContainedIn(newBest.remainingNumbers))
						{
							//Exists a better option
							mutexGLOBAL_Best.unlock();
							return false;
						}
					}
				}
		}
		// new Approx
		++totalGLOBAL;
		GLOBAL_Best[indexN].push_back(newBest);
		if (saveToFile && newBest.grid.totalNumbers <= 4)
		{
			if (PROGRAM_NAME != "./test")
			{
				mutexAPROX.lock();
				newBest.SaveApprox("APROX_" + passwordLevel + ".txt");
				mutexAPROX.unlock();
			}
		}
		mutexGLOBAL_Best.unlock();
		return true;
	}
	mutexGLOBAL_Best.unlock();
	return false;
}

bool operator !=(const LAHC_Node &a, const LAHC_Node &b) {

	if (a.Movelists.size() != b.Movelists.size()) return false;
	if (a.ZerosumSets.size() != b.ZerosumSets.size()) return false;
	if (!a.usedInMovesNumbers.isFullyContainedIn(b.usedInMovesNumbers) || !b.usedInMovesNumbers.isFullyContainedIn(a.usedInMovesNumbers)) return false;
	return true;
}

#ifdef SHUFFLE_CODE
bool Mutate_Shuffle(LAHC_Node& currentBest, LAHC_Node& nodeShuff) {
	ShuffleStrat Shuffle_NMB[MAX_NUMBERS];
	++countShuffles;
	for (int TRIES = 0; TRIES < SHUFFLE_TRIES; ++TRIES)
	{
		nodeShuff.clear();
		nodeShuff.MutateType = currentBest.MutateType;
		for (int i = 0; i < COUNT_NUMBERS; ++i)
			Shuffle_NMB[i].clear();

		vector<int> ML = currentBest.Movelists;
		do_shuffle<vector<int>>(ML);
		for (auto& L : ML)
		{
			auto& ST = currentBest.grid.Plan_NMB[L];
			for (auto& m : ST.linear) {
				nodeShuff.grid.ShuffleMove(m, Shuffle_NMB);
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
#ifndef DISABLE_OPTIONALS
inline bool ELE_Mutate_Genestealers_4(LAHC_Node& currentBest, LAHC_Node& candidate)
{

	if (currentBest.Movelists.size() == 0) return false;

	int IgnoreThis = -1;
	if (currentBest.Movelists.size() > 2 && rnd.NextInt(1000) < 350)
	{ //Remove a completed
		IgnoreThis = currentBest.Movelists[rnd.NextInt((uint32_t)currentBest.Movelists.size())];
	}

	vector<pair<Strategy*, Move>> optionals; //Some Gene + a destination that can be cancelled
	LAHC_Node tmpNode;
	tmpNode.clear();
	//Remove all in completed sets

	for (auto&N : currentBest.ZerosumSets)
		if (N != IgnoreThis)
		{
			if (currentBest.grid.Plan_NMB[N].linear.size() > 1 || rnd.NextInt(1000) < 500)
				for (auto&m : currentBest.grid.Plan_NMB[N].linear)
				{
					tmpNode.grid.Val[m.srcIDX] = 0;
					tmpNode.grid.Val[m.destIDX] = 0;
				}
		}

	vector<Strategy> goodToTry;

	if (rnd.NextInt(1000) < 300)
	{
		for (auto&N : currentBest.ZerosumSets)
			if (N != IgnoreThis)
				currentBest.grid.Plan_NMB[N].enumerateOptionals(tmpNode.grid, optionals);
	}
	else {
		for (int N = 0; N < (int)currentBest.Movelists.size(); ++N)
			if (N != IgnoreThis)
				currentBest.grid.Plan_NMB[N].enumerateOptionals(tmpNode.grid, optionals);
	}

	if (optionals.size() == 0) return false;
	uint32_t mark = rnd.xrandom();
	//Now with incomplete sets
	for (auto&I : currentBest.IncompleteSets) {
		auto& LastMove = currentBest.grid.Plan_NMB[I].linear[currentBest.grid.Plan_NMB[I].linear.size() - 1];

		for (auto& o : optionals)
			if (o.second.srcVal == LastMove.destVal && o.second.destIDX == LastMove.destIDX) //Cancellable!
			{
				Strategy SST = *o.first;
				SST.criticalPath[SST.criticalPath.size() - 1] = o.second;
				SST.criticalPath[SST.criticalPath.size() - 1].destVal = LastMove.destVal;
				SST.linear[SST.linear.size() - 1] = o.second;
				SST.linear[SST.linear.size() - 1].destVal = LastMove.destVal;
				SST.addPlan(currentBest.grid.Plan_NMB[I]);
				goodToTry.push_back(SST);
			}
	}
	vector<int> UnusedNumbers;
	for (int i = 0; i < COUNT_NUMBERS; ++i) {
		auto val = tmpNode.grid.Val[i];
		if (val != 0)
		{
			//Search 
			for (auto& o : optionals)
				if (o.second.srcVal == val && o.second.destIDX == i) //Cancellable!
				{
					Strategy SST = *o.first;
					UnusedNumbers.push_back(SST.criticalPath[SST.criticalPath.size() - 1].destIDX);
					SST.criticalPath[SST.criticalPath.size() - 1] = o.second;
					SST.criticalPath[SST.criticalPath.size() - 1].destVal = val;
					SST.linear[SST.linear.size() - 1] = o.second;
					SST.linear[SST.linear.size() - 1].destVal = val;
					goodToTry.push_back(SST);

				}
		}
	}


	if (goodToTry.size() == 0) return false;

	int chanceCuts = rnd.NextInt(950, 1050);

	do_shuffle<vector<Strategy>>(goodToTry);
	int amountToTry = (int)rnd.NextInt(1, (int)goodToTry.size());
	vector<int> emptyList;
	for (int i = 0; i < amountToTry; ++i)
	{

		if (rnd.NextInt(1000) < 800)
			goodToTry[i].ApplyMoves(candidate.grid, 1050, mark, false, UnusedNumbers);
		else goodToTry[i].ApplyMoves(candidate.grid, 1050, mark, false, emptyList);
	}

	//Now completes
	auto Z = currentBest.ZerosumSets;
	for (auto&i : currentBest.IncompleteSets)
		if (rnd.NextInt(1000) < 800)
		{
			Z.push_back(i);
			//UEL.ApplyMoves(candidate.grid, chanceCuts, mark, withoutEndpoint, UnusedNumbers);
		}
	do_shuffle<vector<int>>(Z);
	do_shuffle<vector<int>>(UnusedNumbers);
	for (auto&z : Z)
		if (z != IgnoreThis)
		{
			currentBest.grid.Plan_NMB[z].ApplyMoves(candidate.grid, chanceCuts, mark, false, UnusedNumbers);
		}
	while (candidate.grid.doRandomMove()) {}
	candidate.MutateType = 4;
	return true;
}

void Inception_Optionals(vector<Strategy*>& opt, Strategy& T)
{
	if (T.optionalCount > 0)
	{
		for (auto& o : T.optionals)
			for (auto& s : o.second)
				Inception_Optionals(opt, s);
	}
	opt.push_back(&T);
}

inline bool Mutate_Genestealers_4(LAHC_Node& currentBest, LAHC_Node& candidate)
{
	vector<Strategy*> optionals;
	for (int N = 0; N < currentBest.Movelists.size(); ++N)
	{
		auto& ST = currentBest.grid.Plan_NMB[N];
		for (auto& u : ST.optionals)
		{
			for (auto& s : u.second)
			{
				Inception_Optionals(optionals, s);
			}
		}
	}

	if (optionals.size() == 0)
		return false;
	do_shuffle<vector<Strategy*>>(optionals);
	//Take some of them
	vector<int> TryMerge;
	uint64_t maskX = 0;
	uint64_t maskY = 0;
	for (int i = 0; i < COUNT_NUMBERS; ++i)
		if (!currentBest.usedInMovesNumbers.get(i))
		{
			TryMerge.push_back(i);
			maskX = 1ULL << (NMB[i].X);
			maskY = 1ULL << (NMB[i].Y);
		}
	/*
	unordered_map<int, vector< tuple<Strategy*, int, bool>>> Optiplus;
	{
	for (auto& I : optionals)
	{
	auto& ST = *I;
	if (ST.linear.size() > 1)
	{

	auto& lastMove = ST.linear[ST.linear.size() - 2];
	bool canAdd = lastMove.srcVal + lastMove.destVal <= W;
	if (Optiplus.find(lastMove.destIDX) == Optiplus.end())
	{
	Optiplus[lastMove.destIDX].reserve(canAdd ? 2 : 1);
	Optiplus[lastMove.destIDX].resize(0);
	}
	tuple<Strategy*, int, bool> P = { I,lastMove.destVal - lastMove.srcVal, false };
	Optiplus[lastMove.destIDX].push_back(P);
	if (canAdd)
	{
	P = { I,lastMove.destVal + lastMove.srcVal, true };
	Optiplus[lastMove.destIDX].push_back(P);
	}
	}
	}
	}

	for (int i=0; i < COUNT_NUMBERS;++i)
	{
	if (!currentBest.usedNumbers.get(i))
	{
	TryMerge.push_back(i);
	//Try to hit some of the Optionals
	vector<Move> posMoves;
	candidate.grid.ExplodeMoves(i, posMoves);
	for (auto& pm : posMoves)
	{
	if (Optiplus.find(pm.destIDX) != Optiplus.end())
	{
	for (auto &oo : Optiplus[pm.destIDX])
	{
	if (get<1>(oo) == pm.srcVal)
	{//JACKPOT
	auto& ST = *get<0>(oo);
	auto& lastMove = ST.linear[ST.linear.size() - 2];
	for (auto&m : ST.linear)
	{
	if (m != lastMove)
	{
	candidate.grid.ApplyMove(m);
	}
	else {
	Move CustMove = lastMove;
	CustMove.sign = (get<2>(oo) ? 1 : 0);
	candidate.grid.ApplyMove(CustMove);
	posMoves.resize(0);
	candidate.grid.ExplodeMoves(i, posMoves);
	for (auto& uj : posMoves)
	{
	if (uj.destIDX == CustMove.destIDX)
	{
	candidate.grid.ApplyMove(uj);
	break;
	}
	}
	break;
	}
	}
	}
	}
	}
	}
	}
	}*/
	int tgtSize = min((int)optionals.size(), (int)rnd.NextInt(1, 3));
	if (optionals.size() > tgtSize)
		optionals.resize(tgtSize);

	unordered_set<size_t> ForbiddenMoves;

	for (auto&o : optionals)
	{
		auto& ST = *o;
		auto& lastMove = ST.linear[ST.linear.size() - 1];
		ForbiddenMoves.emplace(lastMove.CreateHash());
		//for (auto&m : o->linear)
		if (ST.linear.size() > 1)
			for (int xmi = 0; xmi < ST.linear.size() - 1; ++xmi)
			{
				Move m = ST.linear[xmi];
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.srcIDX, ForbiddenMoves);
				if ((ST.linear.size() - 2 == xmi) && (m.srcVal + m.destVal <= W) && (rnd.NextInt(100) < K_ALLOWSUM))
				{
					m.sign = 1;
				}
				else m.sign = 0;
				candidate.grid.ApplyMove(m);
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.destIDX, ForbiddenMoves);
			}
		TryMerge.push_back(lastMove.srcIDX);
		maskX = 1ULL << (NMB[lastMove.srcIDX].X);
		maskY = 1ULL << (NMB[lastMove.srcIDX].Y);
	}
	for (auto&I : currentBest.IncompleteSets)
	{
		auto& NL = currentBest.grid.Plan_NMB[I];
		for (int xmi = 0; xmi < NL.linear.size(); ++xmi)
		{
			Move m = NL.linear[xmi];
			if (NL.linear.size() - 1 == xmi)
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.srcIDX, ForbiddenMoves);
			if ((NL.linear.size() - 1 == xmi) && (rnd.NextInt(100) < K_ALLOWSUM))
			{
				m.sign = 1;
			}
			else m.sign = 0;
			if (ForbiddenMoves.find(m.CreateHash()) == ForbiddenMoves.end())
				candidate.grid.ApplyMove(m);
			if (((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.destIDX, ForbiddenMoves);
		}
	}
	for (auto&I : currentBest.ZerosumSets)
	{
		auto& NL = currentBest.grid.Plan_NMB[I];
		for (int xmi = 0; xmi < NL.linear.size(); ++xmi)
		{
			Move m = NL.linear[xmi];
			if (((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.srcIDX, ForbiddenMoves);
			if (ForbiddenMoves.find(m.CreateHash()) == ForbiddenMoves.end())
				candidate.grid.ApplyMove(m);
			if (((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.destIDX, ForbiddenMoves);
		}
	}
	/*	unordered_map<int, vector< tuple<int, int,bool>>> Incomp; // NMB INDEX -> vector of Strategy ID + CellValue +was (+) or (-)
	{
	vector<int> SHUFI = currentBest.IncompleteSets;
	do_shuffle<vector<int>>(SHUFI);
	for (auto& I : SHUFI) {
	auto& ST = currentBest.grid.Plan_NMB[I];
	auto& lastMove = ST.linear[ST.linear.size() - 1];
	bool canAdd = lastMove.srcVal + lastMove.destVal <= W;
	if (Incomp.find(lastMove.destIDX) == Incomp.end())
	{
	Incomp[lastMove.destIDX].reserve(canAdd?2:1);
	Incomp[lastMove.destIDX].resize(0);
	}
	tuple<int, int, bool> P = { I,lastMove.destVal - lastMove.srcVal, false };
	Incomp[lastMove.destIDX].push_back(P);
	if (canAdd)
	{
	P = { I,lastMove.destVal + lastMove.srcVal, true };
	Incomp[lastMove.destIDX].push_back(P);
	}
	}
	}
	*/


	candidate.MutateType = 4;
	return true;
}
#endif

inline bool Mutate_Opt_4(LAHC_Node& currentBest, LAHC_Node& candidate) {

	int countOpt = rnd.NextInt(1, 3);
	vector<int> MM = currentBest.Movelists;
	do_shuffle<vector<int>>(MM);
	vector<int> TryMerge;
	unordered_set<uint64_t> ForbiddenMoves;
	for (int times = 0; times < countOpt; ++times)
	{
		bool completeOpt = false;
		for (auto& N : MM)
		{
			auto& ST = currentBest.grid.Plan_NMB[N];
			int tam = (int)ST.linear.size();
			if (tam < 2)
				continue;
			int index = rnd.NextInt(tam);
			for (int i = 0; i < tam; ++i)
			{
				++index;
				if (index >= tam)
					index = 0;
				auto& m = ST.linear[index];
				if ((m.sign == 0) && (2 * m.destVal == m.srcVal))
				{ //Optional
					TryMerge.push_back(m.srcIDX);
					ForbiddenMoves.emplace(m.CreateHash());
					if (rnd.NextInt(1000) < 300 && index > 0)
					{ //Also remove previous one
						auto& prev = ST.linear[index - 1];
						if (prev.destIDX == m.srcVal)
						{
							//ForbiddenMoves.emplace(prev.CreateHash());
							TryMerge.push_back(prev.srcIDX);
						}
					}
					completeOpt = true;
					break;
				}
			}
			if (completeOpt)
				break;
		}
	}
	if (TryMerge.size() == 0)
		return false;
	for (auto&I : MM)
	{
		auto& NL = currentBest.grid.Plan_NMB[I];
		for (int xmi = 0; xmi < NL.linear.size(); ++xmi)
		{
			Move m = NL.linear[xmi];
			bool attempt = rnd.NextInt(1000) < 700;
			if (attempt)
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.srcIDX, ForbiddenMoves);

			bool attempt2 = rnd.NextInt(1000) < 100;
			if (attempt2)
			{
				for (auto& IM : TryMerge)
				{
					Move* tm = CACHE_MERGE[IM][m.destIDX];
					if (tm != nullptr && tm->srcVal == candidate.grid.Val[IM])
					{
						Move pega = *tm;
						pega.sign = rnd.NextInt(100) < K_ALLOWSUM ? 1 : 0;
						pega.destVal = candidate.grid.Val[m.destIDX];
						candidate.grid.ApplyMove(pega);
					}
				}
			}
			if (ForbiddenMoves.find(m.CreateHash()) == ForbiddenMoves.end())
				candidate.grid.ApplyMove(m);

			if (attempt)
				for (auto& IM : TryMerge)
					candidate.grid.doMerge(IM, m.destIDX, ForbiddenMoves);
		}
	}
	candidate.MutateType = 4;
	return true;
}


inline void Mutate_Completes_0(LAHC_Node& currentBest, LAHC_Node& candidate)
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
		int chanceCuts = rnd.NextInt(970, 1050);
		uint32_t mark = (uint32_t)rnd.xrandom();

		//Some chance to use incompletes
		for (auto&i : currentBest.IncompleteSets)
			if (rnd.NextInt(1000) < 800)
			{
				Z.push_back(i);
				//UEL.ApplyMoves(candidate.grid, chanceCuts, mark, withoutEndpoint, UnusedNumbers);
			}
		do_shuffle<vector<int>>(Z);


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
			currentBest.grid.Plan_NMB[z].ApplyMoves(candidate.grid, chanceCuts, mark, withoutEndpoint, UnusedNumbers);
		}
	}
}


/* Random move, but from endings*/
int rndLowBound(const vector<size_t>& ListSizes, const int& valPoint)
{
	int selList = (int)(lower_bound(ListSizes.begin(), ListSizes.end(), valPoint, [](auto &a, auto &b) { return a <= b; }) - ListSizes.begin());
	return min(max(0, selList), (int)ListSizes.size() - 1);
}

//Pick a random move
void selectorRandom(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves)
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


void selectorRandomTail(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves)
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
void selectorOptional(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves)
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
void selectorOrphanMoves(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves, const vector<int>& orphanNumbers)
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
void selectorCrossMoves(unordered_set<int>& tmpPointsToRemove, const LAHC_Node& currentBest, const vector<size_t>& ListSizes, const size_t& totalMoves, const vector<int>& pendingNumbers)
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
inline void Mutate_Strings_4(const LAHC_Node& currentBest, LAHC_Node& candidate)
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
				selectorOrphanMoves(tmpPointsToRemove, currentBest, ListSizes, totalMoves, orphanmoves);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL)
				selectorRandomTail(tmpPointsToRemove, currentBest, ListSizes, totalMoves);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL + K_SEL_VAL_OPT)
				selectorOptional(tmpPointsToRemove, currentBest, ListSizes, totalMoves);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL + K_SEL_VAL_OPT + K_SEL_VAL_RANDOM)
				selectorRandom(tmpPointsToRemove, currentBest, ListSizes, totalMoves);
			else if (chance < K_SEL_VAL_ORPHAN + K_SEL_VAL_TAIL + K_SEL_VAL_OPT + K_SEL_VAL_RANDOM + K_SEL_VAL_CROSS)
				selectorCrossMoves(tmpPointsToRemove, currentBest, ListSizes, totalMoves, pendingNumbers);
			else selectorRandom(tmpPointsToRemove, currentBest, ListSizes, totalMoves);
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
			assert(selList < currentBest.Movelists.size());
			auto& ML = currentBest.grid.Plan_NMB[currentBest.Movelists[selList]];
			assert(ML.linear.size() > 0);
			assert(selPoint < ML.linear.size());
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



		do_shuffle<vector<int>>(pendingNumbers);
		vector<Move> tmpMoves;
		//Now repeat all, trying to Merge / explode but without forbidden moves.
		vector<int> Listas = currentBest.IncompleteSets;// currentBest.Movelists;
		if (rnd.NextInt(1000) < 500)
		{
			do_shuffle<vector<int>>(Listas);
			auto HJU = currentBest.ZerosumSets;
			do_shuffle<vector<int>>(HJU);
			for (auto&I : HJU)
				Listas.push_back(I);

		}
		else {
			for (auto&I : currentBest.ZerosumSets)
				Listas.push_back(I);
			do_shuffle<vector<int>>(Listas);
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
									candidate.grid.doMerge(t, uwm.srcIDX);
								}
					}
					bool forced = false;
					if (rnd.NextInt(1000) < K_CHANCE_EXPLODE)
					{
						if (((1ULL << NMB[uwm.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[uwm.srcIDX].Y) & maskY) != 0)
							for (auto& unused : pendingNumbers)
							{
								if (candidate.grid.doForceJoin(unused, uwm.destIDX, 1000)) {
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
									candidate.grid.doMerge(t, uwm.destIDX);
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
		/*	for (auto&L : Listas)
		{
		auto& ML = currentBest.grid.Plan_NMB[L];
		for (auto&uwm : ML.linear)
		if (forbiddenMoves.find(uwm.CreateHash()) == forbiddenMoves.end()
		&& (candidate.grid.Val[uwm.srcIDX] != 0 && candidate.grid.Val[uwm.destIDX] != 0))
		{

		Move m = uwm;
		if (rnd.NextInt(1000) < K_CHANCE_MERGE)
		{
		for (auto& t : endings)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.srcIDX);

		if (((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
		for (auto& t : pendingNumbers)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.srcIDX);
		for (auto& t : cutIDs)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.srcIDX);
		}
		if (rnd.NextInt(1000) < K_CHANCE_EXPLODE)
		{
		for (auto& t : endings)
		if (t != m.srcIDX && t != m.destIDX && candidate.grid.doForceJoin(t, m.destIDX, 1000))
		m.destVal = candidate.grid.Val[m.destIDX];
		if (((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
		for (auto& t : pendingNumbers)
		if (t != m.srcIDX && t != m.destIDX && candidate.grid.doForceJoin(t, m.destIDX, 1000))
		m.destVal = candidate.grid.Val[m.destIDX];
		for (auto& t : cutIDs)
		if (t != m.srcIDX && t != m.destIDX && candidate.grid.doForceJoin(t, m.destIDX, 1000))
		m.destVal = candidate.grid.Val[m.destIDX];
		}
		if (candidate.grid.Val[m.destIDX] == 0 || candidate.grid.Val[m.srcIDX] == 0)
		continue;
		if (candidate.grid.ApplyMove(m)) {
		//	consumed.set(m.destIDX);
		}
		if (rnd.NextInt(1000) < K_CHANCE_MERGE)
		{
		for (auto& t : endings)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.destIDX);
		if (((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
		for (auto& t : pendingNumbers)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.destIDX);
		for (auto& t : cutIDs)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.destIDX);
		}

		// for (auto& t : endings)
		// candidate.grid.doForceJoin(t, m.destIDX, K_CHANCE_EXPLODE);
		// for (auto& t : pendingNumbers)
		// candidate.grid.doForceJoin(t, m.destIDX, K_CHANCE_EXPLODE);
		// for (auto& t : cutIDs)
		// candidate.grid.doForceJoin(t, m.destIDX, K_CHANCE_EXPLODE);


		}
		#ifdef _MSC_VER
		else {
		countRemove += 0;
		}
		#endif



		if (ML.Endpoint == 0)
		{
		endings.push_back(ML.linear.back().destIDX);
		}
		}
		#ifndef fuera
		pendingNumbers.resize(0);
		for (int i = 0; i < COUNT_NUMBERS; ++i)
		if (candidate.grid.Val[i] > 0)
		pendingNumbers.push_back(i);
		do_shuffle<vector<int>>(pendingNumbers);
		if (pendingNumbers.size() > 5)
		pendingNumbers.resize(5);
		//Now repeat all, trying to Merge / explode but without forbidden moves.
		do_shuffle<vector<int>>(Listas);
		for (auto&L : Listas)
		{
		auto& ML = currentBest.grid.Plan_NMB[L];
		for (auto&uwm : ML.linear)
		if (candidate.grid.Val[uwm.srcIDX] != 0 && candidate.grid.Val[uwm.destIDX] != 0)
		{
		Move m = uwm;
		if (candidate.grid.Val[m.srcIDX] == 0)
		continue;
		auto chanceMerge = rnd.NextInt(1000);
		if (chanceMerge < K_CHANCE_MERGE)
		{
		for (auto& t : endings)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.srcIDX);
		if (((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
		for (auto& t : pendingNumbers)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.srcIDX);
		for (auto& t : cutIDs)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.srcIDX);
		}
		auto chanceEXP = rnd.NextInt(1000);
		if (chanceEXP < K_CHANCE_EXPLODE)
		{
		for (auto& t : endings)
		if (t != m.srcIDX && t != m.destIDX &&candidate.grid.doForceJoin(t, m.destIDX, 1000))
		m.destVal = candidate.grid.Val[m.destIDX];
		if (((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
		for (auto& t : pendingNumbers)
		if (t != m.srcIDX && t != m.destIDX && candidate.grid.doForceJoin(t, m.destIDX, 1000))
		m.destVal = candidate.grid.Val[m.destIDX];
		for (auto& t : cutIDs)
		if (t != m.srcIDX && t != m.destIDX &&candidate.grid.doForceJoin(t, m.destIDX, 1000))
		m.destVal = candidate.grid.Val[m.destIDX];
		}

		m.srcVal = candidate.grid.Val[m.srcIDX];
		m.destVal = candidate.grid.Val[m.destIDX];
		if (m.srcVal == 0 || m.destVal == 0)
		continue;
		if (candidate.grid.ApplyMove(m)) {
		//	consumed.set(m.destIDX);
		}
		if (rnd.NextInt(1000) < K_CHANCE_MERGE)
		{
		for (auto& t : endings)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.destIDX);
		if (((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
		for (auto& t : pendingNumbers)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.destIDX);
		for (auto& t : cutIDs)
		if (t != m.srcIDX && t != m.destIDX)
		candidate.grid.doMerge(t, m.destIDX);
		}
		// for (auto& t : endings)
		// candidate.grid.doForceJoin(t, m.destIDX, K_CHANCE_EXPLODE);
		// for (auto& t : pendingNumbers)
		// candidate.grid.doForceJoin(t, m.destIDX, K_CHANCE_EXPLODE);
		// for (auto& t : cutIDs)
		// candidate.grid.doForceJoin(t, m.destIDX, K_CHANCE_EXPLODE);

		}
		}
		#endif
		*/
	}
	candidate.MutateType = 4;
}


inline void Mutate_Points_1(LAHC_Node& currentBest, LAHC_Node& candidate)
{
	candidate.MutateType = 1;
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
	do_shuffle<vector<int>>(K);
	for (auto&k : K)
		Add(linearMoves, currentBest.grid.Plan_NMB[k].linear);
	if (linearMoves.size() > 0)
	{
		uint64_t maskX = 0;
		uint64_t maskY = 0;
		bool tryMerge = /*currentBest.Movelists.size() > 8 &&*/ (currentBest.usedInMovesNumbers.count() > COUNT_NUMBERS * 50 / 100) && rnd.NextInt(100) < 80;
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
						candidate.grid.doMerge(unused, linearMoves[i].srcIDX);
					}
				}
				candidate.grid.ApplyMove(linearMoves[i]);
				BU = (((1ULL << NMB[linearMoves[i].destIDX].X) & maskX) != 0 || ((1ULL << NMB[linearMoves[i].destIDX].Y) & maskY) != 0) && tryMerge && (rnd.NextInt(1000) < 900);
				if (BU)
				{
					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, linearMoves[i].destIDX);
					}
				}

			}
	}

}
inline void Mutate_Lists_2(LAHC_Node& currentBest, LAHC_Node& candidate)
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
						candidate.grid.doMerge(unused, UEL.linear[j].srcIDX);
					}
				}
				candidate.grid.ApplyMove(UEL.linear[j]);
				BU = (((1ULL << NMB[UEL.linear[j].destIDX].X) & maskX) != 0 || ((1ULL << NMB[UEL.linear[j].destIDX].Y) & maskY) != 0)
					&& tryMerge && (rnd.NextInt(1000) < 800);
				if (BU)
				{
					for (auto& unused : UnusedNumbers)
					{
						candidate.grid.doMerge(unused, UEL.linear[j].destIDX);
					}
				}

			}
		}
	}
}
inline void Mutate_Special_3(LAHC_Node& currentBest, LAHC_Node& candidate)
{
	candidate.MutateType = 3;
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
			if (currentBest.grid.Plan_NMB[z].linear.size() > 1 || (rnd.NextInt(1000) < 900)) {
				K.push_back(z);
				for (auto&l : currentBest.grid.Plan_NMB[z].linear)
				{
					//linearMoves.push_back(l);
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
					//linearMoves.push_back(l);
					UsedNumbers.set(l.srcIDX);
					//UsedNumbers.set(l.destIDX);
				}
			}
		}
	do_shuffle<vector<int>>(K);



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
		/*	for (auto&u : UnusedNumbers)
		{
		explotar.push_back(pair<int, int>{(int)rnd.NextInt((uint32_t)linearMoves.size()), u});
		}*/
		do_shuffle< vector<pair<int, int>>>(explotar);
		do_shuffle<vector<int>>(UnusedNumbers);
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
						candidate.grid.doMerge(unused, linearMoves[i].srcIDX);
					}
				}
				if (tryMerge && (rnd.NextInt(1000) < 900))
				{
					for (auto& U : explotar)
						if (U.first > i + 1)
						{
							candidate.grid.doMerge(U.second, linearMoves[i].srcIDX);
						}
				}

				if (CloseIDX.get(linearMoves[i].destIDX)) {
					for (auto& U : explotar)
						if (U.first > i)
						{
							tmpMoves.resize(0);
							candidate.grid.ExplodeMoves(U.second, tmpMoves);
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
						candidate.grid.doMerge(unused, linearMoves[i].destIDX);
					}
				}
				if (tryMerge && (rnd.NextInt(1000) < 900))
				{
					for (auto& U : explotar)
						if (U.first > i)
						{
							candidate.grid.doMerge(U.second, linearMoves[i].destIDX);
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
//Will remove a % of moves affecting a X,Y of remaining numbers.
inline bool Mutate_Cross_4(LAHC_Node& currentBest, LAHC_Node& candidate) {
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
	do_shuffle<vector<size_t>>(OPT);
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
	do_shuffle<vector<int>>(K);
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
					candidate.grid.doMerge(unused, lm.srcIDX);
				}
			}
			candidate.grid.ApplyMove(lm);
			BU = (((1ULL << NMB[lm.destIDX].X) & maskX) != 0 || ((1ULL << NMB[lm.destIDX].Y) & maskY) != 0)
				&& tryMerge && (rnd.NextInt(1000) < 900);
			if (BU)
			{
				for (auto& unused : UnusedNumbers)
				{
					candidate.grid.doMerge(unused, lm.destIDX);
				}
			}

		}
	}
	candidate.MutateType = 4;
	return true;
}
bool Mutate_Shuffle_Incomplete_End_5(LAHC_Node& currentBest, LAHC_Node& candidate)
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
	do_shuffle<vector<int>>(INC);
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
	do_shuffle<vector<int>>(TryMerge);
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
						candidate.grid.doMerge(t, m.srcIDX);
				candidate.grid.ApplyMove(m);
				if (rnd.NextInt(1000) < 500)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.destIDX);
			}
		}
		do_shuffle<vector<Move>>(lastMoves);
		for (auto& m : lastMoves)
		{
			m.sign = rnd.NextInt(100) < K_ALLOWSUM ? 1 : 0;
			ASSERT(NMB[m.destIDX].X == lastX);
			ASSERT(NMB[m.destIDX].Y == lastY);
			m.destVal = candidate.grid.Val[lastIDX];
			if (((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
				if (rnd.NextInt(1000) < K5_MERGE0)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.srcIDX);
			candidate.grid.ApplyMove(m);
			if (((1ULL << NMB[m.srcIDX].X) & maskX) != 0 || ((1ULL << NMB[m.srcIDX].Y) & maskY) != 0)
				if (rnd.NextInt(1000) < K5_MERGE0)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.destIDX);
		}
		//if (candidate.grid.Val[lastIDX] != 0) 
		{
			TryMerge.push_back(lastIDX);
			maskX = 1ULL << (NMB[lastIDX].X);
			maskY = 1ULL << (NMB[lastIDX].Y);
		}
	}
	do_shuffle<vector<int>>(TryMerge);
	vector<int> NMI = currentBest.IncompleteSets;
	if (rnd.NextInt(1000) < 500)
	{
		do_shuffle<vector<int>>(NMI);
		auto HJU = currentBest.ZerosumSets;
		do_shuffle<vector<int>>(HJU);
		for (auto&I : HJU)
			NMI.push_back(I);

	}
	else {
		for (auto&I : currentBest.ZerosumSets)
			NMI.push_back(I);
		do_shuffle<vector<int>>(NMI);
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
						candidate.grid.doMerge(t, m.srcIDX);
				candidate.grid.ApplyMove(m);
				if ((((1ULL << NMB[m.destIDX].X) & maskX) != 0 || ((1ULL << NMB[m.destIDX].Y) & maskY) != 0)
					&& rnd.NextInt(1000) < K5_MERGE1)
					for (auto& t : TryMerge)
						candidate.grid.doMerge(t, m.destIDX);
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
		for (int i = 0; i<LM.linear.size(); ++i)
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
			for (int llm = 0; llm<LM->linear.size(); ++llm)
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
/*
struct EXH {
int indexMoveList=-1;
int destIDX;
int endValue;
vector<int> interValue;
vector<Move> moves;
vector<int> movePos;

void create(LAHC_Node& currentBest,int _destIDX, int moveList = -1 ) {

Strategy* LM =nullptr;
indexMoveList = moveList;
if (moveList < 0)
{
//Search the list
for (auto& M:currentBest.Movelists)
{
for (auto&m : currentBest.grid.Plan_NMB[M].linear)
{
if (m.srcIDX == _destIDX || m.destIDX == _destIDX) {
indexMoveList = M;
break;
}
}
if (indexMoveList >= 0)
break;
}
}

if (indexMoveList >=0)
LM = &currentBest.grid.Plan_NMB[indexMoveList];

destIDX = _destIDX;
endValue = currentBest.grid.Val[destIDX];
//for (auto&m : LM.linear)
moves.resize(0);
interValue.resize(0);
movePos.resize(0);
if (LM != nullptr)
for (int llm = 0; llm<LM->linear.size(); ++llm)
{
auto& m = LM->linear[llm];
if (m.destIDX == destIDX)
{
moves.push_back(m);
movePos.push_back(llm);
}
}
}
};
*/
void permutateEXHSender(const LAHC_Node& currentBest, vector<EXHSender>& src, unordered_map<int, vector<EXHSender>>& dest) {
	dest.clear();
	Grid G;
	vector<Move>tmpMoves;
	Stopwatch muchoTiempo;
	muchoTiempo.Start(4 * 1000 * 1000);
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
			uint64_t signos = (1ULL << min(6, (int)E.movePos.size()));
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
				G.ExplodeMoves(E.senderIDX, tmpMoves);
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
		if (SS.size()>1)
			for (int i = 0; i < SS.size() - 1; ++i)
			{
				for (int j = i + 1; j < SS.size(); ++j) {
					if (SS[i].endValue == SS[j].endValue)
					{
						cerr << "error duplicado" << endl;
					}
				}
			}

	}
}
void permutateEXHMerger(const LAHC_Node& currentBest, vector<EXHMerger>& src, vector<EXHMerger>& dest, int searchMergeValue, int forceENDVALUE = -1)
{
	Stopwatch muchoTiempo;
	muchoTiempo.Start(4 * 1000 * 1000);
	dest.resize(0);
	for (auto&E : src)
	{
		do {
			uint64_t signos = (1ULL << min(6, (int)E.movePos.size()));
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

bool Mutate_Exhaustive_6(const LAHC_Node& currentBest, LAHC_Node& candidate, LAHC_Node& tmpCandidate)
{
	bool mejorado = false;
	vector<EXHSender> incompletes;
	//Untouched numbers
	for (auto&N : NMB)
	{
		if (currentBest.untouchedNumbers.get(N.ID) && currentBest.grid.Val[N.ID] >0) {
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
	permutateEXHSender(currentBest, incompletes, permIncompletes);
	int totalPermutations = 0;

	vector<int> numeros;
	for (auto&I : permIncompletes)
	{
		numeros.push_back(I.first);
		totalPermutations += (int)I.second.size();
		std::sort(I.second.begin(), I.second.end(), [=](auto& a, auto& b) {return b.endValue > a.endValue; });
	}

	if (numeros.size()>1)
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
								tmpCandidate.grid.doMerge(UUm.destIDX, m.srcIDX);
								tmpCandidate.grid.doMerge(UUm.destIDX, m.destIDX);
								tmpCandidate.grid.ApplyMove(m);
								tmpCandidate.grid.doMerge(UUm.destIDX, m.srcIDX);
								tmpCandidate.grid.doMerge(UUm.destIDX, m.destIDX);
							}
						}
					while (tmpCandidate.grid.doRandomMove()) {}
					tmpCandidate.MutateType = 6;
					tmpCandidate.CalcStats();
					if (tmpCandidate.Movelists.size() == 0)
						continue;
					++SUPERCOMBINATOR;
					GLOBAL_Best_Insert(tmpCandidate, false);
					if (tmpCandidate.BestScore > candidate.BestScore)
					{
						swap(candidate, tmpCandidate);
					}
					if (candidate.BestScore > currentBest.BestScore)
					{
						mejorado = true;
						++JACKPOT_OK;
					}
					if (IS_SOLVED(candidate, "JACKPOT Inc-Inc"))
						return true;
				}

			}
		}
	if (currentBest.grid.totalNumbers > K_EXHAUSTIVE_MIN_NUMBERS + 4)
		return false;
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
		for (EXHSender& P : AP.second)
		{
			if (solved)
				return false;
			//**********Self Jackpot************************
			if (P.endValue == 0 && P.indexMoveList >= 0)
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
				while (tmpCandidate.grid.doRandomMove()) {}
				tmpCandidate.MutateType = 6;
				tmpCandidate.CalcStats();
				if (tmpCandidate.Movelists.size() == 0)
					continue;
				GLOBAL_Best_Insert(tmpCandidate, false);
				if (tmpCandidate.BestScore > candidate.BestScore)
				{
					swap(candidate, tmpCandidate);
				}
				if (candidate.BestScore > currentBest.BestScore)
				{
					mejorado = true;
					++JACKPOT_OK;
				}
				if (IS_SOLVED(candidate, "JACKPOT Inc-Self0"))
					return true;
				break;
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
					permutateEXHMerger(currentBest, targetMerge, explodeTargetMerge, P.endValue, expectedEndValue);
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
						for (int i = 0; i<dest.lastPosMove; ++i)
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
								tmpCandidate.grid.doMerge(P.senderIDX, m.srcIDX);
								tmpCandidate.grid.doMerge(P.senderIDX, m.destIDX);
								tmpCandidate.grid.ApplyMove(m);
								tmpCandidate.grid.doMerge(P.senderIDX, m.srcIDX);
								tmpCandidate.grid.doMerge(P.senderIDX, m.destIDX);
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
					while (tmpCandidate.grid.doRandomMove()) {
					}
					tmpCandidate.MutateType = 6;
					tmpCandidate.CalcStats();
					if (tmpCandidate.Movelists.size() == 0)
						continue;
					GLOBAL_Best_Insert(tmpCandidate, false);
					if (tmpCandidate.BestScore > candidate.BestScore)
					{
						swap(candidate, tmpCandidate);
					}
					if (candidate.BestScore > currentBest.BestScore)
					{
						mejorado = true;
						++JACKPOT_OK;
					}
					if (IS_SOLVED(candidate, "JACKPOT Inc-Merge"))
						return true;
				}
			}

		}
	return mejorado;
}


inline bool reArmCandidate(LAHC_Node& currentBest, LAHC_Node& tmpCandidate) {


	vector<int> UnusedNumbers;
	//for (auto& n : NMB)
	for (int n = 0; n < NMB.size(); ++n)
	{
		if (currentBest.grid.Val[n] > 0)
			UnusedNumbers.push_back(n);
	}
	tmpCandidate.clear();
	vector<int> INCop = currentBest.IncompleteSets;
	do_shuffle<vector<int>>(INCop);
	for (auto& I : INCop)
	{
		auto& UEL = currentBest.grid.Plan_NMB[I];
		for (auto&m : UEL.linear)
		{
			bool intentaMerge = rnd.NextInt(1000) < 750;
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.srcIDX);
			tmpCandidate.grid.ApplyMove(m);
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.destIDX);
		}
	}
	INCop.resize(0);
	vector<int> COMP = currentBest.ZerosumSets;
	do_shuffle<vector<int>>(COMP);
	for (auto& I : COMP)
	{
		auto& UEL = currentBest.grid.Plan_NMB[I];
		for (auto&m : UEL.linear)
		{
			bool intentaMerge = rnd.NextInt(1000) < 750;
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.srcIDX);
			tmpCandidate.grid.ApplyMove(m);
			if (intentaMerge)
				for (auto& unused : UnusedNumbers)
					tmpCandidate.grid.doMerge(unused, m.destIDX);

		}
	}
	COMP.resize(0);
	while (tmpCandidate.grid.doRandomMove()) {}
	tmpCandidate.MutateType = currentBest.MutateType;
	tmpCandidate.CalcStats();
	return true;
}

inline void Mutate(LAHC_Node& currentBest, LAHC_Node& candidate, double PROB[])
{
	candidate.clear();
	if (currentBest.Movelists.size() > 0)
	{

		//Mutate_Strings_4(currentBest, candidate);

		double P0 = (currentBest.ZerosumSets.size() > 0 ? PROB[0] : 0.0);
		double P2 = PROB[2] + P0;
		double P3 = PROB[3] + P2;
		double P4 = PROB[4] + P3;
		double P5 = (currentBest.IncompleteSets.size() > 0 ? PROB[5] : 0.0) + P4;
		double PTOTAL = PROB[1] + P5;

		double r = rnd.NextFloat()*PTOTAL;//  (double)rnd.NextInt((int)(PTOTAL*10000.0)) / 10000.0;

		if (r < P0) { Mutate_Completes_0(currentBest, candidate); }
		else if (r < P2) { Mutate_Lists_2(currentBest, candidate); }
		else if (r < P3) { Mutate_Special_3(currentBest, candidate); }
		else if (r < P4) { Mutate_Strings_4(currentBest, candidate); }// { Mutate_Cross_4(currentBest, candidate); }
		else if (r < P5) { Mutate_Shuffle_Incomplete_End_5(currentBest, candidate); }

		if (candidate.MutateType == 9)
			Mutate_Points_1(currentBest, candidate);
		//	else  Mutate_Points_1(currentBest, candidate);

	}
	while (candidate.grid.doRandomMove()) {}
	candidate.CalcStats();
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
		GLOBAL_Best_Insert(R, true);
	}
	Recombinators.resize(0);
}

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

				if (GLOBAL_Best_Insert(newCandidate, false))
					++success;
				newCandidate.clear();
			}
		}
		infile.close();
	}
	mutexSAVE.unlock();
	if (success > 0)
		cerr << " Load Recombinators APROX ok/total:" << success << "/" << total << endl;
	if (firstStart && (success > 4) && (success < total - 10))
	{
		cerr << " Cleaning file " << s << endl;
		//Recreate file
		remove(s.c_str());
		for (int d = 0; d<4; ++d)
			for (auto& H : GLOBAL_Best[d])
				H.SaveApprox(s);
	}
}

void Recombinate(LAHC_Node currentBest, LAHC_Node& candidate) {


	candidate.clear();
	int REN = rnd.NextInt(1000);
	if (REN < 250)
	{
		Mutate_Strings_4(currentBest, candidate);
	}
	else if (REN < 400)
	{
		Mutate_Shuffle_Incomplete_End_5(currentBest, candidate);
	}
	else
	{

#ifndef DISABLE_OPTIONALS
		int Ciclo = rnd.NextInt(9);
#else
		int Ciclo = rnd.NextInt(6);
#endif
		switch (Ciclo)
		{
		case 0:Mutate_Completes_0(currentBest, candidate);	break;
		case 1:Mutate_Points_1(currentBest, candidate);	break;
		case 2:Mutate_Lists_2(currentBest, candidate);	break;
		case 3:Mutate_Special_3(currentBest, candidate);	break;
		case 4:Mutate_Cross_4(currentBest, candidate);	break;
		case 5:Mutate_Shuffle_Incomplete_End_5(currentBest, candidate); break;
#ifndef DISABLE_OPTIONALS
		case 6:Mutate_Opt_4(currentBest, candidate); break;
		case 7:Mutate_Genestealers_4(currentBest, candidate); break;
		case 8:ELE_Mutate_Genestealers_4(currentBest, candidate); break;
		default:Mutate_Genestealers_4(currentBest, candidate);	break;
#endif
		}
	}
	//Mutate_Special_3(currentBest, candidate);
	while (candidate.grid.doRandomMove()) {}
	candidate.CalcStats();
}


void ENDGAME_Analyze() {

}

int ALLOW_RECOMBINATE = 0;

//E.K.Burke and Y.Bykov, . "The Late Acceptance Hill-Climbing Heuristic".European Journal of Operational Research. 
void Worker_LAHC(int ID)
{
	ALIGN deque<LAHC_Node> historyBEST;
	bool RecombinateToo = (ID < ALLOW_RECOMBINATE);
	bool isRecombinateActive = false;
	vector<LAHC_Node> LOCAL_Best[MAX_INSERT_NUMBERS];
	const int TAM_HIS = 25;//15;
	double MUTATE_STATS[10] = { 0 };
	double MUTATE_BESTSTATS[10] = { 0 };
	double PROB[10] = { 10.0,15.0,15.0,10.0,30.0,30.0,0.0,0.0,0.0,0.0 };
	int TOTAL_FLASH = 0;
	int SUCCESS_FLASH = 0;
	const int SMALL_LFA = 150;

	int Flashbacks = -1;
	int EXPAND_LFA;
	int FINAL_LFA;
	int CurrLFA;
	//int Lfa;
	long long localTimeLimit = LIMIT_TIME_IMPROVEMENT;
	int realHistoryLength = 0;
	int totalLOCAL = 0;
	if (SIZE_LFA > 0)
	{
		FINAL_LFA = SIZE_LFA;//rnd.NextInt(SIZE_LFA * 90 / 100, SIZE_LFA * 110 / 100);
	}
	else {
		if (LIMIT_TIME_IMPROVEMENT > 125000)
		{
			FINAL_LFA = rnd.NextInt(20000, 40000);
		}
		else if (LIMIT_TIME_IMPROVEMENT > 35000)
		{
			FINAL_LFA = rnd.NextInt(3000, 7000);
		}
		else FINAL_LFA = rnd.NextInt(400, 900);
	}
	EXPAND_LFA = FINAL_LFA + 20 * CAN_INCREASE_LFA;
	cerr << (RecombinateToo ? "RECOMBINATOR " : "LAHC ") << ID << " Lfa size:" << FINAL_LFA << " EXP:" << EXPAND_LFA << endl;
	LAHC_Node lastAccepted;
	LAHC_Node bestStrategy;

	LAHC_Node tmpCandidate;

	vector<double> Fitness;
	LAHC_Node candidate = lastAccepted;

	{ //Produce an initial solution s
		historyBEST.clear();
		Flashbacks = -1;
		lastAccepted.clear();
		Mutate(lastAccepted, candidate, PROB);
		lastAccepted = candidate;
		Fitness.resize(EXPAND_LFA);
		for (int k = 0; k < FINAL_LFA; ++k)
		{
			Fitness[k] = lastAccepted.BestScore; //For all k € {0...Lfa-1} fk:=C(s)
		}
		CurrLFA = SMALL_LFA; //Fast convergence
	}
	bestStrategy = lastAccepted;
	historyBEST.push_back(bestStrategy);
	++realHistoryLength;

	uint64_t v = 0; //First iteration I=0;
	uint64_t I = 0; //First iteration I=0;
	uint64_t previousI = 0;
	bool STALE = false;
	Stopwatch syncClock;
	syncClock.Start(10 * 1000 * 1000);
	Stopwatch notifyClock;
	notifyClock.Start(30 * 1000 * 1000);
	Stopwatch convergenceClock;
	convergenceClock.Start(10 * 1000 * 1000);
	Stopwatch lastImprovement;
	lastImprovement.Start(0);
	Stopwatch timeFromLastReset;
	timeFromLastReset.Start(0);

	long long maxTimeImprovement = 0;
	bool initRecomb = true;
	while (!solved) //Do until a chosen stopping condition
	{
		//Sync global each 10 secs
		if ((ID == 0 || RecombinateToo) && (initRecomb || (((I % 1000) == 0) && syncClock.Timeout())))
		{
			syncClock.Start(10 * 1000 * 1000);
			initRecomb = false;
			//Reload Recombinators
			totalLOCAL = 0;
			int totalLOCAL4 = 0;
			mutexGLOBAL_Best.lock();
			for (int d = 0; d < MAX_INSERT_NUMBERS; ++d)
			{
				LOCAL_Best[d] = GLOBAL_Best[d];
				totalLOCAL += (int)LOCAL_Best[d].size();
				if (d < 4)
					totalLOCAL4 += (int)LOCAL_Best[d].size();
			}
			mutexGLOBAL_Best.unlock();
			if (!isRecombinateActive && RecombinateToo)
			{
				if (ID == 0)	isRecombinateActive = (totalLOCAL > 0);
				else if (ID == 1) {
					if (MAX_INSERT_NUMBERS == 6)
						isRecombinateActive = (totalLOCAL > 0 && LOCAL_Best[5].size() != totalLOCAL);
					else isRecombinateActive = (totalLOCAL4 > CNT_RECOMB);
				}
				else if (ID == 2) {
					if (MAX_INSERT_NUMBERS >= 5)
						isRecombinateActive = (totalLOCAL > 0 && LOCAL_Best[4].size() + LOCAL_Best[5].size() != totalLOCAL);
					else isRecombinateActive = (totalLOCAL4 > CNT_RECOMB);
				}
				else isRecombinateActive = (totalLOCAL4 > CNT_RECOMB);
			}


			//ID 0 is keeper
			if (ID == 0) {
				syncClock.Start(10 * 1000 * 1000);
				int runTotal = 0;
				int runSuccess = 0;
				candidate.clear();
				candidate.grid.doRandomMove();
				candidate.CalcStats();
				tmpCandidate.clear();
				for (int d = 0; d < MAX_INSERT_NUMBERS; ++d)
				{
					for (int i = 0; i < LOCAL_Best[d].size(); ++i)
					{
						if (!LOCAL_Best[d][i].exhausted)
						{
							++runTotal;
							if (Mutate_Exhaustive_6(LOCAL_Best[d][i], candidate, tmpCandidate))
							{
								++runSuccess;
								if (solved)
									return;
							}
							LOCAL_Best[d][i].exhausted = true;
							//Search on GLOBAL to exhaust it
							mutexGLOBAL_Best.lock();
							for (int j = 0; j < GLOBAL_Best[d].size(); ++j)
							{
								if (!GLOBAL_Best[d][j].exhausted && GLOBAL_Best[d][j].remainingNumbers.Equals(LOCAL_Best[d][i].remainingNumbers))
								{
									GLOBAL_Best[d][j].exhausted = true;
									break;
								}
							}
							mutexGLOBAL_Best.unlock();
						}
					}
				}
				auto elapTime = syncClock.EllapsedMilliseconds();
				if (runTotal>0 && elapTime > 200)
					cerr << "      Keeper Exhaustive:" << runSuccess << "/" << runTotal << " Time:" << elapTime << "ms" << endl;
			}

			if (isRecombinateActive && totalLOCAL >0)
			{
				bool completo = false;
				for (int d = MAX_INSERT_NUMBERS - 1; d >= 0; --d)
				{
					if (GLOBAL_Best[d].size() > 0)
					{
						for (int k = 0; k < FINAL_LFA; ++k)
						{
							Fitness[k] = GLOBAL_Best[d][rnd.NextInt((int)GLOBAL_Best[d].size())].BestScore;
						}
						break;
					}
				}
			}
		}

		if (isRecombinateActive)
		{

			int rndList = -1;
			int rndIndex = -1;
			int bigRecomb = rnd.NextInt(1000);
			// 60% P2
			// 80% P3
			// 100% P4
			if (bigRecomb < 200 && LOCAL_Best[0].size()>0)
			{
				rndList = 0;
				rndIndex = rnd.NextInt((int)LOCAL_Best[rndList].size());
			}
			else if (bigRecomb < 400 && LOCAL_Best[1].size()>0) {
				rndList = 1;
				rndIndex = rnd.NextInt((int)LOCAL_Best[rndList].size());
			}
			else if (bigRecomb < 550 && LOCAL_Best[2].size()>0) {
				rndList = 2;
				rndIndex = rnd.NextInt((int)LOCAL_Best[rndList].size());
			}
			else if (MAX_INSERT_NUMBERS >= 4 && bigRecomb < 700 && LOCAL_Best[3].size()>0) {
				rndList = 3;
				rndIndex = rnd.NextInt((int)LOCAL_Best[rndList].size());
			}
			else if (MAX_INSERT_NUMBERS >= 5 && bigRecomb < 850 && LOCAL_Best[4].size()>0) {
				rndList = 4;
				rndIndex = rnd.NextInt((int)LOCAL_Best[rndList].size());
			}
			else  if (MAX_INSERT_NUMBERS == 6 && LOCAL_Best[5].size() > 0) {
				rndList = 5;
				rndIndex = rnd.NextInt((int)LOCAL_Best[rndList].size());
			}

			if (rndList >= 0 && rndIndex >= 0)
			{
#ifdef SHUFFLE_CODE
				if (rnd.NextInt(1000) < 90 && Mutate_Shuffle(LOCAL_Best[rndList][rndIndex], tmpCandidate))
				{
					if (rnd.NextInt(1000) < 200)
						LOCAL_Best[rndList][rndIndex] = tmpCandidate;
					lastAccepted = tmpCandidate;
				}
				else lastAccepted = LOCAL_Best[rndList][rndIndex];
#else
				lastAccepted = LOCAL_Best[rndList][rndIndex];
#endif
			}
			if (rnd.NextInt(1000) < 600)
				Recombinate(lastAccepted, candidate);
			else Mutate(lastAccepted, candidate, PROB);
			if (candidate.IncompleteSets.size()>0 && (candidate.BestScore < lastAccepted.BestScore || rnd.NextInt(1000) < 100))
			{
				bool bb = reArmCandidate(candidate, tmpCandidate);
				if (bb && tmpCandidate.BestScore >= candidate.BestScore)
				{
					swap(tmpCandidate, candidate);
				}
			}
		}
		else
		{
			Mutate(lastAccepted, candidate, PROB); //Construct a candidate solution s*  and Calculate its cost function C(s*)
												   //uint64_t candHash = candidate.usedNumbers.simpleHash();
			if (candidate.IncompleteSets.size() > 0 //We need incomplete sets to reArm
				&& candidate.MutateType != 5 //Not useful on 5
				&& candidate.MutateType != 0 //Not very useful on 0
				&& ((candidate.BestScore < lastAccepted.BestScore && candidate.BestScore < Fitness[v]) //Candidate will be rejected
					|| rnd.NextInt(1000) < 100) //On Accepted candidates just some random chance
				)
			{
#ifdef LOG_REARM									
				++ATTEMPTRearm;
#endif				
				bool bb = reArmCandidate(candidate, tmpCandidate);
				if (bb && tmpCandidate.BestScore > candidate.BestScore

					)
				{
#ifdef LOG_REARM					
					if ((tmpCandidate.BestScore >= lastAccepted.BestScore || tmpCandidate.BestScore >= Fitness[v]))
						++validRearm;
#endif					
					swap(tmpCandidate, candidate);
				}
			}

		}
		MUTATE_STATS[candidate.MutateType] += 1.0;

		if (candidate.grid.totalPoints != 0
			&& (!(candidate.remainingNumbers.Equals(lastAccepted.remainingNumbers) && candidate.BestScore <= lastAccepted.BestScore))
			&& GLOBAL_Best_Insert(candidate, true))
		{
			++newFound;
		}

		//Rescue dying jobs
		if (((I + 1) % 1000 == 0) && convergenceClock.Timeout())
		{

			if (STALE || lastImprovement.EllapsedMilliseconds() > localTimeLimit - 10)
			{
				if ((candidate.BestScore <= bestStrategy.BestScore) && bestStrategy.grid.totalNumbers <= K_EXHAUSTIVE_MIN_NUMBERS)
				{
					bool savingJob = false;
					//Try last 4 from history
					if (historyBEST.size()>1)
						for (int i = (int)historyBEST.size() - 1; i >= 0 && i >= (int)historyBEST.size() - 6; --i) {
							Mutate_Exhaustive_6(historyBEST[i], candidate, tmpCandidate);
							if (candidate.BestScore > bestStrategy.BestScore)
							{
								savingJob = true;
								break;
							}
						}
					if (savingJob)
						cerr << "    Saved ID:" << ID << " from dying! BestScore:" << bestStrategy.BestScore << "<" << candidate.BestScore << " Numbers:" << bestStrategy.grid.totalNumbers << "<" << candidate.grid.totalNumbers << endl;
					//Try lastAccepted
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
						outfile << (isRecombinateActive ? "RECOMBINATOR " : "LAHC ") << "Stats: Lfa:" << FINAL_LFA << " Time From Reset:" << timeFromLastReset.EllapsedMilliseconds() << " maxTimeImprovement:" << maxTimeImprovement / 1000 << "s Last Improvement:" << lastImprovement.EllapsedMilliseconds() / 1000 << "s" << endl;
						outfile.close();
						mutexSAVE.unlock();
						lastAccepted.SaveToFile(filename, maxTimeImprovement, FINAL_LFA, PROB, " Last Accepted -  Flash:" + to_string(SUCCESS_FLASH) + "/" + to_string(TOTAL_FLASH));
						if (isRecombinateActive)
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
				//NNNNNNN: bestStrategy.clear();
				bestStrategy = candidate;
				++realHistoryLength;
				if (historyBEST.size() >= TAM_HIS) { //To reduce size of memory data...
													 //NNNNNNN: historyBEST[0].clear();
					historyBEST.pop_front();
				}
				historyBEST.push_back(LAHC_Node(candidate));
				maxTimeImprovement = max(maxTimeImprovement, lastImprovement.EllapsedMilliseconds());
				lastImprovement.Start(0);
				//GLOBAL_Best_Insert(bestStrategy, true);
			}
			lastAccepted = candidate;
		}

		Fitness[v] = lastAccepted.BestScore; //Insert the current cost into the fitness array fv:=C(s)
		if (!wasAccepted && !isRecombinateActive && historyBEST.size() > 1 && rnd.NextInt(10000) < 3)
		{ //Truquis, reusar más el ENDGAME
			int indexUse = rnd.NextInt(max(0, (int)historyBEST.size() - 5), max(1, (int)historyBEST.size() - 1));
			//lastAccepted = bestStrategy;
			lastAccepted = historyBEST[indexUse];
		}
		++I;//Increment the iteration number I : = I + 1
		++v;//v := I mod Lfa
		if (v >= CurrLFA)
		{
			v = 0;
		}


		if ((I % 1000 == 0) && convergenceClock.Timeout())
		{
			convergenceClock.Start(10 * 1000 * 1000);
			if (COMPUTER_NAME == "numbershift-a" && !isRecombinateActive
				&& stopwatch.EllapsedMilliseconds() > 21 * 60 * 1000 && CNT_RECOMB > 10)
			{
				CNT_RECOMB = 10;
			}
			if (COMPUTER_NAME == "numbershift-a" && !isRecombinateActive)
			{
				if (ID == 0 && totalGLOBAL >= 7)
				{
					RecombinateToo = true;
					CNT_RECOMB = 6;
				}
				if (ID == 1 && totalGLOBAL >= 15)
				{
					RecombinateToo = true;
					CNT_RECOMB = 6;
				}
				mutexGLOBAL_Best.lock();
				if (GLOBAL_Best[0].size() + GLOBAL_Best[1].size() >= 5)
				{
					RecombinateToo = true;
					CNT_RECOMB = 6;
				}
				mutexGLOBAL_Best.unlock();
			}
			if (COMPUTER_NAME == "numbershift-a"
				&& !RecombinateToo
				//&& ID < THREADS/2
				)
			{
				mutexGLOBAL_Best.lock();
				if ((stopwatch.EllapsedMilliseconds() > 21 * 60 * 1000 && totalGLOBAL > CNT_RECOMB)
					|| (GLOBAL_Best[0].size() + GLOBAL_Best[1].size() >= 15))
				{
					RecombinateToo = true;
					//ALLOW_RECOMBINATE = THREADS/2;
					ALLOW_RECOMBINATE = THREADS;
				}
				mutexGLOBAL_Best.unlock();
			}


#ifdef SHUFFLE_CODE			
			//if (I % 60000 == 0)
			//if (lastAccepted.grid.totalNumbers < 30 && secs > LIMIT_TIME_IMPROVEMENT / 1000 * 85 / 100)
			{
				if (Mutate_Shuffle(lastAccepted, candidate)) //Shuffling idea
				{
					swap(lastAccepted, candidate);
				}
			}
#endif		
			double points = 100.0*(double)(MAP.totalPoints - bestStrategy.grid.totalPoints) / (double)MAP.totalPoints;
			double coverage = 100.0*(double)(COUNT_NUMBERS - bestStrategy.grid.totalNumbers)*INV_COUNT_NUMBERS;
			if ((/*points > 85.0 ||*/ coverage > K_LFA_PERC_MIN))
			{
				int TargetLFA = CurrLFA;
				double media = coverage;//(points + coverage)*0.5;
				if (media > K_LFA_PERC_MAX)
				{
					TargetLFA = FINAL_LFA;
				}
				else {
					TargetLFA = SMALL_LFA + max(0, (int)((FINAL_LFA - SMALL_LFA)*(media - K_LFA_PERC_MIN) / max(K_LFA_PERC_MAX - K_LFA_PERC_MIN, 1.0)));
				}

				if (TargetLFA > CurrLFA)
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
				PROB[0] = max(3.0, (100.0*RAT[0] / SUMRAT));
				PROB[1] = max(3.0, (100.0*RAT[1] / SUMRAT));
				PROB[2] = max(3.0, (100.0*RAT[2] / SUMRAT));
				PROB[3] = max(3.0, (100.0*RAT[3] / SUMRAT));
#ifndef DISABLE_OPTIONALS
				if (bestStrategy.grid.totalNumbers >= MAP.totalNumbers * 35 / 100)
					PROB[3] = 0.0;
#endif
				PROB[4] = max(3.0, (100.0*RAT[4] / SUMRAT));
				PROB[5] = max(3.0, (100.0*RAT[5] / SUMRAT));
			}


			for (int k = 0; k < 10; ++k)//Decay, first moves are less relevant now
			{
				MUTATE_STATS[k] = MUTATE_STATS[k] * 0.90;
				MUTATE_BESTSTATS[k] = MUTATE_BESTSTATS[k] * 0.90;
			}
			//if (Recombinators.size() <= CNT_RECOMB) //Not on recombination
			if (!isRecombinateActive )
				if (STALE || lastImprovement.EllapsedMilliseconds() > localTimeLimit)
				{
					/*	if (file_exists("abort.txt"))
					{
					remove("abort.txt");
					abort();
					}*/
					if (STALE)
					{
						cerr << "   ==STALLED";
						STALE = false;
					}

					lastImprovement.Start(0);

					//Increase limit time
					if (CAN_INCREASE_TIME > 0)
					{
						if (localTimeLimit < 80 * 1000)
						{
							localTimeLimit += CAN_INCREASE_TIME;
						}
					}
					if (CAN_INCREASE_LFA > 0)
					{
						FINAL_LFA += CAN_INCREASE_LFA;
						if (FINAL_LFA > EXPAND_LFA)
							FINAL_LFA = EXPAND_LFA;
					}


					if ((historyBEST.size() > 8) && (Flashbacks != 0) && (bestStrategy.grid.totalNumbers <= K_FLASHBACK_NUMBERS || bestStrategy.grid.totalPoints <= K_FLASHBACK_POINTS))
					{
						++TOTAL_FLASH;
						if (Flashbacks == -1)
						{
							Flashbacks = FLASHBACK_COUNT;
						}
						else
							--Flashbacks;
						//						cerr << "   ===>FLASHBACK LAHC ID:" << ID << " F:" << Flashbacks << ".Shrink:" << realHistoryLength << " to ";
						//int ndI1= max(1,(int)historyBEST.size()*(87+Flashbacks)/100 );
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
						//cerr << realHistoryLength + newIndex - historyBEST.size() << " Deque:" << historyBEST.size() << ". Best:" << bestStrategy.BestScore;
						for (int k = 0; k < FINAL_LFA; ++k)
						{
							int NNe = max(0, rnd.NextInt(lowLimit, newIndex));
							Fitness[k] = historyBEST[NNe].BestScore;
						}
						//cerr << " HistoryB:" << historyBEST.back().BestScore << " Flash:" + to_string(SUCCESS_FLASH) + "/" + to_string(TOTAL_FLASH);


						//						for (int JG = newIndex; JG < historyBEST.size(); ++JG) {
						//							historyBEST[JG].clear();
						//						}

						//cerr << "Resizing from " << historyBEST.size() << " to " << newIndex << endl;
						if (historyBEST.size() > newIndex + 1)
							historyBEST.resize(newIndex);

						if (historyBEST.size() > 0)
						{
							//NNNNNNN: lastAccepted.clear();
							lastAccepted = historyBEST.back();
							//Try to Mutate Shuffle here
#ifdef SHUFFLE_CODE							
							{
								if (Mutate_Shuffle(lastAccepted, candidate)) //Shuffling idea
								{
									swap(lastAccepted, candidate);
									//cerr << " Valid Flashback SHUFFLE!";
								}
							}
#endif							
						}
						//cerr << endl;
					}
					else { //Produce an initial solution s
						cerr << "    **RESET LAHC ID:" << ID << " Lfa:" << FINAL_LFA;
						/*
						if (bestStrategy.grid.totalPoints < 7
						|| bestStrategy.grid.totalNumbers < 3
						|| (COUNT_NUMBERS > 300 && bestStrategy.grid.totalNumbers <= 3 && bestStrategy.grid.totalPoints <= 12))
						{
						cerr << "Save as a good approximation" << endl;
						bestStrategy.SaveToFile("LAHC_" + passwordLevel + ".txt", (int)maxTimeImprovement, FINAL_LFA, PROB, " Flash:" + to_string(SUCCESS_FLASH) + "/" + to_string(TOTAL_FLASH));
						}
						*/
						cerr << " Limit Time Improvement:" << localTimeLimit / 1000 << "s" << endl;
						CurrLFA = SMALL_LFA;
						SUCCESS_FLASH = 0;
						TOTAL_FLASH = 0;
						Flashbacks = -1;
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


						Mutate(lastAccepted, candidate, PROB);
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

		if ((I % 1000 == 0) && notifyClock.Timeout())
		{
			notifyClock.Start(30 * 1000 * 1000);
			//I = 0;
			auto tiempo = lastImprovement.EllapsedMilliseconds();
			auto secs = tiempo / 1000;
			auto ms = tiempo % 1000;

			if (bestStrategy.grid.totalPoints < 9 && bestStrategy.grid.totalNumbers < 4)
			{
				cerr << "++++";

			}
			else {
				if (bestStrategy.grid.totalNumbers < 4)
					cerr << "||";
				else cerr << "  ";
				if (bestStrategy.grid.totalPoints < 9)
					cerr << "==";
				else cerr << "  ";

			}

			auto& show = !isRecombinateActive ? bestStrategy : lastAccepted;
			if (isRecombinateActive)
				cerr << "RC:" << ID;
			else cerr << "ID:" << ID;
			cerr << " Pt:" << setw(3) << show.grid.totalPoints;
			cerr << " Nº:" << setw(2) << (show.grid.totalNumbers);
			/*			cerr << " K:" << (int)show.A << "," << fixed << setprecision(2) << show.B << "," << (int)show.C << "," << (int)show.D;
			cerr << "," << setw(2) << (int)show.ExtraC << "," << setw(2) << (int)show.ExtraD;*/
			cerr << " max:" << setw(2) << (maxTimeImprovement / 1000) << "s";
			cerr << (show.IncompleteSets.size() > 0 ? " I" : "  ");
			cerr << " P:[" << setw(2) << (int)PROB[0] << "," << setw(2) << (int)PROB[1] << "," << setw(2) << (int)PROB[2] << "," << setw(2) << (int)PROB[3] << "," << setw(2) << (int)PROB[4] << "," << setw(2) << (int)PROB[5] << "]";
			if (isRecombinateActive)
			{
				cerr << " LCBEST:";
				for (int lc = 0; lc < MAX_INSERT_NUMBERS; ++lc)
					cerr << LOCAL_Best[lc].size() << ",";
			}

			/*
			cerr << " B:[" << setw(2) << MUTATE_BESTSTATS[0] << "," << setw(2) << MUTATE_BESTSTATS[1] << "," << setw(2) << MUTATE_BESTSTATS[2] << "," << setw(2) << MUTATE_BESTSTATS[3] << "," << setw(2) << MUTATE_BESTSTATS[4] << "," << setw(2) << MUTATE_BESTSTATS[5] << "]";
			cerr << " S:[" << setw(2) << (int)MUTATE_STATS[0] << "," << setw(2) << (int)MUTATE_STATS[1] << "," << setw(2) << (int)MUTATE_STATS[2] << "," << setw(2) << (int)MUTATE_STATS[3] << "," << setw(2) << (int)MUTATE_STATS[4] << "," << setw(2) << (int)MUTATE_STATS[5] << "]";
			*/

			//Also X_CELLS and Y_CELLS
			/*		if (bestStrategy.grid.totalNumbers < 6)
			{
			int xr = 0;
			int yr = 0;
			int x2 = 0;
			int y2 = 0;
			for (int x = 0; x < W; ++x)
			{
			if (bestStrategy.grid.X_CELLS[x] > 0) xr++;
			if (bestStrategy.grid.X_CELLS[x] == 2) x2++;
			}
			for (int y = 0; y < H; ++y)
			{
			if (bestStrategy.grid.Y_CELLS[y] > 0) yr++;
			if (bestStrategy.grid.Y_CELLS[y] == 2) y2++;
			}
			cerr << " XYrows:" << xr << "," << yr << " XY2:" << x2 << "," << y2;
			}*/
			if (show.grid.totalNumbers < 4 || isRecombinateActive)
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
			if (isRecombinateActive)
				cerr << " Recombinators:" << (totalLOCAL) << " GLOBAL:" << totalGLOBAL;


			cerr << " I:" << (I - previousI) << endl;
			previousI = I;
		}
		//++Iterations;
	}
}

void initialExhaustive() {
#if 1
	swapDone = 0;
	//for (auto& G : GLOBAL_Best)
	for (int d = 0; d<4; ++d)
		for (int i = 0; i < GLOBAL_Best[d].size(); ++i)
		{
			auto G = GLOBAL_Best[d][i];
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
			GLOBAL_Best[d][i].exhausted = true;
			if (IS_SOLVED(candidate, "EXHAUSTIVE"))
				return;
			else if (candidate.Movelists.size()>0 && candidate.grid.totalNumbers < G.grid.totalNumbers) {
				cerr << "IMPROVED!!!!!:" << G.grid.totalNumbers << " <->" << GLOBAL_Best[d][i].grid.totalNumbers << " -> " << candidate.grid.totalNumbers << endl;
				++swapDone;
			}
		}
	if (swapDone > 0 /*&& (stopwatch.EllapsedMilliseconds()  > 88*60*1000)*/) //90 minutes or more
	{
		cerr << "Recreating APROX file due to " << swapDone << " swaps..." << endl;
		string ss = "APROX_" + passwordLevel + ".txt";
		remove(ss.c_str());
		for (int d = 0; d<4; ++d)
			for (auto&B : GLOBAL_Best[d])
			{
				B.SaveApprox(ss);
			}
		swapDone = 0;
	}
#endif	
}
void ParallelWork() {
	if (PROGRAM_NAME != "./test")
	{
		LoadRecombinators_APROX("APROX_" + passwordLevel + ".txt", true);
		LoadRecombinators_LAHC("LAHC_" + passwordLevel + ".txt"); //Deprecated
	}
	initialExhaustive();
	vector<thread> threads(max(1, THREADS));

	if (!solved)
		for (int i = 0; i < max(1, THREADS); i++)
		{
			cerr << "Creating Thread " << i << " Type: " << "LAHC" << endl;
			//TODO: Create a vector of thread(Worker_LAHC, i);
		}



	//threads[threads.size() - 1] = thread(Worker_DLX, threads.size() - 1);
	int C = 0;
	while (!solved)
	{
		if (++C > 5 * 150)
		{
			if (ALLOW_RECOMBINATE > 0 && (PROGRAM_NAME != "./test"))
			{
				string command = "rm toImport.txt >/dev/null 2>&1;mv EXTERN_" + passwordLevel + ".txt toImport.txt >/dev/null 2>&1";
				system(command.c_str());
				LoadRecombinators_APROX("toImport.txt", true); //Reload APROX
				system("rm toImport.txt >/dev/null 2>&1");
			}

			if (swapDone > 0 && (PROGRAM_NAME != "./test")) //90 minutes or more
			{
				cerr << "Recreating APROX file due to " << swapDone << " swaps...";
				string ss = "SAFE_" + passwordLevel + ".txt";
				mutexAPROX.lock();
				remove(ss.c_str());
				for (int d = 0; d<4; ++d)
					for (auto&B : GLOBAL_Best[d])
					{
						B.SaveApprox(ss);
					}
				mutexAPROX.unlock();
				swapDone = 0;
				string commando = "mv SAFE_" + passwordLevel + ".txt APROX_" + passwordLevel + ".txt";
				system(commando.c_str());
				cerr << "done.All OK" << endl;
			}


			cerr << "LEVEL:" << level << " P:" << setw(4) << MAP.totalPoints << " Numbers:" << setw(4) << MAP.totalNumbers;
			cerr << " SIM/Sec:" << SimCount / (stopwatch.EllapsedMilliseconds() / 1000);


			if (newFound == 0)
				performance = stopwatch.EllapsedMilliseconds() / 1000 * THREADS;
			else performance = stopwatch.EllapsedMilliseconds() / 1000 * THREADS / newFound;
			cerr << " Performance:" << performance << "s/found(" << newFound << ")";
			cerr << " JACKPOTS:" << JACKPOT_OK << "/" << JACKPOT_TOTAL << " ";
			cerr << "SUPERCOMBINATOR:" << SUPERCOMBINATOR << " ";
			cerr << "MAX_PERMUTATIONS:" << MAX_PERMUTATIONS;

#ifdef LOG_REARM								
			uint64_t ATR = ATTEMPTRearm;
			if (ATR == 0)
				ATTEMPTRearm = 1;
			cerr << " Valid Rearm:" << validRearm - prevRearm << "/" << validRearm << " T:" << ATTEMPTRearm << " : " << (validRearm == 0 ? 0 : validRearm * 100 / ATTEMPTRearm) << "%";
			uint64_t tmp = validRearm;
			prevRearm = tmp;
#endif			
#ifdef SHUFFLE_CODE			
			//	cerr << " Valid Shuffles:" << correctShuffles << "/" << countShuffles << " : " << (countShuffles == 0 ? 0 : correctShuffles * 100 / countShuffles) << "%";
#endif			
			cerr << " GLOBAL_Best:" << GLOBAL_Best[0].size() << "," << GLOBAL_Best[1].size() << "," << GLOBAL_Best[2].size() << "," << GLOBAL_Best[3].size() << "=" << totalGLOBAL;

			/*			cerr << " BestScore:";
			if (totalGLOBAL == 0)
			{
			cerr << "0,0.00,0,0,0,0";
			}
			else {
			LAHC_Node* show2 = nullptr;
			for (int d = 0; d < 4; ++d)
			if (GLOBAL_Best[d].size() > 0)
			{
			show2 = &GLOBAL_Best[d][0];
			}
			if (show2 != nullptr)
			{
			cerr << " K:" << (int)show2->A << "," << fixed << setprecision(2) << show2->B << "," << (int)show2->C << "," << (int)show2->D;
			cerr << "," << setw(2) << (int)show2->ExtraC << "," << setw(2) << (int)show2->ExtraD << " = " << (int)show2->BestScore;
			}

			}*/
			cerr << endl;
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
#ifdef USE_CACHE_MOVES
	//CACHE_MOVES
	Grid searchGRID = MAP;
	vector<Move> moves;
	moves.reserve(4);
	for (int srcID = 0; srcID < COUNT_NUMBERS; ++srcID)
	{
		//Not zero
		for (int srcVal = 1; srcVal < MAX_W + 1; ++srcVal)
		{
			moves.resize(0);
			//ExplodeMoves
			searchGRID.Val[srcID] = srcVal;
			searchGRID.ExplodeMoves(srcID, moves);
			if (moves.size() > 0)
			{
				do_shuffle<vector<Move>>(moves);
				CACHE_MOVES[srcID][srcVal] = new vector<Move>();
				for (auto&m : moves)
					if (m.sign == 0)
					{
						CACHE_MOVES[srcID][srcVal]->push_back(m);
					}
			}
		}
		searchGRID.Val[srcID] = MAP.Val[srcID];
	}
#endif
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
	{
		PROGRAM_NAME = argv[0];
		system("hostname > hostname.txt");
		ifstream f("hostname.txt");
		if (f.good()) { std::getline(f, COMPUTER_NAME); f.close(); }

		if (PROGRAM_NAME == "./test")
		{
			ifstream fn("randomKey.txt");
			uint64_t randomKey = 0;
			if (fn.good()) { fn >> randomKey; fn.close(); }
			cerr << "Random KEY es:" << randomKey << endl;
			rnd = Random(randomKey);

		}
	}


	newFound = 0;
	totalGLOBAL = 0;
#ifdef LOG_REARM									
	prevRearm = 0;
	validRearm = 0;
	ATTEMPTRearm = 0;
#endif
#ifdef SHUFFLE_CODE
	correctShuffles = 0;
	countShuffles = 0;
#endif	
	JACKPOT_TOTAL = 0; JACKPOT_OK = 0; SUPERCOMBINATOR = 0; MAX_PERMUTATIONS = 0;
	stacksize();
	stopwatch.Start(0);

	USE_SA = isPartOf("SA", argv[0]);
	cerr << argv[0] << " USE_SA:" << USE_SA << endl;

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
		ALLOW_RECOMBINATE = atoi(argv[10]);
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
		cerr << "Width overflow" << endl;
		assert(W <= MAX_W);
		assert(H <= MAX_H);
	}
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
	INV_COUNT_NUMBERS = 1.0 / (double)COUNT_NUMBERS;
	INV_GROUPS = 1.0 / (double)GROUPS;
	INV_POINTS = 1.0 / (double)MAP.totalPoints;
	INV_SQUAREPOINTS = 1.0 / (double)MAP.squaredPoints;


	INV_MINIPUNTO = 1.0 / ((double)MAP.totalPoints*(double)GROUPS);

	INV_RESX = 0.5 / (10.0*(double)W);
	INV_RESY = 0.5 / (10.0*(double)H);


	cerr << "LEVEL:" << level << " NUMBERS:" << COUNT_NUMBERS << " SPAWNS:" << spawns << " GROUPS:" << GROUPS << " W,H:" << W << "," << H << endl;

	PopulateCacheMoves();
	CreateNeighbours();
	//testOrphans();
	ParallelWork();

	return 0;
}