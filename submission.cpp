// Highest 125728.708
#define PROGNAME ""
#define hash ___hash

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <ctime>
#include <thread>
#include <cassert>
#include <new>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include "vm/boc.h"
#include "td/utils/base64.h"

#undef hash
using namespace std;
namespace utils {
    td::BufferSlice change_mode_to_compress(td::BufferSlice &data, int mode = 4) {
        td::Ref <vm::Cell> root = vm::std_boc_deserialize(data).move_as_ok();
        return vm::std_boc_serialize(root, mode).move_as_ok();
    }

    td::BufferSlice change_mode_back_to_decompress(td::BufferSlice &data) {
        auto root = vm::std_boc_deserialize(data).move_as_ok();
        return vm::std_boc_serialize(root, 31).move_as_ok();
    }

    td::BufferSlice
    compress_by_function(td::BufferSlice &data, function<vector<uint8_t>(const vector<uint8_t> &)> compress) {
        vector<uint8_t> input(data.size());
        memcpy(input.data(), data.data(), data.size());
        vector<uint8_t> output = compress(input);
        return td::BufferSlice(reinterpret_cast<char *>(output.data()), output.size());
    }

    td::BufferSlice
    decompress_by_function(td::BufferSlice &data, function<vector<uint8_t>(const vector<uint8_t> &)> decompress) {
        vector<uint8_t> compressed(data.size());
        memcpy(compressed.data(), data.data(), data.size());
        vector<uint8_t> output = decompress(compressed);
        return td::BufferSlice(reinterpret_cast<char *>(output.data()), output.size());
    }

    std::vector<unsigned char> append_iterations(const td::BufferSlice &data, int iterations) {
        std::vector<unsigned char> result;
        result.insert(result.end(), data.data(), data.data() + data.size());
        result.push_back(static_cast<unsigned char>(iterations));
        return result;
    }
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <cassert>
#define NDEBUG  // remove for debugging (turns on Array bound checks)

#include <assert.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#endif
#ifdef WINDOWS
#include <windows.h>
#endif
#ifndef DEFAULT_OPTION
#define DEFAULT_OPTION 5
#endif
typedef unsigned char U8;
typedef unsigned short U16;
typedef unsigned int U32;
#ifndef WINDOWS

inline int min(int a, int b) { return a < b ? a : b; }

inline int max(int a, int b) { return a < b ? b : a; }

#endif



#define MEM (0x10000 << 4)


// Used memory and time
class Stats {
    int memused, maxmem;
    clock_t start_time;
public:
    void alloc(int n) {
        memused += n;
        if (memused > maxmem) maxmem = memused;
    }
    Stats():memused(0), maxmem(0) {
        start_time = clock();
    }
    void print() const {
        // printf("Time %1.2f sec, used %d bytes of memory\n",
            //    double(clock() - start_time) / CLOCKS_PER_SEC, maxmem);
    }
} stats;

// Array
template<class T,int ALIGN=0> class Array {
private:
    int n, reserved;
    char *ptr;
    T *data;
    void create(int i);
public:
    explicit Array(int i = 0) { create(i); }
    ~Array();
    T&operator[](int i) { return data[i]; }
    const T&operator[](int i) const { return data[i]; }
    int size() const { return n; }
private:
    Array(const Array&);
    Array&operator=(const Array&);
};
template<class T,int ALIGN> void Array<T,ALIGN>::create(int i) {
    n = reserved = i;
    if (i <= 0) {
        data = 0, ptr = 0;
        return;
    }
    const int sz = ALIGN + n * sizeof(T);
    stats.alloc(sz);
    ptr = (char*)calloc(sz, 1);
    if (!ptr) throw "Out of memory";
    data = (ALIGN ? (T*)(ptr + ALIGN - (((uintptr_t)ptr) & (ALIGN - 1))):(T*)ptr);
}
template<class T,int ALIGN> Array<T,ALIGN>::~Array() {
    stats.alloc(-ALIGN - n * sizeof(T));
    free(ptr);
}

// Random generator
static class Random {
    Array<U32> table;
    int i;
public:
    Random():table(64) {
        table[0] = 123456789;
        table[1] = 987654321;
        for (int j = 0; j < 62; j++) {
            table[j + 2] = table[j + 1] * 11 + table[j] * 23 / 16;
        }
        i = 0;
    }
    U32 operator()() {
        return ++i, table[i & 63] = table[(i - 24) & 63] ^ table[(i - 55) & 63];
    }
} rnd;

// Buffer - array of n last bytes
static int pos;
class Buf {
    Array<U8> b;
public:
    Buf(int i = 0):b(i) {}
    U8& operator[](int i) {
        return b[i & (b.size() - 1)];
    }
    int operator()(int i) const {
        return b[(pos - i) & (b.size() - 1)];
    }
    int size() const{
        return b.size();
    }
};

// Global variables
static FILE* infile, *outfile;
static int y = 0, c0 = 1, bpos = 0;
static U32 c4 = 0;
static Buf buf(MEM * 8);
long len4;
// TOO_SMALL 80'000, 90'000 in original, 94'000 (108'000) breaks 2s. So 80'000 is safe.
#define TOO_SMALL (len4 < 80'000)
// SMALL 130'000, 150'000 in original. MAX 1.6s.
#define SMALL (len4 < 130'000)
// MID 300'000, 400'000 in original. MAX 1.7s.
// TODO: try changin the upper bound, I didn't try, try 260'000 if u wanna be most safe
#define MID (len4 < 300'000)

// Logarithm
static class Ilog {
    Array<U8> t;
public:
    int operator()(U16 x) const { return t[x]; }
    Ilog():t(65536) {
        U32 x = 14155776;
        for (int i = 2; i < 65536; ++i) {
            x += 774541002 / (i * 2 - 1);
            t[i] = x >> 24;
        }
    }
} ilog;

// Precomputed table for pt(i) = 16384 / (i + i + 3)
static class Ptable {
    Array<int> t;
public:
    int operator()(U16 x) const { return t[x]; }
    Ptable():t(1024) {
        for (int i = 0; i < 1024; ++i) {
            t[i] = 16384 / (i + i + 3);
        }
    }
} pt;

// Hash
static U32 hash(U32 a, U32 b, U32 c = 0xffffffff) {
    U32 h = a * 200002979u + b * 30005491u + c * 50004239u + 4114959990u;
    return h ^ h >> 9 ^ a >> 2 ^ b >> 3 ^ c >> 4 ^ 0x4000000;
}


// State table
static const U8 State_table[256][2] = {
        {1,2},{3,5},{4,6},{7,10},{8,12},{9,13},{11,14},{15,19},{16,23},{17,24},
        {18,25},{20,27},{21,28},{22,29},{26,30},{31,33},{32,35},{32,35},{32,35},
        {32,35},{34,37},{34,37},{34,37},{34,37},{34,37},{34,37},{36,39},{36,39},
        {36,39},{36,39},{38,40},{41,43},{42,45},{42,45},{44,47},{44,47},{46,49},
        {46,49},{48,51},{48,51},{50,52},{53,43},{54,57},{54,57},{56,59},{56,59},
        {58,61},{58,61},{60,63},{60,63},{62,65},{62,65},{50,66},{67,55},{68,57},
        {68,57},{70,73},{70,73},{72,75},{72,75},{74,77},{74,77},{76,79},{76,79},
        {62,81},{62,81},{64,82},{83,69},{84,71},{84,71},{86,73},{86,73},{44,59},
        {44,59},{58,61},{58,61},{60,49},{60,49},{76,89},{76,89},{78,91},{78,91},
        {80,92},{93,69},{94,87},{94,87},{96,45},{96,45},{48,99},{48,99},{88,101},
        {88,101},{80,102},{103,69},{104,87},{104,87},{106,57},{106,57},{62,109},
        {62,109},{88,111},{88,111},{80,112},{113,85},{114,87},{114,87},{116,57},
        {116,57},{62,119},{62,119},{88,121},{88,121},{90,122},{123,85},{124,97},
        {124,97},{126,57},{126,57},{62,129},{62,129},{98,131},{98,131},{90,132},
        {133,85},{134,97},{134,97},{136,57},{136,57},{62,139},{62,139},{98,141},
        {98,141},{90,142},{143,95},{144,97},{144,97},{68,57},{68,57},{62,81},{62,81},
        {98,147},{98,147},{100,148},{149,95},{150,107},{150,107},{108,151},{108,151},
        {100,152},{153,95},{154,107},{108,155},{100,156},{157,95},{158,107},{108,159},
        {100,160},{161,105},{162,107},{108,163},{110,164},{165,105},{166,117},
        {118,167},{110,168},{169,105},{170,117},{118,171},{110,172},{173,105},
        {174,117},{118,175},{110,176},{177,105},{178,117},{118,179},{110,180},
        {181,115},{182,117},{118,183},{120,184},{185,115},{186,127},{128,187},
        {120,188},{189,115},{190,127},{128,191},{120,192},{193,115},{194,127},
        {128,195},{120,196},{197,115},{198,127},{128,199},{120,200},{201,115},
        {202,127},{128,203},{120,204},{205,115},{206,127},{128,207},{120,208},
        {209,125},{210,127},{128,211},{130,212},{213,125},{214,137},{138,215},
        {130,216},{217,125},{218,137},{138,219},{130,220},{221,125},{222,137},
        {138,223},{130,224},{225,125},{226,137},{138,227},{130,228},{229,125},
        {230,137},{138,231},{130,232},{233,125},{234,137},{138,235},{130,236},
        {237,125},{238,137},{138,239},{130,240},{241,125},{242,137},{138,243},
        {130,244},{245,135},{246,137},{138,247},{140,248},{249,135},{250,69},{80,251},
        {140,252},{249,135},{250,69},{80,251},{140,252}};

static int squash(int d) {
    static const int t[33] = {
            1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,1546,2047,2549,2994,3348,
            3607,3785,3901,3975,4022,4050,4068,4079,4085,4089,4092,4093,4094};
    if (d > 2047) return 4095;
    if (d < -2047) return 0;
    int w = d & 127;
    d = (d >> 7) + 16;
    return (t[d] * (128 - w) + t[(d + 1)] * w + 64) >> 7;
}

class Stretch {
    Array<short> t;
public:
    int operator()(int x) const { return t[x]; }
    Stretch():t(4096) {
        int j = 0;
        for (int x = -2047; x <= 2047; ++x) {
            int i = squash(x);
            while (j <= i) t[j++] = x;
        }
        t[4095] = 2047;
    }
} stretch;

// #undef __SSE2__
#if defined(__SSE2__)
#include <emmintrin.h>
static int dot_product (const short* const t, const short* const w, int n) {
    __m128i sum = _mm_setzero_si128 ();
    while ((n -= 8) >= 0) {
        __m128i tmp = _mm_madd_epi16 (*(__m128i *) &t[n], *(__m128i *) &w[n]);
        tmp = _mm_srai_epi32 (tmp, 8);
        sum = _mm_add_epi32 (sum, tmp);
    }
    sum = _mm_add_epi32 (sum, _mm_srli_si128 (sum, 8));
    sum = _mm_add_epi32 (sum, _mm_srli_si128 (sum, 4));
    return _mm_cvtsi128_si32 (sum);
}

static void train (const short* const t, short* const w, int n, const int e) {
    if (e) {
        const __m128i one = _mm_set1_epi16 (1);
        const __m128i err = _mm_set1_epi16 (short(e));
        while ((n -= 8) >= 0) {
            __m128i tmp = _mm_adds_epi16 (*(__m128i *) &t[n], *(__m128i *) &t[n]);
            tmp = _mm_mulhi_epi16 (tmp, err);
            tmp = _mm_adds_epi16 (tmp, one);
            tmp = _mm_srai_epi16 (tmp, 1);
            tmp = _mm_adds_epi16 (tmp, *(__m128i *) &w[n]);
            *(__m128i *) &w[n] = tmp;
        }
    }
}
#else
#endif

// Mixer - combines models using neural networkss
class Mixer {
    const int N, M, S;
    Array<short,16> tx, wx;
    Array<int> cxt, pr;
    int ncxt, base, nx;
    Mixer *mp;
public:
    Mixer(int n, int m, int s = 1, int w = 0):
            N((n + 7) & -8), M(m), S(s), tx(N), wx(N * M),
            cxt(S), pr(S), ncxt(0), base(0), nx(0), mp(0) {
        for (int i = 0; i < S; ++i) {
            pr[i] = 2048;
        }
        for (int j = 0; j < N * M; ++j) {
            wx[j] = w;
        }
        if (S > 1) mp = new Mixer(S, 1, 1);
    }
    void update() {
        for (int i = 0; i < ncxt; ++i) {
            int err = ((y << 12) - pr[i]) * 6;
            train(&tx[0], &wx[cxt[i] * N], nx, err);
        }
        nx = base = ncxt = 0;
    }
    void add(int x) { tx[nx++] = x; }
    void set(int cx, int range) {
        cxt[ncxt++] = base + cx;
        base += range;
    }
    int p() {
        while (nx & 7) tx[nx++] = 0;
        if (mp) {
            mp->update();
            for (int i = 0; i < ncxt; ++i) {
                pr[i] = squash(dot_product(&tx[0], &wx[cxt[i] * N], nx) >> 4);
                mp->add(stretch(pr[i]));
            }
            mp->set(0, 1);
            return mp->p();
        } else {
            return pr[0] = squash(dot_product(&tx[0], &wx[0], nx) >> 8);
        }
    }
    ~Mixer() { delete mp; }
};

// APM - stationary map combining a context and an input probability.
class APM {
    int index;
    const int N;
    Array<U16> t;
public:
    APM(int n):index(0), N(n), t(n * 33) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < 33; ++j) {
                t[i * 33 + j] = i == 0 ? squash((j - 16) * 128) * 16 : t[j];
            }
        }
    }
    int p(int pr = 2048, int cxt = 0, int rate = 6) {
        pr = stretch(pr);
        int g = (y << 16) + (y << rate) - y - y;
        t[index] += (g - t[index]) >> rate;
        t[index + 1] += (g - t[index + 1]) >> rate;
        const int w = pr & 127;
        index = ((pr + 2048) >> 7) + cxt * 33;
        return(t[index] * (128 - w) + t[index + 1] * w) >> 11;
    }
};

//  StateMap - maps a context to a probability
class StateMap {
protected:
    const int N;
    int cxt;
    Array<U32> t;
public:
    StateMap(int n = 256):N(n), cxt(0), t(n) {
        for (int i = 0; i < N; ++i) {
            t[i] = 1 << 31;
        }
    }
    int p(int cx) {
        U32 *p = &t[cxt], p0 = p[0];
        int n = p0 & 1023, pr = p0 >> 10;
        if (n < 1023) ++p0;
        else p0 = (p0 & 0xfffffc00) | 1023;
        p0 += (((y << 22) - pr) >> 3) * pt(n) & 0xfffffc00;
        p[0] = p0;
        return t[cxt = cx] >> 20;
    }
};

class SmallStationaryContextMap {
    Array<U16> t;
    int cxt;
    U16 *cp;
public:
    SmallStationaryContextMap(int m):t(m / 2), cxt(0) {
        for (int i = 0; i < t.size(); ++i) {
            t[i] = 32768;
        }
        cp = &t[0];
    }
    void set(U32 cx) { cxt = cx * 256 & (t.size() - 256); }
    void mix(Mixer&m, int rate = 7) {
        *cp += ((y << 16) - (*cp) + (1 << (rate - 1))) >> rate;
        cp = &t[cxt + c0];
        m.add(stretch((*cp) >> 4));
    }
};

class ContextMap {
    const int C;
    class E {
        U16 chk[7];
        U8 last;
    public:
        U8 bh[7][7];
        U8* get(U16 chk);
    };
    Array<E, 64> t;
    Array<U8*> cp, cp0, runp;
    Array<U32> cxt;
    StateMap *sm;
    int cn;
    void update(U32 cx, int c);
public:
    ContextMap(int m, int c = 1);
    ~ContextMap();
    void set(U32 cx, int next = -1);
    int mix(Mixer&m);
};
inline U8*ContextMap::E::get(U16 ch) {
    if (chk[last & 15] == ch) return &bh[last & 15][0];
    int b = 0xffff, bi = 0;
    for (int i = 0; i < 7; ++i) {
        if (chk[i] == ch) return last = last << 4 | i, (U8*) &bh[i][0];
        int pri = bh[i][0];
        if (pri < b && (last & 15) != i && last >> 4 != i) b = pri, bi = i;
    }
    return last = 0xf0 | bi, chk[bi] = ch, (U8*)memset(&bh[bi][0], 0, 7);
}
ContextMap::ContextMap(int m, int c):C(c), t(m>>6), cp(c), cp0(c), runp(c),
                                     cxt(c), cn(0) {
    sm = new StateMap[C];
    for (int i = 0; i < C; ++i) {
        cp0[i] = cp[i] = &t[0].bh[0][0];
        runp[i] = cp[i] + 3;
    }
}
ContextMap::~ContextMap() {
    delete[] sm;
}
inline void ContextMap::set(U32 cx, int next) {
    int i = cn++;
    i &= next;
    cx = cx * 987654323 + i;
    cx = cx << 16 | cx >> 16;
    cxt[i] = cx * 123456791 + i;
}
int ContextMap::mix(Mixer &m) {
    int result = 0;
    for (int i = 0; i < cn; ++i) {
        if (cp[i]) {
            int ns = State_table[*cp[i]][y];
            if (ns >= 204 && rnd() << ((452 - ns) >> 3)) ns -= 4;
            *cp[i] = ns;
        }
        if (bpos > 1 && runp[i][0] == 0) {
            cp[i] = 0;
        } else {
            switch(bpos) {
                case 1: case 3: case 6: cp[i] = cp0[i] + 1 + (c0 & 1); break;
                case 4: case 7: cp[i] = cp0[i] + 3 + (c0 & 3); break;
                case 2: case 5: cp0[i] = cp[i] = t[(cxt[i] + c0) & (t.size() - 1)].get(cxt[i] >> 16); break;
                default: {
                    cp0[i] = cp[i] = t[(cxt[i] + c0) & (t.size() - 1)].get(cxt[i] >> 16);
                    if (cp0[i][3] == 2) {
                        const int c = cp0[i][4] + 256;
                        U8 *p = t[(cxt[i] + (c >> 6)) & (t.size() - 1)].get(cxt[i] >> 16);
                        p[0] = 1 + ((c >> 5) & 1);
                        p[1 + ((c >> 5) & 1)] = 1 + ((c >> 4) & 1);
                        p[3 + ((c >> 4) & 3)] = 1 + ((c >> 3) & 1);
                        p = t[(cxt[i] + (c >> 3)) & (t.size() - 1)].get(cxt[i] >> 16);
                        p[0] = 1 + ((c >> 2) & 1);
                        p[1 + ((c >> 2) & 1)] = 1 + ((c >> 1) & 1);
                        p[3 + ((c >> 1) & 3)] = 1 + (c & 1);
                        cp0[i][6] = 0;
                    }
                    int c1 = buf(1);
                    if (runp[i][0] == 0) {
                        runp[i][0] = 2, runp[i][1] = c1;
                    } else if (runp[i][1] != c1) {
                        runp[i][0] = 3, runp[i][1] = c1;
                    } else if (runp[i][0] < 254) {
                        runp[i][0] += 2;
                    }
                    runp[i] = cp0[i] + 3;
                } break;
            }
        }
        if ((runp[i][1] + 256) >> (8 - bpos) == c0) {
            int rc = runp[i][0];
            int b = (runp[i][1] >> (7 - bpos) & 1) * 2 - 1;
            int c = ilog(rc+1) << (2 + (~rc & 1));
            m.add((b * c << 1) + 75);//40,724 100,763 75,781
        } else {
            m.add(0);
        }
        int p;
        if (cp[i]) {
            result += (*cp[i] > 0);
            p = sm[i].p(*cp[i]);
        } else {
            p = sm[i].p(0);
        }
        m.add(stretch(p));
    }
    if (bpos == 7) cn=0;
    return result;
}

// Match submodel
static int matchModel(Mixer& m) {
    const int MAXLEN = 0xfffe;
    static Array<int> t(MEM);
    static int h = 0, ptr = 0, len = 0, result = 0, posnl = 0;
    static SmallStationaryContextMap scm1(0x20000), scm2(0x20000);
    if (!bpos) {
        h = (h * 997 * 8 + buf(1) + 1) & (t.size() - 1);
        if (len) {
            ++len, ++ptr;
        } else {
            ptr = t[h];
            if (ptr && pos - ptr < buf.size()) {
                while (buf(len + 1) == buf[ptr - len - 1] && len < MAXLEN) ++len;
            }
        }
        t[h] = pos;
        result = len;
        if (SMALL) scm1.set(pos);
        if (buf(1) == 0xff || buf(1) == '\r' || buf(1) == '\n') posnl = pos;
        scm2.set(min(pos - posnl, 255));
    }
    if (len) {
        if (buf(1) == buf[ptr - 1] && c0 == (buf[ptr] + 256) >> (8 - bpos)) {
            if (len > MAXLEN) len = MAXLEN;
            if (buf[ptr] >> (7 - bpos) & 1) {
                m.add(ilog(len) << 2);
                m.add(min(len, 32) << 6);
            } else {
                m.add(-(ilog(len) << 2));
                m.add(-(min(len, 32) << 6));
            }
        } else {
            len=0;
            m.add(0);
            m.add(0);
        }
    } else {
        m.add(0);
        m.add(0);
    }
    if (SMALL) scm1.mix(m);
    scm2.mix(m);
    return result;
}



void sparseModel(Mixer &m) {
    static ContextMap cm(MEM * 4, 8), scm(MEM, 8);
    if (bpos == 0) {
        cm.set(c4 & 0x00ff00ff);
        cm.set(c4 & 0xff0000ff);
        cm.set(buf(1) | buf(5) << 8);
        cm.set(buf(1) | buf(6) << 8);
        cm.set(c4 & 0x00ffff00);
        cm.set(c4 & 0xff00ff00);
        cm.set(buf(3) | buf(6) << 8);
        cm.set(buf(4) | buf(8) << 8);
        for (int i = 0; i < 8; ++i)
            scm.set(buf(i + 1));
    }
    cm.mix(m);
    scm.mix(m);
}

void sparseModel12(Mixer &m) {
    static ContextMap cm(MEM, 4), scm(MEM, 8); // Reduced number of contexts
    if (bpos == 0) {
        cm.set(c4 & 0x00ff00ff);
        cm.set(c4 & 0xff0000ff);
        // cm.set(buf(1) | buf(5) << 8);
        //cm.set(buf(1) | buf(6) << 8);
        cm.set(c4 & 0x00ffff00);
        cm.set(c4 & 0xff00ff00);
        //cm.set(c4 & 0x0000ffff);
        //cm.set(c4 & 0x00ffff00);
        //cm.set(buf(3) | buf(6) << 8);
        //cm.set(buf(4) | buf(8) << 8);
        for (int i = 0; i < 8; ++i)
            scm.set(buf(i + 1));
    }
    cm.mix(m);
    scm.mix(m);
}

void sparseModel6(Mixer &m) {
    static ContextMap cm(MEM, 4), scm(MEM, 4); // Reduced number of contexts
    if (bpos == 0) {
        cm.set(c4 & 0x00ff00ff);
        cm.set(c4 & 0xff00ff00);
        for (int i = 0; i < 4; ++i)
            scm.set(buf(i + 1));
    }
    cm.mix(m);
    scm.mix(m);
}

void sparseModel4(Mixer& m){
    static ContextMap cm(MEM, 4);
    if (bpos == 0) {
        cm.set(c4 & 0x000000ff);
        cm.set(buf(1) | buf(5) << 8);
        cm.set(c4 & 0xff0000ff);
        cm.set(buf(3));
    }
    cm.mix(m);
}

namespace old {
    inline int llog(U32 x) {
        if (x >= 0x1000000)
            return 256 + ilog(x >> 16);
        else if (x >= 0x10000)
            return 128 + ilog(x >> 8);
        else
            return ilog(x);
    }

    void wordModel(Mixer &m) {
        static U32 word0 = 0, word1 = 0, word2 = 0, word3 = 0, word4 = 0;
        static U32 text0 = 0;
        static ContextMap cm(MEM * 32, 14);
        static Array<int> wpos(MEM);
        static int nl1 = -3, nl = -2;
        if (bpos == 0) {
            int c = c4 & 255;
            if (c >= 'A' && c <= 'Z')
                c += 'a' - 'A';
            if (c >= 'a' && c <= 'z' || c >= 128) {
                word0 = word0 * 263 * 4 + c;
                text0 = text0 * 997 * 16 + c;
            } else if (word0) {
                word4 = word3 * 11;
                word3 = word2 * 7;
                word2 = word1 * 5;
                word1 = word0 * 3;
                word0 = 0;
            }
            if (c == 10) nl1 = nl, nl = pos - 1;
            int col = min(255, pos - nl), above = buf[nl1 + col];
            U32 h = word0 * 271 + buf(1);
            cm.set(h);
            cm.set(word0);
            cm.set(h + word1);
            cm.set(word0 + word1 * 17);
            cm.set(h + word2);
            cm.set(h + word1 + word2);
            cm.set(h + word3);
            cm.set(h + word4);
            cm.set(text0 & 0xffff);
            cm.set(text0 & 0xfffff);
            cm.set(col << 8 | above);
            cm.set(col << 8 | buf(1));
            cm.set(buf(1) << 8 | above);
            cm.set(col);
        }
        cm.mix(m);
    }

    void recordModel(Mixer &m) {
        static Array<int> cpos1(256), cpos2(256), cpos3(256), cpos4(256);
        static Array<int> wpos1(0x10000);
        static int rlen = 2, rlen1 = 3, rlen2 = 4;
        static int rcount1 = 0, rcount2 = 0;
        static ContextMap cm(MEM * 4, 7);
        if (!bpos) {
            int c = buf(1);
            int w = c4 & 0xffff;
            int r = pos - cpos1[c];
            if (r > 1 && r == cpos1[c] - cpos2[c]
                && r == cpos2[c] - cpos3[c] && r == cpos3[c] - cpos4[c]
                && (r > 15 || (c == buf(r * 5 + 1)) && c == buf(r * 6 + 1))) {
                if (r == rlen1)++rcount1;
                else if (r == rlen2)++rcount2;
                else if (rcount1 > rcount2) rlen2 = r, rcount2 = 1;
                else rlen1 = r, rcount1 = 1;
            }
            if (rcount1 > 15 && rlen != rlen1) rlen = rlen1, rcount1 = rcount2 = 0;
            if (rcount2 > 15 && rlen != rlen2) rlen = rlen2, rcount1 = rcount2 = 0;
            cm.set(buf(1) << 8 | min(255, pos - cpos1[buf(1)]));
            cm.set(buf(1) << 17 | buf(2) << 9 | llog(pos - wpos1[w]) >> 2);
            int col = pos % rlen;
            cm.set(buf(1) << 8 | buf(rlen));
            cm.set(rlen | buf(rlen) << 10 | buf(rlen * 2) << 18);
            cm.set(rlen | buf(rlen) << 10 | col << 18);
            cm.set(rlen | buf(1) << 10 | col << 18);
            cm.set(col | rlen << 12);
            cpos4[c] = cpos3[c];
            cpos3[c] = cpos2[c];
            cpos2[c] = cpos1[c];
            cpos1[c] = pos;
            wpos1[w] = pos;
        }
        cm.mix(m);
    }

    void indirectModel(Mixer &m) {
        static ContextMap cm(MEM * 4, 5);
        static Array<U32> t1(256);
        static Array<U16> t2(0x10000);
        if (!bpos) {
            U32 &r1 = t1[buf(2)];
            r1 = r1 << 8 | buf(1);
            U16 &r2 = t2[buf(3) << 8 | buf(2)];
            r2 = r2 << 8 | buf(1);
            U32 t = buf(1) | t1[buf(1)] << 8;
            cm.set(t & 0xffff);
            cm.set(t & 0xffffff);
            cm.set(t);
            t = buf(2) << 8 | buf(1) | t2[buf(2) << 8 | buf(1)] << 16;
            cm.set(t & 0xffffff);
            cm.set(t);
        }
        cm.mix(m);
    }

    // DMC submodel
    struct DMCNode{
        unsigned int nx[2];
        U8 state;
        unsigned int c0:12, c1:12;
    };
    static void dmcModel(Mixer& m) {
    static int top = 0, curr = 0;
    static Array<DMCNode> t(MEM * 2);
    static StateMap sm;
    static int threshold = 0x100;
    if (top > 0 && top < t.size()) {
        int next = t[curr].nx[y];
        int n = y ? t[curr].c1 : t[curr].c0;
        int nn = t[next].c0 + t[next].c1;
        if (n >= threshold * 2 && nn - n >= threshold * 3) {
            int r = n * 0x1000 / nn;
            t[next].c0 -= t[top].c0 = t[next].c0 * r >> 12;
            t[next].c1 -= t[top].c1 = t[next].c1 * r >> 12;
            t[top].nx[0] = t[next].nx[0];
            t[top].nx[1] = t[next].nx[1];
            t[top].state = t[next].state;
            t[curr].nx[y] = top++;
            if (top == (t.size() * 4) / 8) {
                threshold = 0x200;
            } else if (top == (t.size() * 6) / 8) {
                threshold = 0x300;
            }
        }
    }
    if (top == t.size() && bpos == 1) top = 0;
    if (top == 0) {
        for (int i = 0; i < 0x100; ++i) {
            for (int j = 0; j < 0x100; ++j) {
                if (i < 0x7f) {
                    t[j * 0x100 + i].nx[0] = j * 0x100 + i * 2 + 1;
                    t[j * 0x100 + i].nx[1] = j * 0x100 + i * 2 + 2;
                } else {
                    t[j * 0x100 + i].nx[0] = (i - 0x7f) * 0x100;
                    t[j * 0x100 + i].nx[1] = (i + 1) * 0x100;
                }
                t[j * 0x100 + i].c0 = 0xc0;
                t[j * 0x100 + i].c1 = 0xc0;
            }
        }
        top = 0x10000;
        curr = 0;
        threshold = 0x100;
    }
    if (y) {
        if (t[curr].c1 < 3800) t[curr].c1 += 0x100;
    } else if (t[curr].c0 < 3800) t[curr].c0 += 0x100;
    t[curr].state = State_table[t[curr].state][y];
    int ns = State_table[t[curr].state][y];
    if (ns >= 204 && rnd() << ((452 - ns) >> 3)) ns -= 4;
    t[curr].state = ns;
    curr = t[curr].nx[y];
    const int pr1 = sm.p(t[curr].state);
    const int n1 = t[curr].c1;
    const int n0 = t[curr].c0;
    const int pr2 = (n1 + 5) * 0x1000 / (n0 + n1 + 10);
    m.add(stretch(pr1));
    m.add(stretch(pr2));
}
}

void sparseModelz(Mixer& m, int seenbefore, int howmany) {
  static ContextMap cm(MEM*8, 40+2);
  if (bpos==0) {
    // cm.set(seenbefore);
    // cm.set(howmany);
    //cm.set(buf(1)|buf(5)<<8);
    cm.set(buf(1)|buf(6)<<8);
    // cm.set(buf(3)|buf(6)<<8);//bad
    // cm.set(buf(4) | buf(8) << 8);
    cm.set(buf(1)|buf(3)<<8|buf(5)<<16);
    // cm.set(buf(2)|buf(4)<<8|buf(6)<<16);//bad
    // cm.set(c4&0x00f0f0ff);//bad
    cm.set(c4&0x00ff00ff);
    // cm.set(c4&0xff0000ff);//bad
    cm.set(c4&0x00f8f8f8);
    // cm.set(c4&0xf8f8f8f8);//bad
    // cm.set(c4&0x00f0f0f0);//bad
    // cm.set(c4&0xf0f0f0f0);//bad
    //cm.set(c4&0x00e0e0e0);//bad
    //cm.set(c4&0xe0e0e0e0);//bad
    //cm.set(c4&0x810000c1);//bad
    //cm.set(c4&0xC3CCC38C);//bad
    cm.set(c4&0x0081CC81);
    //cm.set(c4&0x00c10081);//bad
    for (int i = 1; i< 12; ++i) {
    //   cm.set(seenbefore|buf(i)<<8); // it's bad aslan
      cm.set((buf(i+2)<<8)|buf(i+1));
    //   cm.set((buf(i+3)<<8)|buf(i+1));
    }
  }
  cm.mix(m);
}

void sparseModelzfk(Mixer& m, int seenbefore, int howmany, int ub) {
  static ContextMap cm(MEM*4, 22);
  if (bpos==0) {
    // cm.set(seenbefore);
    // cm.set(howmany);
    //cm.set(buf(1)|buf(5)<<8);
    cm.set(buf(1)|buf(6)<<8);
    // cm.set(buf(3)|buf(6)<<8);//bad
    // cm.set(buf(4) | buf(8) << 8);
    cm.set(buf(1)|buf(3)<<8|buf(5)<<16);
    // cm.set(buf(2)|buf(4)<<8|buf(6)<<16);//bad
    // cm.set(c4&0x00f0f0ff);//bad
    cm.set(c4&0x00ff00ff);
    // cm.set(c4&0xff0000ff);//bad
    cm.set(c4&0x00f8f8f8);
    // cm.set(c4&0xf8f8f8f8);//bad
    // cm.set(c4&0x00f0f0f0);//bad
    // cm.set(c4&0xf0f0f0f0);//bad
    //cm.set(c4&0x00e0e0e0);//bad
    //cm.set(c4&0xe0e0e0e0);//bad
    //cm.set(c4&0x810000c1);//bad
    //cm.set(c4&0xC3CCC38C);//bad
    cm.set(c4&0x0081CC81);
    //cm.set(c4&0x00c10081);//bad
    for (int i = 1; i< ub; ++i) {
    //   cm.set(seenbefore|buf(i)<<8); // it's bad aslan
      cm.set((buf(i+2)<<8)|buf(i+1));
    //   cm.set((buf(i+3)<<8)|buf(i+1));
    }
  }
  cm.mix(m);
}


void distanceModel(Mixer &m) {
    static ContextMap cr(MEM, 3);
    if (bpos == 0) {
        static int pos00 = 0, pos20 = 0, posnl = 0;
        int c = c4 & 0xff;
        if (c == 0x00) pos00 = pos;
        if (c == 0x20) pos20 = pos;
        // if (c == 0xff || c == '\r' || c == '\n') posnl = pos;
        cr.set(min(pos - pos00, 255) | (c << 8));
        cr.set(min(pos - pos20, 255) | (c << 8));
        cr.set(min(pos - posnl, 255) | ((c << 8) + 234567));
    }
    cr.mix(m);
}
void indirectModel(Mixer &m) {
    static ContextMap cm(MEM * 4, 5);
    static U32 t1[256];
    static U16 t2[0x10000];
    static U16 t3[0x8000];

    if (!bpos) {
        U32 d = c4 & 0xffff, c = d & 255, d2 = (buf(1) & 31) + 32 * (buf(2) & 31) + 1024 * (buf(3) & 31);
        U32 &r1 = t1[d >> 8];
        r1 = r1 << 8 | c;
        U16 &r2 = t2[c4 >> 8 & 0xffff];
        r2 = r2 << 8 | c;
        U16 &r3 = t3[(buf(2) & 31) + 32 * (buf(3) & 31) + 1024 * (buf(4) & 31)];
        r3 = r3 << 8 | c;
        const U32 t = c | t1[c] << 8;
        const U32 t0 = d | t2[d] << 16;
        const U32 ta = d2 | t3[d2] << 16;
        cm.set(t);
        cm.set(t0);
        cm.set(ta);
        cm.set(t & 0xff00);
        // cm.set(t0 & 0xff0000);//badorbelow
        // cm.set(ta & 0xff0000);//badorabove
        // cm.set(t & 0xffff); //bad
        // cm.set(t0 & 0xffffff);//bad
        cm.set(ta & 0xffffff);
    }
    cm.mix(m);
}
inline int llog(U32 x) {
    if (x >= 0x1000000)
        return 256 + ilog(x >> 16);
    else if (x >= 0x10000)
        return 128 + ilog(x >> 8);
    else
        return ilog(x);
}
int recordModel(Mixer &m, int rrlen = 0) {
    static int cpos1[256], cpos2[256], cpos3[256], cpos4[256];
    static int wpos1[0x10000]; // buf(1..2) -> last position
    static int rlen = 2, rlen1 = 3, rlen2 = 4, rlen3 = 5, rlenl = 0;  // run length and 2 candidates
    static int rcount1 = 0, rcount2 = 0, rcount3 = 0;  // candidate counts
    static ContextMap cm(32768, 3), cn(32768 / 2, 3), co(32768 * 2, 3), cp(MEM, 3);
    // Find record length
    if (!bpos) {
        int w = c4 & 0xffff, c = w & 255, d = w >> 8;
        int r = pos - cpos1[c];
        if (r > 1) {
            if (rrlen == 0) {
                if ((r == cpos1[c] - cpos2[c] || r == cpos2[c] - cpos3[c] || r == cpos3[c] - cpos4[c]
                    ) && (r > 10 || ((c == buf(r * 5 + 1)) && c == buf(r * 6 + 1)))) {
                    if (r == rlen1) ++rcount1;
                    else if (r == rlen2) ++rcount2;
                    else if (r == rlen3) ++rcount3;
                    else if (rcount1 > rcount2) rlen2 = r, rcount2 = 1;
                    else if (rcount2 > rcount3) rlen3 = r, rcount3 = 1;
                    else rlen1 = r, rcount1 = 1;
                }
            } else rlen = rrlen;
        }

        if (rcount1 > 12 && rlen != rlen1 && rlenl * 2 != rlen1)
            rlenl = rlen = rlen1, rcount1 = rcount2 = rcount3 = 0/*, printf("R1: %d  \n",rlen)*/;
        if (rcount2 > 18 && rlen != rlen2 && rlenl * 2 != rlen2)
            rlenl = rlen, rlen = rlen2, rcount1 = rcount2 = rcount3 = 0/*, printf("R2: %d  \n",rlen)*/;
        if (rcount3 > 24 && rlen != rlen3 && rlenl * 2 != rlen3)
            rlenl = rlen, rlen = rlen3, rcount1 = rcount2 = rcount3 = 0/*, printf("R3: %d  \n",rlen)*/;

        // Set 2 dimensional contexts
        assert(rlen > 0);
        cm.set(c << 8 | (min(255, pos - cpos1[c]) / 4));
        cm.set(w << 9 | llog(pos - wpos1[w]) >> 2);

        cm.set(rlen | buf(rlen) << 10 | buf(rlen * 2) << 18);
        cn.set(w | rlen << 8);
        cn.set(d | rlen << 16);
        cn.set(c | rlen << 8);

        co.set(buf(1) << 8 | min(255, pos - cpos1[buf(1)]));
        //co.set(buf(1) << 17 | buf(2) << 9 | llog(pos - wpos1[w]) >> 2);//bad
        int col = pos % rlen;
        co.set(buf(1) << 8 | buf(rlen));

        //cp.set(w*16);
        //cp.set(d*32);
        //cp.set(c*64);
        cp.set(rlen | buf(rlen) << 10 | col << 18);
        cp.set(rlen | buf(1) << 10 | col << 18);
        cp.set(col | rlen << 12);
        // update last context positions
        cpos4[c] = cpos3[c];
        cpos3[c] = cpos2[c];
        cpos2[c] = cpos1[c];
        cpos1[c] = pos;
        wpos1[w] = pos;
    }
    cm.mix(m);
    cn.mix(m);
    co.mix(m);
    cp.mix(m);
    return rlen;
}
struct DMCNode {  // 12 bytes
    unsigned int nx[2];  // next pointers
    U8 state;  // bit history
    unsigned int c0: 12, c1: 12;  // counts * 256
};
#define nex(state, sel) State_table[state][sel]
void dmcModel(Mixer &m) {
    static int top = 0, curr = 0;  // allocated, current node
    static Array<DMCNode> t(MEM * 2);  // state graph
    static StateMap sm;
    static int threshold = 256;

    // clone next state
    if (top > 0 && top < t.size()) {
        int next = t[curr].nx[y];
        int n = y ? t[curr].c1 : t[curr].c0;
        int nn = t[next].c0 + t[next].c1;
        if (n >= threshold * 2 && nn - n >= threshold * 3) {
            int r = n * 4096 / nn;
            assert(r >= 0 && r <= 4096);
            t[next].c0 -= t[top].c0 = t[next].c0 * r >> 12;
            t[next].c1 -= t[top].c1 = t[next].c1 * r >> 12;
            t[top].nx[0] = t[next].nx[0];
            t[top].nx[1] = t[next].nx[1];
            t[top].state = t[next].state;
            t[curr].nx[y] = top;
            ++top;
            if (top == MEM * 2) threshold = 512;
            if (top == MEM * 3) threshold = 768;
        }
    }

    // Initialize to a bytewise order 1 model at startup or when flushing memory
    if (top == t.size() && bpos == 1) top = 0;
    if (top == 0) {
        assert(t.size() >= 65536);
        for (int i = 0; i < 256; ++i) {
            for (int j = 0; j < 256; ++j) {
                if (i < 127) {
                    t[j * 256 + i].nx[0] = j * 256 + i * 2 + 1;
                    t[j * 256 + i].nx[1] = j * 256 + i * 2 + 2;
                } else {
                    t[j * 256 + i].nx[0] = (i - 127) * 256;
                    t[j * 256 + i].nx[1] = (i + 1) * 256;
                }
                t[j * 256 + i].c0 = 128;
                t[j * 256 + i].c1 = 128;
            }
        }
        top = 65536;
        curr = 0;
        threshold = 256;
    }

    // update count, state
    if (y) {
        if (t[curr].c1 < 3800) t[curr].c1 += 256;
    } else if (t[curr].c0 < 3800) t[curr].c0 += 256;
    t[curr].state = nex(t[curr].state, y);
    curr = t[curr].nx[y];

    // predict
    const int pr1 = sm.p(t[curr].state);
    const int n1 = t[curr].c1;
    const int n0 = t[curr].c0;
    const int pr2 = (n1 + 5) * 4096 / (n0 + n1 + 10);
    m.add(stretch(pr1));
    m.add(stretch(pr2));
}


static APM a1(0x100), a2(0x100), a3(0x100), a4(0x100);

static int predictNextTooSmall(){
  static ContextMap cm(MEM * 64, 30);
  static Mixer mixer(300, 1160, 5);
//   static APM a1(0x100), a2(0x100), a3(0x100);
  static U32 cxt[15], t1[0x100];
  static U16 t2[0x10000];
  static U32 mask = 0, mask2 = 0, word0 = 0, word1 = 0;

  c0 += c0 + y;
  if (c0 >= 256) {
    buf[pos++] = c0;
    c4 = (c4 << 8) + c0 - 256;
    c0 = 1;
  }   
  bpos = (bpos + 1) & 7;  
  int c1 = c4 & 0xff, c2 = (c4 & 0xff00) >> 8, c3 = (c4 & 0xff0000) >> 16;

  mixer.update();
  int ismatch = ilog(matchModel(mixer));
  dmcModel(mixer);

  if (bpos == 0) {
    // for (int i = 14; i > 0; --i) {
    //   cxt[i] = hash(cxt[i - 1], c1);
    // }
    cm.set(0);
    cm.set(c1);
    cm.set(c4 & 0x0000ffff);
    cm.set(c4 & 0x00ffffff);
    cm.set(c4);
    // for (int i = 0; i < 7; ++i) cm.set(cxt[i]);
    // cm.set(cxt[8]);        
    // cm.set(cxt[14]);
    /* Main was:
    cm.set(cxt[5]);
    cm.set(cxt[6]);
    cm.set(cxt[14]);
     */
    cm.set(c4 & 0xf8f8c0ff);
    cm.set(c4 & 0x00e0e0e0);
    cm.set(c4 & 0xffc0ff80);
    cm.set(ismatch | (c4 & 0xffff0000));
    cm.set(ismatch | (c4 & 0x0000ff00));
    cm.set(ismatch | (c4 & 0x00ff0000));
    mask = (mask << 3) | (!c1 ? 0 : (c1 == 255) ? 4 : (c1 < 16) ? 5 : (c1 < 64) ? 6 : 7);
    cm.set(mask);
    // 0.181 total score increase
    //mask2 = (mask2 << 3) | ((mask >> 27) & 7);
    //cm.set(hash(mask << 5, mask2 << 2));

    // by refernce change
    U32& ic1r = t1[c2];
    ic1r = ic1r << 8 | c1;
    U16& ic2r = t2[(buf(3) << 8) | c2];
    ic2r = ic2r << 8 | c1;

    const U32 ic1 = c1 | t1[c1] << 8;
    const U32 ic2 = ((c2 << 8) | c1) | t2[(c2 << 8) | c1] << 16;
    cm.set((ic1 >> 8) & ((1 << 16) - 1));
    cm.set((ic2 >> 16) & ((1 << 8) - 1));
    cm.set(ic1 & ((1 << 16) - 1));
    cm.set(ic2 & ((1 << 24) - 1));  
  }

  int o = cm.mix(mixer);

  sparseModelzfk(mixer, ismatch, o, 12);
  indirectModel(mixer);
  recordModel(mixer);
  distanceModel(mixer);
  //if (len4 < 40'000)
  //  old::wordModel(mixer);


  mixer.set(c1 + 8, 264);
  mixer.set(c0, 256);
  mixer.set(o + ((c1 > 32) << 4) + ((bpos == 0) << 5) + ((c1 == c2) << 6), 128);
  mixer.set(c2, 256);
  mixer.set(ismatch, 256);  
     
    int pr0 = mixer.p();
  return (a1.p(pr0, c0) * 11
         + a2.p(pr0, c3) * 10
         + a3.p(pr0, c1) * 6
         + a4.p(pr0, c2, 5) * 5
         + 16) >> 5; // Probability adjusted with 3 APMs                             
}

static int predictNextSmall(){
  static ContextMap cm(MEM * 32, 30);
  static Mixer mixer(300, 1160, 5);
//   static APM a1(0x100), a2(0x100), a3(0x100);
  static U32 cxt[15], t1[0x100];
  static U16 t2[0x10000];
  static U32 mask = 0, mask2 = 0, word0 = 0, word1 = 0;

  c0 += c0 + y;
  if (c0 >= 256) {
    buf[pos++] = c0;
    c4 = (c4 << 8) + c0 - 256;
    c0 = 1;
  }   
  bpos = (bpos + 1) & 7;  
  int c1 = c4 & 0xff, c2 = (c4 & 0xff00) >> 8, c3 = (c4 & 0xff0000) >> 16;

  mixer.update();
  int ismatch = ilog(matchModel(mixer));
  dmcModel(mixer);

  if (bpos == 0) {
    // for (int i = 14; i > 0; --i) {
    //   cxt[i] = hash(cxt[i - 1], c1);
    // }
    cm.set(0);
    cm.set(c1);
    cm.set(c4 & 0x0000ffff);
    cm.set(c4 & 0x00ffffff);
    cm.set(c4);
    // for (int i = 0; i < 7; ++i) cm.set(cxt[i]);
    // cm.set(cxt[8]);        
    // cm.set(cxt[14]);
    // cm.set(cxt[5]);
    // cm.set(cxt[6]);
    // cm.set(cxt[14]);
    cm.set(c4 & 0xf8f8c0ff);
    cm.set(c4 & 0x00e0e0e0);
    cm.set(c4 & 0xffc0ff80);
    cm.set(ismatch | (c4 & 0xffff0000));
    cm.set(ismatch | (c4 & 0x0000ff00));
    cm.set(ismatch | (c4 & 0x00ff0000));
    mask = (mask << 3) | (!c1 ? 0 : (c1 == 255) ? 4 : (c1 < 16) ? 5 : (c1 < 64) ? 6 : 7);
    cm.set(mask);
    // //0.181 total score increase
    // mask2 = (mask2 << 3) | ((mask >> 27) & 7);
    // cm.set(hash(mask << 5, mask2 << 2));

    // by refernce change
    U32& ic1r = t1[c2];
    ic1r = ic1r << 8 | c1;
    U16& ic2r = t2[(buf(3) << 8) | c2];
    ic2r = ic2r << 8 | c1;

    const U32 ic1 = c1 | t1[c1] << 8;
    const U32 ic2 = ((c2 << 8) | c1) | t2[(c2 << 8) | c1] << 16;
    cm.set((ic1 >> 8) & ((1 << 16) - 1));
    cm.set((ic2 >> 16) & ((1 << 8) - 1));
    cm.set(ic1 & ((1 << 16) - 1));
    cm.set(ic2 & ((1 << 24) - 1));  
  }
  
  int o = cm.mix(mixer);

  sparseModelz(mixer, ismatch, o);
//   distanceModel(mixer);
//    if (len4 < 110'000) recordModel(mixer); // adds 8 total points, up to 1.78s
//     // old::recordModel(mixer);
//    else 
   indirectModel(mixer);    

  mixer.set(c1 + 8, 264);
  mixer.set(c0, 256);
  mixer.set(o + ((c1 > 32) << 4) + ((bpos == 0) << 5) + ((c1 == c2) << 6), 128);
  mixer.set(c2, 256);
  mixer.set(ismatch, 256);  
     
  int pr0 = mixer.p();
  return (a1.p(pr0, c0) * 11
         + a2.p(pr0, c3) * 10
         + a3.p(pr0, c1) * 6
         + a4.p(pr0, c2, 5) * 5
         + 16) >> 5; // Probability adjusted with 3 APMs                                          
}

static int predictNextMid(){
  static ContextMap cm(MEM * 64, 30);
  static Mixer mixer(300, 1160, 5);
//   static APM a1(0x100), a2(0x100), a3(0x100), a4(0x100);
  static U32 cxt[15], t1[0x100];
  static U16 t2[0x10000];
  static U32 mask = 0, mask2 = 0, word0 = 0, word1 = 0;

  c0 += c0 + y;
  if (c0 >= 256) {
    buf[pos++] = c0;
    c4 = (c4 << 8) + c0 - 256;
    c0 = 1;
  }   
  bpos = (bpos + 1) & 7;  
  int c1 = c4 & 0xff, c2 = (c4 & 0xff00) >> 8, c3 = (c4 & 0xff0000) >> 16;

  mixer.update();
int ismatch = ilog(matchModel(mixer));
//   dmcModel(mixer);

  if (bpos == 0) {
    // for (int i = 14; i > 0; --i) {
    //   cxt[i] = hash(cxt[i - 1], c1);
    // }
    cm.set(0);
    cm.set(c1);
    cm.set(c4 & 0x0000ffff);
    cm.set(c4 & 0x00ffffff);
    cm.set(c4);
    // for (int i = 0; i < 7; ++i) cm.set(cxt[i]);
    // cm.set(cxt[8]);        
    // cm.set(cxt[14]);
    // cm.set(cxt[5]);
    // cm.set(cxt[6]);
    // cm.set(cxt[14]);
    cm.set(c4 & 0xf8f8c0ff);
    cm.set(c4 & 0x00e0e0e0);
    cm.set(c4 & 0xffc0ff80);
     cm.set(ismatch | (c4 & 0xffff0000));
     cm.set(ismatch | (c4 & 0x0000ff00));
     cm.set(ismatch | (c4 & 0x00ff0000));
    mask = (mask << 3) | (!c1 ? 0 : (c1 == 255) ? 4 : (c1 < 16) ? 5 : (c1 < 64) ? 6 : 7);
    cm.set(mask);
    // //0.181 total score increase
    // mask2 = (mask2 << 3) | ((mask >> 27) & 7);
    // cm.set(hash(mask << 5, mask2 << 2));

    // by refernce change
    U32& ic1r = t1[c2];
    ic1r = ic1r << 8 | c1;
    U16& ic2r = t2[(buf(3) << 8) | c2];
    ic2r = ic2r << 8 | c1;

    const U32 ic1 = c1 | t1[c1] << 8;
    const U32 ic2 = ((c2 << 8) | c1) | t2[(c2 << 8) | c1] << 16;
    cm.set((ic1 >> 8) & ((1 << 16) - 1));
    cm.set((ic2 >> 16) & ((1 << 8) - 1));
    cm.set(ic1 & ((1 << 16) - 1));
    cm.set(ic2 & ((1 << 24) - 1));  
  }
  
  int o = cm.mix(mixer);
//   sparseModelz(mixer, ismatch, o);//1296,timelimit
    // sparseModel(mixer);//1291,timelimit
    // sparseModel12(mixer);//1290,timelimit
//    recordModel(mixer);//1289,timelimit
    //  distanceModel(mixer);//1288
//   indirectModel(mixer);//1287
    // sparseModel4(mixer);//1287
    // sparseModel6(mixer);//1287

   if (len4 < 160'000)
    sparseModelzfk(mixer, ismatch, o, 8);
   else if (len4 < 200'000)
    sparseModelzfk(mixer, ismatch, o, 4);
//     distanceModel(mixer);



  mixer.set(c1 + 8, 264);
  mixer.set(c0, 256);
  mixer.set(o + ((c1 > 32) << 4) + ((bpos == 0) << 5) + ((c1 == c2) << 6), 128);
  mixer.set(c2, 256);
  mixer.set(ismatch, 256);  
     
  int pr0 = mixer.p();//1288.097
  return (a1.p(pr0, c0) * 11
         + a2.p(pr0, c3) * 10
         + a3.p(pr0, c1) * 6
         + a4.p(pr0, c2, 5) * 5
         + 16) >> 5; // Probability adjusted with 3 APMs                                                            
}

static int predictNextBig(){
  static ContextMap cm(MEM * 64, 30);
  static Mixer mixer(300, 1160, 5);
//   static APM a1(0x100), a2(0x100), a3(0x100);
  static U32 cxt[6], t1[0x100];
  static U16 t2[0x10000];
  static U32 mask = 0, mask2 = 0, word0 = 0, word1 = 0;

  c0 += c0 + y;
  if (c0 >= 256) {
    buf[pos++] = c0;
    c4 = (c4 << 8) + c0 - 256;
    c0 = 1;
  }   
  bpos = (bpos + 1) & 7;  
  int c1 = c4 & 0xff, c2 = (c4 & 0xff00) >> 8, c3 = (c4 & 0xff0000) >> 16;

  mixer.update();
//   dmcModel(mixer);

  int ismatch = ilog(matchModel(mixer));
  if (bpos == 0) {
    cm.set(0);
    cm.set(c1);
    // cm.set(c4 & 0x0000ffff);
    cm.set(c4 & 0x00ffffff);
    // cm.set(c4);
    //cm.set(cxt[5]);
    // cm.set(cxt[6]);
    // cm.set(cxt[14]);
    // cm.set(c4 & 0xf8f8c0ff);
    // cm.set(c4 & 0x00e0e0e0);
    cm.set(c4 & 0xffc0ff80);
    // cm.set(ismatch | (c4 & 0xffff0000));
    cm.set(ismatch | (c4 & 0x0000ff00));
    // cm.set(ismatch | (c4 & 0x00ff0000));
    mask = (mask << 3) | (!c1 ? 0 : (c1 == 255) ? 4 : (c1 < 16) ? 5 : (c1 < 64) ? 6 : 7);
    cm.set(mask);

    // by refernce change
    U32& ic1r = t1[c2];
    ic1r = ic1r << 8 | c1;
    U16& ic2r = t2[(buf(3) << 8) | c2];
    ic2r = ic2r << 8 | c1;

    const U32 ic1 = c1 | t1[c1] << 8;
    const U32 ic2 = ((c2 << 8) | c1) | t2[(c2 << 8) | c1] << 16;
    cm.set((ic1 >> 8) & ((1 << 16) - 1));
    cm.set((ic2 >> 16) & ((1 << 8) - 1));
    // cm.set(ic1 & ((1 << 16) - 1));
    cm.set(ic2 & ((1 << 24) - 1));  
  }
  
  int o = cm.mix(mixer);

  mixer.set(c1 + 8, 264);
  mixer.set(c0, 256);
  mixer.set(o + ((c1 > 32) << 4) + ((bpos == 0) << 5) + ((c1 == c2) << 6), 128);
  mixer.set(c2, 256);
  mixer.set(ismatch, 256);  
     
  int pr0 = mixer.p();
  return (a1.p(pr0, c0) * 11
         + a2.p(pr0, c3) * 10
         + a3.p(pr0, c1) * 6
         + a4.p(pr0, c2, 5) * 5
         + 16) >> 5; // Probability adjusted with 3 APMs                           
}

// Main model - predicts next bit probability from previous data
static int predictNext() {
    if (TOO_SMALL) return predictNextTooSmall();
    if (SMALL) return predictNextSmall();
    if (MID) return predictNextMid();
    else return predictNextBig();
}


// Encoder - Arithmetic coding
class Encoder {
private:
    const bool mode;
    U32 x, x1, x2;
    int p;
    int code(int i = 0) {
        p += p < 2048;
        U32 xmid = x1 +((x2 - x1) >> 12) * p + (((x2 - x1) & 0xfff) * p >> 12);
        if (!mode) y = x <= xmid; else y = i;
        y ? (x2 = xmid) : (x1 = xmid+1);
        // if (SMALL) p = old::predictNext();
        // else 
        p = predictNext(); // Update models and predict next bit probability
        while (((x1 ^ x2) & 0xff000000) == 0) {
            if (mode) putc(x2 >> 24, outfile);
            x1 <<= 8;
            x2 = (x2 << 8) + 255;
            if (!mode) x = (x << 8) + (getc(infile) & 255);
        }
        return y;
    }
public:
    Encoder(bool m): mode(m), x(0), x1(0), x2(0xffffffff), p(2048) {
        if (!mode) {
            for (int i = 0; i < 4; ++i) {
                x = (x << 8) + (getc(infile) & 255);
            }
        }
    }
    void flush() {
        putc(x1 >> 24, outfile);
    }
    void compress(int c) {
        for (int i = 7; i >= 0; --i) {
            code((c >> i) & 1);
        }
        // cm.set(c);
    }
    int decompress() {
        int c = 0;
        for (int i = 0; i < 8; ++i) {
            c += c + code();
        }
        // cm.set(c);
        return c;
    }
};


void paq_compress8() {
    if (((infile = fopen("i", "rb")) == NULL) || ((outfile = fopen("o", "wb")) == NULL)) {
        // printf("Cannot read/write file\n");
        exit(1);
    }
    Encoder en(true);
    fseek(infile, 0, SEEK_END);
    long len = ftell(infile);
    fseek(infile, 0, SEEK_SET);
    en.compress(len >> 24);
    en.compress(len >> 16);
    en.compress(len >> 8);
    en.compress(len);
    len4 = len;
    for (long i = 0; i < len; i++) {
        en.compress(getc(infile));
    }
    en.flush();
    fclose(infile);
    fclose(outfile);
}


void paq_decompress8() {
    if (((infile = fopen("o", "rb")) == NULL) || ((outfile = fopen("i", "wb")) == NULL)) {
        // printf("Cannot read/write file\n");
        exit(1);
    }
    Encoder en(false);
    long len;
    len = en.decompress() << 24;
    len |= en.decompress() << 16;
    len |= en.decompress() << 8;
    len |= en.decompress();
    len4 = len;
    for (long i = 0; i < len; i++) {
        putc(en.decompress(), outfile);
    }
    fclose(infile);
    fclose(outfile);
}


int main() {
    // ios_base::sync_with_stdio(false);
    // cin.tie(NULL);
    // cout.tie(NULL);
    std::string mode;
    std::cin >> mode;
    CHECK(mode == "compress" || mode == "decompress");
    std::string base64_data;
    std::cin >> base64_data;
    CHECK(!base64_data.empty());
    td::BufferSlice data(td::base64_decode(base64_data).move_as_ok());
    if (mode == "compress") {
        int sz = data.size();
        auto data4 = utils::change_mode_to_compress(data, 4);
        data = data4.clone();
        std::vector<unsigned char> compressed_data;
        std::ofstream out("i", std::ios::binary);
        out.write(reinterpret_cast<const char *>(data.data()), data.size());
        out.close();
        paq_compress8();
        std::ifstream in("o", std::ios::binary);
        compressed_data = std::vector<unsigned char>((std::istreambuf_iterator<char>(in)),
                                                     std::istreambuf_iterator<char>());
        in.close();
        data = td::BufferSlice(reinterpret_cast<char *>(compressed_data.data()), compressed_data.size());
    } else {
        std::ofstream out("o", std::ios::binary);
        out.write(reinterpret_cast<const char *>(data.data()), data.size());
        out.close();
        paq_decompress8();
        std::ifstream in("i", std::ios::binary);
        std::vector<unsigned char> decompressed_data = std::vector<unsigned char>(
                (std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        in.close();
        data = td::BufferSlice(reinterpret_cast<char *>(decompressed_data.data()), decompressed_data.size());
        data = utils::change_mode_back_to_decompress(data);
    }
    std::cout << td::base64_encode(data) << std::endl;
}//try900onesfjlefile//kepolddmc2//stable???3-2ndsol3//5//2?