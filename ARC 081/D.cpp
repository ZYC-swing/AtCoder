#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cmath>
#include <string>
#include <cstring>
#include <map>
#include <set>
#include <bitset>
#include <sstream>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <climits>
#include <ctype.h>
using namespace std;

#define PI acos(-1.0)
#define all(x) (x).begin(), (x).end()
#define pb push_back
#define fi first
#define se second

#define eps 1e-8
#define mod 1000000007

typedef long long ll;
typedef pair<int, int> pii;
typedef vector<vector<int>> vvi;

struct Matrix {
    vvi data;
    int r, c;
    Matrix(int row, int col, bool identity = false) : r(row), c(col) {
        data.assign(row, vector<int>(col, 0));
        if (identity) {
            for (int i = 0; i < r; ++i) {
                data[i][i] = 1;
            }
        }
    }
    Matrix operator * (Matrix& other) {
        int m = r, n = c, p = other.c;
        Matrix res(m, p);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < p; ++j) {
                for (int k = 0; k < n; ++k) {
                    res.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return res;
    }
};

// 扩展欧几里得算法 res[1] * a + res[2] * b = res[0] = gcd(a, b)
vector<int> extendGcd(int a, int b) {
    if (b == 0) {
        return {a, 1, 0};
    } else {
        vector<int> tmp = extendGcd(b, a % b);
        return {tmp[0], tmp[2], tmp[1] - (a / b) * tmp[2]};
    }
}

// 矩阵快速幂
Matrix matrix_power(Matrix base, ll exp) {
    int n = base.r;
    Matrix res(n, n, true);
    while (exp) {
        if (exp & 1) {
            res = res * base;
        }
        base = base * base;
        exp >>= 1;
    }
    return res;
}

// 带模快速幂
ll power_mod(ll base, int exp) {
    ll res = 1;
    while (exp) {
        if (exp & 1) res = res * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return res;
}

// 快速幂
ll power(ll base, int exp) {
    ll res = 1;
    while (exp) {
        if (exp & 1) res *= base;
        base *= base;
        exp >>= 1;
    }
    return res;
}

// 求逆元
ll inv(ll a) {
    return power_mod(a, mod - 2);
}

#define MAXFAC 100010

ll fac[MAXFAC];

void initFac() {
    fac[0] = 1;
    for (int i = 1; i < MAXFAC; ++i) {
        fac[i] = i * fac[i - 1] % mod;
    }
}

ll Combine(ll a, ll b) {
    return (fac[a] * inv(fac[b]) % mod) * inv(fac[a - b]) % mod;
}

/******************************** template ********************************/
#define MAXN 2010
int H, W;
string grid[MAXN];
int foo[MAXN][MAXN], h[MAXN], l[MAXN], r[MAXN], lbound[MAXN], rbound[MAXN];

int main() {
    cin >> H >> W;
    for (int i = 0; i < H; ++i) {
        cin >> grid[i];
    }

    memset(foo, 0, sizeof(foo));
    for (int i = 0; i < H - 1; ++i) {
        for (int j = 0; j < W - 1; ++j) {
            int tmp = 0;
            for (int x = i; x <= i + 1; ++x) {
                for (int y = j; y <= j + 1; ++y) {
                    tmp += grid[x][y] == '#';
                }
            }
            if (tmp % 2) foo[i + 1][j + 1] = 1;
        }
    }

//    cout << endl;
//    for (int i = 0; i < H + 1; ++i) {
//        for (int j = 0; j < W + 1; ++j) {
//            printf("%d%c", foo[i][j], " \n"[j == W]);
//        }
//    }
    memset(h, 0, sizeof(h));
    memset(l, 0, sizeof(l));
    int ind = 0;
    for (int i = 0; i <= W; ++i) {
        lbound[i] = ind;
        if (foo[0][i] == 1) {
            ind = i;
        }
    }
    ind = W;
    for (int i = W; i >= 0; --i) {
        r[i] = W;
        rbound[i] = ind;
        if (foo[0][i] == 1) {
            ind = i;
        }
    }

    int res = W;
    for (int i = 1; i <= H; ++i) {
        for (int j = 0; j <= W; ++j) {
            if (foo[i - 1][j] == 1) {
                h[j] = 1;
                l[j] = 0;
                r[j] = W;
            } else {
                ++h[j];
                l[j] = max(l[j], lbound[j]);
                r[j] = min(r[j], rbound[j]);
                res = max(res, h[j] * (r[j] - l[j]));
            }
        }

        int ind = 0;
        for (int j = 0; j <= W; ++j) {
            lbound[j] = ind;
            if (foo[i][j] == 1) {
                ind = j;
            }
        }
        ind = W;
        for (int j = W; j >= 0; --j) {
            rbound[j] = ind;
            if (foo[i][j] == 1) {
                ind = j;
            }
        }
    }

    cout << res;

    return 0;
}