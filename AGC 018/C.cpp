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
#define MAXN 100010

int x, y, z, n;
int g[MAXN], s[MAXN], foo[MAXN];

ll sliver[MAXN], gold[MAXN];
priority_queue<int, vector<int>, greater<int>> pq;

bool cmp(int& a, int& b) {
    return (ll)g[a] - s[a] < (ll)g[b] - s[b];
}

int main() {
    cin >> x >> y >> z;
    n = x + y + z;

    ll sum = 0;
    for (int i = 0; i < n; ++i) {
        int b;
        scanf("%d%d%d", g + i, s + i, &b);
        sum += b;
        g[i] -= b;
        s[i] -= b;
        foo[i] = i;
    }

    sort(foo, foo + n, cmp);

    for (int i = 0; i < y; ++i) {
        pq.push(s[foo[i]]);
        sliver[y - 1] += s[foo[i]];
    }

    for (int i = y; i < n - x; ++i) {
        pq.push(s[foo[i]]);
        sliver[i] = sliver[i - 1] + s[foo[i]] - pq.top();
        pq.pop();
    }

    while (!pq.empty()) pq.pop();
    for (int i = n - 1; i >= n - x; --i) {
        pq.push(g[foo[i]]);
        gold[n - x] += g[foo[i]];
    }
    for (int i = n - x - 1; i >= y; --i) {
        pq.push(g[foo[i]]);
        gold[i] = gold[i + 1] + g[foo[i]] - pq.top();
        pq.pop();
    }

    ll res = LLONG_MIN;
    for (int i = y - 1; i < n - x; ++i) {
        res = max(res, sliver[i] + gold[i + 1]);
    }

    cout << sum + res;

    return 0;
}