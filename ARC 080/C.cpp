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
#define MAXN 200010
#define L (id << 1)
#define R ((id << 1) | 1)

int n;
int p[MAXN], q[MAXN];

int oddmin[MAXN << 2], evenmin[MAXN << 2];

void build(int id, int l, int r) {
    if (l == r) {
        oddmin[id] = l;
        evenmin[id] = -1;
        return;
    }

    int mid = (l + r) >> 1;
    build(L, l, mid);
    build(R, mid + 1, r);

    if ((mid - l) % 2) {
        if (p[oddmin[L]] < p[oddmin[R]]) {
            oddmin[id] = oddmin[L];
        } else {
            oddmin[id] = oddmin[R];
        }

        if (evenmin[L] == -1) {
            evenmin[id] = evenmin[R];
        } else {
            if (evenmin[R] == -1 || p[evenmin[L]] < p[evenmin[R]]) {
                evenmin[id] = evenmin[L];
            } else {
                evenmin[id] = evenmin[R];
            }
        }
    } else {
        if (evenmin[R] == -1 || p[oddmin[L]] < p[evenmin[R]]) {
            oddmin[id] = oddmin[L];
        } else {
            oddmin[id] = evenmin[R];
        }

        if (evenmin[L] == -1 || p[oddmin[R]] < p[evenmin[L]]) {
            evenmin[id] = oddmin[R];
        } else {
            evenmin[id] = evenmin[L];
        }
    }
}

int query(int id, int l, int r, int ql, int qr) {
    if (l > qr || r < ql) return -1;

    bool tmp = ((l - ql) % 2 == 0);

    if (l >= ql && r <= qr) {
        return tmp ? oddmin[id] : evenmin[id];
    }

    int mid = (l + r) >> 1;

    int res1 = query(L, l, mid, ql, qr);
    int res2 = query(R, mid + 1, r, ql, qr);

    if (res1 == -1) return res2;
    if (res2 == -1) return res1;
    if (p[res1] < p[res2]) {
        return res1;
    } else {
        return res2;
    }
}

int main() {
    cin >> n;
    for (int i = 0; i < n; ++i) {
        scanf("%d", p + i);
    }

    build(1, 0, n - 1);

    // left, right, oddminimum
    auto cmp = [](vector<int>& a, vector<int>& b) {
        return p[a[2]] > p[b[2]];
    };
    priority_queue<vector<int>, vvi, decltype(cmp)> pq(cmp);
    int cur = 0;

    pq.push({0, n - 1, query(1, 0, n - 1, 0, n - 1)});
    while (!pq.empty()) {
        vector<int> vec = pq.top();
        pq.pop();
        int end = query(1, 0, n - 1, vec[2] + 1, vec[1]);
        if (vec[0] != vec[2]) {
            pq.push({vec[0], vec[2] - 1, query(1, 0, n - 1, vec[0], vec[2] - 1)});
        }
        if (vec[2] + 1 != end) {
            pq.push({vec[2] + 1, end - 1, query(1, 0, n - 1, vec[2] + 1, end - 1)});
        }
        if (end != vec[1]) {
            pq.push({end + 1, vec[1], query(1, 0, n - 1, end + 1, vec[1])});
        }
        q[cur++] = p[vec[2]];
        q[cur++] = p[end];
    }

    for (int i = 0; i < n; ++i) {
        printf("%d%c", q[i], " \n"[i == n - 1]);
    }

    return 0;
}