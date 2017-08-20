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
string s;
int dp[MAXN][26];
vector<int> pos[26];

int main() {
    memset(dp, 0, sizeof(dp));
    cin >> s;
    dp[s.length() - 1][s.back() - 'a'] = 1;
    for (int i = s.length() - 2; i >= 0; --i) {
        int tmp = MAXN;
        for (int j = 0; j < 26; ++j) {
            tmp = min(tmp, dp[i + 1][j]);
        }
        for (int j = 0; j < 26; ++j) {
            dp[i][j] = dp[i + 1][j];
        }
        dp[i][s[i] - 'a'] = 1 + tmp;
        pos[s[i] - 'a'].pb(i);
    }
    for (int i = 0; i < 26; ++i) {
        reverse(all(pos[i]));
    }

    int mini = MAXN;
    for (int i = 0; i < 26; ++i) {
        mini = min(mini, dp[0][i]);
    }

    for (int i = 0; i < 26; ++i) {
        if (dp[0][i] == mini) {
            if (mini == 0) {
                cout << (char)('a' + i);
                return 0;
            }
            string tmp;
            char ch = 'a' + i;
            int cur = -1;
            for (int j = 0; j < mini; ++j) {
                int next = *upper_bound(all(pos[ch - 'a']), cur);
                tmp.pb(ch);
                cur = next;
                for (int k = 0; k < 26; ++k) {
                    if (dp[next + 1][k] == dp[next][ch - 'a'] - 1) {
                        ch = 'a' + k;
                        break;
                    }
                }
            }
            tmp.pb(ch);
            cout << tmp;
            return 0;
        }
    }

    return 0;
}