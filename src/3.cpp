#include <bits/stdc++.h>
using namespace std;

/*** ---------- Solver (same objective as problem) ---------- ***/
class Solution {
public:
    double computeMinVariance(const vector<vector<double>>& grids,
                              const vector<double>& blocks) {
        vector<double> G;
        for (const auto& row : grids) for (double x : row) G.push_back(x);
        const int n = (int)G.size();
        if (n <= 1) return 0.0; // sample variance for 1 point -> 0

        sort(G.begin(), G.end());                   // ascending
        vector<double> B = blocks;
        sort(B.begin(), B.end(), greater<double>());// descending

        // Precompute for const part sum (g_i - t)^2 = sum_g2 - 2t*sum_g + n*t^2
        double sum_g = 0.0, sum_g2 = 0.0;
        double gmin = 1e300, gmax = -1e300;
        for (double g : G) {
            sum_g += g; sum_g2 += g * g;
            gmin = min(gmin, g); gmax = max(gmax, g);
        }
        double bmin = 1e300, bmax = -1e300;
        for (double b : B) { bmin = min(bmin, b); bmax = max(bmax, b); }

        // Bound t inside [min possible y, max possible y]
        double lo = gmin + bmin, hi = gmax + bmax;

        auto eval_at_t = [&](double t) -> double {
            // const term from grids about t
            double const_term = sum_g2 - 2.0 * t * sum_g + n * t * t;

            // weights for the ordered selection (k-th selected gets weight 2*(g[k]-t))
            vector<double> w(n);
            for (int k = 0; k < n; ++k) w[k] = 2.0 * (G[k] - t);

            // DP over blocks: dp[k] = min cost selecting exactly k items so far
            const double INF = 1e300;
            vector<double> dp(n + 1, INF), next(n + 1, INF);
            dp[0] = 0.0;

            for (double b : B) {
                // skip
                for (int k = 0; k <= n; ++k) next[k] = dp[k];
                // take -> becomes the (k+1)-th largest among selected (since B is desc)
                for (int k = 0; k < n; ++k) {
                    if (dp[k] >= INF/2) continue;
                    double cand = dp[k] + b * b + w[k] * b;
                    if (cand < next[k + 1]) next[k + 1] = cand;
                }
                dp.swap(next);
            }

            double best = dp[n]; // minimal extra cost from choosing n blocks
            return (const_term + best) / n; // sample variance
        };

        // Ternary search on convex function
        for (int it = 0; it < 80; ++it) {
            double m1 = (2.0 * lo + hi) / 3.0;
            double m2 = (lo + 2.0 * hi) / 3.0;
            double f1 = eval_at_t(m1);
            double f2 = eval_at_t(m2);
            if (f1 < f2) hi = m2; else lo = m1;
        }
        return eval_at_t((lo + hi) * 0.5);
    }
};

/*** ---------- Utilities for testcase printing ---------- ***/
static string fmt_scientific(double x) {
    ostringstream oss;
    oss.setf(std::ios::scientific);
    oss << setprecision(4) << x;
    return oss.str();
}
static void print_grids(const vector<vector<double>>& g) {
    cout << "grids = [";
    for (size_t i = 0; i < g.size(); ++i) {
        if (i) cout << ", ";
        cout << "[";
        for (size_t j = 0; j < g[i].size(); ++j) {
            if (j) cout << ", ";
            cout << std::fixed << setprecision(3) << g[i][j];
        }
        cout << "]";
    }
    cout << "]\n";
}
static void print_blocks(const vector<double>& b) {
    cout << "blocks = [";
    for (size_t i = 0; i < b.size(); ++i) {
        if (i) cout << ", ";
        cout << std::fixed << setprecision(3) << b[i];
    }
    cout << "]\n";
}

/*** ---------- Random generator ---------- ***/
static double rnd_uniform(mt19937_64& rng, double a, double b) {
    uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}
static vector<vector<double>> make_grids(int M, int N, mt19937_64& rng, double lo, double hi) {
    vector<vector<double>> g(M, vector<double>(N));
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            g[i][j] = rnd_uniform(rng, lo, hi);
    return g;
}
static vector<double> make_blocks(int K, mt19937_64& rng, double lo, double hi) {
    vector<double> b(K);
    for (int i = 0; i < K; ++i) b[i] = rnd_uniform(rng, lo, hi);
    return b;
}

/*** ---------- Built-in cases ---------- ***/
static void built_in_cases() {
    Solution sol;

    auto run_and_print = [&](const vector<vector<double>>& grids,
                             const vector<double>& blocks,
                             const string& title) {
        cout << "==== " << title << " ====\n";
        print_grids(grids);
        print_blocks(blocks);
        double ans = sol.computeMinVariance(grids, blocks);
        cout << "expected = " << fmt_scientific(ans) << "\n\n";
    };

    // 1) Problem sample
    run_and_print({{-0.03, 0.02}}, {-0.01, 0.00, 0.01}, "Sample (M=1,N=2,K=3)");

    // 2) All zeros -> variance 0
    run_and_print({{0.0, 0.0},{0.0, 0.0}}, {0.0,0.0,0.0,0.0}, "All zeros (M=2,N=2,K=4)");

    // 3) K==MN with perfect cancel pairs -> can reach 0
    run_and_print({{0.010, -0.010}}, {-0.010, 0.010}, "Perfect cancel (M=1,N=2,K=2)");

    // 4) K>MN, mixed values
    run_and_print({{-0.020, 0.015},{0.005, -0.012}},
                  {-0.030, -0.005, 0.000, 0.008, 0.020},
                  "Mixed small values (M=2,N=2,K=5)");

    // 5) Rectangular 2x3, near-zero noise
    run_and_print({{-0.004, 0.006, -0.002},{0.001, -0.003, 0.002}},
                  {-0.007, -0.001, 0.000, 0.001, 0.003, 0.006, 0.009},
                  "Near-zero noise (M=2,N=3,K=7)");

    // 6) 3x3, wider range
    run_and_print({{-0.030, 0.025, -0.015},
                   { 0.010,-0.008,  0.005},
                   { 0.000, 0.012, -0.020}},
                  {-0.040, -0.020, -0.010, -0.005, 0.000, 0.006, 0.015, 0.022, 0.030, 0.040, 0.050},
                  "Wide range (M=3,N=3,K=11)");
}

/*** ---------- Custom random cases from stdin ---------- ***/
static void custom_cases_from_stdin() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    if (!(cin >> T)) {
        // If no T given, fall back to built-in
        built_in_cases();
        return;
    }
    Solution sol;
    for (int tc = 1; tc <= T; ++tc) {
        int M, N, K;
        double gmin, gmax, bmin, bmax;
        uint64_t seed;
        cin >> M >> N >> K >> gmin >> gmax >> bmin >> bmax >> seed;
        if (K < M * N) K = M * N;

        mt19937_64 rng(seed);
        auto grids = make_grids(M, N, rng, gmin, gmax);
        auto blocks = make_blocks(K, rng, bmin, bmax);

        cout << "==== Case #" << tc << " ====\n";
        print_grids(grids);
        print_blocks(blocks);
        double ans = sol.computeMinVariance(grids, blocks);
        cout << "expected = " << fmt_scientific(ans) << "\n\n";
    }
}

int main() {
    // If stdin has content starting with integer T, use custom mode; else use built-ins.
    // We peek one character: if it's a digit or '-', assume custom mode and pass to parser.
    int c = cin.peek();
    if (c == EOF) { built_in_cases(); return 0; }
    if (isdigit(c) || c=='-' || c=='+') { custom_cases_from_stdin(); }
    else { built_in_cases(); }
    return 0;
}
