#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
using namespace std;

class Solution {
public:
    double computeMinVariance(const vector<vector<double>>& grids,
                              const vector<double>& blocks) {
        // Flatten grids to G and gather constants
        vector<double> G;
        for (const auto& row : grids) for (double x : row) G.push_back(x);
        const int n = (int)G.size();
        if (n <= 1) return 0.0; // sample variance with one point -> 0 here

        sort(G.begin(), G.end());                       // ascending
        vector<double> B = blocks;
        sort(B.begin(), B.end(), greater<double>());    // descending

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
            return (const_term + best) / (n - 1); // sample variance
        };

        // Ternary search on convex function
        for (int it = 0; it < 80; ++it) {
            double m1 = (2.0 * lo + hi) / 3.0;
            double m2 = (lo + 2.0 * hi) / 3.0;
            double f1 = eval_at_t(m1);
            double f2 = eval_at_t(m2);
            if (f1 < f2) hi = m2; else lo = m1;
        }
        double ans = eval_at_t((lo + hi) * 0.5);
        return ans;
    }
};
