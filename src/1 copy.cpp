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
        if (n <= 1) return 0.0; 

        sort(G.begin(), G.end());                   
        vector<double> B = blocks;
        sort(B.begin(), B.end(), greater<double>());

        
        double sum_g = 0.0, sum_g2 = 0.0;
        double gmin = 1e300, gmax = -1e300;
        for (double g : G) {
            sum_g += g; sum_g2 += g * g;
            gmin = min(gmin, g); gmax = max(gmax, g);
        }
        double bmin = 1e300, bmax = -1e300;
        for (double b : B) { bmin = min(bmin, b); bmax = max(bmax, b); }

        
        double lo = gmin + bmin, hi = gmax + bmax;

        auto eval_at_t = [&](double t) -> double {
            
            double const_term = sum_g2 - 2.0 * t * sum_g + n * t * t;

            
            vector<double> w(n);
            for (int k = 0; k < n; ++k) w[k] = 2.0 * (G[k] - t);

            
            const double INF = 1e300;
            vector<double> dp(n + 1, INF), next(n + 1, INF);
            dp[0] = 0.0;

            for (double b : B) {
                
                for (int k = 0; k <= n; ++k) next[k] = dp[k];
                
                for (int k = 0; k < n; ++k) {
                    if (dp[k] >= INF/2) continue;
                    double cand = dp[k] + b * b + w[k] * b;
                    if (cand < next[k + 1]) next[k + 1] = cand;
                }
                dp.swap(next);
            }

            double best = dp[n]; 
            return (const_term + best) / n; 
        };

        
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

