#include <vector>
#include <unordered_map>
#include <cmath>
#include <string>
#include <algorithm>

class Solution {
private:
    
    std::unordered_map<std::string, double> memo;
    
    
    double solve(std::vector<int>& stack1, std::vector<int>& stack2, 
                const std::vector<int>& sequence, int currentIndex) {
        if (currentIndex == sequence.size()) {
            return 0.0; 
        }
        
        
        std::string state = "";
        for (int x : stack1) state += std::to_string(x) + ",";
        state += "|";
        for (int x : stack2) state += std::to_string(x) + ",";
        state += "|" + std::to_string(currentIndex);
        
        
        if (memo.count(state)) return memo[state];
        
        int target = sequence[currentIndex];
        double result = 0.0;
if (!stack1.empty() && stack1.back() == target) {
    stack1.pop_back();
    result = 1.0 + solve(stack1, stack2, sequence, currentIndex + 1); 
    stack1.push_back(target);
}

else if (!stack2.empty() && stack2.back() == target) {
    stack2.pop_back();
    result = 1.0 + solve(stack1, stack2, sequence, currentIndex + 1); 
    stack2.push_back(target);
}
else {
    std::vector<int> temp1 = stack1, temp2 = stack2;
    double moves1 = 0.0, moves2 = 0.0;
    bool found1 = false, found2 = false;

    
    for (int i = (int)temp1.size() - 1; i >= 0; --i) {
        if (temp1[i] == target) {
            found1 = true;
            int movesCount = (int)temp1.size() - i - 1; 
            
            while ((int)temp1.size() - 1 > i) {
                temp2.push_back(temp1.back());
                temp1.pop_back();
            }
            
            temp1.pop_back();
            moves1 = movesCount + 1.0 + solve(temp1, temp2, sequence, currentIndex + 1);
            break;
        }
    }

    temp1 = stack1; temp2 = stack2;

    
    for (int i = (int)temp2.size() - 1; i >= 0; --i) {
        if (temp2[i] == target) {
            found2 = true;
            int movesCount = (int)temp2.size() - i - 1;
            while ((int)temp2.size() - 1 > i) {
                temp1.push_back(temp2.back());
                temp2.pop_back();
            }
            temp2.pop_back(); 
            moves2 = movesCount + 1.0 + solve(temp1, temp2, sequence, currentIndex + 1);
            break;
        }
    }

    
    if (found1 && found2) result = std::min(moves1, moves2);
    else if (found1)      result = moves1;
    else if (found2)      result = moves2;
}
        memo[state] = result;
        return result;
    }
    
public:
    double expectedMoveSteps(int n, const std::vector<int>& sequence) {
        double total = 0.0;
        int totalConfigurations = 1 << n; 
        
        
        for (int mask = 0; mask < totalConfigurations; mask++) {
            std::vector<int> stack1, stack2;
            
            
            for (int i = 0; i < n; i++) {
                if (mask & (1 << i)) {
                    stack1.push_back(i + 1);
                } else {
                    stack2.push_back(i + 1);
                }
            }
            
            memo.clear();
            double current = solve(stack1, stack2, sequence, 0);
            total += current;
        }
        
        
        return total / totalConfigurations;
    }
};