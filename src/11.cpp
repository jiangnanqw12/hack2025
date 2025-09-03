#include <vector>
#include <unordered_map>
#include <cmath>
#include <string>
#include <algorithm>

class Solution {
private:
    // 使用哈希表来记忆化递归结果，避免重复计算
    std::unordered_map<std::string, double> memo;
    
    // 辅助函数：计算给定配置下的期望移动次数
    double solve(std::vector<int>& stack1, std::vector<int>& stack2, 
                const std::vector<int>& sequence, int currentIndex) {
        if (currentIndex == sequence.size()) {
            return 0.0; // 所有数据都处理完毕
        }
        
        // 创建状态字符串作为哈希表的键
        std::string state = "";
        for (int x : stack1) state += std::to_string(x) + ",";
        state += "|";
        for (int x : stack2) state += std::to_string(x) + ",";
        state += "|" + std::to_string(currentIndex);
        
        // 如果已经计算过这个状态，直接返回记忆化的结果
        if (memo.count(state)) return memo[state];
        
        int target = sequence[currentIndex];
        double result = 0.0;
        
        // 检查目标数字是否在stack1的顶部
        if (!stack1.empty() && stack1.back() == target) {
            stack1.pop_back();
            result = solve(stack1, stack2, sequence, currentIndex + 1);
            stack1.push_back(target);
        }
        // 检查目标数字是否在stack2的顶部
        else if (!stack2.empty() && stack2.back() == target) {
            stack2.pop_back();
            result = solve(stack1, stack2, sequence, currentIndex + 1);
            stack2.push_back(target);
        }
        // 如果目标数字在stack1中但不在顶部
        else {
            std::vector<int> temp1 = stack1;
            std::vector<int> temp2 = stack2;
            double moves1 = 0.0, moves2 = 0.0;
            bool found1 = false, found2 = false;
            
            // 在stack1中查找并计算移动次数
            for (int i = temp1.size() - 1; i >= 0; i--) {
                if (temp1[i] == target) {
                    found1 = true;
                    // 移动上面的数字到stack2
                    int movesCount = temp1.size() - i - 1;
                    while (temp1.size() > i) {
                        temp2.push_back(temp1.back());
                        temp1.pop_back();
                    }
                    moves1 = movesCount + solve(temp1, temp2, sequence, currentIndex + 1);
                    break;
                }
            }
            
            temp1 = stack1;
            temp2 = stack2;
            
            // 在stack2中查找并计算移动次数
            for (int i = temp2.size() - 1; i >= 0; i--) {
                if (temp2[i] == target) {
                    found2 = true;
                    // 移动上面的数字到stack1
                    int movesCount = temp2.size() - i - 1;
                    while (temp2.size() > i) {
                        temp1.push_back(temp2.back());
                        temp2.pop_back();
                    }
                    moves2 = movesCount + solve(temp1, temp2, sequence, currentIndex + 1);
                    break;
                }
            }
            
            if (found1 && found2) {
                result = std::min(moves1, moves2);
            } else if (found1) {
                result = moves1;
            } else if (found2) {
                result = moves2;
            }
        }
        
        memo[state] = result;
        return result;
    }
    
public:
    double expectedMoveSteps(int n, const std::vector<int>& sequence) {
        double total = 0.0;
        int totalConfigurations = 1 << n; // 2^n种可能的配置
        
        // 枚举所有可能的初始配置
        for (int mask = 0; mask < totalConfigurations; mask++) {
            std::vector<int> stack1, stack2;
            
            // 根据二进制位分配数字到两个栈中
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
        
        // 返回期望值（所有配置的平均值）
        return total / totalConfigurations;
    }
};