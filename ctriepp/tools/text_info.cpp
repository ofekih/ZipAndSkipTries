#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <cmath>

#define LOG(message) std::cerr << "[DEBUG] " << message << std::endl
#define TRACE(message)                                     \
    std::cerr << "[TRACE] " << __FILE__ << "/" << __LINE__ \
              << "/ " #message " : " << message << std::endl
#define COUT(message) std::cout << message << std::endl

int main(int argc, char** argv) {
    if (argc < 2) {
        exit(1);
    }

    std::string text_path = argv[1];
    COUT(text_path);

    std::ifstream file(argv[1]);
    if (!file.is_open()) {
        throw std::runtime_error("Couldn't read " + text_path + ".");
    }

    // std::vector<std::string> lines;
    std::string line;
    std::set<char> char_set;
    std::set<std::string> string_set;
    long num_string = 0;
    long total_length = 0;
    long max_length = 0;
    long total_lcpLength = 0;
    long max_lcpLength = 0;
    std::string prev;

    std::map<long, long> length2count;
    std::map<long, long> lcpLength2count;

    while (std::getline(file, line)) {
        total_length += line.size();
        num_string++;
        string_set.insert(line);
        long length = line.size();
        max_length = std::max(length, max_length);
        for (size_t i = 0; i < line.size(); i++) {
            char_set.insert(line[i]);
        }

        {
            int pow = (length == 0)? 0 : std::ceil(std::log2(length)) + 1;
            if(length2count.count(pow) == 0){
                length2count[pow] = 1;
            }else{
                ++length2count[pow];
            }
        }

        if(!prev.empty()) {
            long lcpLength = 0;
            long minLength = std::min(prev.size(), line.size());
            for (int i = 0; i < minLength; ++i){
               if (prev[i] == line[i]){
                    ++lcpLength;
                } else {
                    break;
                }
            }
            total_lcpLength+= lcpLength;
            max_lcpLength = std::max(max_lcpLength, lcpLength);

            int pow = (lcpLength == 0)? 0 : std::ceil(std::log2(lcpLength)) + 1;
            if(lcpLength2count.count(pow) == 0){
                lcpLength2count[pow] = 1;
            }else{
                ++lcpLength2count[pow];
            }
        }
        prev = line;
    }

    TRACE(total_length);
    TRACE(num_string);
    TRACE((double)total_length / num_string);
    TRACE(max_length);
    TRACE(char_set.size());
    TRACE(string_set.size());
    TRACE(total_lcpLength);
    TRACE(max_lcpLength);
    TRACE((double)total_lcpLength / num_string);

    for (int i = 0; i < length2count.size(); i++) {
            std::cout << (i == 0 ? 0 : (1 << (i - 1))) << " " << length2count[i] << std::endl;
    }

    for (int i = 0; i < lcpLength2count.size(); i++) {
            std::cout << (i == 0 ? 0 : (1 << (i - 1))) << " " << lcpLength2count[i] << std::endl;
    }
}
