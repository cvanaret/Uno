#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>

#define RED    "\x1B[31m"
#define GREEN   "\x1B[32m"
#define YELLOW  "\x1B[33m"
#define BLUE    "\x1B[34m"
#define MAGENTA "\x1B[35m"
#define CYAN    "\x1B[36m"
#define WHITE   "\x1B[37m"
#define RESET "\x1B[0m"

enum Level {
    ERROR = 0, WARNING, INFO, DEBUG
};

class Logger {
public:
    static Level logger_level;
};

template <typename T>
const Level& operator<<(const Level& level, T& element) {
    if (level <= Logger::logger_level) {
        std::cout << element;
    }
    return level;
}

template <typename T>
const Level& operator<<(const Level& level, const T& element) {
    if (level <= Logger::logger_level) {
        std::cout << element;
    }
    return level;
}

#endif // LOGGER_H
