#ifndef UNO_RANGE_H
#define UNO_RANGE_H

// https://en.wikipedia.org/wiki/Generator_(computer_programming)#C++
class Range {
public:
   explicit Range(size_t end_index);
   Range(size_t start_index, size_t end_index);

   // iterable functions
   [[nodiscard]] const Range& begin() const;
   [[nodiscard]] const Range& end() const;

   // iterator functions
   [[nodiscard]] bool operator!=(const Range&) const;
   void operator++();
   [[nodiscard]] size_t operator*() const;

protected:
   size_t end_index{};
   size_t current_index{};
};

inline Range::Range(size_t end_index): end_index(end_index), current_index(0) {
}

inline Range::Range(size_t start_index, size_t end_index): end_index(end_index), current_index(start_index) {
}

inline const Range& Range::begin() const {
   return *this;
}

inline const Range& Range::end() const {
   return *this;
}

inline bool Range::operator!=(const Range&) const {
   return (this->current_index < this->end_index);
}

inline void Range::operator++() {
   ++this->current_index;
}

inline size_t Range::operator*() const {
   return this->current_index;
}

#endif // UNO_RANGE_H