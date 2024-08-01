#include <city.h>
#include "Fast.hpp"

class BitString {
    // static const Uchar UCHAR_MASK = 0xFF;
    // static const Ulong ALL_MASK = 0xFFFFFFFFFFFFFFFFLL;
    // std::vector<Ulong> bits_;

   public:
    const std::string *text_;
    Uint size_;

    inline BitString() : text_(nullptr), size_(0) {}

    inline BitString(const std::string &text) : text_(&text) {
        size_ = text.size() == 0 ? 0 : text.size() * CHAR_SIZE;
    }

    inline BitString(const BitString &bitString)
        : text_(bitString.text_), size_(bitString.size()) {}

    inline BitString(const BitString &bitString, Uint size)
        : text_(bitString.text_), size_(size) {
        assert(size_ <= bitString.text_->size() * CHAR_SIZE);
    }

    inline Uint size() const { return size_; }

    inline bool at(const Uint &index) const {
        assert(text_ != nullptr);
        assert(index < size_);
        return (text_->at(index / 8) & (1 << (index % 8))) != 0;
    }

    inline static Uint getLCPLength(const BitString &a, const BitString &b) {
        if (a.text_ == nullptr || b.text_ == nullptr) {
            return 0;
        }

        Uint result = 0;
        Uint minLength = std::min(a.size(), b.size());
        for (size_t i = 0; i < minLength / 8; ++i) {
            char a_char = a.text_->at(i);
            char b_char = b.text_->at(i);
            if (a_char != b_char) {
                break;
            }
            result += 8;
        }
        if (result < minLength) {
            char a_char = a.text_->at(result / 8);
            char b_char = b.text_->at(result / 8);
            if (a_char == b_char) {
                result += 8;
            } else {
                result += __builtin_ffs(a_char ^ b_char) - 1;
            }
            if (minLength < result) {
                result = minLength;
            }
        }

        // TRACE(result);
        assert(result <= minLength);

        return result;
    }

    // is b prefix of a.
    inline static bool isPrefix(const BitString &a, const BitString &b) {
        return getLCPLength(a, b) == b.size();
    }

    // inline BitString *SubBitString(Uint size) const {
    //     BitString *result = new BitString(text_);
    //     result->size_ = size;
    //     return result;
    // }

    inline std::string toString() const {
        if (text_ == nullptr) {
            return std::string();
        }

        std::string result = text_->substr(0, size_ / 8);

        if ((size_ / 8) * 8 < size_) {
            std::string tail(9, '\0');
            for (size_t i = (size_ / 8) * 8; i < size_; i++) {
                tail[i % 8] = at(i) ? '1' : '0';
            }
            result += " + " + tail;
        }
        // TRACE(text_);

        return result;
    }

    inline std::string toBinaryString() const {
        if (text_ == nullptr) {
            return std::string();
        }

        std::string result(size_ + 1, '\0');
        for (size_t i = 0; i < size_; i++) {
            result[i] = at(i) ? '1' : '0';
        }
        return result;
    }

    BitString &operator=(const BitString &a) {
        text_ = a.text_;
        size_ = a.size_;

        return (*this);
    }

    inline bool operator==(const BitString &a) const {
        if (size() != a.size()) {
            return false;
        }
        return getLCPLength(*this, a) == size();
    }

    inline bool operator!=(const BitString &a) const { return !(*this == a); }

    struct Hash {
        inline Ulong operator()(const BitString &key) const {
            if (key.text_ == nullptr) {
                return 0;
            }

            assert(key.text_ != nullptr);
            // assert(key.text_->size() != 0);
            // assert(key.size() != 0);
            Ulong result = 0;
            if (0 < key.size() / 8) {
                result = CityHash64(key.text_->data(), key.size() / 8);
            }
            if (0 < key.size() % 8) {
                char tail = key.text_->at(key.size() / 8)
                            << (8 - (key.size() % 8));
                result ^= CityHash64(&tail, 1);
            }

            return result;
        }
    };

    struct Equal {
        inline bool operator()(const BitString &a, const BitString &b) {
            return a == b;
        }
    };
};

namespace std {
ostream &operator<<(ostream &os, const BitString &bitString) {
    os << bitString.toString();
    return os;
}
}  // namespace std
