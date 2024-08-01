#include <chrono>
#include <fstream>
#include <iostream>

#undef NDEBUG

#include "CTriePP.hpp"

using namespace ctriepp;

void testTrie(const std::vector<std::string>& texts, bool is_inverse = false) {
    CTriePP<bool, false> tree;

    for (size_t i = 0; i < texts.size(); i++) {
        tree.insert(&(texts[i]), true);
        tree.print();
    }
    tree.print();

    for (size_t i = 0; i < texts.size(); i++) {
        assert(tree.contains(texts[i]));
    }

    for (size_t i = 0; i < texts.size(); i++) {
        for (size_t j = 0; j < texts[i].size(); j++) {
            std::string prefix = texts[i].substr(0, j);
            assert(tree.containsPrefix(prefix));
        }
    }

    for (size_t i = 0; i < texts.size(); i++) {
        for (size_t j = 0; j < texts[i].size(); j++) {
            std::string prefix = texts[i].substr(0, j) + "X";
            assert(!tree.contains(prefix));
            assert(!tree.containsPrefix(prefix));
        }
    }

    LOG(texts.size());

    for (size_t i = 0; i < texts.size(); i++) {
        size_t index;
        if (is_inverse) {
            index = texts.size() - i - 1;
        } else {
            index = i;
        }
        tree.erase(texts[index]);
        assert(!tree.contains(texts[index]));
        tree.print();
    }
}

// key=\"journals/corr/cs-CY-0109091\">
int main(int argc, char** argv) {
    assert(LongString::charSize(0) == 0);
    assert(LongString::charSize(1) == 1);
    char char_m1 = -1;
    Ulong long_char_m1 = 0x00000000000000FFLL;
    assert(((Ulong)(Uchar)char_m1) == long_char_m1);

    assert(LongString::getCharLCPLength(LongString::toLong("abc"),
                                        LongString::toLong("abc")) == 3);
    assert(LongString::getCharLCPLength(LongString::toLong("abc"),
                                        LongString::toLong("abcx")) == 3);
    assert(LongString::getCharLCPLength(LongString::toLong("abcx"),
                                        LongString::toLong("abc")) == 3);
    assert(LongString::getCharLCPLength(LongString::toLong("abcd"),
                                        LongString::toLong("abcx")) == 3);
    assert(LongString::getCharLCPLength(LongString::toLong("abcdefgh"),
                                        LongString::toLong("abcdefgx")) == 7);
    assert(LongString::getCharLCPLength(
               LongString::toSubLong(LongString::toLong("00000001"), 0, 7),
               LongString::toLong("00000000")) == 7);

    assert(LongString::toSubLong(1LL | (2LL << 8) | (3LL << 16), 0, 3) ==
           (1LL | (2LL << 8) | (3LL << 16)));
    assert(LongString::toSubLong(1LL | (2LL << 8) | (3LL << 16), 1, 3) ==
           (2LL | (3LL << 8)));
    assert(LongString::toSubLong(1LL | (2LL << 8) | (3LL << 16), 0, 2) ==
           (1LL | (2LL << 8)));
    assert(LongString::toSubLong(1LL | (2LL << 8) | (3LL << 16), 1, 2) == 2LL);
    assert(LongString::toSubLong(3544385890265608240LL, 0, 8) ==
           3544385890265608240LL);

    assert(Fast::mostSignificantBit(0) == -1);

    assert(Fast::twoFattest(0, 8) == 8);
    assert(Fast::twoFattest(1, 8) == 8);
    assert(Fast::twoFattest(0, 9) == 8);
    assert(Fast::twoFattest(0, 4) == 4);
    assert(Fast::twoFattest(0, 7) == 4);
    assert(Fast::twoFattest(4, 7) == 6);
    assert(Fast::twoFattest(-1, 8) == 0);
    assert(Fast::twoFattest(8, 8) == 0);

    Factory<bool> factory;
    assert(factory.size() == 0);
    assert(factory.make() == 0);
    assert(factory.size() == 1);
    factory.erase(0);
    assert(factory.size() == 0);
    assert(factory.make() == 0);
    assert(factory.size() == 1);
    assert(factory.make() == 1);
    assert(factory.size() == 2);
    assert(factory.make() == 2);
    assert(factory.size() == 3);
    factory.erase(1);
    assert(factory.size() == 2);
    assert(factory.make() == 1);
    assert(factory.size() == 3);

    //     AlphabetAwareZFastTrie<Uint> tFastTrie;
    //     tFastTrie.add(1LL, 1);
    //     tFastTrie.add(1LL | (1LL << 8), 2);
    //     tFastTrie.add(1LL | (2LL << 8) | (3LL << 16), 3);
    //     tFastTrie.add(1LL | (2LL << 8) | (4LL << 16), 4);
    //     assert(tFastTrie.containsPrefix(1LL | (1LL << 8)));
    //     assert(tFastTrie.containsPrefix(1LL | (2LL << 8)));
    //     assert(tFastTrie.containsPrefix(1LL | (2LL << 8) | (4LL << 16)));
    //     assert(tFastTrie.containsPrefix(1LL | (2LL << 8) | (3LL << 16)));
    //     tFastTrie.add(3544385890265608240LL);

    std::vector<std::string> texts;
    texts.push_back("A");
    texts.push_back("AA");
    texts.push_back("AAA");
    texts.push_back("AAAA");
    texts.push_back("AAAAA");
    texts.push_back("AAAAAA");
    texts.push_back("AAAAAAA");
    texts.push_back("AAAAAAAA");
    texts.push_back("B");
    testTrie(texts);
    testTrie(texts, true);

    texts.clear();
    texts.push_back("0000");
    texts.push_back("0000000100000011");
    texts.push_back("00000002");
    texts.push_back("0000000100000012");
    texts.push_back("000000010000001100000021");
    testTrie(texts);
    testTrie(texts, true);

    texts.clear();
    texts.push_back("0000000100000011");
    texts.push_back("00000001");
    texts.push_back("000000010000");
    texts.push_back("00000");
    testTrie(texts);
    testTrie(texts, true);

    texts.clear();
    texts.push_back("0a");
    texts.push_back("0b");
    texts.push_back("0c");
    texts.push_back("0d");
    texts.push_back("0e");
    texts.push_back("0f");
    testTrie(texts);
    testTrie(texts, true);

    texts.clear();
    texts.push_back("aa");
    texts.push_back("a");
    texts.push_back("b");
    texts.push_back("ba");
    testTrie(texts);
    testTrie(texts, true);

    LOG("Test Success");
    return 0;
}
