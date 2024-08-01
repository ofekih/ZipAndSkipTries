#include <chrono>
#include <fstream>
#include <iostream>

#undef NDEBUG

#include "ZFastTrie.hpp"

void testTrie(const std::vector<std::string>& texts, bool is_inverse = false) {
    ZFastTrie<bool, false> tree;

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
    assert(Fast::mostSignificantBit(0) == -1);

    assert(Fast::twoFattest(0, 0) == 0);
    assert(Fast::twoFattest(0, 8) == 8);
    assert(Fast::twoFattest(1, 8) == 8);
    assert(Fast::twoFattest(0, 9) == 8);
    assert(Fast::twoFattest(0, 4) == 4);
    assert(Fast::twoFattest(0, 7) == 4);
    assert(Fast::twoFattest(4, 7) == 6);
    assert(Fast::twoFattest(-1, 8) == 0);
    assert(Fast::twoFattest(8, 8) == 0);

    //     assert(BitString::getLCPLength(BitString("a"), BitString("e")) ==
    //     2);

    TRACE(BitString("1").toBinaryString());
    TRACE(BitString("9").toBinaryString());
    assert(BitString::getLCPLength(BitString("1"), BitString("9")) == 3);

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
    texts.push_back("0");
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
