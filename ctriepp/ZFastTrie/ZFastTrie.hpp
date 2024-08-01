#include "BitString.hpp"

template <typename Value, const Value EMPTY_VALUE = -1>
class ZFastTrie {
    class Node {
       public:
        // typedef rigtorp::HashMap<Uchar, Node *, Fast::UcharHash,
        //                          Fast::UcharEqual>
        //     Char2NodeMap;
        // typedef std::unordered_map<Uchar, Node *> Char2NodeMap;
        // typedef std::map<Uchar, Node *> Char2NodeMap;

        Value value_;
        BitString extent_;
        Int nameLength_;
        Node *leftChild_;
        Node *rightChild_;

        inline Node(const Value &value, const BitString &extent)
            : value_(value),
              extent_(extent),
              nameLength_(1),
              leftChild_(nullptr),
              rightChild_(nullptr) {}

        inline Node(const Value &value, const BitString &extent,
                    const Int &nameLength)
            : value_(value),
              extent_(extent),
              nameLength_(nameLength),
              leftChild_(nullptr),
              rightChild_(nullptr) {}

        inline void set(const Value &value, const BitString &extent,
                        const Int &nameLength) {
            value_ = value;
            extent_ = extent;
            nameLength_ = nameLength;
        }

        inline void setExtent(const BitString &extent) {
            if (extent.size() == 0) {
                nameLength_ = 0;
            }
            assert(nameLength_ <= extent.size());
            extent_ = extent;
        }
        inline Uint handleLength() const {
            return Fast::twoFattest(nameLength_ - 1, extentLength());
        }
        inline BitString handle() const {
            return extentLength() <= 0 ? BitString()
                                       : BitString(extent_, handleLength());
        }
        inline BitString extent() const { return extent_; }

        inline Uint extentLength() const { return extent_.size(); }

        inline bool isLeaf() const {
            return leftChild_ == nullptr && rightChild_ == nullptr;
        }

        inline bool key() const { return extent_.at(nameLength_ - 1); }

        inline Uint sizeChildren() const {
            Uint size = 0;
            if (leftChild_ != nullptr) {
                ++size;
            }
            if (rightChild_ != nullptr) {
                ++size;
            }
            return size;
        }

        inline void insertChild(Node *child, Uint lcpLength) {
            // Uint lcpLength = BitString::getLCPLength(extent(),
            // child->extent());
            if (child->extent_.at(lcpLength)) {
                rightChild_ = child;
            } else {
                leftChild_ = child;
            }
        }

        inline void eraseChild(bool key) {
            if (key) {
                rightChild_ = nullptr;
            } else {
                leftChild_ = nullptr;
            }
        }

        inline Node *getChild() {
            if (leftChild_ != nullptr) {
                return leftChild_;
            }
            if (rightChild_ != nullptr) {
                return rightChild_;
            }
            return nullptr;
        }

        inline void print(const std::string indent = "") const {
            INFO(indent << "ZFastTrie::Node \"" << extent().toString()
                        << "\", value: " << value_
                        << ", nameLength: " << nameLength_);
            if (leftChild_ != nullptr) {
                INFO(indent << "  children Left");
                assert(extentLength() + 1 == leftChild_->nameLength_);
                leftChild_->print(indent + "    ");
            }
            if (rightChild_ != nullptr) {
                INFO(indent << "  children Right");
                assert(extentLength() + 1 == rightChild_->nameLength_);
                rightChild_->print(indent + "    ");
            }
        }
    };

    typedef rigtorp::HashMap<BitString, Node *, BitString::Hash,
                             BitString::Equal>
        Handle2NodeMap;

    Int size_;

    Node *root_;

    Handle2NodeMap handle2NodeMap_;

   public:
    inline ZFastTrie()
        : size_(0), root_(nullptr), handle2NodeMap_(16, BitString()) {}

    inline void insert(const std::string *newText, const Value &value) {
        insert(BitString(*newText), value);
    }
    inline void insert(const BitString &newText, const Value &value) {
        LOG("ZFastTrie::insert(" << newText << ", value: " << value << ")");
        // assert(!contains(newText));
        assert(value != EMPTY_VALUE);

        if (size_ == 0) {
            root_ = new Node(value, newText);
        } else {
            Node *exitNode = getExitNode(newText);
            assert(exitNode != nullptr);
            Int lcpLength =
                BitString::getLCPLength(exitNode->extent(), newText);

            if (lcpLength < exitNode->extentLength()) {
                const BitString &exitNode_extent = exitNode->extent();

                const BitString &newExtent =
                    (lcpLength == 0) ? BitString()
                                     : BitString(exitNode->extent(), lcpLength);
                if (exitNode->isLeaf()) {
                    exitNode->setExtent(newExtent);
                    insertHandle2NodeMap(exitNode);
                } else {
                    const uint &newHandleLength =
                        Fast::twoFattest(lcpLength, newExtent.size());
                    bool isChangeHandle =
                        (exitNode->handleLength() != newHandleLength);
                    if (isChangeHandle) {
                        eraseHandle2NodeMap(exitNode->handle());
                        exitNode->setExtent(newExtent);
                        insertHandle2NodeMap(exitNode);
                    } else {
                        exitNode->setExtent(newExtent);
                    }
                }

                // make new internal node
                // swap new internal node and exitNode
                Node *newNode =
                    new Node(exitNode->value_, exitNode_extent, lcpLength + 1);
                swapChildren(exitNode, newNode);

                exitNode->insertChild(newNode, lcpLength);

                if (!newNode->isLeaf()) {
                    insertHandle2NodeMap(newNode);
                }

                if (lcpLength < newText.size()) {
                    // make new leaf node
                    Node *newTextNode = new Node(value, newText, lcpLength + 1);
                    exitNode->value_ = EMPTY_VALUE;
                    exitNode->insertChild(newTextNode, lcpLength);
                    newTextNode->value_ = value;
                } else {
                    exitNode->value_ = value;
                }
            } else {
                if (lcpLength == newText.size()) {
                    if (exitNode->value_ == EMPTY_VALUE) {
                        exitNode->value_ = value;
                    } else {
                        WARN("new text already exist.");
                    }
                } else {
                    if (exitNode->isLeaf()) {
                        insertHandle2NodeMap(exitNode);
                    }
                    Node *newTextNode = new Node(value, newText, lcpLength + 1);
                    exitNode->insertChild(newTextNode, lcpLength);
                    newTextNode->value_ = value;
                }
            }
        }

        ++size_;
        if (!contains(newText)) {
            print();
        }
        assert(handle2NodeMap_.size() < size_);
        assert(contains(newText));
        assert(containsPrefix(newText));
    }

    inline void erase(const std::string &targetText) {
        erase(BitString(targetText));
    }

    inline void erase(const BitString &targetText) {
        LOG("ZFastTrie::erase " << targetText.toString());
        // print();
        assert(contains(targetText));
        Node *targetNode = getExitNode(targetText);
        assert(targetNode != nullptr);
        if (size_ == 0) {
            root_ = nullptr;
        } else if (size_ <= 1) {
            root_ = nullptr;
            eraseHandle2NodeMap(targetNode->handle());
        } else {
            Int lcpLength =
                BitString::getLCPLength(targetNode->extent(), targetText);
            assert(lcpLength == targetText.size());
            assert(targetText.size() <= targetNode->extentLength());
            if (targetNode->isLeaf()) {
                if (targetNode == root_) {
                    eraseHandle2NodeMap(targetNode->handle());
                } else {
                    Node *parentNode;
                    if (targetNode->nameLength_ <= 1) {
                        parentNode = root_;
                    } else {
                        parentNode = getExitNode(BitString(
                            targetNode->extent(), targetNode->nameLength_ - 1));
                    }

                    parentNode->eraseChild(targetNode->key());
                    eraseHandle2NodeMap(targetNode->handle());

                    if (parentNode->isLeaf()) {
                        eraseHandle2NodeMap(parentNode->handle());
                    } else if (parentNode->sizeChildren() == 1 &&
                               parentNode->value_ == EMPTY_VALUE) {
                        LOG("swap parent and child node");
                        // swap parent and child node
                        eraseHandle2NodeMap(parentNode->handle());
                        Node *childNode = parentNode->getChild();

                        parentNode->eraseChild(childNode->key());
                        swapChildren(parentNode, childNode);

                        parentNode->set(childNode->value_, childNode->extent(),
                                        parentNode->nameLength_);

                        if (!parentNode->isLeaf()) {
                            eraseHandle2NodeMap(childNode->handle());
                            insertHandle2NodeMap(parentNode);
                        }
                        // delete parentNode
                    }
                }
            } else if (targetNode->sizeChildren() == 1) {
                // delete internal node
                Node *childNode = targetNode->getChild();
                targetNode->eraseChild(childNode->key());
                eraseHandle2NodeMap(targetNode->handle());
                targetNode->set(childNode->value_, childNode->extent(),
                                targetNode->nameLength_);

                swapChildren(targetNode, childNode);
                if (!targetNode->isLeaf()) {
                    eraseHandle2NodeMap(childNode->handle());
                    insertHandle2NodeMap(targetNode);
                }
            } else {
                targetNode->value_ = EMPTY_VALUE;
            }
        }
        --size_;
        if (contains(targetText)) {
            print();
        }
        assert(!contains(targetText));
    }

    inline void update(const BitString &targetText, Value value) {
        Node *exitNode = getExitNode(targetText);
        if (exitNode == nullptr) {
            WARN("ZFastTrie::update");
            WARN("Not Found key");
        } else {
            Value &exitNodeValue = exitNode->value_;
            if (exitNodeValue == EMPTY_VALUE) {
                WARN("ZFastTrie::update");
                WARN("Not Found key");
            } else {
                exitNodeValue = value;
            }
        }
    }

    inline void insertHandle2NodeMap(Node *node) {
        if (node->extentLength() == 0) {
            return;
        }
        if (node->handle().size() == 0) {
            return;
        }
        const BitString &handle = node->handle();
        assert(handle2NodeMap_.count(handle) == 0);
        handle2NodeMap_[handle] = node;
        assert(handle2NodeMap_.count(node->handle()) != 0);
    }

    inline void eraseHandle2NodeMap(const BitString &handle) {
        if (handle.size() == 0) {
            return;
        }
        // assert(handle2NodeMap_.count(hangit dle) != 0);
        handle2NodeMap_.erase(handle);
        // assert(handle2NodeMap_.count(handle) == 0);
    }

    inline bool contains(const std::string &pattern) const {
        return contains(BitString(pattern));
    }

    inline bool contains(const BitString &pattern) const {
        Node *exitNode = getExitNode(pattern);
        if (exitNode == nullptr) {
            return false;
        } else {
            return (exitNode->extent() == pattern) &&
                   (exitNode->value_ != EMPTY_VALUE);
        }
    }

    inline bool containsPrefix(const std::string &pattern) const {
        return containsPrefix(BitString(pattern));
    }

    inline bool containsPrefix(const BitString &pattern) const {
        Node *node = getExitNode(pattern);
        return node == nullptr ? false
                               : BitString::isPrefix(node->extent(), pattern);
    }

    inline Node *getExitNode(const BitString &pattern) const {
        // LOG("getExitNode of " << pattern.toString());

        Int pattern_length = pattern.size();
        Int a = 0;
        Int b = pattern_length;
        Int f;
        Node *node = nullptr;
        Node *result = root_;

        while (0 < b - a) {
            f = Fast::twoFattest(a, b);
            // TRACE(BitString::Hash()(node->handle()));
            node = getNode(BitString(pattern, f));

            if (node != nullptr) {
                a = node->extentLength();
                assert(f <= a);
                result = node;
            } else {
                b = f - 1;
            }
        }

        if (result != nullptr) {
            Uint lcpLength = BitString::getLCPLength(result->extent(), pattern);
            if (lcpLength == result->extentLength() &&
                lcpLength < pattern_length) {
                Node *next = pattern.at(lcpLength) ? result->rightChild_
                                                   : result->leftChild_;
                if (next != nullptr) {
                    result = next;
                }
            }
        }

        return result;
    }

    inline Node *getNode(const BitString &handle) const {
        auto iter = handle2NodeMap_.find(handle);
        return iter == handle2NodeMap_.end() ? nullptr : iter->second;
    }

    inline static void swapChildren(Node *a, Node *b) {
        std::swap(a->leftChild_, b->leftChild_);
        std::swap(a->rightChild_, b->rightChild_);
    }

    inline void print(const std::string indent = "") const {
        INFO(indent << "print ZFastTrie");
        auto iter = handle2NodeMap_.begin();
        for (; iter != handle2NodeMap_.end(); ++iter) {
            INFO(indent << "    \"" << iter->first << "\"\t-> \""
                        << iter->second->extent() << "\", "
                        << iter->first.size() << " -> "
                        << iter->second->extent().size());
        }
        if (root_ != nullptr) {
            root_->print(indent);
        }
    }
};
