// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MxStruct.hpp"

namespace uno {

    mxArray* MxStruct::operator[](const std::string& key) const {
        auto it = this->data.find(key);
        return (it != this->data.end()) ? it->second : nullptr;
    }

    void MxStruct::insert(const std::string& key, mxArray* value) {
        this->data[key] = value;
    }

    bool MxStruct::contains(const std::string& key) const {
        return this->data.find(key) != this->data.end();
    }

} // namespace 