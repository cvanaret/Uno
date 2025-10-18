// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MXSTRUCT_H
#define UNO_MXSTRUCT_H

#include <string>
#include <map>
#include "mex.h"

namespace uno {

    class MxStruct {
        
    public:
        mxArray* operator[](const std::string& key) const;
        void insert(const std::string& key, mxArray* value);
        bool contains(const std::string& key) const;
        auto begin() const { return this->data.begin(); }
        auto end() const { return this->data.end(); }

    private:
        std::map<std::string, mxArray*> data;
    };
    
}; // namespace


#endif // UNO_MXSTRUCT_H