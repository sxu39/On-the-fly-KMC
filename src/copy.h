#ifndef COPY_H
#define COPY_H

#include <string>

using namespace std;

class Copy
{
    public:
        void copy_file(string src_path, string dst_path);
        void move_file(string src_path, string dst_path);
};
#endif
