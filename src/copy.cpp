#include "copy.h"

using namespace std;

void Copy::copy_file(string src_path, string dst_path)
{
    string cmd_copy = "cp -rf " + src_path + "   " + dst_path;
    (void)system(cmd_copy.c_str());
}

void Copy::move_file(string src_path, string dst_path)
{
    string cmd_mv = "mv " + src_path + "   " + dst_path;
    (void)system(cmd_mv.c_str());
}
