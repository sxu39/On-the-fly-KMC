#ifndef PARA_H
#define PARA_H

class Para
{
    public:
        Para();
        ~Para();
    public:
         int yes_count;
	 double upper_pos;
	 double upper_devi;
    public:
         int get_yes_count();
	 double get_upper_pos();
	 double get_upper_devi();
};
#endif
