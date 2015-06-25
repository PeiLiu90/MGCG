#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

int mod(const int & a, const int & b)
{
    if(a%b<0)
    {
        return a%b+b;
    }
    else
    {
        return a%b;
    }
}

#endif // TOOLS_H_INCLUDED
