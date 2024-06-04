#include "MyPoint.h"

ostream& operator<<(ostream& os, const myPoint& pt)
{
    os << "("<<pt.px<<", "<<pt.py<<")";
    return os;
}

ostream& operator<<(ostream& os, const myRealPoint& pt)
{
    os << "("<<pt.px<<", "<<pt.py<<")";
    return os;
}
