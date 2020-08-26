#ifndef LOKI_CPP_MATRIX2PICTURE_H
#define LOKI_CPP_MATRIX2PICTURE_H

#include "LoKI-B/LinearAlgebra.h"
#include <iostream>
#include <fstream>

namespace loki {

/** \todo See if there something else that can be created
 *        this easily and more widely supported, instead
 *        of xpm. 
 */
inline void writeXPM(const Matrix& m, std::ostream& os)
{
    os << "/* XPM */" << std::endl;
    os << "char* const name[] = {" << std::endl;
    // height, width, number of colors, number of chars per pixel
    os << '"' << m.rows() << ' ' << m.cols() << " 3 1\"" << std::endl;
    // declare color code characters w (white), b (blue) anr r (red)
    os << "\"w c #ffffff\"" << std::endl;
    os << "\"b c #0000ff\"" << std::endl;
    os << "\"r c #ff0000\"" << std::endl;
    for (Matrix::Index r=0; r!= m.rows(); ++r)
    {
        os << '"';
        for (Matrix::Index c=0; c!= m.cols(); ++c)
        {
            // <0: red
            // =0: white
            // >0: blue
            const double v = m(r,c);
            os << (v>0 ? 'b' : v<0 ? 'r' : 'w');
        }
        os << '"';
        if (r!=m.rows()-1)
            os << ',';
        os << '\n';
    }
    os << "};\n";
}

inline void writeXPM(const Matrix& m, std::string f)
{
    std::ofstream os(f);
    if (!os)
    {
        throw std::runtime_error("Could not open file '" + f + "' for writing.");
    }
    writeXPM(m,os);
}

} // namespace loki

#endif // LOKI_CPP_MATRIX2PICTURE_H
