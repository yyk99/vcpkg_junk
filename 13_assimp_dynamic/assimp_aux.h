#pragma once

namespace {
std::ostream &operator<<(std::ostream &ss, aiVector3D const &v) {
    ss << "{" << std::setprecision(9) << v.x << "," << std::setprecision(9)
       << v.y << "," << std::setprecision(9) << v.z << "}";
    return ss;
};

std::ostream &operator<<(std::ostream &ss, aiAABB const &bb) {
    ss << "{" << bb.mMin << "," << bb.mMax << "}";
    return ss;
};

std::ostream &operator<<(std::ostream &ss, aiMatrix4x4 const &m) {
    ss << "[[" << m.a1 << ", " << m.a2 << ", " << m.a3 << ", " << m.a4 << "], "
       << "[" << m.b1 << ", " << m.b2 << ", " << m.b3 << ", " << m.b4 << "], "
       << "[" << m.c1 << ", " << m.c2 << ", " << m.c3 << ", " << m.c4 << "], "
       << "[" << m.d1 << ", " << m.d2 << ", " << m.d3 << ", " << m.d4 << "]]";
    return ss;
}

std::ostream &operator<<(std::ostream &ss, aiString const &s)
{
    ss << s.C_Str();
    return ss;
}

} // namespace

