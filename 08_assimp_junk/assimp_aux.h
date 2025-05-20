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
} // namespace

