#pragma once

inline float min(float a, float b)
{
    return (a < b) ? a : b;
}

inline float max(float a, float b)
{
    return (a > b) ? a : b;
}

inline float abs(const float a)
{
    return (a < 0.0f) ? -a : a;
}
