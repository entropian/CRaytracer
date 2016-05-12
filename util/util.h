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

bool getNextTokenInFile(char buffer[], FILE* fp)
{
    char c = fgetc(fp);
    while((c ==  ' ' || c == '\n') && c != EOF)
    {
        c = fgetc(fp);
    }
    if(c != EOF)
    {
        int i;
        for(i = 0 ; c != ' ' && c != '\n' && c != EOF; i++)
        {
            buffer[i] = c;
            c = fgetc(fp);
        }
        buffer[i] = '\0';
        return true;
    }else
    {
        return false;
    }
}
