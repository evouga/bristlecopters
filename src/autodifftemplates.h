#ifndef AUTODIFFTEMPLATES_H
#define AUTODIFFTEMPLATES_H


template<typename T> void diff(const T *vec1, const T *vec2, T *difference)
{
    for(int i=0; i<3; i++)
        difference[i] = vec1[i]-vec2[i];
}

template<typename T> T norm(const T *vec)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += vec[i]*vec[i];
    return sqrt(result);
}

template<typename T> void cross(const T *vec1, const T *vec2, T *result)
{
    result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

template<typename T> void normalize(T *vec)
{
    T n = norm(vec);
    for(int i=0; i<3; i++)
        vec[i] /= n;
}

template<typename T> T dot(const T *vec1, const T *vec2)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += vec1[i]*vec2[i];
    return result;
}

template<typename T> T triarea(T len1, T len2, T len3)
{
    T semi = (len1+len2+len3)/2.0;
    return sqrt(semi*(semi-len1)*(semi-len2)*(semi-len3));
}

template<typename T> T sintri(T side1, T side2, T oppside)
{
    // s1 s2 sin theta = 2A
    T area = triarea(side1, side2, oppside);
    return (area*2.0)/side1/side2;
}

template<typename T> T costri(T side1, T side2, T oppside)
{
    // law of cosines
    return (side1*side1+side2*side2-oppside*oppside)/(2.0*side1*side2);
}

template<typename T> T cottri(T side1, T side2, T oppside)
{
    return costri(side1, side2, oppside)/sintri(side1, side2, oppside);
}

template<typename T> T Ledge(T valcenter, T valnb, T midlen, T side1, T side2, T opp1, T opp2)
{
    T cot1 = cottri(side1, opp1, midlen);
    T cot2 = cottri(side2, opp2, midlen);
    T cotweight = (cot1+cot2)/2.0;
    T diffval = valcenter-valnb;
    return cotweight*diffval;
}

template<typename T> T L(T valcenter, int numnbs, T *valnbs, T *midlens, T *rightopplens)
{
    T result = 0;
    for(int i=0; i<numnbs; i++)
    {
        T valnb = valnbs[i];
        T midlen = midlens[i];
        int previdx = (numnbs+i-1) % numnbs;
        T side1 = midlens[previdx];
        T opp1 = rightopplens[previdx];
        int nextidx = (i+1) % numnbs;
        T side2 = midlens[nextidx];
        T opp2 = rightopplens[i];

        result += Ledge(valcenter, valnb, midlen, side1, side2, opp1, opp2);
    }
}

template<typename T> T dualareatri(T side1, T side2, T oppside)
{
    T cot1 = cottri(side1, oppside, side2);
    T cot2 = cottri(side2, oppside, side1);
    return (side2*side2*cot1+side1*side1*cot2)/8.0;
}

template<typename T> T dualarea(int numnbs, T *lens, T *rightopplens)
{
    T result = 0;
    for(int i=0; i<numnbs; i++)
    {
        int nextidx = (i+1)%numnbs;
        T side1 = lens[i];
        T side2 = lens[nextidx];
        T oppside = rightopplens[i];
        result += dualareatri(side1, side2, oppside);
    }
}

#endif // AUTODIFFTEMPLATES_H
