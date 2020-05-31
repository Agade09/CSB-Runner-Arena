#pragma once
#include <bits/stdc++.h>
using namespace std;

struct vec{
    float x,y;
    inline bool operator!=(const vec &a)const noexcept{
        return x!=a.x || y!=a.y;
    }
    inline float norm2()const noexcept{
        return x*x+y*y;
    }
    inline float norm()const noexcept{
        return sqrt(norm2());
    }
    inline vec normalised()const noexcept{
    	const float n{norm()};
    	return vec{x/n,y/n};
    }
    inline vec normalised_safe()const noexcept{
    	const float n{norm()};
    	if(n==0){
    		return vec{1,0};
    	}
    	return vec{x/n,y/n};
    }
    inline void operator/=(const float a)noexcept{
        x/=a;
        y/=a;
    }
    inline void normalise()noexcept{
        *this/=norm();
    }
    inline float dot(const vec &a)const noexcept{
    	return x*a.x+y*a.y;
    }
    inline float cross(const vec &a)const noexcept{
    	return -x*a.y+y*a.x;
    }
    inline vec operator-(const vec &a)const noexcept{
    	return vec{x-a.x,y-a.y};
    }
    inline vec operator-()const noexcept{
    	return vec{-x,-y};
    }
    inline vec operator+(const vec &a)const noexcept{
    	return vec{x+a.x,y+a.y};
    }
    template<typename T>inline vec operator*(const T a)const noexcept{
    	return vec{x*a,y*a};
    }
    inline float operator*(const vec &a)const noexcept{
        return x*a.x+y*a.y;
    }
    inline void operator+=(const vec &a)noexcept{
    	x+=a.x;
    	y+=a.y;
    }
    inline void operator-=(const vec &a)noexcept{
        x-=a.x;
        y-=a.y;
    }
    inline void operator*=(const float a)noexcept{
        x*=a;
        y*=a;
    }
    inline vec operator/(const float a)const noexcept{
        return vec{x/a,y/a};
    }
    inline void round_vec(){
        x=round(x);
        y=round(y);
    }
};

inline ostream& operator<<(ostream &os,const vec &r){
	os << r.x << " " << r.y;
	return os;
}

inline istream& operator>>(istream &is,vec &r){
    is >> r.x >> r.y;
    return is;
}