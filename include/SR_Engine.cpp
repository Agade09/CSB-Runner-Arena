#pragma once
#include <bits/stdc++.h>
#include "vec.cpp"
#include "assert.cpp"
using namespace std;
constexpr float MaxSpeed{1133.3333333333333333333+200},DegToRad{M_PI/180.0};
constexpr int CP_Radius{600},Pod_Radius{400},N_Laps{3};
constexpr int Timeout_Length{50};

struct pod{
    vec r,v;
    float angle;//Degrees
    int next,lap;
};

struct action{
    int thrust;
    float delta_angle;//Degrees
};

ostream& operator<<(ostream &os,const action &mv){
    os << mv.thrust << " " << mv.delta_angle;
    return os;
}

istream& operator>>(istream &is,action &mv){
    is >> mv.thrust >> mv.delta_angle;
    return is;
}

typedef vector<vec> Map;

struct state{
	pod p;
	int turn{0},timeout;
	Map C;
	inline bool finished_race()const{
		if(p.lap==3){
			return true;
		}
		return false;
	}
	inline bool game_over()const{
		return finished_race() || timeout<=0;
	}
	inline int remaining_CPs(const pod &p)const{
		if(p.lap==3){
			return 0;
		}
		return max(0,2-p.lap)*C.size()+(p.next==0?1:C.size()-p.next+1); 
	}
};

struct collision{
    double t;
};

const vector<vector<vec>> Maps = {
    {{1000, 4500}, {2500, 3905}, {4000, 5095}, {5500, 3905}, {7000, 5095}, {8500, 3905}, {10000, 5095}, {11500, 3905}},
    {{1000, 1000}, {12000, 1000}, {12500, 2500}, {13000, 4000}, {12500, 5500}, {12000, 7000}},
    {{8179, 7909}, {11727, 5704}, {11009, 3026}, {10111, 1169}, {5835, 7503}, {1380, 2538}, {4716, 1269}, {4025, 5146}},
    {{6910, 1656}, {14908, 1849}, {2485, 3249}, {5533, 6258}, {12561, 1063}, {1589, 6883}, {13542, 2666}, {13967, 6917}},
    {{4912, 4817}, {9882, 5377}, {3692, 3080}, {3562, 1207}, {4231, 7534}, {14823, 6471}, {10974, 1853}, {9374, 3740}},
    {{13332, 4114}, {5874, 7746}, {7491, 4801}, {14268, 6672}, {2796, 1751}, {1039, 2272}, {6600, 1874}, {13467, 2208}},
    {{5674, 4795}, {9623, 7597}, {12512, 6231}, {4927, 3377}, {8358, 6630}, {4459, 7216}, {10301, 2326}, {2145, 3943}},
    {{12271, 7160}, {14203, 4266}, {3186, 5112}, {8012, 5958}, {2554, 6642}, {5870, 4648}, {11089, 2403}, {9144, 2389}},
    {{14086, 1366}, {1779, 2501}, {5391, 2200}, {13348, 4290}, {6144, 4176}, {11687, 5637}, {14990, 3490}, {3569, 7566}},
    {{9302, 5289}, {6419, 7692}, {2099, 4297}, {13329, 3186}, {13870, 7169}, {13469, 1115}, {5176, 5061}, {1260, 7235}},
    {{6752, 5734}, {10177, 7892}, {5146, 7584}, {11531, 1216}, {1596, 5797}, {8306, 3554}, {5814, 2529}, {9471, 5505}},
    {{9476, 3253}, {10312, 1696}, {2902, 6897}, {5072, 7852}, {5918, 1004}, {3176, 2282}, {14227, 2261}, {9986, 5567}},
    {{10353, 1986}, {2757, 4659}, {3358, 2838}},
    {{13048, 3493}, {9614, 3949}, {6999, 2367}, {12067, 4880}, {8525, 5705}, {6759, 5582}, {14646, 5876}, {4158, 3179}},
    {{8565, 6690}, {10713, 3220}, {7321, 7928}, {1578, 3893}, {6882, 2145}, {8878, 3844}, {1025, 7671}, {3637, 6578}},
    {{7483, 7350}, {14154, 4505}, {3917, 7630}, {9957, 6899}, {8070, 3272}, {1884, 1763}, {3155, 3640}, {10140, 1152}},
    {{4295, 7416}, {12579, 6780}, {6585, 5187}, {2804, 6546}, {5038, 4810}, {1702, 1007}, {10114, 1658}, {8425, 6507}},
    {{8333, 6039}, {13617, 5740}, {13457, 3465}, {6659, 7011}, {12132, 6914}, {10277, 1624}, {3740, 2896}, {9054, 7429}},
    {{13853, 1419}, {12855, 3432}, {2453, 7829}, {8173, 7778}, {1428, 4878}, {10194, 3223}, {2814, 2394}, {11452, 1809}},
    {{3144, 4250}, {6352, 4948}, {14725, 5968}, {4864, 7961}, {8442, 1307}, {14501, 3206}, {12630, 7105}, {1767, 6800}},
    {{7355, 6865}, {11967, 7228}, {4501, 7146}, {2977, 5349}, {9592, 4217}, {11713, 4176}, {10485, 2979}, {6139, 1981}},
    {{2589, 1765}, {7918, 4590}, {5921, 4279}, {10590, 2077}, {9780, 6425}, {5945, 6701}, {14440, 3369}, {4988, 2966}},
    {{6732, 4875}, {2541, 1997}, {13969, 3703}, {11421, 5223}, {7687, 7371}, {2560, 4311}, {3857, 5771}, {14273, 1692}},
    {{11141, 4590}, {3431, 6328}, {4284, 2801}},
    {{5125, 6049}, {6292, 4792}, {7679, 7898}, {3140, 5406}, {4676, 4325}, {8348, 3287}, {10258, 2927}, {1620, 1867}},
    {{3416, 2572}, {3994, 6091}, {6110, 7235}, {1493, 4089}, {1537, 7029}, {1594, 2079}, {8993, 5700}, {13129, 7028}},
    {{1613, 4944}, {13640, 7073}, {2072, 1872}, {14854, 3078}, {4484, 2083}, {10084, 5389}, {7002, 1561}, {8127, 6064}},
    {{13389, 4971}, {5988, 4410}, {8092, 1250}, {14259, 3041}, {13657, 6973}, {7445, 6601}, {4240, 5278}, {1662, 1335}},
    {{4290, 7419}, {9451, 1669}, {12371, 7190}, {3974, 3907}, {11155, 4435}, {3274, 2000}, {14666, 4335}, {6054, 1285}},
    {{3486, 5885}, {10938, 1106}, {6113, 2576}, {10667, 4140}, {13926, 1263}, {1638, 2764}, {7838, 7775}, {14491, 2824}},
    {{10449, 5356}, {12806, 1444}, {2108, 2201}, {6362, 3893}, {10672, 3222}, {12535, 6285}, {6657, 2160}, {13184, 4759}},
    {{3899, 5119}, {11225, 5109}, {14442, 2223}, {6176, 4313}, {7409, 5908}, {10162, 1057}, {8179, 7625}, {3765, 1798}},
    {{9663, 4767}, {9850, 7663}, {8480, 3611}, {9745, 2178}, {7110, 1185}, {11651, 3625}, {2446, 7330}, {12226, 7367}},
    {{14377, 7927}, {10282, 6710}, {4391, 4202}, {6951, 1538}, {12324, 6293}, {13854, 5681}, {1738, 6745}, {8578, 4429}},
    {{1515, 4019}, {2544, 5921}, {4463, 1963}, {6827, 5569}, {11309, 6057}, {14948, 4812}, {13964, 2231}, {7365, 3517}},
    {{4045, 6842}, {1465, 1643}, {6862, 6261}, {1491, 5253}, {14727, 3179}, {6082, 3183}, {13219, 4669}, {13310, 7720}},
    {{8957, 1732}, {1578, 6406}, {11732, 5976}, {11436, 2425}, {14054, 2637}, {2551, 3170}, {9647, 5917}, {5920, 3971}},
    {{12894, 4304}, {10314, 6249}, {4368, 6750}, {2885, 5776}, {9302, 4546}, {3233, 2294}, {14572, 4689}, {7955, 5693}},
    {{7237, 7984}, {4874, 2369}, {5346, 5225}, {6399, 3637}, {3165, 1448}, {9308, 1323}, {12931, 4556}, {2178, 5950}},
    {{3610, 6881}, {14942, 4738}, {10731, 1634}, {3880, 4129}, {3262, 1431}, {7435, 6793}, {8388, 4702}, {13784, 1301}},
    {{14682, 5560}, {12208, 5614}, {10579, 1951}, {9412, 3808}, {11739, 7467}, {2559, 5223}, {6841, 4150}, {6838, 6961}},
    {{10542, 1448}, {13016, 2024}, {12094, 4489}, {6045, 6111}, {3079, 6468}, {9520, 4308}, {6688, 1488}, {12845, 6967}},
    {{9212, 2566}, {8596, 4908}, {10021, 1067}, {6747, 5088}, {2067, 4245}, {10389, 6367}, {4507, 4352}, {5075, 1492}},
    {{8150, 7321}, {4285, 3513}, {12095, 3060}, {2383, 4893}, {9220, 1508}, {10207, 6519}, {5204, 7579}, {13766, 1956}},
    {{3627, 5131}, {2915, 7241}, {12639, 1171}, {6549, 3529}, {13500, 3687}, {1746, 5642}, {8351, 4522}, {1839, 3045}},
    {{12000, 1000}, {12500, 2500}, {12500, 5500}, {12000, 7000}, {8000, 7000}, {7500, 5500}, {7500, 2500}, {8000, 1000}},
    {{1000, 1000}, {15000, 8000}, {1000, 8000}, {15000, 1000}, {1000, 4500}, {15000, 4500}},
    {{1711, 3942}, {10892, 5399}, {4058, 1092}, {6112, 2872}, {1961, 6027}, {7148, 4594}, {7994, 1062}},
    {{12603, 1090}, {1043, 1446}, {10158, 1241}, {13789, 7502}, {7456, 3627}, {6218, 1993}, {7117, 6546}, {5163, 7350}},
    {{9214, 6145}, {1271, 7171}, {14407, 3329}, {10949, 2136}, {2443, 4165}, {5665, 6432}, {3079, 1942}, {4019, 5141}}};
    
inline double Angular_Distance(const double a,const double b)noexcept{
    const double diff{abs(a-b)};
    return min(diff,360-diff);
}

inline float Angle_Sum(const float a,const float b)noexcept{
    float sum{a+b};
    return sum+(sum>360?-360:sum<0?360:0);
}

inline double Closest_Angle(const double a,const double t){
    if(Angular_Distance(a,t)<=18){
        return t;
    }
    else{
        const double a1{Angle_Sum(a,18.0)},a2{Angle_Sum(a,-18.0)};
        return Angular_Distance(a1,t)<Angular_Distance(a2,t)?a1:a2;
    }
}

inline double Angle(const vec &d){
	if(d.x==0 && d.y==0){
		return 0;//Avoid domain error
	}
	const double rad{atan2(d.y,d.x)};// atan2() returns radians [-Pi,Pi]
    return (rad>0?rad:2*M_PI+rad)/DegToRad;//[0,360]
}

inline void Passes_Checkpoint(const pod &p,const vec &CP,list<collision> &C,const double T,const float Simulation_Time)noexcept{
    vec R=p.r-CP;
    double R2{R.norm2()},RV{R*p.v},V2{p.v.norm2()},det{RV*RV-V2*(R2-pow(CP_Radius,2))};
    if(RV<0.0 && det>=0.0){
        double t{T-(RV+sqrt(det))/V2};
        if(t<Simulation_Time){
            C.push_back(collision{t});
        }
    }
}

inline void ScanForAllCollisions(list<collision> &C,const state &S,const float Simulation_Time)noexcept{
    Passes_Checkpoint(S.p,S.C[S.p.next],C,0,Simulation_Time);
}

inline void ScanForNewCollisions(list<collision> &C,const state &S,const double T,const float Simulation_Time)noexcept{
    Passes_Checkpoint(S.p,S.C[S.p.next],C,T,Simulation_Time);
}

template <bool is_coll> void Simulation_Step(const collision &c,state &S,list<collision> &C,double &T,const float Simulation_Time){
    S.p.r+=S.p.v*(c.t-T);
    T=c.t;
    if(is_coll){
        C.clear();
        pod &p=S.p;
         ++p.next;
        if(p.next==S.C.size()){
            p.next=0;
        }
        else if(p.next==1){
            ++p.lap;
        }
        S.timeout=Timeout_Length;
        ScanForNewCollisions(C,S,T,Simulation_Time);
    }
}

void Simulate_Pod_Move(pod &p,const action &mv){
    p.angle=Angle_Sum(p.angle,mv.delta_angle);
    vec dir;
    sincosf(p.angle*DegToRad,&dir.y,&dir.x);
    p.v+=dir*mv.thrust;
}

void advance_time(state &S,const float Simulation_Time){
    double T{0};
    list<collision> C;
    ScanForAllCollisions(C,S,Simulation_Time);
    //int CP_Passage{0};
    while(!C.empty()){
        collision first=*min_element(C.begin(),C.end(),[](const collision &a,const collision &b){return a.t<b.t;});
        Simulation_Step<true>(first,S,C,T,Simulation_Time);
        //++CP_Passage;
    }
    /*if(CP_Passage>=2){
        for(const vec &CP:S.C){
            cerr << CP << " ";
        }
        cerr << endl;
    }*/
    Simulation_Step<false>(collision{Simulation_Time},S,C,T,Simulation_Time);
    pod &p{S.p};
    p.r=vec{round(p.r.x),round(p.r.y)};
    p.v=vec{trunc(0.85f*p.v.x),trunc(0.85f*p.v.y)};
}

void Simulate(state &S,const action &mv){
    Simulate_Pod_Move(S.p,mv);
    ++S.turn;
    --S.timeout;
    advance_time(S,1);
}