#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>
#include <sstream>
#include <unistd.h>
#include <sys/wait.h>
#include <ext/stdio_filebuf.h>
#include <poll.h>
#include <omp.h>
using namespace std;
using namespace std::chrono;

struct vec{
    int x,y;
    inline bool operator!=(const vec &a)const noexcept{
        return x!=a.x || y!=a.y;
    }
    inline bool operator==(const int a)const noexcept{
        return x==a && y==a;
    }
    inline double norm2()const noexcept{
        return static_cast<double>(x)*static_cast<double>(x)+static_cast<double>(y)*static_cast<double>(y);
    }
    inline double norm()const noexcept{
        return sqrt(norm2());
    }
    inline vec operator-(const vec &a)const noexcept{
    	return vec{x-a.x,y-a.y};
    }
};

struct vecf{
    double x,y;
    inline double norm2()const noexcept{
        return x*x+y*y;
    }
    inline double norm()const noexcept{
        return sqrt(norm2());
    }
    inline void operator/=(const double a)noexcept{
        x/=a;
        y/=a;
    }
    inline void normalize()noexcept{
        *this/=norm();
    }
};

struct play{
	vec t;
	string thrust;
};

typedef vector<vec> Map;

constexpr int PIPE_READ{0},PIPE_WRITE{1};
constexpr int N_L{3},N{1};
constexpr double DegToRad{M_PI/180.0};
constexpr bool Debug_AI{false},Timeout{false};
const vector<Map> Maps{{{12443,1345},{10527,5973},{3601,5164},{13557,7577}},{{7495,6936},{6021,5375},{11281,2838}},{{5436,2813},{10342,3362},{11183,5397},{7263,6669}},{{6103,2192},{3011,5214},{6280,7779},{14130,7732},{13855,1202},{10238,4941}},{{13910,1964},{8047,3252},{2694,6993},{10063,5943}},{{10543,5039},{13114,2337},{4569,2164},{7375,4967},{3324,7235},{14588,7687}},{{11185,5450},{7275,6640},{5399,2831},{10324,3337}},{{11217,5429},{7267,6636},{5445,2818},{10327,3364}},{{10517,6007},{3605,5165},{13564,7580},{12430,1357}},{{5975,5356},{11330,2792},{7470,6967}},{{13056,1891},{6545,7863},{7504,1379},{12687,7085},{4053,4677}},{{7982,7911},{13271,5511},{9588,1411},{3631,4408}},{{6537,7841},{7480,1388},{12725,7122},{4040,4665},{13026,1899}},{{4123,7390},{13530,2333},{12965,7230},{5624,2560}},{{8670,7466},{7189,2140},{3581,5262},{13853,5060},{10663,2253}},{{12919,7244},{5663,2575},{4073,7436},{13489,2365}},{{8713,7474},{7218,2184},{3625,5283},{13861,5105},{10664,2252}},{{7796,853},{7634,5978},{3134,7547},{9549,4370},{14522,7750},{6344,4278}},{{7481,1337},{12719,7126},{4083,4657},{13053,1894},{6557,7820}},{{6019,5341},{11321,2843},{7526,6933}},{{13480,2354},{12964,7233},{5665,2581},{4126,7427}},{{5439,2855},{10352,3388},{11186,5418},{7261,6668}},{{13892,1935},{7991,3278},{2681,7041},{10043,5985}},{{11484,6057},{9075,1848},{5021,5245}},{{5964,4225},{14648,1424},{3435,7217},{9450,7218}},{{13118,2339},{4536,2174},{7325,4933},{3338,7201},{14556,7712},{10533,5081}},{{3581,5180},{13562,7615},{12489,1332},{10512,5959}},{{3569,5200},{13573,7613},{12461,1343},{10539,5966}},{{13305,5551},{9589,1392},{3648,4422},{8003,7872}},{{7323,4951},{3310,7225},{14610,7682},{10563,5060},{13085,2305},{4559,2164}},{{7511,6954},{6015,5340},{11326,2835}},{{6294,4263},{7779,833},{7658,5965},{3131,7511},{9500,4405},{14548,7762}},{{11293,2806},{7475,6916},{6011,5387}},{{3649,4444},{7977,7878},{13275,5556},{9576,1422}},{{8693,7480},{7193,2138},{3600,5308},{13824,5102},{10692,2303}},{{9540,4393},{14498,7791},{6329,4309},{7787,871},{7689,5949},{3144,7534}},{{13873,1225},{10251,4932},{6124,2188},{3041,5165},{6297,7734},{14091,7731}},{{7689,5950},{3135,7511},{9511,4388},{14539,7757},{6316,4302},{7812,854}},{{11480,6053},{9125,1824},{5005,5236}},{{12960,7232},{5636,2583},{4099,7403},{13492,2322}},{{13508,2341},{12939,7242},{5629,2567},{4109,7407}},{{8028,3249},{2650,7031},{10030,5946},{13948,1934}},{{8026,3253},{2664,7027},{10063,5986},{13944,1934}},{{8677,7452},{7171,2173},{3594,5288},{13856,5068},{10680,2287}},{{13021,1928},{6550,7818},{7460,1352},{12670,7110},{4083,4637}},{{6304,4304},{7804,838},{7662,5948},{3144,7567},{9515,4390},{14540,7773}},{{5423,2839},{10296,3382},{11214,5411},{7253,6660}},{{13838,5092},{10701,2262},{8676,7483},{7209,2163},{3625,5265}},{{11312,2842},{7477,6926},{5977,5370}},{{13856,5076},{10707,2297},{8694,7446},{7222,2186},{3614,5268}},{{13473,2340},{12910,7200},{5666,2573},{4106,7438}},{{12466,1348},{10522,5985},{3566,5206},{13581,7593}},{{13890,1947},{8034,3271},{2664,7044},{10059,5993}},{{13488,2356},{12940,7225},{5644,2598},{4074,7418}},{{11307,2831},{7509,6954},{5974,5346}},{{3326,7200},{14562,7690},{10560,5086},{13110,2346},{4574,2200},{7324,4943}},{{10668,2262},{8726,7470},{7228,2175},{3601,5257},{13812,5077}},{{12487,1373},{10514,5986},{3570,5167},{13592,7624}},{{4094,7391},{13479,2342},{12922,7200},{5636,2587}},{{3571,5165},{13555,7573},{12456,1372},{10546,5963}},{{9511,4369},{14518,7759},{6332,4313},{7816,846},{7669,5969},{3124,7540}},{{7484,6925},{5982,5383},{11273,2820}},{{5619,2557},{4094,7421},{13526,2339},{12967,7203}},{{6589,7846},{7506,1337},{12670,7104},{4047,4652},{13048,1900}},{{4993,5285},{11500,6098},{9080,1867}},{{9571,1423},{3622,4436},{8011,7917},{13326,5548}},{{10305,3369},{11230,5402},{7247,6675},{5414,2825}},{{10044,5962},{13913,1914},{8009,3254},{2662,7034}},{{10337,3353},{11205,5400},{7243,6668},{5430,2826}},{{5630,2554},{4070,7406},{13522,2332},{12950,7244}},{{5412,2838},{10304,3385},{11190,5436},{7288,6681}},{{13510,2334},{12933,7205},{5654,2585},{4097,7398}},{{12920,7229},{5637,2567},{4110,7440},{13522,2311}},{{5453,2836},{10315,3392},{11190,5410},{7256,6667}},{{3569,5160},{13565,7573},{12465,1338},{10538,5953}},{{3612,5292},{13838,5093},{10679,2303},{8711,7474},{7192,2145}},{{7230,6665},{5396,2841},{10298,3380},{11191,5425}},{{10583,5067},{13107,2302},{4580,2160},{7377,4955},{3293,7212},{14597,7700}},{{6024,5388},{11298,2846},{7486,6936}},{{14516,7770},{6307,4298},{7777,872},{7663,5998},{3148,7553},{9546,4364}},{{11477,6051},{9120,1826},{5038,5237}},{{10518,6004},{3573,5169},{13599,7604},{12467,1366}},{{3610,5193},{13582,7602},{12485,1359},{10550,5984}},{{13867,5100},{10654,2253},{8708,7475},{7207,2187},{3580,5268}},{{13329,5565},{9556,1399},{3642,4431},{7993,7887}},{{14503,7799},{6304,4306},{7816,871},{7683,5976},{3121,7549},{9499,4376}},{{3349,7214},{14565,7720},{10583,5080},{13095,2332},{4553,2203},{7359,4915}},{{11476,6069},{9103,1857},{4981,5255}},{{12457,1364},{10548,5957},{3588,5201},{13597,7585}},{{5972,4258},{14684,1391},{3479,7207},{9440,7260}},{{7332,4937},{3329,7257},{14603,7715},{10539,5032},{13117,2293},{4561,2197}},{{2647,7009},{10015,5953},{13898,1916},{8033,3251}},{{6092,2194},{3034,5207},{6264,7749},{14118,7734},{13908,1242},{10222,4916}},{{10318,3394},{11224,5436},{7276,6650},{5440,2859}},{{3346,7201},{14591,7705},{10558,5050},{13097,2307},{4565,2165},{7354,4933}},{{8018,7896},{13272,5534},{9540,1395},{3619,4395}},{{11232,5445},{7241,6681},{5443,2866},{10293,3392}},{{9446,7252},{5995,4223},{14673,1385},{3450,7206}},{{7262,6627},{5405,2848},{10315,3377},{11173,5431}},{{9100,1851},{5014,5235},{11494,6060}}};

struct pod{
    vec r,v;
    int angle,next,lap;
    inline bool operator!=(const pod &a)const noexcept{
        return r!=a.r || v!=a.v || angle!=a.angle || next!=a.next || lap!=a.lap;
    }
};

struct AI{
	int id,pid,outPipe,errPipe;
	string name;
	ostream *in;
	__gnu_cxx::stdio_filebuf<char> inBuff;
	pod p;
	bool dead;
	AI()=default;
	inline AI(const AI &a):id(a.id),name(a.name),p(a.p),dead(a.dead){}
	inline ~AI(){
		kill(pid,SIGKILL);
		close(errPipe);
		close(outPipe);
		delete in;
	}
};

void StartProcess(AI &Bot){
	int StdinPipe[2];
  	int StdoutPipe[2];
  	int StderrPipe[2];
  	if(pipe(StdinPipe)<0){
	    perror("allocating pipe for child input redirect");
  	}
  	if(pipe(StdoutPipe)<0){
	    close(StdinPipe[PIPE_READ]);
	    close(StdinPipe[PIPE_WRITE]);
	    perror("allocating pipe for child output redirect");
  	}
  	if(pipe(StderrPipe)<0){
  		close(StderrPipe[PIPE_READ]);
  		close(StderrPipe[PIPE_WRITE]);
  		perror("allocating pipe for child stderr redirect failed");
  	}
  	int nchild{fork()};
  	//cerr << "nchild: " << nchild << " " << Bot.name << endl;
  	if(nchild==0){//Child process
	    if(dup2(StdinPipe[PIPE_READ],STDIN_FILENO)==-1){// redirect stdin
			perror("redirecting stdin");
			return;
	    }
	    if(dup2(StdoutPipe[PIPE_WRITE],STDOUT_FILENO)==-1){// redirect stdout
			perror("redirecting stdout");
			return;
	    }
	    if(dup2(StderrPipe[PIPE_WRITE],STDERR_FILENO)==-1){// redirect stderr
			perror("redirecting stderr");
			return;
	    }
	    close(StdinPipe[PIPE_READ]);
	    close(StdinPipe[PIPE_WRITE]);
	    close(StdoutPipe[PIPE_READ]);
	    close(StdoutPipe[PIPE_WRITE]);
	    close(StderrPipe[PIPE_READ]);
	    close(StderrPipe[PIPE_WRITE]);
	    execl(Bot.name.c_str(),Bot.name.c_str(),(char*)0);//(char*)0 is really important
	    //If you get past the previous line its an error
	    perror("exec of the child process");
  	}
  	else if(nchild>0){//Parent process
  		close(StdinPipe[PIPE_READ]);//Parent does not read from stdin of child
    	close(StdoutPipe[PIPE_WRITE]);//Parent does not write to stdout of child
    	close(StderrPipe[PIPE_WRITE]);//Parent does not write to stderr of child
    	Bot.inBuff=__gnu_cxx::stdio_filebuf<char>(StdinPipe[PIPE_WRITE], std::ios::out);
    	Bot.in=new ostream(&Bot.inBuff);
		//Bot.outBuff=__gnu_cxx::stdio_filebuf<char>(StdoutPipe[PIPE_READ], std::ios::in);
		//Bot.out = new istream(&Bot.outBuff);
		Bot.outPipe=StdoutPipe[PIPE_READ];
    	Bot.errPipe=StderrPipe[PIPE_READ];
    	Bot.pid=nchild;
  	}
  	else{//failed to create child
  		close(StdinPipe[PIPE_READ]);
	    close(StdinPipe[PIPE_WRITE]);
	    close(StdoutPipe[PIPE_READ]);
	    close(StdoutPipe[PIPE_WRITE]);
	    perror("Failed to create child process");
  	}
}

int Duplicates(const vector<Map> &Maps){
	int duplicates{0};
	for(int k=0;k<Maps.size();++k){
		for(int j=k+1;j<Maps.size();++j){
			const Map &m{Maps[k]},&m2{Maps[j]};
			bool same{m.size()==m2.size()};
			if(same){
				for(int i=0;i<m.size();++i){
					if(m[i]!=m2[i]){
						same=false;
						break;
					}
				}
			}
			if(same){
				++duplicates;
			}
		}
	}
	return duplicates;
}

inline int Angular_Distance(const int a,const int b){
    int diff{abs(a-b)};
    return min(diff,360-diff);
}

inline double Angle_Sum(const double a,const double b){
    double sum{a+b};
    return sum+(sum>360?-360:sum<0?360:0);
}

inline double Angle(const vec &d){
    double a{round(acos(static_cast<double>(d.x)/d.norm())*180/M_PI)};//Should I keep this round()?
    if(d.y<0){
        a=360-a;
    }
    return a;
}

inline int Closest_Angle(const int a,const int t){
    if(Angular_Distance(a,t)<=18){
        return t;
    }
    else{
        double a1{Angle_Sum(a,18.0)},a2{Angle_Sum(a,-18.0)};
        return Angular_Distance(a1,t)<Angular_Distance(a2,t)?a1:a2;
    }
}

inline bool Passes_Checkpoint(const vec &r,const vecf &v,const vec &t){
    double Rx{static_cast<double>(r.x-t.x)},Ry{static_cast<double>(r.y-t.y)},R2{Rx*Rx+Ry*Ry};
    double RV{Rx*v.x+Ry*v.y},V2{v.x*v.x+v.y*v.y},det{RV*RV-V2*(R2-600*600)};
    return RV<0.0 && det>=0.0 && -(RV+sqrt(det))/V2<1;
}

bool Make_Move(const string &MoveStr,AI &Bot,const Map &C,const int turn){
	stringstream ss(MoveStr);
	play Move;
	ss >> Move.t.x >> Move.t.y >> Move.thrust;
	pod &p=Bot.p;
	if(Move.thrust=="SHIELD"){
		cerr << "Found SHIELD" << endl;
	}
	int thrust=Move.thrust=="BOOST"?650:Move.thrust=="SHIELD"?0:stoi(Move.thrust);
	double desired_angle{Angle(Move.t-p.r)};
	p.angle=turn==1?desired_angle:Closest_Angle(p.angle,desired_angle);
    vecf new_v{p.v.x+cos(p.angle*DegToRad)*thrust,p.v.y+sin(p.angle*DegToRad)*thrust};
    if(Passes_Checkpoint(p.r,new_v,C[p.next])){
        ++p.next;
        if(p.next==C.size()){
            p.next=0;
        }
        else if(p.next==1){
            ++p.lap;
        }
    }
    p.r.x+=round(new_v.x);
    p.r.y+=round(new_v.y);
    p.v.x=0.85*new_v.x;
    p.v.y=0.85*new_v.y;
    return true;
}

inline string EmptyPipe(const int f){
	char out[5000];
	fill(out,out+5000,0);
	int ptr{0};
	pollfd outpoll{f,POLLIN};
	while(poll(&outpoll,1,0)){
		ssize_t bytes_read{read(f,out+(ptr++),1)};
		if(bytes_read<1){
			cerr <<"Pipe is empty" << endl;
			exit(0);
		}
	}
	return string(out);
}

bool IsValidMove(const string &Move){
	stringstream ss(Move);
	for(int i=0;i<2;++i){
		play PodMove;
		if(!(ss >> PodMove.t.x >> PodMove.t.y >> PodMove.thrust)){
			return false;
		}
	}
	return true;
}

string GetMove(AI &Bot,const int turn){
	pollfd outpoll{Bot.outPipe,POLLIN|POLLPRI};
	time_point<system_clock> Start_Time{system_clock::now()};
	string Move;
	while( (!Timeout || static_cast<duration<double>>(system_clock::now()-Start_Time).count()<(turn==1?1:0.15)) && !IsValidMove(Move)){
		double TimeLeft{(turn==1?1:0.15)-static_cast<duration<double>>(system_clock::now()-Start_Time).count()};
		if(poll(&outpoll,1,Timeout?TimeLeft*1000:-1)){
			Move+=EmptyPipe(Bot.outPipe);
		}
	}
	return IsValidMove(Move)?Move:"";
}

int Play_Game(vector <string> &Bot_Names,const Map &C){
	vector<AI> Bot(Bot_Names.size());
	for(int i=0;i<Bot_Names.size();++i){
		Bot[i].id=i;
		Bot[i].name=Bot_Names[i];
		StartProcess(Bot[i]);
		Bot[i].dead=false;
		Bot[i].p.r=C[0];
		Bot[i].p.v=vec{0,0};
		Bot[i].p.lap=0;
		Bot[i].p.next=1;
		Bot[i].p.angle=0;//Could do better?
	}
	int turn{0};
	//Feed first turn inputs
	for(const AI &b:Bot){
		*b.in << N_L << " " << C.size() << endl;
		for(const vec &cp:C){
			*b.in << cp.x << " " << cp.y << endl;
		}
	}
	while(++turn>0){
		for(int i=0;i<Bot.size();++i){
			if(kill(Bot[i].pid,0)!=ESRCH){//Check if process is still running
				//Feed turn inputs
				for(int j=0;j<2;++j){//Duplicate pod info
					*Bot[i].in << Bot[i].p.r.x << " " << Bot[i].p.r.y << " " << Bot[i].p.v.x << " " << Bot[i].p.v.y << " " << Bot[i].p.angle << " " << Bot[i].p.next << endl;
				}
				for(int j=0;j<2;++j){//Feed fake far away opponent position
					*Bot[i].in << -10000 << " " << -10000 << " " << 0 << " " << 0 << " " << 0 << " " << 1 << endl;
				}
				string Move{GetMove(Bot[i],turn)};
				if(Move!=""){
					if(!Make_Move(Move,Bot[i],C,turn)){
						return -1;
					}
					if(Bot[i].p.lap==N_L){
						return turn;
					}
				}
				else{
					cerr << "Loss by timeout of AI " << i << endl;
					return -1;
				}
				string err_str{EmptyPipe(Bot[i].errPipe)};
				if(Debug_AI){
					ofstream err_out("log.txt",ios::app);
					err_out << err_str << endl;
				}
			}
		}
	}
	return -2;
}

int main(int argc,char **argv){
	time_point<system_clock> Start_Time{system_clock::now()};
	cerr << Maps.size() << " different maps" << endl;
	int dupes{Duplicates(Maps)};
	if(dupes>0){
		cerr << "Warning: found " << dupes << " duplicate maps" << endl;
	}
	if(argc<2){
		cerr << "Program takes 1 input, the names of the AI to test" << endl;
		return 0;
	}
	vector<string> Bot_Names(N);
	for(int i=0;i<N;++i){
		Bot_Names[i]=argv[i+1];
	}
	cout << "Testing AI " << Bot_Names[0] << endl;
	for(int i=0;i<Bot_Names.size();++i){
		ifstream Test{Bot_Names[i].c_str()};
		if(!Test){
			cerr << Bot_Names[i] << " couldn't be found" << endl;
			return 0;
		}
		Test.close();
	}
	int total_turns{0},maps_done{0};
	#pragma omp parallel for num_threads(1) shared(total_turns,maps_done,Bot_Names) //schedule(dynamic)
	for(int i=0;i<Maps.size();++i){
		int turns{Play_Game(Bot_Names,Maps[i])};
		if(turns>0){
			#pragma omp atomic
			total_turns+=turns;
		}
		else{
			cerr << "Bot didn't complete map" << endl;
		}
		#pragma omp atomic
		++maps_done;
		cerr << 100.0*maps_done/Maps.size() << "% done \r";
	}
	cerr << "Total turns: " << total_turns << endl;
	cerr << static_cast<duration<double>>(system_clock::now()-Start_Time).count() << " s" << endl;
}