#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <cmath>
#include <array>
#include <fstream>
#include <random>
#include <chrono>
#include <sstream>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/ioctl.h>
#include <poll.h>
#include <omp.h>
#include <algorithm>
#include <thread>
using namespace std;
using namespace std::chrono;

struct vec{
    double x,y;
    inline bool operator!=(const vec &a)const noexcept{
        return x!=a.x || y!=a.y;
    }
    inline void operator+=(const vec &a)noexcept{
    	x+=a.x;
    	y+=a.y;
    }
    inline double norm2()const noexcept{
        return pow(x,2)+pow(y,2);
    }
    inline double norm()const noexcept{
        return sqrt(norm2());
    }
    inline vec operator+(const vec &a)const noexcept{
    	return vec{x+a.x,y+a.y};
    }
    inline vec operator-(const vec &a)const noexcept{
    	return vec{x-a.x,y-a.y};
    }
    inline vec operator*(const double &a)const noexcept{
    	return vec{x*a,y*a};
    }
};

inline ostream& operator<<(ostream &os,const vec &r)noexcept{
    os << r.x << " " << r.y;
    return os;
}

inline istream& operator>>(istream &is,vec &r)noexcept{
    is >> r.x >> r.y;
    return is;
}

struct action{
	vec t;
	string thrust;
};

typedef vector<vec> map;

constexpr int PIPE_READ{0},PIPE_WRITE{1};
constexpr int N_L{3},N{1};
constexpr double DegToRad{M_PI/180.0};
constexpr bool Debug_AI{true},Timeout{false};
constexpr double FirstTurnTime{1*(Timeout?1:10)},TimeLimit{0.15*(Timeout?1:10)};
const vector<map> Maps{{{12460,1350},{10540,5980},{3580,5180},{13580,7600}},{{3600,5280},{13840,5080},{10680,2280},{8700,7460},{7200,2160}},{{4560,2180},{7350,4940},{3320,7230},{14580,7700},{10560,5060},{13100,2320}},{{5010,5260},{11480,6080},{9100,1840}},{{14660,1410},{3450,7220},{9420,7240},{5970,4240}},{{3640,4420},{8000,7900},{13300,5540},{9560,1400}},{{4100,7420},{13500,2340},{12940,7220},{5640,2580}},{{14520,7780},{6320,4290},{7800,860},{7660,5970},{3140,7540},{9520,4380}},{{10040,5970},{13920,1940},{8020,3260},{2670,7020}},{{7500,6940},{6000,5360},{11300,2820}},{{4060,4660},{13040,1900},{6560,7840},{7480,1360},{12700,7100}},{{3020,5190},{6280,7760},{14100,7760},{13880,1220},{10240,4920},{6100,2200}},{{10323,3366},{11203,5425},{7259,6656},{5425,2838}}};

bool stop{false};//Global flag to stop all arena threads when SIGTERM is received

struct pod{
    vec r,v;
    int next,lap,time;
    double angle;
};

struct AI{
	int id,pid,outPipe,errPipe,inPipe;
	string name;
	pod p;
	inline void stop(){
        if(alive()){
            kill(pid,SIGTERM);
            int status;
            waitpid(pid,&status,0);//It is necessary to read the exit code for the process to stop
            if(!WIFEXITED(status)){//If not exited normally try to "kill -9" the process
                kill(pid,SIGKILL);
            }
        }
    }
	inline bool alive()const{
		return kill(pid,0)!=-1;//Check if process is still running
	}
	inline void Feed_Inputs(const string &inputs){
        if(write(inPipe,&inputs[0],inputs.size())!=inputs.size()){
            throw(5);
        }
    }
	inline ~AI(){
		close(errPipe);
		close(outPipe);
		close(inPipe);
		stop();
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
        execl(Bot.name.c_str(),Bot.name.c_str(),(char*)NULL);//(char*)Null is really important
        //If you get past the previous line its an error
        perror("exec of the child process");
    }
    else if(nchild>0){//Parent process
        close(StdinPipe[PIPE_READ]);//Parent does not read from stdin of child
        close(StdoutPipe[PIPE_WRITE]);//Parent does not write to stdout of child
        close(StderrPipe[PIPE_WRITE]);//Parent does not write to stderr of child
        Bot.inPipe=StdinPipe[PIPE_WRITE];
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

inline double Angular_Distance(const double a,const double b)noexcept{
    double diff{abs(a-b)};
    return min(diff,360-diff);
}

inline double Angle_Sum(const double a,const double b)noexcept{
    double sum{a+b};
    return sum+(sum>360?-360:sum<0?360:0);
}

inline double Angle(const vec &d)noexcept{
	double n{d.norm()};
	double a{acos(d.x/n)*180/M_PI};
    if(d.y<0){
        a=360-a;
    }
    return n>0?a:0;
}

inline double Closest_Angle(const double a,const double t){
    if(Angular_Distance(a,t)<=18){
        return t;
    }
    else{
        double a1{Angle_Sum(a,18.0)},a2{Angle_Sum(a,-18.0)};
        return Angular_Distance(a1,t)<Angular_Distance(a2,t)?a1:a2;
    }
}

inline bool Passes_Checkpoint(const vec &r,const vec &v,const vec &t){
    double Rx{static_cast<double>(r.x-t.x)},Ry{static_cast<double>(r.y-t.y)},R2{Rx*Rx+Ry*Ry};
    double RV{Rx*v.x+Ry*v.y},V2{v.x*v.x+v.y*v.y},det{RV*RV-V2*(R2-600*600)};
    return RV<0.0 && det>=0.0 && -(RV+sqrt(det))/V2<1;
}

void Make_Move(const string &MoveStr,AI &Bot,const map &C,const int turn){
	stringstream ss(MoveStr);
	action Move;
	ss >> Move.t >> Move.thrust;
	pod &p=Bot.p;
	if(Move.thrust=="SHIELD"){
		cerr << "Found SHIELD" << endl;
	}
	const int thrust{Move.thrust=="BOOST"?650:Move.thrust=="SHIELD"?0:stoi(Move.thrust)};
	const double desired_angle{Angle(Move.t-p.r)};
	p.angle=turn==1?desired_angle:Closest_Angle(p.angle,desired_angle);
    vec new_v{p.v+vec{cos(p.angle*DegToRad),sin(p.angle*DegToRad)}*thrust};
    if(Passes_Checkpoint(p.r,new_v,C[p.next])){
    	p.time=100;
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
    p.v=new_v*0.85;
    --p.time;
}

inline string EmptyPipe(const int fd){
    int nbytes;
    if(ioctl(fd,FIONREAD,&nbytes)<0){
        throw(4);
    }
    string out;
    out.resize(nbytes);
    if(read(fd,&out[0],nbytes)<0){
        throw(4);
    }
    return out;
}

bool IsValidMove(const string &Move){
	stringstream ss(Move);
	for(int i=0;i<2;++i){
		action PodMove;
		if(!(ss >> PodMove.t >> PodMove.thrust)){
			return false;
		}
	}
	return true;
}

string GetMove(AI &Bot,const int turn){
    pollfd outpoll{Bot.outPipe,POLLIN};
    time_point<system_clock> Start_Time{system_clock::now()};
    string out;
    while(static_cast<duration<double>>(system_clock::now()-Start_Time).count()<(turn==1?FirstTurnTime:TimeLimit)){
        double TimeLeft{(turn==1?FirstTurnTime:TimeLimit)-static_cast<duration<double>>(system_clock::now()-Start_Time).count()};
        if(poll(&outpoll,1,TimeLeft)){
            out+=EmptyPipe(Bot.outPipe);
            if(IsValidMove(out)){
            	return out;
            }
        }
    }
    throw(1);
}

inline map Generate_Map()noexcept{
	default_random_engine generator(system_clock::now().time_since_epoch().count());
	uniform_int_distribution<int> Map_Distrib(0,Maps.size()-1);
	map C{Maps[Map_Distrib(generator)]};
	uniform_int_distribution<int> Rotate_Distrib(0,C.size()-1),Delta_Distrib(-30,30);
	rotate(C.begin(),C.begin()+Rotate_Distrib(generator),C.end());
	for(vec &cp:C){
		cp+=vec{static_cast<double>(Delta_Distrib(generator)),static_cast<double>(Delta_Distrib(generator))};
	}
	return C;
}

int Play_Game(array <string,N> &Bot_Names){
	map C=Generate_Map();
	vector<AI> Bot(Bot_Names.size());
	for(int i=0;i<Bot_Names.size();++i){
		Bot[i].id=i;
		Bot[i].name=Bot_Names[i];
		StartProcess(Bot[i]);
		Bot[i].p.r=C[0];
		Bot[i].p.v=vec{0,0};
		Bot[i].p.lap=0;
		Bot[i].p.next=1;
		Bot[i].p.angle=Angle(C[1]-C[0]);//Could do better?
		Bot[i].p.time=100;
	}
	int turn{0};
	//Feed first turn inputs
	for(AI &b:Bot){
		stringstream ss;
		ss << N_L << endl << C.size() << endl;
		for(const vec &cp:C){
			ss << cp << endl;
		}
		b.Feed_Inputs(ss.str());
	}
	while(++turn>0 && !stop){
		for(int i=0;i<Bot.size();++i){
			if(Bot[i].alive()){
				//Feed turn inputs
				stringstream ss;
				for(int j=0;j<2;++j){//Duplicate pod info
					ss << Bot[i].p.r << " " << round(Bot[i].p.v.x) << " " << round(Bot[i].p.v.y) << " " << round(Bot[i].p.angle) << " " << Bot[i].p.next << endl;
				}
				for(int j=0;j<2;++j){//Feed fake far away opponent position
					ss << -10000 << " " << -10000 << " " << 0 << " " << 0 << " " << 0 << " " << 1 << endl;
				}
				try{
					Bot[i].Feed_Inputs(ss.str());
					string Move{GetMove(Bot[i],turn)};
					string err_str{EmptyPipe(Bot[i].errPipe)};
					if(Debug_AI){
						ofstream err_out("log.txt",ios::app);
						err_out << err_str << endl;
					}
					Make_Move(Move,Bot[i],C,turn);
					if(Bot[i].p.lap==N_L){
						return turn;
					}
					else if(Bot[i].p.time==0){
						return -2;
					}
				}
				catch(...){
					cerr << "Loss by timeout" << endl;
					return -1;
				}
			}
			else{
				return -2;
			}
		}
	}
	return -2;
}

void StopArena(const int signum){
    stop=true;
}

int main(int argc,char **argv){
	if(argc<2){
		cerr << "Program takes 1 input, the names of the AI to test" << endl;
		return 0;
	}
	int N_Threads{1};
	if(argc>=3){//Optional N_Threads parameter
		N_Threads=min(2*omp_get_num_procs(),max(1,atoi(argv[2])));
		cerr << "Running " << N_Threads << " arena threads" << endl;
	}
	array<string,N> Bot_Names;
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
	signal(SIGTERM,StopArena);//Register SIGTERM signal handler so the arena can cleanup when you kill it
    signal(SIGPIPE,SIG_IGN);//Ignore SIGPIPE to avoid the arena crashing when an AI crashes
	int maps_done{0},maps_failed{0};
	double average_turns{0};
	#pragma omp parallel num_threads(N_Threads) shared(average_turns,maps_done,Bot_Names)
	while(!stop){
		int turns{Play_Game(Bot_Names)};
		if(turns>0){
			#pragma omp atomic
			average_turns+=(turns-average_turns)/(++maps_done);
		}
		else{
            ++maps_failed;
			//cerr << "Bot didn't complete map" << endl;
		}
		#pragma omp critical
		cerr << "Average turns: " << setprecision(5) << average_turns << " with " << maps_done << " maps done with " << setprecision(2) << 100.0*maps_failed/(maps_failed+maps_done) << "% maps failed" << endl;
	}
}